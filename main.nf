#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
nextflow.enable.strict = true
nextflow.preview.output = true

// Default parameter values
params.help = false

params.samples = ""
params.multiqc_config = "${projectDir}/multiqc_config.yaml"

params.downsample = params.downsample ?: [:]
params.downsample.enabled = params.downsample.enabled ?: false
params.downsample.sample_size = params.downsample.sample_size ?: 1000
params.downsample.seed = params.downsample.seed ?: 1

params.merge = params.merge ?: [:]
params.merge.enabled = params.merge.enabled ?: false
params.merge.pear_args = params.merge.pear_args ?: ""

params.trim = params.trim ?: [:]
params.trim.enabled = params.trim.enabled ?: false
params.trim.cutadapt_args = params.trim.cutadapt_args ?: ""

params.quantify = params.quantify ?: [:]
params.quantify.enabled = params.quantify.enabled ?: false
params.quantify.libspec = params.quantify.libspec ?: ""
params.quantify.dnacomb_args = params.quantify.dnacomb_args ?:  ""

// Include processes
include {
    fastqc as fastqc_input; fastqc as fastqc_downsampled; fastqc as fastqc_filtered;
    fastqc as fastqc_merged; fastqc as fastqc_unmerged; fastqc as fastqc_discarded;
    fastqc as fastqc_trimmed; fastqc as fastqc_untrimmed;
    multiqc;
    seqkit_stats;
    seqtk;
    pear;
    cutadapt;
    dnacomb;
} from './src/pipelines.nf'

process qc_counts {
    cpus 1
    memory '8 GB'
    queue 'normal'

    input:
    path dnacomb_output
    path library
    val roots

    output:
    path "count_qc.pdf", emit: counts

    script:
    """
    qc_counts.R --library ${library} --roots ${roots}
    """

    stub:
    """
    touch count_qc.pdf
    """
}

workflow {
    main:
    // Print help
    if ( params.help ) {
        help = """
Nextflow pipeline processing generic sequencing read data in a configurable style. Allows you to do standard processing steps (SeqKit filtering, SeqTK downsampling, PEAR merging, CutAdapt trimming, DNAComb counting) with options to include/exclude and configure each stage and automatic FastQC/MultiQC output. Written in a simple portable style in one file to make it each to copy and customise for other projects.

    Call signature:
        nextflow run bin/process_reads.nf -profile lsf -config meta/kite_hap1_pe.config

    Arguments:
        --help       : Print this message
        -config      : Path to config file of parameters (Nextflow core argument)
        -profile     : Nextflow exectution profile (Nextflow core argument)

    Required parameters (these can also be passed as --arg):
        - samples   : Path to the sample sheet CSV file (name,read1,read2)
    """.stripMargin()

        println(help)
        exit(0)
    }

    if ( params.samples == "" ) {
        error "No sample sheet specified, exiting"
        exit 1
    }

    // Load sample sheet
    // Read channels follow NF-Core with structure [meta, reads] where
    // meta is [id: str, single_end: bool]
    input_reads = Channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .map { row ->
            def read1 = file(row.read1)
            def read2 = row.read2 ? file(row.read2) : null
            def reads = [read1, read2].findAll()
            def meta = [id: row.name, single_end: read2 == null]
            tuple(meta, reads)
        }

    // Run FastQC on input reads
    fastqc_input(input_reads, channel.value("raw"))
    fastqc_zips = fastqc_input.out.zip.collect{it[1]}.ifEmpty([])
    fastqc_htmls = fastqc_input.out.html.collect{it[1]}.ifEmpty([])

    // SeqTK Downsample
    if ( params.downsample.enabled ) {
        // Run SeqTK
        seqtk(
            input_reads,
            channel.value(params.downsample.sample_size),
            channel.value(params.downsample.seed)
        )
        downsampled_reads = seqtk.out.reads
        seqtk_out = seqtk.out.reads

        // Run FastQC
        fastqc_downsampled(downsampled_reads, channel.value("downsampled"))
        fastqc_zips = fastqc_zips.mix(fastqc_downsampled.out.zip.collect{it[1]}.ifEmpty([]))
        fastqc_htmls = fastqc_htmls.mix(fastqc_downsampled.out.html.collect{it[1]}.ifEmpty([]))
    } else {
        downsampled_reads = input_reads
        seqtk_out = Channel.empty()
    }

    // PEAR Merging
    if ( params.merge.enabled ) {
        // Run PEAR
        pear(downsampled_reads, channel.value(params.merge.pear_args))
        merged_reads = pear.out.assembled.map{ meta, reads -> tuple([id: meta.id, single_end: true], reads)}
        pear_discarded = pear.out.discarded.map{ meta, reads -> tuple([id: meta.id, single_end: true], reads)}

        pear_out = pear.out.assembled.mix(pear.out.unassembled, pear.out.discarded)

        // Run FastQC
        fastqc_merged(merged_reads, channel.value("merged"))
        fastqc_unmerged(pear.out.unassembled, channel.value("unmerged_unassembled"))
        fastqc_discarded(pear_discarded, channel.value("unmerged_discarded"))

        fastqc_zips = fastqc_zips.mix(
            fastqc_merged.out.zip.collect{it[1]}.ifEmpty([]),
            fastqc_unmerged.out.zip.collect{it[1]}.ifEmpty([]),
            fastqc_discarded.out.zip.collect{it[1]}.ifEmpty([]),
        )
        fastqc_htmls = fastqc_htmls.mix(
            fastqc_merged.out.html.collect{it[1]}.ifEmpty([]),
            fastqc_unmerged.out.html.collect{it[1]}.ifEmpty([]),
            fastqc_discarded.out.html.collect{it[1]}.ifEmpty([]),
        )
    } else {
        merged_reads = downsampled_reads
        pear_out = Channel.empty()
    }

    // CutAdapt Trimming
    if ( params.trim.enabled ) {
        // Run CutAdapt
        cutadapt(merged_reads, channel.value(params.trim.cutadapt_args))
        trimmed_reads = cutadapt.out.reads

        cutadapt_out = cutadapt.out.reads.mix(
            cutadapt.out.untrimmed_reads, cutadapt.out.json, cutadapt.out.log
        )

        // Run FastQC
        fastqc_trimmed(cutadapt.out.reads, channel.value("trimmed"))
        fastqc_untrimmed(cutadapt.out.untrimmed_reads, channel.value("untrimmed"))
        fastqc_zips = fastqc_zips.mix(
            fastqc_trimmed.out.zip.collect{it[1]}.ifEmpty([]),
            fastqc_untrimmed.out.zip.collect{it[1]}.ifEmpty([])
        )
        fastqc_htmls = fastqc_htmls.mix(
            fastqc_trimmed.out.html.collect{it[1]}.ifEmpty([]),
            fastqc_untrimmed.out.html.collect{it[1]}.ifEmpty([])
        )
    } else {
        trimmed_reads = merged_reads
        cutadapt_out = Channel.empty()
    }

    // DNAComb Quantification
    if ( params.quantify.enabled ) {
        // Read and parse libspec JSON
        def libspec_path = file(params.quantify.libspec)
        def libspec_json = new groovy.json.JsonSlurper().parse(libspec_path)
        def library_path = file(libspec_json.library)

        // Optionally check existence
        if (!library_path.exists()) {
            error "Library file '${libspec_json.library}' (from ${libspec_path}) not found"
        }

        dnacomb(
            trimmed_reads,
            channel.fromPath(libspec_path).first(),
            channel.fromPath(library_path).first(),
            channel.value(params.quantify.dnacomb_args)
        )

        dnacomb_out = dnacomb.out.counts.mix(
            dnacomb.out.library_counts, dnacomb.out.log, dnacomb.out.summary
        )

        qc_in = dnacomb_out.map{x->x[1]}.collect()
        roots = dnacomb.out.counts.map{x->x[1].baseName.replaceFirst(".counts", "")}.reduce{a,b-> a + " " + b}
        qc_counts(qc_in, channel.fromPath(library_path).first(), roots)
    } else {
        dnacomb_out = Channel.empty()
    }

    // MultiQC QC Summary
    multiqc(
        fastqc_zips.mix(cutadapt.out.json.collect{it[1]}.ifEmpty([])).collect(),
        channel.fromPath( params.multiqc_config ).first()
    )

    // SeqKit Stats summary of all files
    seqkit_stat_files = input_reads.mix(
        downsampled_reads,
        merged_reads,
        trimmed_reads
    ).flatMap { i -> i[1] }
     .unique()
     .collect()

    seqkit_stats(seqkit_stat_files)

    publish:
    fastqc = fastqc_htmls.mix(fastqc_zips)
    multiqc = multiqc.out.report.mix(multiqc.out.zip, multiqc.out.data, multiqc.out.plots)
    seqkit_stats = seqkit_stats.out.tsv
    seqtk = seqtk_out
    pear = pear_out
    cutadapt = cutadapt_out
    dnacomb = dnacomb_out
    count_qc = qc_counts.out
}

output {
    fastqc {
        path "fastqc"
    }

    multiqc {
        path "."
    }

    seqkit_stats {
        path "."
    }

    seqtk {
        path "downsampled"
    }

    pear {
        path "merged"
    }

    cutadapt {
        path "trimmed"
    }

    dnacomb {
        path "counts"
    }

    count_qc {
        path "counts"
    }
}

