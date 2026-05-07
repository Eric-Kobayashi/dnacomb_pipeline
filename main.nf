#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
nextflow.enable.strict = true

// Default parameter values
params.help = false

params.samples = ""
params.multiqc_config = "${projectDir}/multiqc_config.yaml"
params.dnacomb_qc_rmd = "${projectDir}/bin/dnacomb_qc_report.Rmd"

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
params.quantify.library = params.quantify.library ?: []
params.quantify.dnacomb_args = params.quantify.dnacomb_args ?:  ""

params.qc = params.qc ?: [:]
params.qc.seqkit = params.qc.seqkit ?: false

// Include processes
include {
    fastqc as fastqc_input; fastqc as fastqc_downsampled;
    fastqc as fastqc_merged; fastqc as fastqc_unmerged; fastqc as fastqc_discarded;
    fastqc as fastqc_trimmed; fastqc as fastqc_untrimmed;
    multiqc;
    seqkit_stats;
    seqtk;
    pear;
    cutadapt;
    dnacomb;
    qc_counts;
    count_records;
    combine_record_counts;
} from './src/pipelines.nf'

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
    input_reads = channel
        .fromPath(params.samples)
        .splitCsv(header: true)
        .map { row ->
            def read1 = file(row.read1)
            def read2 = row.read2 ? file(row.read2) : null
            def reads = [read1, read2].findAll()
            def meta = [id: row.name, single_end: read2 == null, label: "raw", stage: 0]
            tuple(meta, reads)
        }

    // Run FastQC on input reads
    fastqc_input(input_reads)
    fastqc_zips = fastqc_input.out.zip.collect{x->x[1]}.ifEmpty([])
    fastqc_htmls = fastqc_input.out.html.collect{x->x[1]}.ifEmpty([])

    // Establish channel for record counting
    seq_files = channel.empty().mix(input_reads)

    // SeqTK Downsample
    if ( params.downsample.enabled ) {
        // Run SeqTK
        seqtk(
            input_reads,
            channel.value(params.downsample.sample_size),
            channel.value(params.downsample.seed)
        )
        seqtk_out = seqtk.out.reads.map { meta, reads ->
            tuple(meta + [label: "downsampled", stage: meta.stage + 1], reads)
        }

        downsampled_reads = seqtk_out

        seq_files = seq_files.mix(seqtk_out)

        // Run FastQC
        fastqc_downsampled(downsampled_reads)
        fastqc_zips = fastqc_zips.mix(fastqc_downsampled.out.zip.collect{x->x[1]}.ifEmpty([]))
        fastqc_htmls = fastqc_htmls.mix(fastqc_downsampled.out.html.collect{x->x[1]}.ifEmpty([]))
    } else {
        downsampled_reads = input_reads
        seqtk_out = channel.empty()
    }

    // PEAR Merging
    if ( params.merge.enabled ) {
        // Run PEAR
        pear(downsampled_reads, channel.value(params.merge.pear_args))

        merged_reads = pear.out.assembled.map { meta, reads ->
            tuple(meta + [label: "merged", stage: meta.stage + 1, single_end: true], reads)
        }

        unmerged_reads = pear.out.unassembled.map { meta, reads ->
            tuple(meta + [label: "unmerged", stage: meta.stage + 1], reads)
        }

        discarded_reads = pear.out.discarded.map { meta, reads ->
            tuple(meta + [label: "discarded", stage: meta.stage + 1, single_end: true], reads)
        }

        pear_out = merged_reads.mix(unmerged_reads, discarded_reads)

        seq_files = seq_files.mix(pear_out)

        // Run FastQC
        fastqc_merged(merged_reads)
        fastqc_unmerged(unmerged_reads)
        fastqc_discarded(discarded_reads)

        fastqc_zips = fastqc_zips.mix(
            fastqc_merged.out.zip.collect{x->x[1]}.ifEmpty([]),
            fastqc_unmerged.out.zip.collect{x->x[1]}.ifEmpty([]),
            fastqc_discarded.out.zip.collect{x->x[1]}.ifEmpty([]),
        )
        fastqc_htmls = fastqc_htmls.mix(
            fastqc_merged.out.html.collect{x->x[1]}.ifEmpty([]),
            fastqc_unmerged.out.html.collect{x->x[1]}.ifEmpty([]),
            fastqc_discarded.out.html.collect{x->x[1]}.ifEmpty([]),
        )
    } else {
        merged_reads = downsampled_reads
        pear_out = channel.empty()
    }

    // CutAdapt Trimming
    if ( params.trim.enabled ) {
        // Run CutAdapt
        cutadapt(merged_reads, channel.value(params.trim.cutadapt_args))

        trimmed_reads = cutadapt.out.reads.map { meta, reads ->
            tuple(meta + [label: "trimmed", stage: meta.stage + 1], reads)
        }

        untrimmed_reads = cutadapt.out.untrimmed_reads.map { meta, reads ->
            tuple(meta + [label: "untrimmed", stage: meta.stage + 1], reads)
        }

        trim_json = cutadapt.out.json.map { meta, json ->
            tuple(meta + [label: "trimmed", stage: meta.stage + 1], json)
        }

        trim_log = cutadapt.out.log.map { meta, log ->
            tuple(meta + [label: "trimmed", stage: meta.stage + 1], log)
        }

        cutadapt_out = trimmed_reads.mix(
            untrimmed_reads, trim_json, trim_log
        )

        seq_files = seq_files.mix(trimmed_reads, untrimmed_reads)

        // Run FastQC
        fastqc_trimmed(trimmed_reads)
        fastqc_untrimmed(untrimmed_reads)
        fastqc_zips = fastqc_zips.mix(
            fastqc_trimmed.out.zip.collect{x->x[1]}.ifEmpty([]),
            fastqc_untrimmed.out.zip.collect{x->x[1]}.ifEmpty([])
        )
        fastqc_htmls = fastqc_htmls.mix(
            fastqc_trimmed.out.html.collect{x->x[1]}.ifEmpty([]),
            fastqc_untrimmed.out.html.collect{x->x[1]}.ifEmpty([])
        )
    } else {
        trimmed_reads = merged_reads
        cutadapt_out = channel.empty()
        trim_json = channel.empty()
    }

    // Count records in each sequence file produced
    count_records(seq_files)
    combine_record_counts(count_records.out.rows.collect())

    // DNAComb Quantification
    if ( params.quantify.enabled ) {
        // Read and parse libspec JSON
        def libspec_path = file(params.quantify.libspec)

        // Collect and validate all library files
        def library_files = params.quantify.library.collect { lib -> file(lib) }
        library_files.each { lib ->
            if (!lib.exists()) {
                error "Library file '${lib}' not found"
            }
        }

        dnacomb(
            trimmed_reads,
            channel.fromPath(libspec_path).first(),
            channel.fromPath(library_files).collect(),
            channel.value(params.quantify.dnacomb_args)
        )

        dnacomb_out = dnacomb.out.counts.mix(
            dnacomb.out.library_counts, dnacomb.out.log, dnacomb.out.summary, dnacomb.out.filtered
        )

        qc_in = dnacomb_out.map{x->x[1]}.collect()
        roots = dnacomb.out.counts.map{x->x[1].baseName.replaceFirst(/\.counts/, "")}.reduce{a,b-> a + " " + b}

        qc_counts(
            qc_in,
            channel.fromPath(library_files).collect(),
            combine_record_counts.out.tsv.collect(),
            channel.fromPath( params.dnacomb_qc_rmd ).first(),
            roots
        )

        qc_counts_out = qc_counts.out.qc_html
    } else {
        dnacomb_out = channel.empty()
        qc_counts_out = channel.empty()
    }

    // MultiQC QC Summary
    multiqc(
        fastqc_zips.mix(trim_json.collect{x->x[1]}.ifEmpty([])).collect(),
        channel.fromPath( params.multiqc_config ).first()
    )

    if ( params.qc.seqkit ) {
        // SeqKit Stats summary of all files
        seqkit_stat_files = input_reads.mix(
            downsampled_reads,
            merged_reads,
            trimmed_reads
        ).flatMap { i -> i[1] }
        .unique()
        .collect()

        seqkit_stats(seqkit_stat_files)
        seqkit_out = seqkit_stats.out.tsv
    } else {
        seqkit_out = channel.empty()
    }


    publish:
    fastqc = fastqc_htmls.mix(fastqc_zips)
    multiqc = multiqc.out.report.mix(multiqc.out.zip, multiqc.out.data, multiqc.out.plots)
    seqkit_stats = seqkit_out
    seqtk = seqtk_out
    pear = pear_out
    cutadapt = cutadapt_out
    dnacomb = dnacomb_out
    count_qc = qc_counts_out
    record_counts = combine_record_counts.out.tsv
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

    record_counts {
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
        path "."
    }
}

