#!/usr/bin/env nextflow
// General processes and workflows to share across NF pipelines
nextflow.enable.dsl = 2

// All rules taking/emitting reads assume an input channel using the NF-Core structure.
// It contains a tuple [meta, reads] where meta is [id: str, single_end: bool] and
// reads is either [r1, r2] or [r1]

process fastqc {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)
    val(suffix)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip

    script:
    if (meta.single_end) {
        """
        [ ! -f ${meta.id}_${suffix}.fastq.gz ] && ln -s $reads ${meta.id}_${suffix}.fastq.gz
        fastqc --threads $task.cpus ${meta.id}_${suffix}.fastq.gz
        """
    } else {
        """
        [ ! -f ${meta.id}_f_${suffix}.fastq.gz ] && ln -s ${reads[0]} ${meta.id}_f_${suffix}.fastq.gz
        [ ! -f ${meta.id}__r_${suffix}.fastq.gz ] && ln -s ${reads[1]} ${meta.id}_r_${suffix}.fastq.gz
        fastqc --threads $task.cpus ${meta.id}_f_${suffix}.fastq.gz ${meta.id}_r_${suffix}.fastq.gz
        """
    }

    stub:
    if (meta.single_end) {
        """
        touch ${meta.id}_${suffix}.zip ${meta.id}_${suffix}.html
        """
    } else {
        """
        touch ${meta.id}_f_${suffix}.zip ${meta.id}_r_${suffix}.zip
        touch ${meta.id}_f_${suffix}.html ${meta.id}_r_${suffix}.html
        """
    }
}

process multiqc {
    tag "multiqc"

    input:
    path zip_files
    path config

    output:
    path "multiqc_report.html"       , emit: report
    path "multiqc_report.zip"        , emit: zip
    path "multiqc_report_data/*"     , emit: data
    path "multiqc_report_plots/*"    , optional:true, emit: plots

    script:
    """
    multiqc -c ${config} .
    zip -r multiqc_report.zip multiqc_report.html multiqc_report_data multiqc_report_plots
    """

    stub:
    """
    mkdir multiqc_report_data
    touch multiqc_report.html multiqc_report_data/multiqc_data.json multiqc_report.zip
    """
}

process seqkit_stats {
    tag "seqkit stats"

    input:
    path reads

    output:
    path "seqkit_stats.tsv", emit: tsv

    script:
    """
    seqkit stats * > seqkit_stats.tsv
    """

    stub:
    """
    touch seqkit_stats.tsv
    """
}

process seqtk {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)
    val sample_size
    val seed

    output:
    tuple val(meta), path("*_downsampled.fq.gz"), emit: reads

    script:
    if (meta.single_end) {
        """
        seqtk sample -s${seed} ${reads[0]} ${sample_size} | gzip --no-name > ${meta.id}_downsampled.fq.gz
        """
    } else {
        """
        seqtk sample -s${seed} ${reads[0]} ${sample_size} | gzip --no-name > ${meta.id}_f_downsampled.fq.gz
        seqtk sample -s${seed} ${reads[1]} ${sample_size} | gzip --no-name > ${meta.id}_r_downsampled.fq.gz
        """
    }

    stub:
    if (meta.single_end) {
        """
        touch ${meta.id}_downsampled.fq.gz
        """
    } else {
        """
        touch ${meta.id}_f_downsampled.fq.gz ${meta.id}_r_downsampled.fq.gz
        """
    }
}

process pear {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)
    val args

    output:
    tuple val(meta), path("*_merged.assembled.fastq.gz"), emit: assembled
    tuple val(meta), path("*_merged.unassembled.*.fastq.gz"), emit: unassembled
    tuple val(meta), path("*_merged.discarded.fastq.gz"), emit: discarded

    script:
    if (meta.single_end) {
        error("Trying to merge single end reads")
    } else {
        """
        gunzip -f ${reads[0]}
        gunzip -f ${reads[1]}
        pear -f ${reads[0].baseName} -r ${reads[1].baseName} -o ${meta.id}_merged -j $task.cpus $args
        gzip -f ${meta.id}_merged.assembled.fastq
        gzip -f ${meta.id}_merged.unassembled.forward.fastq
        gzip -f ${meta.id}_merged.unassembled.reverse.fastq
        gzip -f ${meta.id}_merged.discarded.fastq
        """
    }

    stub:
    if (meta.single_end) {
        error("Trying to merge single end reads")
    } else {
        """
        touch \\
           ${meta.id}_merged.assembled.fastq.gz \\
           ${meta.id}_merged.unassembled.forward.fastq.gz \\
           ${meta.id}_merged.unassembled.reverse.fastq.gz \\
           ${meta.id}_merged.discarded.fastq.gz
        """
    }
}

process cutadapt {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)
    val args

    output:
    tuple val(meta), path('*{,_f,_r}_trimmed.fastq.gz')  , emit: reads
    tuple val(meta), path('*{,_f,_r}_untrimmed.fastq.gz'), emit: untrimmed_reads, optional: true
    tuple val(meta), path('*.log')                       , emit: log
    tuple val(meta), path('*.json')                      , emit: json

    script:
    if (meta.single_end) {
        """
        cutadapt \\
          --cores $task.cpus \\
          --json=${meta.id}.cutadapt.json \\
          $args \\
          -o ${meta.id}_trimmed.fastq.gz \\
          --untrimmed-output ${meta.id}_untrimmed.fastq.gz \\
          $reads \\
        > ${meta.id}.cutadapt.log
        """
    } else {
        """
        cutadapt \\
          --cores $task.cpus \\
          --json=${meta.id}.cutadapt.json \\
          $args \\
          -o ${meta.id}_f_trimmed.fastq.gz \\
          -p ${meta.id}_r_trimmed.fastq.gz \\
          --untrimmed-output ${meta.id}_f_untrimmed.fastq.gz \\
          --untrimmed-paired-output ${meta.id}_r_untrimmed.fastq.gz \\
          $reads \\
        > ${meta.id}.cutadapt.log
        """
    }

    stub:
    if (meta.single_end) {
        """
        touch \\
          ${meta.id}_trimmed.fastq.gz \\
          ${meta.id}_untrimmed.fastq.gz \\
          ${meta.id}.cutadapt.log \\
          ${meta.id}.cutadapt.json
        """
    } else {
        """
        touch \\
          ${meta.id}_f_trimmed.fastq.gz \\
          ${meta.id}_r_trimmed.fastq.gz \\
          ${meta.id}_f_untrimmed.fastq.gz \\
          ${meta.id}_r_untrimmed.fastq.gz \\
          ${meta.id}.cutadapt.log \\
          ${meta.id}.cutadapt.json
        """
    }
}

process dnacomb {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)
    path libspec
    path library
    val args

    output:
    tuple val(meta), path("*.counts.tsv"), emit: counts
    tuple val(meta), path("*.library_counts.tsv"), emit: library_counts, optional: true
    tuple val(meta), path("*.summary.tsv"), emit: summary
    tuple val(meta), path("*.log"), emit: log

    script:
    ls = libspec ? "--library-spec nf_patched_libspec.json" : ""
    r2 = !meta.single_end ? reads[1] : ""
    cmdargs = "--verbose ${ls} --output ${meta.id} ${args}"
    """
    jq '.library = "${library}"' ${libspec} > nf_patched_libspec.json
    jq empty nf_patched_libspec.json || { echo "Patched JSON is invalid"; exit 1; }
    dnacomb ${cmdargs} ${reads[0]} ${r2} > ${meta.id}.log 2>&1
    """

    stub:
    """
    jq '.library = "${library}"' ${libspec} > nf_patched_libspec.json
    jq empty nf_patched_libspec.json || { echo "Patched JSON is invalid"; exit 1; }
    touch ${meta.id}.counts.tsv ${meta.id}.library_counts.tsv ${meta.id}.summary.tsv ${meta.id}.log
    """
}
