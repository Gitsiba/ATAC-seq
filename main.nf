#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.samplesheet = null
params.reference_fasta = null
params.genome_size = 'hs'
params.min_mapq = 30
params.outdir = 'results'

process FASTQC_RAW {
    tag "${meta.id}"
    cpus 2
    conda 'bioconda::fastqc=0.12.1'
    publishDir "${params.outdir}/qc/raw", mode: 'copy'
    input:
    tuple val(meta), path(read1), path(read2)
    output:
    path("*_fastqc.*"), emit: reports
    script:
    """
    fastqc --threads ${task.cpus} ${read1} ${read2}
    """
}

process TRIM_READS {
    tag "${meta.id}"
    cpus 2
    conda 'bioconda::fastp=0.23.4'
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    input:
    tuple val(meta), path(read1), path(read2)
    output:
    tuple val(meta), path("${meta.id}_R1.trimmed.fastq.gz"), path("${meta.id}_R2.trimmed.fastq.gz"), emit: reads
    tuple path("${meta.id}.fastp.html"), path("${meta.id}.fastp.json"), emit: qc
    script:
    """
    fastp --thread ${task.cpus} \
      --in1 ${read1} --in2 ${read2} \
      --out1 ${meta.id}_R1.trimmed.fastq.gz \
      --out2 ${meta.id}_R2.trimmed.fastq.gz \
      --html ${meta.id}.fastp.html \
      --json ${meta.id}.fastp.json
    """
}

process FASTQC_TRIMMED {
    tag "${meta.id}"
    cpus 2
    conda 'bioconda::fastqc=0.12.1'
    publishDir "${params.outdir}/qc/trimmed", mode: 'copy'
    input:
    tuple val(meta), path(read1), path(read2)
    output:
    path("*_fastqc.*"), emit: reports
    script:
    """
    fastqc --threads ${task.cpus} ${read1} ${read2}
    """
}

process BUILD_INDEX {
    tag 'reference'
    cpus 4
    conda 'bioconda::bowtie2=2.5.4'
    publishDir "${params.outdir}/reference_index", mode: 'copy'
    input:
    path reference_fasta
    output:
    path("reference.*.bt2*"), emit: index_files
    script:
    """
    bowtie2-build --threads ${task.cpus} ${reference_fasta} reference
    """
}

process ALIGN_AND_DEDUP {
    tag "${meta.id}"
    cpus 4
    conda 'bioconda::bowtie2=2.5.4 bioconda::samtools=1.20'
    publishDir "${params.outdir}/alignment", mode: 'copy'
    input:
    tuple val(meta), path(read1), path(read2)
    path index_files
    output:
    tuple val(meta), path("${meta.id}.dedup.bam"), path("${meta.id}.dedup.bam.bai"), emit: bam
    path("${meta.id}.flagstat.txt"), emit: flagstat
    script:
    """
    bowtie2 --very-sensitive --no-mixed --no-discordant -X 2000 \
      -x reference -1 ${read1} -2 ${read2} -p ${task.cpus} 2> ${meta.id}.bowtie2.log \
      | samtools view -u -f 2 -F 3844 -q ${params.min_mapq} - \
      | samtools sort -n -@ ${task.cpus} -o ${meta.id}.name.bam -

    samtools fixmate -m ${meta.id}.name.bam ${meta.id}.fixmate.bam
    samtools sort -@ ${task.cpus} -o ${meta.id}.coordinate.bam ${meta.id}.fixmate.bam
    samtools markdup -r -@ ${task.cpus} ${meta.id}.coordinate.bam ${meta.id}.dedup.bam
    samtools index -@ ${task.cpus} ${meta.id}.dedup.bam
    samtools flagstat ${meta.id}.dedup.bam > ${meta.id}.flagstat.txt
    """
}

process CALL_PEAKS {
    tag "${meta.id}"
    cpus 2
    conda 'bioconda::macs2=2.2.9.1'
    publishDir "${params.outdir}/peaks", mode: 'copy'
    input:
    tuple val(meta), path(bam), path(bai)
    output:
    path("${meta.id}_peaks.narrowPeak"), emit: narrowpeak
    path("${meta.id}_summits.bed"), emit: summits
    path("${meta.id}_peaks.xls"), emit: peak_table
    script:
    """
    macs2 callpeak -t ${bam} -f BAMPE -g '${params.genome_size}' \
      --keep-dup all -n '${meta.id}' --outdir . 2> ${meta.id}.macs2.log
    """
}

process MULTIQC {
    tag 'multiqc'
    cpus 1
    conda 'bioconda::multiqc=1.25.2'
    publishDir "${params.outdir}/qc/multiqc", mode: 'copy'
    input:
    path reports
    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_report_data"), emit: data
    script:
    """
    multiqc . --force --filename multiqc_report.html
    """
}

workflow {
    if (!params.samplesheet) error "Missing --samplesheet (CSV columns: sample,fastq_1,fastq_2)"
    if (!params.reference_fasta) error "Missing --reference_fasta (reference FASTA)"
    if (!file(params.samplesheet).exists()) error "Samplesheet not found: ${params.samplesheet}"
    if (!file(params.reference_fasta).exists()) error "Reference FASTA not found: ${params.reference_fasta}"

    read_pairs = channel.fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample || !row.fastq_1 || !row.fastq_2)
                error "Each samplesheet row must contain sample, fastq_1, and fastq_2"
            def sample = row.sample.toString().trim()
            if (!(sample ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/))
                error "Invalid sample ID '${sample}'. Use letters, numbers, dot, underscore, or hyphen."
            tuple([id: sample],
                  file(row.fastq_1.toString().trim(), checkIfExists: true),
                  file(row.fastq_2.toString().trim(), checkIfExists: true))
        }

    raw_qc = FASTQC_RAW(read_pairs)
    trimmed = TRIM_READS(read_pairs)
    trimmed_qc = FASTQC_TRIMMED(trimmed.reads)
    index = BUILD_INDEX(file(params.reference_fasta))
    aligned = ALIGN_AND_DEDUP(trimmed.reads, index.index_files)
    CALL_PEAKS(aligned.bam)

    qc_inputs = raw_qc.reports
        .flatten()
        .mix(trimmed_qc.reports.flatten())
        .mix(trimmed.qc.flatMap { item -> item })
        .mix(aligned.flagstat)
        .collect()
    MULTIQC(qc_inputs)
}
