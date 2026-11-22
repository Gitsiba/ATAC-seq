# Reproducible tool image for the ATAC-seq Nextflow workflow.
# Build and push this image to Amazon ECR before submitting AWS Batch jobs.
FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="ATAC-seq workflow tools"
LABEL org.opencontainers.image.description="Pinned command-line tools used by Gitsiba/ATAC-seq"

RUN micromamba install --yes --name base --channel conda-forge --channel bioconda \
    fastqc=0.12.1 \
    fastp=0.23.4 \
    bowtie2=2.5.4 \
    samtools=1.20 \
    macs2=2.2.9.1 \
    multiqc=1.25.2 \
    && micromamba clean --all --yes

# The image is intended for non-interactive Nextflow tasks.
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
