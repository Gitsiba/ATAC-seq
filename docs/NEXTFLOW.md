# Nextflow ATAC-seq workflow

This is a paired-end, single-sample peak-calling starter workflow for demonstrating reproducible workflow engineering. It has **not** been validated as a production clinical or ENCODE-equivalent pipeline.

## Workflow stages

1. FastQC on raw paired-end reads.
2. fastp adapter and quality trimming.
3. FastQC on trimmed reads.
4. Bowtie 2 index construction from the supplied reference FASTA.
5. Paired-end alignment with Bowtie 2.
6. SAMtools name sorting, fixmate, coordinate sorting, duplicate removal, BAM indexing, and flagstat.
7. MACS2 peak calling in paired-end mode (`BAMPE`).
8. MultiQC aggregation of available FastQC, fastp, and flagstat reports.

Nextflow execution timeline, report, trace, and DAG are enabled in `nextflow.config`.

## Requirements

- Nextflow 24.04 or newer and a compatible Java runtime.
- Conda/Mamba for `-profile conda` mode.
- Paired-end FASTQ files and a reference FASTA; use the same genome build for all samples.
- Adequate disk space for reference indexes and intermediate BAMs.

## Input sheet

Copy `samplesheet.example.csv` and replace paths. The required CSV header is exactly:

```csv
sample,fastq_1,fastq_2
sample_A,data/sample_A_R1.fastq.gz,data/sample_A_R2.fastq.gz
```

Paths are resolved relative to the directory from which Nextflow is launched. Sample IDs must contain letters, numbers, dot, underscore, or hyphen. Use paired-end reads; single-end libraries are not supported.

## Run

```bash
nextflow run main.nf \
  -profile conda \
  --samplesheet samplesheet.csv \
  --reference_fasta /path/to/reference.fa \
  --genome_size hs \
  --outdir results
```

Use `--genome_size mm` for the MACS2 mouse preset, or provide an appropriate effective genome size for another organism. The FASTA must match the genome size and any annotation used downstream.

## Main outputs

- Raw and trimmed FastQC HTML/ZIP reports.
- fastp HTML and JSON reports.
- Deduplicated coordinate-sorted BAM and BAI files.
- SAMtools flagstat summaries.
- MACS2 narrowPeak, summits, and peak-table outputs.
- MultiQC HTML report and supporting data.
- Nextflow execution report, timeline, trace, and DAG under `results/pipeline_info/`.

Nextflow work directories are retained for resume support. Do not commit private sequencing data or large FASTQ/BAM outputs.

## CI validation

GitHub Actions runs `nextflow lint` against `main.nf` and `nextflow.config`, then checks that Nextflow can resolve the configuration. This validates syntax/configuration only; it does **not** prove the workflow completes on sequencing data. Check the Actions tab and run a small end-to-end test before treating the pipeline as validated.

## Limitations and interpretation

- Calls peaks separately per sample; it does not create a consensus peak set or perform differential accessibility.
- Does not calculate TSS enrichment or FRiP, apply a genome blacklist, or automatically choose a reference annotation.
- Filters for proper pairs with SAMtools, applies a configurable MAPQ threshold (`--min_mapq`, default 30), excludes unmapped/secondary/QC-fail/duplicate-flagged/supplementary records, then removes duplicates with `samtools markdup -r`. Confirm this policy is suitable for the library and study design.
- MACS2 uses `--keep-dup all` because duplicates have already been removed upstream.
- Does not accept a control/input BAM.
- The default genome size is `hs`; explicitly set it for the target organism.
- Must be tested on a small dataset and reviewed biologically before drawing research conclusions.

## Next engineering milestones

1. Add a tiny end-to-end fixture with a small reference and assert expected outputs.
2. Add blacklist filtering and tested TSS/FRiP metrics.
3. Add replicate metadata, consensus peaks, and a validated differential-accessibility module.
4. Add pinned container images, parameter schema, and software bill of materials.
5. Compare outputs and parameters against [nf-core/atacseq](https://github.com/nf-core/atacseq) and the [ENCODE ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline).

## References

- [Nextflow documentation](https://www.nextflow.io/docs/latest/)
- [nf-core/atacseq](https://github.com/nf-core/atacseq)
- [ENCODE-DCC ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline)
