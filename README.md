# ATAC-seq analysis toolkit

A small, modular repository for exploring **ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing)** analysis. The project currently contains plotting and R analysis scripts plus a Conda environment specification. It is a developing toolkit—not yet a fully automated end-to-end pipeline.

[![Smoke tests](https://github.com/Gitsiba/ATAC-seq/actions/workflows/smoke-tests.yml/badge.svg)](https://github.com/Gitsiba/ATAC-seq/actions/workflows/smoke-tests.yml)

## Project status

| Component | Current status |
|---|---|
| Environment | Conda environment specification is provided in `environment.yml` |
| QC plots | `07_qc_metrics_and_plots.py` creates insert-size and peak-length plots from supplied data; demo data are synthetic and must not be interpreted as biological results |
| Peak annotation | `annotate_peaks.R` annotates a MACS2 narrowPeak file against the hg38 known-gene annotation |
| Differential accessibility | `Diffrential_accessibility.R` analyzes a user-supplied peak count matrix and sample metadata using DESeq2 |
| FASTQ QC, trimming, alignment, peak calling | Tools are listed in the environment, but a fully orchestrated end-to-end workflow is not yet implemented in this repository |

## Workflow overview

The intended analysis stages are:

1. Inspect raw reads with FastQC and MultiQC (MultiQC report generation is not currently automated here).
2. Trim adapters/low-quality bases with fastp or an explicitly documented alternative.
3. Align reads to a matching reference genome with Bowtie 2.
4. Sort/index alignments and assess mapping/duplication metrics with SAMtools and Picard.
5. Call accessible regions with MACS2 using parameters appropriate to the library layout and experimental design.
6. Review QC metrics and annotate peaks.
7. For replicated experiments, build a consistent peak count matrix and test differential accessibility with a suitable statistical design.

**Important:** genome build, chromosome naming, paired/single-end layout, controls, biological replicates, blacklist filtering, duplicate policy, and MACS2 settings must be selected for the actual experiment. This repository does not currently automate all of these decisions.

## Quick start

### 1. Clone and create the environment

```bash
git clone https://github.com/Gitsiba/ATAC-seq.git
cd ATAC-seq
conda env create -f environment.yml
conda activate atac
```

### 2. Generate QC plots

The plotter accepts a tab-separated file containing one fragment length per row and a standard MACS2 narrowPeak file.

```bash
python 07_qc_metrics_and_plots.py \
  --insert-sizes examples/insert_sizes.tsv \
  --peaks examples/peaks.narrowPeak \
  --outdir results/qc
```

The insert-size file should contain a single numeric column (a header is optional). The peak file must be BED/narrowPeak-like, with chromosome, start, and end in columns 1–3. Outputs are saved as PNG files in the chosen directory.

To generate clearly labeled synthetic demonstration plots instead:

```bash
python 07_qc_metrics_and_plots.py --demo --outdir results/demo
```

Synthetic demo plots are software examples only; do not report them as experimental QC.

### 3. Annotate peaks (R)

Install the required Bioconductor packages in an R environment, then provide a real MACS2 narrowPeak file:

```bash
Rscript annotate_peaks.R results/peaks/sample_peaks.narrowPeak results/annotation
```

The script writes a tabular annotation and a pie plot to the output directory. This annotation uses the hg38 known-gene database; use a matching genome build and annotation for other organisms/builds.

### 4. Differential accessibility (R / DESeq2)

Prepare a tab-separated **raw integer peak-count matrix**: rows are consensus peaks, columns are samples, and the first column is a unique peak identifier. Prepare a metadata TSV with sample IDs and a condition column; row names must match the count matrix columns.

```bash
Rscript Diffrential_accessibility.R counts.tsv metadata.tsv results/differential_accessibility
```

Use biological replicates and a design appropriate to the experiment. DESeq2 should receive raw counts, not normalized values, TPM, or FPKM. Do not use this analysis for a single sample per condition without an appropriate statistical design and replication.

## Inputs and outputs

- **Inputs:** FASTQ files for upstream analysis; narrowPeak/BED files for peak-level analysis; raw peak counts and sample metadata for differential accessibility.
- **Outputs:** QC plots, peak annotation table/plot, and DESeq2 result table/MA plot.
- **Reference:** the examples above assume hg38 for peak annotation. The count matrix must be based on a consistent consensus peak set across samples.

## Reproducibility checklist

- Record reference genome build, annotation release, tool versions, and command-line parameters.
- Keep sample IDs consistent between count matrices and metadata.
- Preserve raw counts for differential testing.
- Inspect alignment, library complexity, fragment-size distribution, TSS enrichment, FRiP, and replicate concordance where applicable; these metrics are not all computed by the current scripts.
- Keep private or unpublished sequencing data out of public repositories.
- Store full FASTQ/BAM files outside Git; commit only small, shareable test fixtures and documentation.

## Repository layout

```text
.
├── 07_qc_metrics_and_plots.py
├── annotate_peaks.R
├── Diffrential_accessibility.R
├── environment.yml
├── examples/
└── tests/
```

## Testing

Run the lightweight tests with:

```bash
python -m unittest discover -s tests -v
```

The tests validate script behavior on tiny synthetic fixtures; they do not establish biological validity or replace validation on real sequencing data.

## References

- [ENCODE ATAC-seq standards and processing overview](https://www.encodeproject.org/atac-seq/)
- [ENCODE-DCC ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline)
- [nf-core/atacseq](https://github.com/nf-core/atacseq)

These established pipelines are useful references for a future fully automated, production-style workflow. This repository is an independent learning/toolkit project and should not be represented as ENCODE-validated.

## License

See [LICENSE](LICENSE).
