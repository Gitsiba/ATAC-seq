# ATAC-seq Pipeline

Modular and reproducible pipeline for **ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing)** analysis.  
This workflow takes raw FASTQ files through preprocessing, alignment, peak calling, QC, and visualization, producing publication‑ready outputs.

---
---

## 🚀 Features
- Preprocessing of raw FASTQ files
- Alignment to **hg38** reference genome
- Peak calling with MACS2
- QC metrics (duplication, insert size, peak length distribution)
- Publication‑ready plots
- Modular scripts in Python and R
- Reproducible environment via Conda

---

## ⚙️ Installation
Clone the repository:
```bash
git clone https://github.com/Gitsiba/ATAC-seq.git
cd ATAC-seq
conda env create -f env/environment.yml
conda activate atacseq
Usage

Run the pipeline step by step or via batch runner:
# Example: alignment + peak calling
python scripts/07_qc_metrics_and_plots.py --input data/sam_R1.fastq --output results/
Rscript scripts/annotate_peaks.R results/peaks/sample_peaks.narrowPeak
