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
conda env create -f environment.yml
conda activate atac
```

---

## 📋 Usage

Run the pipeline step by step or via batch runner:
```bash
# Example: QC metrics and plots
python 07_qc_metrics_and_plots.py

# Example: annotate peaks
Rscript annotate_peaks.R

# Example: differential accessibility analysis
Rscript Differential_accessibility.R
```
