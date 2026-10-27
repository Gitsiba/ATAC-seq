#!/usr/bin/env Rscript
# Annotate MACS2 narrowPeak intervals against the hg38 knownGene transcript database.
# Usage: Rscript annotate_peaks.R input.narrowPeak outdir
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || length(args) > 2) {
  stop("Usage: Rscript annotate_peaks.R input.narrowPeak [outdir]", call. = FALSE)
}
peak_file <- args[[1]]
outdir <- if (length(args) == 2) args[[2]] else "results/annotation"
if (!file.exists(peak_file)) stop("Peak file not found: ", peak_file, call. = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

peaks <- readPeakFile(peak_file)
if (length(peaks) == 0) stop("No peaks were read from the input file.", call. = FALSE)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
annotation <- annotatePeak(peaks, TxDb = txdb, tssRegion = c(-3000, 3000))
write.table(as.data.frame(annotation), file = file.path(outdir, "peak_annotation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
png(file.path(outdir, "peak_annotation_pie.png"), width = 1000, height = 800, res = 140)
plotAnnoPie(annotation)
dev.off()
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
