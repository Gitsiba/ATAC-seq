#!/usr/bin/env Rscript
# Differential accessibility from a raw integer peak-count matrix and sample metadata.
# Usage: Rscript Diffrential_accessibility.R counts.tsv metadata.tsv outdir
suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript Diffrential_accessibility.R counts.tsv metadata.tsv outdir", call. = FALSE)
}
counts_path <- args[[1]]
metadata_path <- args[[2]]
outdir <- args[[3]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

counts_df <- read.delim(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(counts_df) < 3) {
  stop("Count matrix must contain a peak ID column and at least two sample columns.", call. = FALSE)
}
peak_ids <- counts_df[[1]]
if (anyDuplicated(peak_ids) || any(is.na(peak_ids)) || any(peak_ids == "")) {
  stop("First column must contain unique, non-empty peak IDs.", call. = FALSE)
}
count_matrix <- as.matrix(counts_df[, -1, drop = FALSE])
rownames(count_matrix) <- peak_ids
suppressWarnings(storage.mode(count_matrix) <- "numeric")
if (anyNA(count_matrix) || any(!is.finite(count_matrix)) || any(count_matrix < 0) ||
    any(count_matrix != floor(count_matrix))) {
  stop("Counts must be finite, non-negative integers.", call. = FALSE)
}
storage.mode(count_matrix) <- "integer"

metadata <- read.delim(metadata_path, row.names = 1, check.names = FALSE,
                       stringsAsFactors = FALSE)
if (!"condition" %in% colnames(metadata)) {
  stop("Metadata must contain a column named 'condition' and sample IDs as row names.", call. = FALSE)
}
if (!setequal(colnames(count_matrix), rownames(metadata))) {
  stop("Sample IDs in metadata must exactly match count-matrix columns.", call. = FALSE)
}
metadata <- metadata[colnames(count_matrix), , drop = FALSE]
metadata$condition <- factor(metadata$condition)
if (nlevels(metadata$condition) < 2) {
  stop("At least two condition levels are required.", call. = FALSE)
}
if (any(table(metadata$condition) < 2)) {
  warning("At least one condition has fewer than two samples; inference may be unreliable.")
}
if (nlevels(metadata$condition) != 2) {
  stop("This starter script currently supports exactly two conditions.", call. = FALSE)
}
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = metadata, design = ~ condition)
dds <- dds[rowSums(counts(dds)) > 0, ]
if (nrow(dds) == 0) stop("No peaks with non-zero counts remain.", call. = FALSE)
dds <- DESeq(dds)
reference_level <- levels(metadata$condition)[1]
contrast_level <- levels(metadata$condition)[2]
res <- results(dds, contrast = c("condition", contrast_level, reference_level))
res_df <- as.data.frame(res)
res_df$peak_id <- rownames(res_df)
res_df <- res_df[, c("peak_id", setdiff(colnames(res_df), "peak_id"))]
write.table(res_df, file = file.path(outdir, "differential_accessibility.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
png(file.path(outdir, "ma_plot.png"), width = 1200, height = 900, res = 150)
plotMA(res, main = paste(contrast_level, "vs", reference_level), ylim = c(-5, 5))
dev.off()
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
