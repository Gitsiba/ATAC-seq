library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

peak <- readPeakFile("results/peaks/sample_peaks.narrowPeak")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
anno <- annotatePeak(peak, TxDb=txdb, tssRegion=c(-3000, 3000))
plotAnnoPie(anno)
R
