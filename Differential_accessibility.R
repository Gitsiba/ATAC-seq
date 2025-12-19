library(DESeq2)

# Example: count matrix from peaks
counts <- matrix(rpois(200, lambda=20), ncol=2)
colnames(counts) <- c("ConditionA", "ConditionB")
dds <- DESeqDataSetFromMatrix(countData=counts,
                              colData=data.frame(condition=c("A","B")),
                              design=~condition)
dds <- DESeq(dds)
res <- results(dds)
plotMA(res)
