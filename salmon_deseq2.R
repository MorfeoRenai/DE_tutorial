# load packages ----
library(DESeq2)
library(tximport)
library(GenomicFeatures)
library(readr)

# import sample Metadata ----
sampleMetadata <- read.csv("sampleMetadata.csv", row.names = "Rownames")
sampleMetadata

# quant.sf file paths ----
files <- file.path("salmon_results", sampleMetadata$sample, "quant.sf")
files
names(files) <- sampleMetadata$sample

# create tx2gene object ----
txdb <- GenomicFeatures::makeTxDbFromGFF("Homo_sapiens.GRCh38.100.chr.gtf")
k <- keys(txdb, keytype = "GENEID")
tx2gene <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")

# reorder columns of tx2gene ----
tx2gene <- tx2gene[, c("TXNAME", "GENEID")]

# check tx2gene ----
head(tx2gene)

# import salmon quantification data ----
txi.salmon <- tximport(files = files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
# check imported data
head(txi.salmon$counts)
# check if sample names match between objects
identical(x = rownames(sampleMetadata), y = colnames(txi.salmon$counts))

# DESeq2 pipeline ----
dds <- DESeqDataSetFromTximport(txi = txi.salmon, colData = sampleMetadata, design = ~group)
dds$group <- relevel(x = dds$group, ref = "Mock")
dds <- DESeq(dds)

# save the results of DE analysis ----
res <- results(dds)

# look at the results ----
summary(res)

# filter significant genes ----
res_005 <- subset(x = res, padj < 0.05) # adjusted p-value
res_01  <- subset(x = res, padj < 0.1 )
res_005_up <- subset(x = res_005, log2FoldChange > 0)

# export data to file ----
write.table(x = res, file = "covid_results.txt", sep = "\t", col.names = NA)

# Visualization of results: MA plot ----
DESeq2::plotMA(res, alpha = 0.05)

