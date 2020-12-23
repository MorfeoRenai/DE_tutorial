# RNA-seq_exercise
## Introduction

Transcriptional response to SARS-CoV-2 infection [link](https://www.biorxiv.org/content/10.1101/2020.03.24.004655v1.full)
The raw data referenced by the paper is present in the SRA (ID: SRP253951). The GEO page [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507)

## Trimming
trimgalore_script.sh is used to trim the sequences. It uses 'trim_galore' . The code is shown below.
mkdir trimgalore_results

for file in ./fastq_files/*
do 
    trim_galore \
        --quality 25 \
        --stringency 5 \
        --length 50 \
        --output_dir ./trimgalore_results \
        --fastqc \
        ${file}
done

The script created the directory fastqc_results and used it as the output directory. It shall contain the trimming logs and also the fastqc results, since we used the --fastqc flag.

The 0.5% of all reads is discarded in average. We can verify it by reading the trimgalore log files.
Quantification

salmon is a python tool used for a “wicked fast” transcript indexing and quantification from RNA-seq data. Since we were already supplied with the index hg38_index, we created another bash file salmon_script.sh with another for-loop in order to quantify the data. The code is shown below.

GEO=$(cat ../geo_accessions.txt)

mkdir ../salmon_results

for i in ${GEO}
do
    SRR=$(grep ${i} ../SraRunTable.txt | cut -d "," -f 1)
    SRR=$(echo ${SRR} | sed "s/ /_trimmed.fq.gz /g")
    SRR=${SRR}_trimmed.fq.gz
    salmon quant -i ../hg38_index --libType A -o ../salmon_results/${i} -r ${SRR}
done

Since each RNA-seq sample is composed by 4 runs with the same GEO accessions ID, we had to use salmon quant with four files as its arguments for each GEO ID.

The results are stored in the salmon_results directory.
MultiQC Report

The multiqc python tool a report generator perfect for summarizing the output form numerous bioinformatics tools in a single HTML file. It’s suffiecient to run multiqc ., it will return generate said HTML file and a data directory.
Differential Expression Analysis

We analyzed the differential expression thanks to R packages tximport, GenomicFeatures and DESeq2.

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

The number of genes differentially expressed at the threshold of adjusted p-value (padj) 0.05 is 325. If we change the padj threshold to 0.1 the number of significant genes will increase to 464. This can be demonstrated by running dim(res_005) and dim(res_01); the number of rows is the number of genes.

The number of upregulated genes at the threshold of 0.05 is instead 221. This can be demonstrated by running dim(res_005_up); the number of rows is the number of genes.

Finally the results will be stored in the covid_results.txt file. Also the plot can be exported; we can find the png file in the Google Drive directory.
