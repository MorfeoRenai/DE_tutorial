# RNA-seq_exercise

## Introduction

This will be a very simple, undergrad level, step-by-step tutorial to differential analysis in RNA-seq. It's based on the *Blanco-Melo et al* 2020 [paper](www.biorxiv.org/content/10.1101/2020.03.24.004655v1.full) on SARS-CoV-2 transcriptional signature. The raw data referenced by the paper is present in the Sequence Read Archive (or SRA) and it can be found with this ID SRP253951. [Here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507) you can find the Gene Expression Omnibus (or GEO) webpage.

## Conda Environment

SET UP CONDA ENV
INSTALL PACKAGES

## Downloading the data

From the GEO webpage we shall access the SRA Run Selector (Fig.1). Select two selection filters, `Cell_line:nhbe` and `Organism:homo sapiens`, 24 runs will be selected.

![Fig.1](img/fig1.jpg)

Now click on the **Metadata** button (Fig.2) the CSV file `SraRunTable.txt` will be downloaded. It contains all the metadata from the selected sequencing runs.

![Fig.2](img/fig2.jpg)

Before going any further, we must note (see Fig.3) that each RNA-seq sample is actually composed by 4 runs, beacause each library was splitted in 4 different lanes during the sequencing. You can see that the Runs from 1 to 4 and from 5 to 8 are associated to only one Experiment ID and one GEO_Accession ID.

![Fig.3](img/fig3.jpg)

Now we can easily download every selected run using the sra-toolkit, in particular `prefetch` and `fastq-dump`.

```sh
VAR=$(cut -d ',' -f 1 SraRunTable.txt | tail -n +2) # select the first field, containg the IDs needed for the sra-toolkit, and eliminate the header "Run"

for i in ${VAR}
    do
        echo "Download SRA sample: ${i}"
        prefetch ${i}                               # for each run ID we prefetch the data...
        fastq-dump --gzip --defline-qual '+' ${i}   # ...and then download it already decrompressed 
    done
```
The `--defline-qual '+'` is needed because without this argument `fastq-dump` for some reasons eliminates the third row (second header) of the fastq read.

## Trimming

WHAT IS TRIMMING?

```sh
mkdir trimgalore_results

for file in ./fastq_files/*
do 
    trim_galore \
        --quality 25 \
        --stringency 5 \
        --length 50 \
        --output_dir ./trimgalore_results \
        --fastqc \    # output dir will contains both the trimming logs and fastqc results
        ${file}
done
```

## Quantification

```sh
GEO=$(cat ../geo_accessions.txt)

mkdir ../salmon_results

for i in ${GEO}
do
    SRR=$(grep ${i} ../SraRunTable.txt | cut -d "," -f 1)
    SRR=$(echo ${SRR} | sed "s/ /_trimmed.fq.gz /g")
    SRR=${SRR}_trimmed.fq.gz
    salmon quant -i ../hg38_index --libType A -o ../salmon_results/${i} -r ${SRR}
done
```

Since each RNA-seq sample is composed by 4 runs with the same GEO accessions ID, we had to use `salmon quant` with four files as its arguments, one for each GEO ID. The results are stored in the salmon_results directory.

## MultiQC Report

The `multiqc` python tool a report generator, perfect for summarizing the output form numerous bioinformatics tools in a single HTML file. It’s suffiecient to run `multiqc .`, it will generate said HTML file and a data directory.

MULTIQC BREAKDOWN

## Differential Expression Analysis

```r
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
```

Finally the results will be stored in the covid_results.txt file. The plot can be exported too.
