# BiocManager install
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

# install Bioconductor packages
packages <- c("tximport", "GenomicFeatures", "DESeq2")
BiocManager::install(packages)

# install CRAN package
install.packages("readr")
