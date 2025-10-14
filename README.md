# Macaque Metadata

This repository contains all scripts required for organizing the metadata associated with the Cayo Santiago rhesus macaque social and -omics data. The scripts are organized by -omic type and each script is numbered in the order for which they should be run.

## Packages Required

The following snippets can be used to install all the required packages:

#### In R:

```         
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("bsseq","GenomicRanges","GenomicFeatures","comethyl", "GENESIS"))

install.packages(c("tidyverse", "here"))
```

#### In Python after cloning the repo:

```         
python -m pip install -r requirements.txt
```
