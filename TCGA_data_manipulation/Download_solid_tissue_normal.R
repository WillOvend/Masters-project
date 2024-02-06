library(tidyverse)
library(SummarizedExperiment)
library(DESeq2)
library(TCGAbiolinks)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(DESeq2)
library(dplyr)
library(tidyr)

###get counts data for liver and hepatocellular carcinoma patients 
query_lihc_all = GDCquery(                       
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor",
  access = "open")

#create summarized experiment variable 
tcga_lihc_data <- GDCprepare(query_lihc_all, summarizedExperiment = TRUE)
#Save the meta data to a dataframe 
gene_metadata <- as.data.frame(rowData(tcga_lihc_data)) 

#create a matrix with transcriptomic data for each case 
lihc_matrix <- assay(tcga_lihc_data, "unstranded")