#load required packages ####

library(tidyverse)
library(SummarizedExperiment)
library(DESeq2)
library(TCGAbiolinks)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(DESeq2)
library(dplyr)
library(tidyr)


#First step: Create a dataframe with genes as row names and a single column 
#which contains median expression values from a single TCGA project.####


###get Primary tumour counts data 
query = GDCquery(                       
  project = "ACC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = "Primary Tumor",
  access = "open")

#download data 
GDCdownload(query)

#create summarized experiment variable 
se <- GDCprepare(query, summarizedExperiment = TRUE)

#Save the meta data to a dataframe 
gene_metadata <- as.data.frame(rowData(se)) 

#create a matrix with transcriptomic data for each case 
count_matrix <- assay(se, "unstranded")

#convert expression matrix into dataframe
data <- as.data.frame(count_matrix)

#save dataframe with gene names mapping to ensembl gene names 
gene_name_metadata <- gene_metadata[7]

# create DESeq2 dataset ####

coldata_pt <- as.data.frame(colData(se))

#retain wanted columns for sample info 

sample_info <- coldata_pt[,c(1,5)]

#save sample metadata as csv 
write.table(sample_info, file = "Sample_info.csv", sep = ',',
            row.names = T, col.names = T)

#read in coldata
coldata <- read.csv("Sample_info.csv")

#Create DESeq dataset using coldata and count data 
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design = ~ 1)

#Now count data is in a DESeq2 object, the DESeq2 package can be used to perform 
#median of ratios normalisation - I have added this step to adjust for library 
#composition and size. These normalisation steps are important for comparing 
#between samples and sine i am going to find a median expression value across 
#all samples, this normalisation is necassary 

#estimate size factors 
dds <- estimateSizeFactors(dds)

#normalise data by median of ratios method 
normalised_counts <- as.data.frame(counts(dds, normalized = TRUE))

#create a column with a median value for each gene 
normalised_counts$median <- apply(normalised_counts, 1, median, na.rm = TRUE)

#create data frame with genes as row names and one column containing median 
#expression values for the ACC project 
median_expression <- as.data.frame(normalised_counts[length(normalised_counts)])
colnames(median_expression) <- "ACC"

#Step 2: Now that the median expression dataframe is created the data from all other TCGA #### 
#projects can be added, which will produce a data frame with genes as row names 
#and TCGA projects as column headers. The data in the matrix will be the median 
#expression values for a gene in a given project/ in a given cancer

#Create a list with all the TCGA projects i wish to mine data from 
TCGA_codes = c("BLCA", "LGG", "BRCA", "CESC", "CHOL", "COAD", 
               "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", 
               "LUAD", "LUSC", "DLBC", "MESO", "OV", "PAAD", "PCPG", "PRAD",
               "READ", "SARC", "SKCM", "STAD", "TGCT", "THYM", "THCA", "UCS", "UCEC", 
               "UVM")

#This for loop is a repetition of the code above so that each iteration a new 
#TCGA projects data median values are added to the median_expression data frame
for (TCGA_code in TCGA_codes) {
    #download data ####
  
  ###get Primary tumour counts data 
  query = GDCquery(                       
    project = paste("TCGA", TCGA_code, sep = "-"),
    data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
    experimental.strategy = "RNA-Seq",
    workflow.type = "STAR - Counts",
    data.type = "Gene Expression Quantification",
    sample.type = "Primary Tumor",
    access = "open")
  
  #if the inputted TCGA code retrives no data, the loop stops here and skips 
  #to the next iteration 
  if (length(query) ==0) {
    next
  }
  #download data 
  GDCdownload(query)
  
  #create summarized experiment variable 
  se <- GDCprepare(query, summarizedExperiment = TRUE)
  
  #Save the meta data to a dataframe 
  gene_metadata <- as.data.frame(rowData(se)) 
  
  #create a matrix with transcriptomic data for each case 
  count_matrix <- assay(se, "unstranded")
  
  #convert expression matrix into dataframe
  data <- as.data.frame(count_matrix)
  
  #save dataframe with gene names mapping to ensembl gene names 
  gene_name_metadata <- gene_metadata[7]
  
  # create DESeq2 dataset ####
  coldata_pt <- as.data.frame(colData(se))
  
  #retain wanted columns for sample info 
  sample_info <- coldata_pt[,c(1,5)]
  
  #save sample metadata as csv 
  write.table(sample_info, file = "Sample_info.csv", sep = ',',
              row.names = T, col.names = T)
  
  #read in coldata
  coldata <- read.csv("Sample_info.csv")
  
  #Create DESeq dataset using coldata and count data 
  dds <- DESeqDataSetFromMatrix(countData = data,
                                colData = coldata,
                                design = ~ 1)
  
  #normalisation as above 
  dds <- estimateSizeFactors(dds)
  normalised_counts <- as.data.frame(counts(dds, normalized = TRUE))
  
  #generate median values column 
  normalised_counts$median <- apply(normalised_counts, 1, median, na.rm = TRUE)
  
  #add median values to median_expression dataframe 
  median_expression[, paste(TCGA_code)] <- normalised_counts$median
  
}

#save median expression values dataframe to a csv file 
write.table(median_expression, "median_expression.csv", col.names = T,
            row.names = T, sep = ",")

#read in to test it has been saved correctly
median_read_in_test <- read.csv("median_expression.csv")
