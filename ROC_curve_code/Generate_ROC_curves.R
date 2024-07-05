library(tidyverse)
library(SummarizedExperiment)
library(DESeq2)
library(TCGAbiolinks)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(DESeq2)
library(dplyr)
library(tidyr)
library(rlist)
install.packages("rlist")
library(pROC)

#Define list of TCGA codes
TCGA_codes <- c("BLCA","BRCA", "CESC", "CHOL", "COAD", 
                "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", 
                "LUAD", "LUSC", "PAAD", "PCPG", "PRAD",
                "READ", "SARC", "SKCM", "STAD", "THYM", "THCA", "UCEC")

#load in gene name metadata
gene_name_metadata <- read.csv("gene_name_metadata.csv")

#Loop through TCGA codes 
for (TCGA_code in TCGA_codes) {
  ###get solid tissue normal counts data 
  query_stn = GDCquery(                       
    project = paste("TCGA", TCGA_code, sep = "-"),
    data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
    experimental.strategy = "RNA-Seq",
    workflow.type = "STAR - Counts",
    data.type = "Gene Expression Quantification",
    sample.type = "Solid Tissue Normal",
    access = "open")
  
  
  #create summarized experiment variable 
  tcga_stn <- GDCprepare(query_stn, summarizedExperiment = TRUE)
  #Save the meta data to a dataframe 
  gene_metadata <- as.data.frame(rowData(tcga_stn)) 
  
  #create a matrix with transcriptomic data for each case 
  stn_matrix <- assay(tcga_stn, "unstranded")
  
  #convert expression matrix into dataframe
  stn_df <- as.data.frame(stn_matrix)
  
  ###get Primary tumour counts data 
  query_pt = GDCquery(                       
    project = paste("TCGA", TCGA_code, sep = "-"),
    data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
    experimental.strategy = "RNA-Seq",
    workflow.type = "STAR - Counts",
    data.type = "Gene Expression Quantification",
    sample.type = "Primary Tumor",
    access = "open")

  
  #create summarized experiment variable 
  tcga_pt <- GDCprepare(query_pt, summarizedExperiment = TRUE)
  
  #create a matrix with transcriptomic data for each case 
  pt_matrix <- assay(tcga_pt, "unstranded")
  
  #convert expression matrix into dataframe
  pt_df <- as.data.frame(pt_matrix)
  
  #merge STN and cancer data and subset lncRNAs ####
  
  #merge solid tissue normal and cancer tissue data raw counts 
  count_data <- merge(stn_df, pt_df, by.x = 0, by.y = 0)
  
  #reset index to gene names 
  row.names(count_data) <- count_data$Row.names
  
  #remove gene names column
  count_data <- count_data[,-1]
  
  # create DESeq2 dataset ####
  
  #check order of samples 
  sample_order <- c(colnames(stn_df), colnames(pt_df))
  
  #save the column data 
  coldata_stn <- as.data.frame(colData(tcga_stn))  
  coldata_pt <- as.data.frame(colData(tcga_pt))
  
  #retain wanted columns for sample info 
  sample_info_stn <- coldata_stn[,c(1,5)]
  sample_info_pt <- coldata_pt[,c(1,5)]
  
  #append the sample info dataframes together 
  sample_info <- bind_rows(sample_info_stn, sample_info_pt)
  colnames(sample_info) <- c("Barcode", "Tissue_type")
  
  #create deseq2 dataset
  dds <- DESeqDataSetFromMatrix(countData = count_data,
                                colData = sample_info,
                                design = ~ Tissue_type)
  #estimate size factors for normalisation
  dds <- estimateSizeFactors(dds)
  
  #normalise data by median of ratios method 
  normalised_expression <-as.data.frame(counts(dds, normalized = TRUE))
  
  #create directory in which to store ROC curves 
  dir.create(paste(TCGA_code, "_RESULTS/", TCGA_code, "_ROC_curves", sep = ""))
  
  #load in markers data 
  markers <- read.csv(paste(TCGA_code, "_RESULTS/", TCGA_code, "_biomarkers.csv", sep = ""))
  
  #convert markers names into list 
  markers_list <- c(markers$gene_name)
  
  #loop through list of markers for cancer project in given iteration of loop
  for(marker in markers_list) {
    
    #subset data to only include a single lncRNA - the marker for this iteration
    one_lnc_data <- normalised_expression[rownames(gene_name_metadata)[which(gene_name_metadata$gene_name == 
                                                                               paste(marker))],]%>%
      gather(key = 'case_id', value = 'counts')
    
    #create confusion matrix dataframe - merge sample info, which has tissue type data and one_lnc_data
    #which has normalised counts for one lncRNA.
    confusion_matrix <- merge(one_lnc_data, sample_info, by.x = "case_id", by.y = "Barcode")
    
    #create true diagnosis column based on tissue type 
    confusion_matrix$true_diagnosis <- ifelse(confusion_matrix$Tissue_type == "Solid Tissue Normal",
                                              0, 1)
    
    #pass true diagnosis and counts to roc function to determine roc curve and AUC 
    roc_obj <- roc(confusion_matrix$true_diagnosis, confusion_matrix$counts)
    par(pty = "s")
    
    #round AUC to 3 decimal places 
    AUC <- round(roc_obj$auc, 3)
    
    #plot roc curve using ggroc
    roc_curve <- ggroc(roc_obj, legacy.axes = T, color = "blue")+
      geom_abline(linetype = "dashed")+.           #dashed line represents a random classifier 
      annotate("text", x = 0.6, y = 0.25, label = paste("AUC = ",AUC, sep = ""))+
      xlab("False positive rate")+
      ylab("True positive rate")+
      theme_classic()
      
    #save plot to ROC curve directory
    ggsave(paste(TCGA_code, "_RESULTS/", TCGA_code, "_ROC_curves/",  
                 marker, "_ROC_curve.png", sep = ""), roc_curve,
           device = "png", 
           height = 5, width = 6, units = "in")
  }
}
 
