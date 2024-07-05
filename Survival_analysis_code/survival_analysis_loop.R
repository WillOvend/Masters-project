library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)
library(tidyr)
library(ggplot2)
library(hdf5r)
library(HDF5Array)

#define list of TCGA codes 
TCGA_codes <- c("BLCA","BRCA", "CESC", "CHOL", "COAD", 
                "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", 
                "LUAD", "LUSC", "PAAD", "PCPG", "PRAD",
                "READ", "SARC", "SKCM", "STAD", "THYM", "THCA", "UCEC")
TCGA_code ="BRCA"

for (TCGA_code in TCGA_codes){
  clinical_data  <- GDCquery_clinic(paste("TCGA-",TCGA_code, sep = ""))
  
  clinical_data$deceased <- ifelse(clinical_data$vital_status == "Alive", FALSE, TRUE) #recoded the data so that 
  #those that are alive are given a value of FALSE and those that are dead are given a value TRUE - this is 
  #boolean code to express the death status of the individual 
  
  clinical_data$overall_survival <- ifelse(clinical_data$vital_status =="Alive",   #created a new column where 
                                           clinical_data$days_to_last_follow_up,  #if alive - assigned days to 
                                           clinical_data$days_to_death) 
  
  #last follow up value and if dead 
  #assigned days to death value 
  
  #save clinical survival data
  survival_data <- clinical_data[, c("submitter_id", "deceased", "overall_survival")]
  write.table(survival_data, paste(TCGA_code, "_RESULTS/clinical_survival_data.csv", sep = ""), row.names = T, col.names = T,
              sep = ",")
  
  #now query for transcriptomic data - 
  query = GDCquery(                       
    project = paste("TCGA-", TCGA_code, sep = ""),
    data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
    experimental.strategy = "RNA-Seq",
    workflow.type = "STAR - Counts",
    data.type = "Gene Expression Quantification",
    sample.type = "Primary Tumor",
    access = "open")
  
  #set 
  transcriptomic_data <- GDCprepare(query, summarizedExperiment = TRUE)
  
  transcriptomic_matrix <- assay(transcriptomic_data, "unstranded")
  
  #saving metadata from summarized experiment object - this is data regarding the individual transcripts being counted
  #including what typr of coding sequence they are attributed to 
  
  gene_metadata <- as.data.frame(rowData(transcriptomic_data))  
  
  #saving column data - this is data regarding the sample, including patient ID, barcode etc
  
  coldata <- as.data.frame(colData(transcriptomic_data))  
  
  #the counts data we have downloaded now needs to be variance stabilised - this is so that we can later 
  #make a cut off for high and low expression, allowing the split of our patients into groups 
  
  #DESeq2 is a package that deals with differential expression of transriptomic data, in this case we give the 
  #function our countdata, contained within the lihc matrix and our column data, contained within the df coldata 
  
  dds <- DESeqDataSetFromMatrix(countData = transcriptomic_matrix,
                                colData = coldata,
                                design = ~ 1)
  
  #variance stabilising transformation 
  vsd <- vst(dds, blind=FALSE)
  transcriptomic_data <- assay(vsd)
  
  #convert to data frame 
  transcriptomic_data <- as.data.frame(transcriptomic_data)
  
  #save transcriptomic data
  write.table(transcriptomic_data, paste(TCGA_code, "_RESULTS/transcriptomic_data.csv", sep = ""), row.names = T,
              col.names = T, sep = ",")
  
}


#read in gene name metadata
gene_name_metadata <- read.csv("gene_name_metadata.csv")
ESCA <- read.csv("ESCA_RESULTS/transcriptomic_data.csv")
coad_markers <- read.csv("ESCA_RESULTS/ESCA_bi")

for (TCGA_code in TCGA_codes) {
  
  #read in survival data
  survival_data <- read.csv(paste(TCGA_code, "_RESULTS/clinical_survival_data.csv", sep = ""))
  
  #read in transcriptomic data 
  transcriptomic_data <- read.csv(paste(TCGA_code, "_RESULTS/transcriptomic_data.csv", sep = "") ,check.names = F)
  
  #read in biomarkers 
  markers <- read.csv(paste(TCGA_code, "_RESULTS/", TCGA_code, "_biomarkers.csv", sep = ""))
  
  #convert markers names into list 
  markers_list <- c(markers$gene_name)
  
  #create results directory
  dir.create(paste(TCGA_code, "_RESULTS/", TCGA_code, "_survival_plots", sep = ""))

  for (marker in markers_list) {
    #select only the row of the transcriptomic data concerning the desired RNA
    #and convert the data to long format with a column with patient case IDs 
    #and a column with counts 
    one_lnc_data <- transcriptomic_data[rownames(gene_name_metadata)[which(gene_name_metadata$gene_name == 
                                                                             paste(marker))],] %>%
    gather(key = 'case_id', value = 'counts')
    
    #get median count value - in the condition that the median value is the minimum counts value 
    #the median is increased. by 0.1. This means that there will be a divide into high and low
    #expression groups even when the median is the minimum value - otherwise all values would be 
    #in the high expression group 
    if (median(one_lnc_data$counts) == min(one_lnc_data$counts)) {
      median_value <- min(one_lnc_data$counts) + 0.1
    } else {
      median_value <- median(one_lnc_data$counts)
    }
   
    
    #use median value to denote high and low expression
    one_lnc_data$strata <- ifelse(one_lnc_data$counts >= median_value, "HIGH", "LOW")
    
    #Add clinical information to lihc_RRN3P2
    #first we need to change the case-id names so they match the submitter-id's from the clinical dataset 
    one_lnc_data$case_id <- gsub('-01.*','', one_lnc_data$case_id) #this says sub everything after the '01' with nothing - gets rid of the end of rht ID 
    one_lnc_data <- merge(one_lnc_data, survival_data, by.x = 'case_id', by.y = 'submitter_id')
    
    #Now we compute the survival curve 
    fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = one_lnc_data)
    
    #and finally plot the survival curve
    survival_curve_plot <- ggsurvplot(fit,
                                      xlab = "Days post diagnosis",
                                      legend.title = "Expression level",
                                      data = one_lnc_data,
                                      pval = T)
    
    #save plot 
    ggsave(paste(TCGA_code, "_RESULTS/", TCGA_code, "_survival_plots/",  
                 marker, "_survival_curve.png", sep = ""), survival_curve_plot$plot,
           device = "png", 
           height = 10, width = 10, units = "in")
    
    
  }
}









