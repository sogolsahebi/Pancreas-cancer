# Format_ICGC.R

# This script is for creating the SE_TCGA.rds file: 39 samples and 20215 genes

# Input Data Specifications for tcga SE obj :
# - icgc_clin: Initial dimensions of 265 x 11 (clinical data)
# - icgc_expr: Initial dimensions of 21701 x 265 (expression data)
# - separate_clin_icgc: Initial dimensions of 217 x 14 (icgc_clinical.csv)
# - separate_clin_quantile: Initial dimensions of 257  16 (train_tune_with_quantile.csv)

# Load necessary libraries
library(SummarizedExperiment)

# Load the data
load("~/BHK lab/Pancreas-cancer/files/panc.rda")

# Access the TCGA SummarizedExperiment object
icgc_se <- pancreasData[[1]]$ICGCMICRO_SumExp

# Access clinical, expression, and annotation data respectively
icgc_clin <- data.frame(colData(icgc_se))   # dim 265 x 11
icgc_expr <- assays(icgc_se)[["exprs"]]     # dim 21701 x 265

# Load and preprocess external two clinical data
separate_clin_icgc <- read.csv("files/icgc_clinical.csv") # dimension 217 x 14
separate_clin_quantile <- read.csv("files/train_tune_with_quantile.csv") 

# Subsetthe quantile dataframe to icgc
separate_clin_quantile <- separate_clin_quantile[separate_clin_quantile$cohort == "icgc",]  #dimension 50 x 16
separate_clin_quantile$slide_id <- gsub("^([^_]+_[^_]+).*", "\\1", separate_clin_quantile$slide_id)

# Add a column to slide_id separate_clin_icgc to match with quantile slide_id
separate_clin_icgc$slide_id <- paste(separate_clin_icgc$patient_id, separate_clin_icgc$other_id, sep = "_")

# merged "separate_clin_icgc" and "separate_clin_quantile" based on slide_id
merged_clin <- merge(separate_clin_icgc, separate_clin_quantile, by = "slide_id")   # dim 51 x 30

# Identical age column, so remove the extra one
merged_clin$sex.y <- NULL
merged_clin$grade.y <- NULL
colnames(merged_clin)[colnames(merged_clin) == "sex.x"] <- "sex"
colnames(merged_clin)[colnames(merged_clin) == "grade.x"] <- "grade"

# merged with 'icgc_clin' we end up having 39 patients.
merged_clin <- merge(icgc_clin,merged_clin, by.x= "unique_patient_ID", by.y = "icgc_id")
merged_clin <- merged_clin[!duplicated(merged_clin$unique_patient_ID), ] # dim 39 x 38

# Identical age column, so remove the extra one
merged_clin$sex.y <- NULL
merged_clin$age.y <- NULL
colnames(merged_clin)[colnames(merged_clin) == "sex.x"] <- "sex"
colnames(merged_clin)[colnames(merged_clin) == "age.x"] <- "age"

#Rename columns since garde columns are not identical
colnames(merged_clin)[colnames(merged_clin) == "grade.x"] <- "grade"
colnames(merged_clin)[colnames(merged_clin) == "grade.y"] <- "grade_pathology"

# expression is log2-transformed and quantile normalized as expression data
range(icgc_expr) #1.022189 14.092020

# Filter columns of tcga_expr based on merged_clin$unique_patient_ID
icgc_expr <- icgc_expr[, colnames(icgc_expr) %in% merged_clin$unique_patient_ID] # dim 21701 x 39

# Gencode19 version si used for annotation                  
load("~/BHK lab/Annotation/Gencode.v19.annotation.RData")
features_df <- features_gene

# Remove Duplicated Gene Names
features_df <- features_df[!duplicated(features_df$gene_name), ]

# Filter and Order Assay Data
assay <- icgc_expr[rownames(icgc_expr) %in% features_df$gene_name, ]   # dim 20215 x 39 
assay <- assay[order(rownames(assay)), ]

# Prepare Row Data for SummarizedExperiment
assay_genes <- features_df[features_df$gene_name %in% rownames(assay), ]
assay_genes$gene_id_no_ver <- gsub("\\..*$", "", assay_genes$gene_id)
assay_genes <- assay_genes[!is.na(assay_genes$start), ]
rownames(assay_genes) <- assay_genes$gene_name
assay_genes <- assay_genes[order(rownames(assay_genes)), ]    # dim 20215 x 15

# Create a SummarizedExperiment object for TCGA
icgc <- SummarizedExperiment(assays = list("gene_expression" = assay), colData = merged_clin, rowData = assay_genes)

# Save the SummarizedExperiment object
saveRDS(icgc, "output data/SE_ICGC.rds")

