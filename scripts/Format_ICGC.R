# Format_ICGC.R

# This script is for creating the SE_TCGA.rds file: 75 samples and 20215 genes

# Input Data Specifications for tcga SE obj :
# - tcga_clin: Initial dimensions of 265 x 11 (clinical data)
# - tcga_expr: Initial dimensions of 21701 x 265 (expression data)
# - separate_clin: Initial dimensions of 217 x 14 (icgc_clinical.csv)

# Load necessary libraries
library(SummarizedExperiment)

# Load the data
load("~/BHK lab/Pancreas-cancer/files/panc.rda")

# Access the TCGA SummarizedExperiment object
icgc_se <- pancreasData[[1]]$ICGCMICRO_SumExp # TODO: Double-check if its not ICGCSEQ_SumExp

# Access clinical, expression, and annotation data respectively
icgc_clin <- data.frame(colData(icgc_se))   # dim 265 x 11
icgc_expr <- assays(icgc_se)[["exprs"]]     # dim 21701 x 265

# Load and preprocess external clinical data
separate_clin <- read.csv("files/icgc_clinical.csv") # dimension 217 x 14

# Merge the clinical data
merged_clin <- merge(icgc_clin, separate_clin, by.x = "unique_patient_ID", by.y = "icgc_id") # dimension 137 x 24

# Identical age column, so remove the extra one
merged_clin$sex.y <- NULL
merged_clin$grade.y <- NULL
colnames(merged_clin)[colnames(merged_clin) == "sex.x"] <- "sex"
colnames(merged_clin)[colnames(merged_clin) == "grade.x"] <- "grade"

# Filter columns of icgc_expr based on merged_clin$unique_patient_ID
icgc_expr <- icgc_expr[, colnames(icgc_expr) %in% merged_clin$unique_patient_ID] # dim 21701 x 75

# Check the range of expression values
range(icgc_expr, na.rm = TRUE)   #1.022189 14.092020

# Gencode19 version si used for annotation                  
load("~/BHK lab/Annotation/Gencode.v19.annotation.RData")
features_df <- features_gene

# Remove Duplicated Gene Names
features_df <- features_df[!duplicated(features_df$gene_name), ]

# Filter and Order Assay Data
assay <- icgc_expr[rownames(icgc_expr) %in% features_df$gene_name, ]   # dimension 20215 x 75 
assay <- assay[order(rownames(assay)), ]

# Prepare Row Data for SummarizedExperiment
assay_genes <- features_df[features_df$gene_name %in% rownames(assay), ]
assay_genes$gene_id_no_ver <- gsub("\\..*$", "", assay_genes$gene_id)
assay_genes <- assay_genes[!is.na(assay_genes$start), ]
rownames(assay_genes) <- assay_genes$gene_name
assay_genes <- assay_genes[order(rownames(assay_genes)), ]    # dimension 20215 x 15

# Keep clin slide IDs that match with expr column names
intersect_ids <- intersect(merged_clin$unique_patient_ID, colnames(assay))
clin <- merged_clin[merged_clin$unique_patient_ID %in% intersect_ids, ]
clin <- clin[!duplicated(clin$unique_patient_ID), ]

# Create a SummarizedExperiment object for TCGA
icgca_se <- SummarizedExperiment(assays = list("gene_expression" = assay), colData = clin, rowData = assay_genes)

# Save the SummarizedExperiment object
saveRDS(icgca_se , "output data/SE_ICGCA.rds")

