# Format_TCGA.R

# This script is for creating the SE_TCGA.rds file: 120 samples and 18416 genes

# Input Data Specifications for tcga SE obj :
# - tcga_clin: Initial dimensions of 146 x 11 (clinical data)
# - tcga_expr: Initial dimensions of 20,501 x 146 (expression data)
# - separate_clin_tcga: Initial dimensions of 148 x 16 (train_tune_with_quantile.csv)

# Load necessary libraries
library(SummarizedExperiment)

# Load the data
load("~/BHK lab/Pancreas-cancer/files/panc.rda")

# Access the TCGA SummarizedExperiment object
tcga_se <- pancreasData[[1]]$TCGA_SumExp

# Access clinical, expression, and annotation data respectively
tcga_clin <- data.frame(colData(tcga_se))
tcga_expr <- assays(tcga_se)[["exprs"]]

# Load and preprocess external clinical data
separate_clin <- read.csv("files/train_tune_with_quantile.csv")

# Subset the TCGA only
separate_clin_tcga <- separate_clin[separate_clin$cohort == "tcga",] # dim 148 x 16

# Rename the same column names for accuracy
colnames(separate_clin_tcga)[colnames(separate_clin_tcga) %in% c("grade", "age", "sex")] <- c("grade_pathology", "age_pathology", "sex_pathology")

# Replace hyphens with underscores in slide_id
separate_clin_tcga$slide_id <- gsub("-", "_", separate_clin_tcga$slide_id)

# Remove an underscore followed by any three characters at the end of slide_id
separate_clin_tcga$slide_id <- gsub("_.{3}$", "", separate_clin_tcga$slide_id)

# Merge the clinical data
merged_clin <- merge(tcga_clin, separate_clin_tcga, by.x = "unique_patient_ID", by.y = "slide_id") # dim 120 x 26

# Identical age column, so remove the extra one
merged_clin$age_pathology <- NULL

# Note: The stage, stage_pathology are differnt 
# View(merged_clin[, c("grade", "grade_pathology")])

# Columns 'sex' and 'sex_pathology' are identical, remove 'sex2'
merged_clin$sex_pathology <- NULL

# Harmonize column names format in tcga_expr to match merged_clin IDs
colnames(tcga_expr) <- gsub("\\.", "_", colnames(tcga_expr))

# Filter columns of tcga_expr based on merged_clin$unique_patient_ID
tcga_expr <- tcga_expr[, colnames(tcga_expr) %in% merged_clin$unique_patient_ID]

# Remove rows with all NA values
tcga_expr_cleaned <- tcga_expr[!apply(tcga_expr, 1, function(x) all(is.na(x))), ] #13 rows removed

# Double-check the range of expression values
range(tcga_expr, na.rm = TRUE) # 0.00000 to 22.62807

# Gencode19 version si used for annotation                  
load("~/BHK lab/Annotation/Gencode.v19.annotation.RData")
features_df <- features_gene

# Remove Duplicated Gene Names
features_df <- features_df[!duplicated(features_df$gene_name), ]

# Filter and Order Assay Data
assay <- tcga_expr[rownames(tcga_expr) %in% features_df$gene_name, ]   # dimension 18416 x 120 
assay <- assay[order(rownames(assay)), ]
 
# Prepare Row Data for SummarizedExperiment
assay_genes <- features_df[features_df$gene_name %in% rownames(assay), ]
assay_genes$gene_id_no_ver <- gsub("\\..*$", "", assay_genes$gene_id)
assay_genes <- assay_genes[!is.na(assay_genes$start), ]
rownames(assay_genes) <- assay_genes$gene_name
assay_genes <- assay_genes[order(rownames(assay_genes)), ]    # dimension 18416 x 120

# Create a SummarizedExperiment object for TCGA
tcga_se <- SummarizedExperiment(assays = list("gene_expression" = assay), colData = merged_clin, rowData = assay_genes)

# Save the SummarizedExperiment object
saveRDS(tcga_se, "output data/SE_TCGA.rds")

