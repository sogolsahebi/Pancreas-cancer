# Format_TCGA.R

# This script is for creating the SE_TCGA.rds file.

# Purpose: Create the SE_TCGA.rds file, integrating clinical, expression, and annotation data.

# Input Data Specifications:
# - tcga_clin: Initial dimensions of 146 x 11 (clinical data)
# - tcga_expr: Initial dimensions of 20,501 x 146 (expression data)
# - t_annot: Initial dimensions of 20,501 x 4 (annotation data)

# Output Data Specifications:
# - SE_TCGA: A SummarizedExperiment object containing:
#   - Clinical data with dimensions of 120 x 24
#   - Expression data with dimensions of 20,501 x 120
#   - Annotation data with dimensions of 20,501 x 4

# Load necessary libraries
library(SummarizedExperiment)

# Load the data
load("~/BHK lab/Pancreas-cancer/files/panc.rda")

# Access the TCGA SummarizedExperiment object
tcga_se <- pancreasData[[1]]$TCGA_SumExp

# Access clinical, expression, and annotation data respectively
tcga_clin <- data.frame(colData(tcga_se))
tcga_expr <- assays(tcga_se)[["exprs"]]
tcga_annot <- data.frame(rowData(tcga_se))

# Load and preprocess external clinical data
separate_clin <- read.csv("files/train_tune_with_quantile.csv")

# Subset the TCGA only
separate_clin_tcga <- separate_clin[separate_clin$cohort == "tcga",] # dim 148 x 16

# Rename the same column names for accuracy
colnames(separate_clin_tcga)[colnames(separate_clin_tcga) %in% c("grade", "age", "sex")] <- c("grade2", "age2", "sex2")

# Replace hyphens with underscores in slide_id
separate_clin_tcga$slide_id <- gsub("-", "_", separate_clin_tcga$slide_id)

# Remove an underscore followed by any three characters at the end of slide_id
separate_clin_tcga$slide_id <- gsub("_.{3}$", "", separate_clin_tcga$slide_id)

# Merge the clinical data
merged_clin <- merge(tcga_clin, separate_clin_tcga, by.x = "unique_patient_ID", by.y = "slide_id") # dim 120 x 26

# Identical age column, so remove the extra one
if(identical(merged_clin$age, merged_clin$age2)) merged_clin$age2 <- NULL

# TODO: The stage, stage2 are not the same; mismatch happening.
View(merged_clin[, c("grade", "grade2")])

# Replace "female" with "F" and "male" with "M" in sex2 column
merged_clin$sex2[merged_clin$sex2 == "female"] <- "F"
merged_clin$sex2[merged_clin$sex2 == "male"] <- "M"

# Convert sex2 to factor with specified levels
merged_clin$sex2 <- factor(merged_clin$sex2, levels = c("F", "M"))

# If 'sex' and 'sex2' are identical, remove 'sex2'
if(identical(merged_clin$sex, merged_clin$sex2)) merged_clin$sex2 <- NULL

# Harmonize column names format in tcga_expr to match merged_clin IDs
colnames(tcga_expr) <- gsub("\\.", "_", colnames(tcga_expr))

# Filter columns of tcga_expr based on merged_clin$unique_patient_ID
tcga_expr <- tcga_expr[, colnames(tcga_expr) %in% merged_clin$unique_patient_ID]

# TODO: There are 1560 NA values in expression data. Consider addressing these.
sum(is.na(tcga_expr))

# Double-check the range of expression values
range(tcga_expr, na.rm = TRUE) # 0.00000 to 22.62807

# TODO: I use Gencode19 version results in the loss of 18416 genes. 
gencode_file <- "~/BHK lab/Annotation/Gencode.v19.annotation.RData"
load(gencode_file)
features_df <- features_gene

# Create a SummarizedExperiment object for TCGA
tcga_se_result <- SummarizedExperiment(assays = list(exprs = tcga_expr), colData = merged_clin, rowData = DataFrame(tcga_annot))

# Save the SummarizedExperiment object
saveRDS(tcga_se_result, "output data/SE_TCGA.rds")

