# Format_ICGC.R

# This script is for creating the SE_TCGA.rds file: 19 samples and 15051 genes

# Input Data Specifications for tcga SE obj :
# - icgc_clin: Initial dimensions of 95 x 11 (clinical data)
# - icgc_expr: Initial dimensions of 15965 x 95 (expression data)
# - separate_clin_icgc: Initial dimensions of 217 x 14 (icgc_clinical.csv)
# - separate_clin_quantile: Initial dimensions of 257  16 (train_tune_with_quantile.csv)

# Load necessary libraries
library(SummarizedExperiment)

# Load the data
load("~/BHK lab/Pancreas-cancer/files/panc.rda")

# Access the TCGA SummarizedExperiment object
icgc_se <- pancreasData[[1]]$ICGCSEQ_SumExp

# Access clinical, expression, and annotation data respectively
icgc_clin <- data.frame(colData(icgc_se))   # dim 95 x 11
icgc_expr <- assays(icgc_se)[["exprs"]]     # dim 15965 x 95

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

length(intersect(merged_clin$icgc_id, colnames(icgc_expr))) # 4 patient matching!

#### only 4 patient matching ####

