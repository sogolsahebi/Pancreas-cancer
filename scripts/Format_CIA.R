# Format_CIA.R

# This script is for creating SE_CIA.rds: 50 samples and 51627 genes

# Initial input:
# - clin has dimensions 59 x 16 (clinical data)
# - expr has dimensions 60483 x 1472 (expression data)

# load libraries 
library(SummarizedExperiment)
library(data.table)

# Read the clinical CSV file for 3 studies
separate_clin <- read.csv("files/train_tune_with_quantile.csv")

# Subset separate_clin to ones that have CIA only
clin <- separate_clin[separate_clin$cohort == "cia",]

# Clean slide_id in clin
clin$slide_id_copy <- gsub("-[^-]+$", "", clin$slide_id)

# Optionally, convert the data.table to a data.frame if needed
expr <- as.data.frame(fread("files/ALL_RNA-Seq_Expr_WashU_FPKM.tsv")) # dim 60483 x 1472

# Only matching samples with '-T' at the end. 
colnames(expr) <- sub("-T$", "", colnames(expr))

# Remove duplicate gene_names
expr<- expr[!duplicated(expr$gene_name), ] # dim 58387 x 1472

# Store the gene_names for later
gene_row <- expr$gene_name

# Filter columns of expr based clin$slide_id
expr <- expr[, colnames(expr) %in% clin$slide_id_copy] # dim 58387 x 50

# Convert to log(expr + 1)
expr <- apply(apply(expr,2,as.character),2,as.numeric)
expr <- log2(expr + 1)

# Set expr rownames
rownames(expr) <- gene_row

# Use Gencode version 19 to add annot 
load("~/BHK lab/Annotation/Gencode.v19.annotation.RData")
features_df <- features_gene

# Remove duplicated gene names
features_df <- features_df[!duplicated(features_df$gene_name), ]

# Filter and order assay data based on gene IDs
assay <- expr[rownames(expr) %in% features_df$gene_name, ]
assay <- assay[order(rownames(assay)), ]    # dim 51627 x 50

# Prepare Row Data for SummarizedExperiment
assay_genes <- features_df[features_df$gene_name %in% rownames(assay), ]
assay_genes$gene_id_no_ver <- gsub("\\..*$", "", assay_genes$gene_id)
assay_genes <- assay_genes[!is.na(assay_genes$start), ]
rownames(assay_genes) <- assay_genes$gene_name
assay_genes <- assay_genes[order(rownames(assay_genes)), ]  # Dimension 51627 x 15

# Keep clin slide IDs that match with expr column names
clin <- clin[clin$slide_id_copy %in% colnames(expr), ]

rownames(clin) = clin$slide_id_copy
patient = intersect( colnames(expr) , rownames(clin) )
clin = clin[ patient , ]
expr = expr[ , patient ]

# Create a SummarizedExperiment object for CIA
cia_se <- SummarizedExperiment(assays = list("gene_expression" = assay), colData = clin, rowData = assay_genes)

# Save the SummarizedExperiment object
saveRDS(cia_se, "output data/SE_CIA.rds")
