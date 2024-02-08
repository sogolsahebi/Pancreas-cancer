# Format_CIA.R

# This scrip is for making SE_TCGA.rds file which include.
# - clin is dimesnion 59 x 16.
# - expr is dimesion 60483 x 1472
# - annot is dimesnion 60483 x 2


# This scrip is for making SE_CIA.rds : 50 samples and 60483 genes

# load libraries 
# read the clin csv file for 3 studies
seprate_clin <- read.csv("files/train_tune_with_quantile.csv")

# subest seprate clin to ones that have CIA only 
clin <- seprate_clin[seprate_clin$cohort == "cia",]

# Assuming 'expr' is your dataframe
clin$slide_id <- gsub("-[^-]+$", "", clin$slide_id )

# read expr file of of cia study 
expr <- data.frame(read.table("files/ALL_RNA-Seq_Expr_WashU_FPKM.tsv", header = TRUE, sep = "\t"))

#TODO if we use Gencode 40 we have 31300 genes , and v19 we have 33008 genes.
annot <- expr[, c("gene_id", "gene_name")]
gene_row <- expr$gene_id

# Convert to log(expr + 1)
expr <- apply(apply(expr,2,as.character),2,as.numeric)
expr <- log2(expr + 1)

#set rownames
rownames(expr) <- gene_row
  
store_row <- rownames(expr)
# reformat colnames of expr
colnames(expr) <- gsub("^X", "", colnames(expr))
colnames(expr) <- gsub("\\.T$", "", colnames(expr))
colnames(expr) <- gsub("\\.", "-", colnames(expr))

# Filter columns of expr based clin$slide_id
expr <- expr[, colnames(expr) %in% clin$slide_id] # dim 60483 x 50


# use Gencode version 19 to add annot
annot <- annot[, annot$gene_id %in% rownames(expr)]

# TODO: load gencode V19 will remove 27475 genes , so I dont use Gencode.
# gencode_file <- "~/BHK lab/Annotation/Gencode.v19.annotation.RData"
# load(gencode_file)
# features_df <- features_gene
# merge_annot <- merge(annot, features_df, by = "gene_id" )

clin <- clin[clin$slide_id %in% colnames(expr), ]

# Create a SE object for CIA
cia_se <- SummarizedExperiment(assays = list(exprs = expr), colData = clin, rowData = DataFrame(annot))

# Save the SummarizedExperiment object
saveRDS(cia_se, "output data/SE_CIA.rds")
