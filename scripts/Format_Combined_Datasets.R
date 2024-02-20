# Format_Combined_Datasets.R

# This script is for creating SE_Combined.rds, combining three datas sets : 209 samples and 15376 genes

# > dim(clin)
# [1] 209   2
# > dim(expr)
# [1] 15376   209

# Initial input descriptions:
# - SE_CIA.rds: Summarizes the 'CIA' study with 50 samples and 51,627 genes
# - SE_TCGA.rds: Summarizes the 'TCGA' study with 120 samples and 18,416 genes
# - SE_ICGC.rds: Summarizes the 'ICGC' study with 39 samples and 20,215 genes

# Set directory 
dir <- "~/BHK lab/Pancreas-cancer/output data/"

# Reading three dataset SE objects
SE_TCGA <- readRDS(file.path(dir, "SE_TCGA.rds"))
SE_CIA <- readRDS(file.path(dir, "SE_CIA.rds"))
SE_ICGC <- readRDS(file.path(dir, "SE_ICGC.rds"))

# Extract RNA (gene_expression), clinical, and annotation data
# TCGA dataset
expr_tcga <- assay(SE_TCGA)
clin_tcga <- data.frame(colData(SE_TCGA))
annot_tcga <- data.frame(rowData(SE_TCGA))

# CIA dataset
expr_cia <- assay(SE_CIA)
clin_cia <- data.frame(colData(SE_CIA))
annot_cia <- data.frame(rowData(SE_CIA))

# ICGC dataset
expr_icgc <- assay(SE_ICGC) 
clin_icgc <- data.frame(colData(SE_ICGC))
annot_icgc <- data.frame(rowData(SE_ICGC))

# Remove 13 NAs from TCGA expression data
gene_na <- apply(expr_tcga, 1, function(x) any(is.na(x)))
expr_tcga <- expr_tcga[!gene_na, ]

# Remove zero/low expressed genes (FPKM) for CIA and TCGA
# For TCGA
r <- as.numeric(apply(expr_tcga, 1, function(i) sum(round(i, 6) == round(log10(1), 6))))
remove <- which(r > dim(expr_tcga)[2] * 0.5)  # 50% of samples    
expr_tcga <- expr_tcga[-remove, ]  # dim 16260 x 120

# For CIA
r <- as.numeric(apply(expr_cia, 1, function(i) sum(round(i, 6) == round(log2(1), 6))))
remove <- which(r > dim(expr_cia)[2] * 0.5) 
expr_cia <- expr_cia[-remove, ] # dim 30816 x 50

# Find common genes across datasets
int <- intersect(intersect(rownames(expr_tcga), rownames(expr_cia)), rownames(expr_icgc)) # 15376 genes

# Filter expression data to include only common genes
expr_tcga <- expr_tcga[rownames(expr_tcga) %in% int, ]
expr_cia <- expr_cia[rownames(expr_cia) %in% int, ]
expr_icgc <- expr_icgc[rownames(expr_icgc) %in% int, ]

# Scale data for PCA analysis
scaled_tcga <- t(scale(t(expr_tcga)))
scaled_cia <- t(scale(t(expr_cia)))
scaled_icgc <- t(scale(t(expr_icgc)))

# Combine expression data from three cohorts for PCA
combined_expr <- cbind(scaled_cia[order(rownames(scaled_cia)), ], 
                       scaled_tcga[order(rownames(scaled_tcga)), ], 
                       scaled_icgc[order(rownames(scaled_icgc)), ])   #dim 15376 x 209

# Combine annotation data
annot_cia <- annot_cia[rownames(annot_cia) %in% rownames(scaled_cia),]
annot_tcga <- annot_tcga[rownames(annot_tcga) %in% rownames(scaled_tcga),]
annot_icgc <- annot_icgc[rownames(annot_icgc) %in% rownames(scaled_icgc),]

combined_annot <- cbind(annot_cia[order(rownames(annot_cia)), ],    # dim 15376 x 45
                        annot_tcga[order(rownames(annot_tcga)), ], 
                        annot_icgc[order(rownames(annot_icgc)), ])

# Combine clinical data
combined_clin <- data.frame(study = c(rep("CIA", nrow(clin_cia)),
                                      rep("TCGA", nrow(clin_tcga)),
                                      rep("ICGC", nrow(clin_icgc))),
                            group = c(clin_cia[, "low_high"], 
                                      clin_tcga[, "low_high"], 
                                      clin_icgc[, "low_high"]))     # 209 samples

group <- factor(combined_clin$study)

# PCA plot settings
cols <- c("#9970ab", "#74add1", "#636363") # Colors for each study

# Perform PCA
pca.results <- prcomp(t(combined_expr))
var.res <- pca.results$sdev^2 / sum(pca.results$sdev^2)

#PCA plot
plot(pca.results$x[,1], pca.results$x[,2], pch=19, xlab='PC1 (11%)', ylab='PC2 (8%)', main='PCA: TCGA, CIA & ICGC')
legend('topright', levels(group), pch=19, col=cols, cex=0.8)
dev.off() # 1

# Save merged and processed data (TCGA, CIA, & ICGC)
SE_Combined <- SummarizedExperiment(assays= list("gene_expression"= combined_expr),
                             rowData= combined_annot,
                             colData= combined_clin)

saveRDS(SE_Combined , file="~/BHK lab/Pancreas-cancer/output data/SE_Combined.rds")
