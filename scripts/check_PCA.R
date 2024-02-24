# add libraries 
library(sva)
library(limma)
library(stats)

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

# # Scale data for PCA analysis
# scaled_tcga <- t(scale(t(expr_tcga)))
# scaled_cia <- t(scale(t(expr_cia)))
# scaled_icgc <- t(scale(t(expr_icgc)))

# # Combine expression data from three cohorts for PCA
# combined_expr <- cbind(scaled_cia[order(rownames(scaled_cia)), ], 
#                        scaled_tcga[order(rownames(scaled_tcga)), ], 
#                        scaled_icgc[order(rownames(scaled_icgc)), ])   #dim 15376 x 209

# Combine expression data from three cohorts for PCA
combined_expr <- cbind(expr_cia[order(rownames(expr_cia)), ], 
                       expr_tcga[order(rownames(expr_tcga)), ], 
                       expr_icgc[order(rownames(expr_icgc)), ])   

# Combine annotation data
annot_cia <- annot_cia[rownames(annot_cia) %in% rownames(expr_cia),]
annot_tcga <- annot_tcga[rownames(annot_tcga) %in% rownames(expr_tcga),]
annot_icgc <- annot_icgc[rownames(annot_icgc) %in% rownames(expr_icgc),]


combined_annot <- cbind(annot_cia[order(rownames(annot_cia)), ],    # dim 15376 x 45
                        annot_tcga[order(rownames(annot_tcga)), ], 
                        annot_icgc[order(rownames(annot_icgc)), ])

# Combine clinical data
combined_clin <- data.frame(study = c(rep("CIA", nrow(clin_cia)),
                                      rep("TCGA", nrow(clin_tcga)),
                                      rep("ICGC", nrow(clin_icgc))),
                            group = c(clin_cia[, "low_high"], 
                                      clin_tcga[, "low_high"], 
                                      clin_icgc[, "low_high"]))   # 209 samples

group <- factor(combined_clin$study)

# PCA plot settings
cols <- c("#9970ab", "#74add1", "#636363") # Colors for each study

# Perform PCA
pca.results <- prcomp(t(combined_expr))
var.res <- pca.results$sdev^2 / sum(pca.results$sdev^2)* 100
var.res <- round(var.res, 2) # Round to 2 decimal places

#PCA plot
xlab_text <- paste("PC1 (", var.res[1], "%)", sep = "")
ylab_text <- paste("PC2 (", var.res[2], "%)", sep = "")

plot(pca.results$x[,1], pca.results$x[,2], pch=19, xlab= xlab_text , ylab= ylab_text, main='PCA: TCGA, CIA & ICGC')
legend('topright', levels(group), pch=19, col=cols, cex=0.8)
ind= group == levels(group)[1]
points(pca.results$x[ind,1],pca.results$x[ind,2], pch=19,col=cols[1])
ind= group == levels(group)[2]
points(pca.results$x[ind,1],pca.results$x[ind,2], pch=19,col=cols[2])
ind= group == levels(group)[3]
points(pca.results$x[ind,1],pca.results$x[ind,2], pch=19,col=cols[3])
legend('topright',paste(levels(group)), pch=19,col=cols, cex=0.8)

# dev.off()

#--- Removing the Batch Effect with sva Package ---

# Defining model matrices
mod = model.matrix(~as.factor(combined_clin$group), data = combined_clin)
mod0 = model.matrix(~1, data = combined_clin)

# Estimating the number of surrogate variables
#n.sv = num.sv(combined_expr, mod, method = "leek")
n.sv = num.sv(combined_expr, mod, method = "be")
n.sv

var.res[1:4]
# Estimating surrogate variables
svobj = sva(combined_expr, mod, mod0, n.sv = n.sv)

# method1: surrogate variables using 'f.pvalue ' function.

# Initial p-value calculation without surrogate variable correction
pValues = f.pvalue(combined_expr, mod, mod0)
qValues = p.adjust(pValues, method = "BH")

# Adjusting the linear model for surrogate variables
modSv = cbind(mod, svobj$sv)
mod0Sv = cbind(mod0, svobj$sv)

# P-values and q-values adjusted with surrogate variable correction
pValuesSv = f.pvalue(combined_expr, modSv, mod0Sv)
q3ValuesSv = p.adjust(pValuesSv, method = "BH")

# method2: surrogate variables using limma package

# Creating an sva-adjusted version of data for 'limma' model fit

fit = lmFit(combined_expr, modSv)

# Constructing the Contrast Matrix
numCovs <- ncol(mod)
numSVs <- svobj$n.sv

# Formulating a standard-sized contrast matrix for differential analysis
contrast.matrix <- matrix(c(-1, 0, 1, -1, rep(0, (numCovs + numSVs) - 2)), 
                          nrow = numCovs + numSVs, ncol = 2, byrow = FALSE)

# Enhancing linear fit to embody output deconvoluted from sva with established contrasts
fitContrasts = contrasts.fit(fit, contrast.matrix)

# Acquisition of Empirical Bayes statics used in meta-comparison
eb = eBayes(fitContrasts)

# Inducting a hypostasis test with F-adjustment towards display of places
topTableF(eb, adjust = "BH")

# Output example for the demonstration
tT <- topTable(eb, adjust = "BH", sort.by = "F")
head(tT)



