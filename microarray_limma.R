# microarray_limma.R
# Differential expression analysis using limma (example: GSE15852)

# Load libraries
library(GEOquery)
library(limma)
library(pheatmap)

# 1. Download dataset
gset <- getGEO("GSE15852", GSEMatrix = TRUE)[[1]]

# 2. Extract expression matrix
exprs_mat <- exprs(gset)

# 3. Build design matrix (Normal vs Tumor)
group <- factor(pData(gset)$characteristics_ch1.1)  # adjust depending on metadata
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# 4. Fit linear model
fit <- lmFit(exprs_mat, design)
contrast.matrix <- makeContrasts(Tumor-Normal, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 5. Extract top DE genes
top_genes <- topTable(fit2, adjust="fdr", number=50)
write.csv(top_genes, "results/microarray_top_genes.csv")

# 6. Volcano plot
plot(fit2$coefficients[,1], -log10(fit2$p.value[,1]),
     pch=20, main="Volcano Plot",
     xlab="Log2 Fold Change", ylab="-log10 p-value",
     col=ifelse(fit2$p.value[,1] < 0.05, "red", "black"))
abline(h=-log10(0.05), col="blue", lty=2)

# 7. Heatmap of top DE genes
top_ids <- rownames(top_genes)
pheatmap(exprs_mat[top_ids,], scale="row",
         annotation_col=data.frame(Group=group))
