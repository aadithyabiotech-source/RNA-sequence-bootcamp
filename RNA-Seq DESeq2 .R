# rnaseq_deseq2.R
# RNA-Seq differential expression analysis using DESeq2 (airway dataset)

# Load libraries
library(airway)
library(DESeq2)
library(pheatmap)

# 1. Load dataset
data("airway")
airway

# 2. Build DESeq2 object
dds <- DESeqDataSet(airway, design = ~ cell + dex)

# 3. Run DESeq2 pipeline
dds <- DESeq(dds)
res <- results(dds)

# 4. Save results
write.csv(as.data.frame(res), "results/rnaseq_deseq2_results.csv")

# 5. Volcano plot
plot(res$log2FoldChange, -log10(res$pvalue),
     pch=20, main="Volcano Plot",
     xlab="Log2 Fold Change", ylab="-log10 p-value",
     col=ifelse(res$pvalue < 0.05, "red", "black"))
abline(h=-log10(0.05), col="blue", lty=2)

# 6. Variance Stabilizing Transformation
vsd <- vst(dds, blind=FALSE)

# 7. PCA plot
plotPCA(vsd, intgroup=c("dex","cell"))

# 8. Heatmap of top variable genes
mat <- assay(vsd)
mat <- mat[complete.cases(mat), ]
topVarGenes <- head(order(rowVars(mat), decreasing=TRUE), 50)
mat_top <- mat[topVarGenes, ]
annotation_col <- as.data.frame(colData(dds))
pheatmap(mat_top, scale="row", annotation_col=annotation_col)
