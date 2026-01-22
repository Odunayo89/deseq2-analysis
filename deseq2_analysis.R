############################################################
# DESeq2 RNA-seq Differential Expression Analysis
# Author: Adekunle Ajiboye
# Description:
#   End-to-end DESeq2 workflow including:
#   - Data loading
#   - Filtering
#   - Differential expression analysis
#   - MA plot
#   - Volcano plot
#   - Heatmap of top genes
############################################################

library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)

cnt <- read.csv("counts.csv", row.names = 1, check.names = FALSE)
met <- read.csv("metadata2.csv", row.names = 1, check.names = FALSE)

stopifnot(all(colnames(cnt) == rownames(met)))

dds <- DESeqDataSetFromMatrix(
  countData = cnt,
  colData   = met,
  design    = ~ dexamethasone
)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds$dexamethasone <- factor(dds$dexamethasone)
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

deg <- DESeq(dds)

res <- results(deg)
write.csv(as.data.frame(res), "outcome.csv")

res0.05 <- results(deg, alpha = 0.05)

plotMA(res, ylim = c(-5, 5))

res_df <- as.data.frame(res)
res_df <- res_df[!is.na(res_df$padj), ]

best_genes <- res_df %>% arrange(padj) %>% head(10)
best_genes$gene <- rownames(best_genes)
write.csv(best_genes, "best_genes.csv", row.names = FALSE)

vol <- res_df %>% mutate(neglog10padj = -log10(padj))

png("volcano_plot_labeled.png", width = 1000, height = 800)
print(
  ggplot(vol, aes(
    x = log2FoldChange,
    y = neglog10padj,
    color = padj < 0.05 & abs(log2FoldChange) > 1
  )) +
    geom_point(alpha = 0.7) +
    geom_text(
      data = best_genes,
      aes(x = log2FoldChange, y = -log10(padj), label = gene),
      hjust = -0.2, vjust = 0.5, size = 3
    )
)
dev.off()

vsd <- vst(deg, blind = FALSE)
top30 <- rownames(res_df[order(res_df$padj), ])[1:30]
mat <- assay(vsd)[top30, ]
mat_z <- t(scale(t(mat)))

ann <- as.data.frame(colData(deg)[, "dexamethasone", drop = FALSE])
colnames(ann) <- "condition"
ann <- ann[colnames(mat_z), , drop = FALSE]

png("heatmap_top30.png", width = 1100, height = 900)
pheatmap(mat_z, annotation_col = ann)
dev.off()

save.image("deseq2_checkpoint_day3.RData")
