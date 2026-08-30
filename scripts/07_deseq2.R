# Differential expression analysis with DESeq2
# Drosophila melanogaster RNA-seq
#
# Comparison:
# RNAi vs Untreated
#
# Input:
#   results/counts/all_samples_gene_counts.txt
#   results/metadata.csv
#
# Output:
#   results/DEG/
#   results/plots/

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

# ------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------

project_dir <- "/mnt/d/RNAseq"

count_file <- file.path(
  project_dir,
  "results/counts/all_samples_gene_counts.txt"
)

metadata_file <- file.path(
  project_dir,
  "results/metadata.csv"
)

result_dir <- file.path(
  project_dir,
  "results/DEG"
)

plot_dir <- file.path(
  project_dir,
  "results/plots"
)

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Read count matrix
# ------------------------------------------------------------

counts_raw <- read.delim(
  count_file,
  comment.char = "#",
  check.names = FALSE
)

# featureCounts contains annotation columns before sample columns
sample_names <- c(
  "SRR031708",
  "SRR031712",
  "SRR031713",
  "SRR031720",
  "SRR031721",
  "SRR031722"
)

counts <- counts_raw[, c("Geneid", sample_names)]

rownames(counts) <- counts$Geneid
counts$Geneid <- NULL

counts <- as.matrix(counts)

# Make sure counts are numeric/integer
counts <- round(counts)

# ------------------------------------------------------------
# 3. Read metadata
# ------------------------------------------------------------

metadata <- read.csv(
  metadata_file,
  row.names = 1,
  stringsAsFactors = FALSE
)

metadata$condition <- factor(
  metadata$condition,
  levels = c("Untreated", "RNAi")
)

# Make sure sample order is identical
metadata <- metadata[colnames(counts), , drop = FALSE]

stopifnot(
  identical(
    rownames(metadata),
    colnames(counts)
  )
)

# ------------------------------------------------------------
# 4. Create DESeq2 dataset
# ------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ condition
)

# ------------------------------------------------------------
# 5. Low-count filtering
# ------------------------------------------------------------

keep <- rowSums(counts(dds) >= 10) >= 3

dds <- dds[keep, ]

cat("Genes before filtering:", nrow(counts), "\n")
cat("Genes after filtering:", nrow(dds), "\n")

# ------------------------------------------------------------
# 6. Run DESeq2
# ------------------------------------------------------------

dds <- DESeq(dds)

# ------------------------------------------------------------
# 7. Differential expression results
# ------------------------------------------------------------

res <- results(
  dds,
  contrast = c("condition", "RNAi", "Untreated"),
  alpha = 0.1
)

res_df <- as.data.frame(res)

res_df$gene <- rownames(res_df)

# Reorder columns
res_df <- res_df[, c(
  "gene",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)]

# Save complete DESeq2 results
write.csv(
  res_df,
  file.path(result_dir, "DESeq2_all_results.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Define DEGs
# ------------------------------------------------------------

DEG <- subset(
  res_df,
  !is.na(padj) &
    padj < 0.1 &
    abs(log2FoldChange) > 1
)

DEG$direction <- ifelse(
  DEG$log2FoldChange > 0,
  "Up",
  "Down"
)

# ------------------------------------------------------------
# 9. Gene annotation
# ------------------------------------------------------------

gene_map_raw <- read.delim(
  file.path(
    project_dir,
    "reference/BDGP5.25.62/Drosophila_melanogaster.BDGP5.25.62.gtf"
  ),
  header = FALSE,
  comment.char = "#",
  sep = "\t",
  stringsAsFactors = FALSE
)

# Extract gene_id and gene_name from GTF attributes
extract_attribute <- function(x, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  m <- regexec(pattern, x)
  regmatches(x, m)[, 2]
}

gene_ids <- extract_attribute(
  gene_map_raw$V9,
  "gene_id"
)

gene_names <- extract_attribute(
  gene_map_raw$V9,
  "gene_name"
)

gene_map <- data.frame(
  gene = gene_ids,
  gene_name = gene_names,
  stringsAsFactors = FALSE
)

gene_map <- gene_map[
  !duplicated(gene_map$gene),
]

# Add annotation
res_df_annotated <- merge(
  res_df,
  gene_map,
  by = "gene",
  all.x = TRUE
)

DEG_annotated <- merge(
  DEG,
  gene_map,
  by = "gene",
  all.x = TRUE
)

# ------------------------------------------------------------
# 10. Save annotated results
# ------------------------------------------------------------

write.csv(
  res_df_annotated,
  file.path(
    result_dir,
    "DESeq2_all_results_annotated.csv"
  ),
  row.names = FALSE
)

write.csv(
  DEG_annotated,
  file.path(
    result_dir,
    "DEG_final.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. Print summary
# ------------------------------------------------------------

cat("\n===== DESeq2 SUMMARY =====\n")
cat("Genes tested:", nrow(res_df), "\n")
cat("DEGs:", nrow(DEG), "\n")
cat("Upregulated:", sum(DEG$direction == "Up"), "\n")
cat("Downregulated:", sum(DEG$direction == "Down"), "\n")

# ------------------------------------------------------------
# 12. PCA
# ------------------------------------------------------------

vsd <- vst(dds, blind = FALSE)

pcaData <- plotPCA(
  vsd,
  intgroup = "condition",
  returnData = TRUE
)

percentVar <- round(
  100 * attr(pcaData, "percentVar")
)

pca_plot <- ggplot(
  pcaData,
  aes(
    x = PC1,
    y = PC2,
    color = condition,
    label = name
  )
) +
  geom_point(size = 4) +
  geom_text(
    vjust = -0.8,
    size = 3
  ) +
  xlab(
    paste0("PC1: ", percentVar[1], "% variance")
  ) +
  ylab(
    paste0("PC2: ", percentVar[2], "% variance")
  ) +
  theme_classic()

ggsave(
  file.path(plot_dir, "PCA_RNAi_vs_Untreated.png"),
  pca_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 13. MA plot
# ------------------------------------------------------------

png(
  file.path(
    plot_dir,
    "MA_RNAi_vs_Untreated.png"
  ),
  width = 2000,
  height = 1600,
  res = 300
)

plotMA(
  res,
  ylim = c(-5, 5),
  alpha = 0.1,
  main = "RNAi vs Untreated"
)

dev.off()

# ------------------------------------------------------------
# 14. Volcano plot
# ------------------------------------------------------------

volcano_df <- res_df_annotated

volcano_df$significance <- "NS"

volcano_df$significance[
  !is.na(volcano_df$padj) &
    volcano_df$padj < 0.1 &
    volcano_df$log2FoldChange > 1
] <- "Up"

volcano_df$significance[
  !is.na(volcano_df$padj) &
    volcano_df$padj < 0.1 &
    volcano_df$log2FoldChange < -1
] <- "Down"

volcano_df$neglog10padj <- -log10(
  pmax(volcano_df$padj, .Machine$double.xmin)
)

volcano_plot <- ggplot(
  volcano_df,
  aes(
    x = log2FoldChange,
    y = neglog10padj
  )
) +
  geom_point(
    aes(color = significance),
    alpha = 0.6,
    size = 1.5
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.1),
    linetype = "dashed"
  ) +
  labs(
    title = "RNAi vs Untreated",
    x = "log2 fold change",
    y = "-log10 adjusted p-value"
  ) +
  theme_classic()

ggsave(
  file.path(
    plot_dir,
    "Volcano_RNAi_vs_Untreated.png"
  ),
  volcano_plot,
  width = 7,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 15. Heatmap of DEGs
# ------------------------------------------------------------

deg_ids <- DEG$gene

if (length(deg_ids) > 1) {

  heatmap_mat <- assay(vsd)[
    intersect(deg_ids, rownames(vsd)),
    ,
    drop = FALSE
  ]

  # Limit visualization to the strongest DEGs
  heatmap_genes <- head(
    rownames(
      DEG[
        order(DEG$padj),
        ]
    ),
    50
  )

  heatmap_genes <- intersect(
    heatmap_genes,
    rownames(heatmap_mat)
  )

  heatmap_mat <- heatmap_mat[
    heatmap_genes,
    ,
    drop = FALSE
  ]

  annotation_col <- data.frame(
    condition = metadata$condition
  )

  rownames(annotation_col) <- rownames(metadata)

  png(
    file.path(
      plot_dir,
      "DEG_heatmap.png"
    ),
    width = 2200,
    height = 2400,
    res = 300
  )

  pheatmap(
    heatmap_mat,
    scale = "row",
    annotation_col = annotation_col,
    show_rownames = TRUE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize = 8,
    main = "Top DEGs"
  )

  dev.off()
}

cat("\nAnalysis completed successfully.\n")
cat("Results:", result_dir, "\n")
cat("Plots:", plot_dir, "\n")
