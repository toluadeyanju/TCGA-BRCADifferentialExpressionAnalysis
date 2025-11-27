###############################################################################
# Project: TCGA-BRCADifferentialExpressionAnalysis
# Author: Toluwanimi Adeyanju (toluwanimiadeyanju@gmail.com; ORCID:0000-0002-4433-0496)
# Date: 2025-11-27
# R Version: 4.5.2
#
# Description:
#   Processes publicly available TCGA-BRCA data:
#     1. Download / prepare counts & clinical metadata (TCGAbiolinks)
#     2. Differential expression analysis (DESeq2)
#     3. GO and KEGG enrichment (clusterProfiler)
###############################################################################
###############################################
# PACKAGES
###############################################
#############################
# 0. Package installation   #
#############################

# Bioconductor installation
install.packages("BiocManager")

# Core Bioc packages
BiocManager::install(c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "DESeq2",
  "limma",
  "edgeR",
  "org.Hs.eg.db",
  "clusterProfiler"
))

# CRAN visualization packages
install.packages(c("ggplot2", "pheatmap", "RColorBrewer", "dplyr", "stringr"))

# EnhancedVolcano (Bioc or GitHub if needed)
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  BiocManager::install("EnhancedVolcano")
}

#############################
# 1. Load libraries         #
#############################

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(limma)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(EnhancedVolcano)
library(dplyr)
library(stringr)

#############################
# 2. Download / prepare data#
#############################

# Query TCGA-BRCA HTSeq counts
query <- GDCquery(
  project      = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts"
)

GDCdownload(query)
data <- GDCprepare(query)   # SummarizedExperiment

# (Optional) Save object for reuse
dir.create("data", showWarnings = FALSE)
save(data, file = "data/TCGA_BRCA_SE.rda")

# To load later:
# load("data/TCGA_BRCA_SE.rda")

#############################
# 3. Inspect data           #
#############################

# Available assays (e.g., raw counts)
assayNames(data)

# Clinical / metadata snapshot
colData(data)[1:5, 1:10]

#############################
# 4. Prepare DESeq2 object  #
#############################

# Check subtype column
colnames(colData(data))

# Filter samples with PAM50 subtype annotation
subtype_col <- "paper_BRCA_Subtype_PAM50"

brca_filtered <- data[, !is.na(colData(data)[[subtype_col]])]

# Set factor levels for PAM50
colData(brca_filtered)[[subtype_col]] <-
  factor(colData(brca_filtered)[[subtype_col]],
         levels = c("LumA", "LumB", "Her2", "Basal", "Normal"))

# Build DESeqDataSet
dds <- DESeqDataSet(brca_filtered, design = ~ paper_BRCA_Subtype_PAM50)

# Filter low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2
dds <- DESeq(dds)

#############################
# 5. Differential expression#
#############################

# Contrast: Basal vs LumA (Basal - LumA)
res <- results(dds, contrast = c("paper_BRCA_Subtype_PAM50", "Basal", "LumA"))
res <- res[order(res$padj), ]

head(res)

#############################
# 6. Gene annotation        #
#############################

# Clean Ensembl IDs (remove version suffix like ".1")
clean_ids <- gsub("\\..*", "", rownames(res))

# Add gene symbols
res$symbol <- mapIds(org.Hs.eg.db,
                     keys     = clean_ids,
                     column   = "SYMBOL",
                     keytype  = "ENSEMBL",
                     multiVals = "first")

head(res[, c("symbol", "log2FoldChange", "padj")])

#############################
# 7. Volcano plot           #
#############################

p <- EnhancedVolcano(
  res,
  lab              = res$symbol,
  x                = "log2FoldChange",
  y                = "padj",
  xlab             = bquote(~Log[2]~ "fold change"),
  ylab             = bquote(~-Log[10]~ "adjusted p-value"),
  title            = "Basal vs LumA (PAM50)",
  subtitle         = "DESeq2 results",
  pCutoff          = 0.05,
  FCcutoff         = 1.5,
  pointSize        = 2.0,
  labSize          = 4.0,
  colAlpha         = 1
)

print(p)

ggplot2::ggsave("volcano_tall_selectlabels2.png",
                plot = p,
                width = 18,
                height = 14,   # increase height to spread y axis
                dpi = 300)

#############################
# 8. Heatmap of top genes   #
#############################
library(SummarizedExperiment)

# Top 50 DE genes
topGenes <- head(order(res$padj), 50)

# Variance-stabilizing transformation for counts
vsd <- vst(dds, blind = FALSE)
expr_matrix <- assay(vsd)[topGenes, ]

# Sample annotation (PAM50 subtype)
annotation_col <- data.frame(
  Subtype = colData(brca_filtered)$paper_BRCA_Subtype_PAM50
)
rownames(annotation_col) <- colnames(expr_matrix)

# Add gene symbols to rows
gene_ids <- gsub("\\..*", "", rownames(expr_matrix))
symbols <- mapIds(org.Hs.eg.db,
                  keys     = gene_ids,
                  column   = "SYMBOL",
                  keytype  = "ENSEMBL",
                  multiVals = "first")

rownames(expr_matrix) <- ifelse(!is.na(symbols), symbols, rownames(expr_matrix))

# Plot heatmap
pheatmap(
  expr_matrix,
  annotation_col          = annotation_col,
  show_rownames           = TRUE,
  fontsize_row            = 6,
  fontsize_col            = 8,
  scale                   = "row",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main                    = "Top 50 Differentially Expressed Genes (Basal vs LumA)"
)

#############################
# 9. Functional enrichment  #
#############################

# Significant genes for enrichment
sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

# Convert Ensembl → Entrez
gene_ids_sig <- gsub("\\..*", "", rownames(sig_genes))
entrez <- mapIds(org.Hs.eg.db,
                 keys     = gene_ids_sig,
                 column   = "ENTREZID",
                 keytype  = "ENSEMBL",
                 multiVals = "first")
entrez <- na.omit(entrez)

## 9a. GO Biological Process
ego <- enrichGO(
  gene          = entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

dotplot(ego, showCategory = 15, title = "GO Biological Processes: Basal vs LumA")

## 9b. KEGG pathways
ekegg <- enrichKEGG(
  gene         = entrez,
  organism     = "hsa",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.5
)

p <- dotplot(ekegg, showCategory = 15, title = "KEGG Pathways: Basal vs LumA") +
  theme(
    axis.text.y = element_text(size = 4, lineheight = 1, margin = margin(r = 15)),
    plot.title  = element_text(size = 10, face = "plain")
  )

print(p)

ggplot2::ggsave("KEGG_dotplot_spaced.png", p,
                width = 10, height = 14, dpi = 300)

#############################
# 10. Combined GO+KEGG plot #
#############################

# Convert to data frames
ego_df   <- as.data.frame(ego)   %>% mutate(Source = "GO_BP")
ekegg_df <- as.data.frame(ekegg) %>% mutate(Source = "KEGG")

# Combine
combined <- bind_rows(ego_df, ekegg_df) %>%
  select(ID, Description, GeneRatio, Count, p.adjust, Source)

# Convert GeneRatio "a/b" to numeric
combined$GeneRatio <- sapply(
  strsplit(as.character(combined$GeneRatio), "/"),
  function(x) as.numeric(x[1]) / as.numeric(x[2])
)

# Top 10 per source
combined_top <- combined %>%
  group_by(Source) %>%
  slice_min(order_by = p.adjust, n = 10) %>%
  ungroup() %>%
  arrange(p.adjust)

# Wrap long descriptions
combined_top$Description <- str_wrap(combined_top$Description, width = 40)

# Bubble plot
ggplot(combined_top, aes(
  x     = GeneRatio,
  y     = reorder(Description, GeneRatio),
  color = Source,
  size  = Count
)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("GO_BP" = "#1f77b4", "KEGG" = "#ff7f0e")) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme_minimal(base_size = 13) +
  labs(
    title = "GO and KEGG Enrichment Summary (Basal vs LumA)",
    x     = "Gene Ratio",
    y     = NULL,
    color = "Database"
  ) +
  theme(
    plot.title  = element_text(hjust = 0.5, size = 12, face = "bold"),
    axis.text.y = element_text(size = 5, lineheight = 1, margin = margin(r = 20)),
    axis.text.x = element_text(size = 10),
    legend.position = "right"
  )

