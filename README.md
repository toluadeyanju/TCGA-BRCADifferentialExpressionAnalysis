TCGA-BRCA Differential Expression Analysis (PAM50 Subtypes)
Exploring differences between Basal and Luminal A subtypes

This repository contains the R script used to perform differential genetic expression analysis and functional enrichment analysis using publicly available data from the The Cancer Genome Atlas Project.

The primary script tcga_brca_de_analysis.R executes the following steps
1. Data Acquisition: Using the TCGABiolinks package to query and download gene counts and associated clinical metadata for TCGA BRCA dataset
2. Data Preparation: Filters samples based on the availability of the paper_BRCA_Subtype_PAM50 annotation and removes genes with low overall counts.
3. Differential Expression: Performs Differential Genetic Expression Analysis (DGEA) using the DESeq2 package, comparing the Basal-like subtype against the Luminal A (LumA) subtype.
4. Gene Annotation: Maps Ensembl Gene IDs to HUGO Gene Symbols using the org.Hs.eg.db annotation package.
5. Visualization: Generates a Volcano Plot and a Heatmap to visualize the differentially expressed genes (DEGs).
6. Functional Enrichment: Performs Gene Ontology (GO) Biological Process and KEGG Pathway enrichment analysis on the significant DEGs using clusterProfiler.

Prerequisites and Installation
The script requires a combination of Bioconductor and CRAN packages.

R Version
Tested with R version 4.5.2.

To save time and prevent redownload of the dataset anytime the dataset is run, the TCGA data was saved locally on the device to enable faster loading on next analysis using this: (save(data, file = "data/TCGA_BRCA_SE.rda") 

Expected Outputs
The script generates several high-quality plots essential for interpreting the differential expression results. These plots are saved in the plots/ directory.

1. Volcano Plot
Visualizes the significance ($\text{-Log}_{10}$ adjusted p-value) versus the magnitude ($\text{Log}_2$ Fold Change) of all expressed genes in the Basal vs. LumA contrast.
File: plots/volcano_plot_basal_vs_luma.png
Cutoffs: |Log2FC| > 1.5 and p-adjusted < 0.05.

2. Heatmap
Displays the expression of the top 50 differentially expressed genes across samples, clustered by gene expression similarity. The top annotation bar shows the PAM50 subtype for each sample.
File: plots/heatmap_top50_basal_vs_luma.png
Data: Variance-stabilizing transformed (vst) counts, scaled by row.

3. Functional Enrichment Plots

  A. GO Biological Process Dot Plot

  File: plots/GO_dotplot_basal_vs_luma.png

  Shows the top 15 enriched Gene Ontology Biological Process terms.

  B. KEGG Pathway Dot Plot

  File: plots/KEGG_dotplot_basal_vs_luma.png

  Shows the top 15 enriched KEGG pathways.

  C. Combined Summary Bubble Plot (GO & KEGG)

  File: plots/combined_enrichment_basal_vs_luma.png

  Combines the top 10 most significant terms from both GO-BP and KEGG into a single bubble plot for easy comparison. The bubble size represents the number of genes found       (Count), and the position represents the Gene Ratio (proportion of DEGs in the pathway).
