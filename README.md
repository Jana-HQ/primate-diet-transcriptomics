# Primate Diet Transcriptomics

## Overview
This repository contains an end-to-end bioinformatics pipeline analyzing the relationship between primate diets and gene expression profiles. By processing transcriptomic data, this project identifies statistically significant over-expressed and under-expressed genes associated with different dietary habits using differential expression analysis and dimensionality reduction techniques.

## Pipeline Architecture
The analysis is broken down into three modular scripts, simulating a standard reproducible bioinformatics workflow:

1. **Preprocessing & EDA (`01_preprocess_and_eda.R`):** Fetches the raw matrix, applies log-transformation, visualizes data distributions, and establishes baseline sample relationships via hierarchical clustering.
2. **Differential Expression (`02_differential_expression.R`):** Conducts two-sample T-tests and calculates Fold Changes across dietary categories. Applies stringent biological thresholds (Fold Change >= 1.0, P-value <= 0.05) to isolate significant genes, outputting volcano plots and dendrogram-clustered heatmaps.
3. **Dimensionality Reduction (`03_pca_analysis.R`):** Performs Principal Component Analysis (PCA) to cluster samples by diet and extracts the highest-loading genes driving biological variance.

## Biological Insights & Enrichment Analysis
Downstream enrichment analysis of the Differentially Expressed Genes (DEGs) revealed significant shifts in metabolic pathways based on diet:
* **Human Cafeteria Diet:** Showed a significant downregulation in **Long-Chain Fatty Acid Metabolism**, specifically reducing the liver's ability to oxidize long-chain fatty acids in the mitochondria.
* **General Pathway Impacts:** Significant variations were observed in **Lipid & Sterol Biosynthesis**, **Acyl-CoA metabolic processes**, and **Small molecule biosynthetic processes** (e.g., GO:0044283).
* **Presentation:** For a deep dive into the biological interpretation, PCA clustering, and pathway analysis, please view the [Project Presentation (PDF)](docs/Presentation_DGE_Mice_Diets.pdf) included in this repository.

## Repository Structure
```text
primate-diet-transcriptomics/
├── data/                   # Generated intermediate data (.rds)
├── outputs/                # Generated visual outputs (PDFs, CSVs)
├── docs/                   # Project report and presentation materials
├── scripts/
│   ├── 01_preprocess_and_eda.R
│   ├── 02_differential_expression.R
│   └── 03_pca_analysis.R
└── README.md
```

## Sources
* **Primary Study:** Weng K, Hu H, Xu AG, Khaitovich P, Somel M. Mechanisms of dietary response in mice and primates: a role for EGR1 in regulating the reaction to human-specific nutritional content. PLoS One. 2012;7(8):e43915. doi: 10.1371/journal.pone.0043915. Epub 2012 Aug 24. PMID: 22937124; PMCID: PMC3427207.
* **Dataset:** Moustafa, A. Primates Diet Gene Expression Dataset. Retrieved from [GitHub](https://github.com/ahmedmoustafa/gene-expression-datasets/blob/main/datasets/primates_diet/primates_diet.tsv).