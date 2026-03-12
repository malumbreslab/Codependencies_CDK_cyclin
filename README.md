# Codependencies_CDK_Cyclin

## Description

This github repository contains the analysis code used in "A reference map of CDK-cyclin dependencies across human cancers" published in bioRxiv:

Link:

## Table of contents

- [Data](#Data)
- [Scripts](#Scripts)


## Data

Data including metadata, transcriptomics, gene effect, copy number and mutations were downloaded from DepMap 25Q3(Broad Institute, 2025) and PRISM screening data from from DepMap 19Q4 (Corsello et al., 2020). All data needed to evaluate the conclusions in the paper are present in the paper and/or the Supplementary Materials.

## Scripts

A. Codependencies and clustering

- [A_Codependencies_and_clustering.R](/scripts/A_Codependencies_and_clustering.R): Analysis of codependencies using gene effect data, clustering of cell lines and characterization of each cluster, extracting expression markers.

B. PRISM screening
   
- [B_PRISM_screening.R](/scripts/B_PRISM_screening.R): Analysis of CDKi profiles in cell lines.
  
C. Efficacy and selectivity
- [C_Efficacy_and_selectivity.R](/scripts/C_Efficacy_and_selectivity.R): Using previous studies, analysis of selectivity and efficacy of each geens in each cluster o group of interest.
  
D. Comparisons
- [D_Comparisons.R](/scripts/D_Comparisons.R): Comparison and plotting between groups of interest and clusters.

E. Funtional analysis
- [E_GSEA_analysis_only.ipynb](/scripts/E_GSEA_analysis_only.ipynb): Funtional analysis using GSEApy.

F. TCGA validation
- [F_TCGA_validation.ipynb](/scripts/F_TCGA_validation.ipynb): Validation of pancancer markers from different clusters in TCGA data.

G. GSE validation
- [G_GSE_validation.ipynb](/scripts/G_GSE_validation.ipynb): Validation of breast cancer markers from different clusters in BRCA in published data.

## Authors

All uploaded scripts have been done by:

- Marcos Malumbres
- Agustín Sánchez-Belmonte
