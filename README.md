DNA Methylation Analysis in Plasmablastic Lymphoma



Overview:

Plasmablastic lymphoma (PBL) is a rare and highly aggressive non-Hodgkin lymphoma derived from B cells in a plasmablast state, frequently associated with immunodeficiency, Epstein–Barr virus infection, and alterations in the MYC gene. Despite advances in its genetic characterization, the epigenetic mechanisms underlying its biological heterogeneity remain poorly understood.
This project was developed as part of a Bachelor's thesis at the IDIBAPS Epigenetics group and focuses on the analysis of DNA methylation profiles in PBL using R-based bioinformatics approaches.



Objective:

The main objective of this study is to characterize epigenetic alterations in plasmablastic lymphoma and evaluate whether DNA methylation patterns can define distinct biological subgroups.



Data and Methods:

DNA methylation data generated using Illumina arrays from PBL patient samples were analyzed using computational and statistical approaches in R.
The workflow includes:
- Data preprocessing and quality control
- Identification of differentially methylated CpG sites
- Unsupervised analysis of methylation profiles
- Estimation of proliferative history using the epiCMIT mitotic clock
- Comparative analysis with normal B-cell differentiation stages



Key Analyses:

- Unsupervised clustering to explore epigenetic heterogeneity
- Differential methylation analysis
- Variance-based feature selection
- Heatmap visualization of methylation patterns
- Integration with biological context of B-cell differentiation



Main Findings:

The analysis reveals a high degree of epigenetic heterogeneity among PBL samples, largely associated with differences in accumulated proliferative history.
Key observations include:
- Predominantly stochastic DNA methylation changes
- Lack of strong enrichment in canonical regulatory regions (e.g., promoters or CpG islands)
- Overlap between tumor methylation patterns and physiological B-cell differentiation trajectories
- Evidence that proliferative history plays a major role in shaping the PBL methylome



Institution: IDIBAPS – Epigenetics Research Group

Author: Sofia Lombera Tolrà (Bachelor’s Thesis in Bioinformatics and Epigenetics)

Notes: This repository contains scripts for reproducible analysis of DNA methylation data, including unsupervised analysis, differential methylation, and visualization workflows.c
