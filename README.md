Indoor Microbiome and Lung Health Study
This repository contains the R code and associated input/output files used in the analysis for the paper:

"Indoor Airborne Bacterial Exposure and Lung Health in Adults: A Cross-Sectional Study"

Contents
scripts/: Full R pipeline for data processing, analysis, and figure generation.

data/: Participant metadata and small supporting datasets.

figures/: Final plots and heatmaps generated for the manuscript.

results/: Output tables from multivariable models and Mantel tests.

Reproducibility
This repository includes:

Bacterial diversity and load analysis (Shannon diversity, qPCR, endotoxin)

CLR transformation and scaling by interquartile range (IQR)

Multivariable linear regression models stratified by sex

Beta-diversity association testing (Mantel tests)

Visualizations including heatmaps and association plots

Note: Large microbiome data files (e.g., full OTU tables) were too large to upload to GitHub and are available separately through Zenodo. See "Data Availability" below.

Data Availability
Large microbiome datasets are available from Zenodo:

🔗 Access full datasets via Zenodo

Small files such as metadata are provided directly in the data/ folder.

Citation
Please cite the paper if using this code or analysis workflow in your research:

Indoor Airborne Bacterial Exposure and Lung Health in Adults: A Cross-Sectional Study
Author names, Year.
