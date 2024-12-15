# Spatial Head and Neck Squamous Cell Carcinoma Paper

This repository contains the data and code for our paper:

> Dor Simkin, Thomas F Barrett ... Sidharth V Puram, Itay Tirosh,\
> (2024). *Spatial analysis of head and neck cancer reveals two ecosystems with distinct modes of epithelial-to-mesenchymal transition*. NOT PUBLISHED.

## Contents

This **repository** contains:

-   [:file_folder: analysis_scripts](/analysis_scripts):
    -   [:file_folder: functions](/analysis_scripts/functions):\
        Functions containing R scripts. These functions are sourced as needed in each of the analysis scripts.
    -   Individual analysis scripts. The flow of the analysis is as follows:
        1.  **MP_Generation.R** - Run Leiden and NMF clustering, pool them and recluster into metaprograms (MPs).
        2.  **Spot_Assignment.R** - Label spots by their dominant MP identity and their composite MP proportions (by scoring-based deconvolution), and compute chromosomal copy-number aberrations (CNAs).
        3.  **State_Proximity_and_Zonation.R** - Define per-sample malignant zones, compute MP enrichment over each zone and summarize their per-sample and overall spatial distribution. Characterize tumor's growth pattern as measured by tumor coherence (the degree of spatial clustering of malignant spots).
        4.  **Overall_Survival_Analysis.R** - Calculate the 5-year overall survival associated with each MP expression using Cox proportional hazards regression.
        5.  **TMA_Analysis.R** - Assess the agreement of p-EMT patterns across different TMA cores belonging to the same tumor (patient).
        6.  **CellLines_RNAseq_Analysis.R** - Process, normalize, and run differential gene expression analysis of bulk RNAseq from two cell lines (JHU006, 93VU-147T), each treated with several different conditions (OSM, TGFB, IL1B, SPP1).
        7.  **Figures_Generation.R** - Generate the remaining paper figures using the results files.
-   [:file_folder: config_scripts](/config_scripts): Configuration scripts:
    1.  **get_data.R** - Fetch spatial transcriptomic (10X Visium) data from GEO, download, and arrange it in a proper folder structure for downstream analysis.
    2.  **get_data_REVIEWERS.R** - Temporary script for review purposes ([see instructions below](#for-reviewers))
    3.  **get_supp_files.R** - Download supplementary files generated in the analysis steps and organize them in their correct folder structure.
    4.  **restore_r_env_func.R** - Helper function to restore the R environment provided, suited for downstream analysis. This is used to automatically install all required packages.
-   [:file_folder: metadata](/metadata): Contains:
    -   Sample-level metadata (`samples_metadata.rds`).
    -   Spot-level metadata (`merged_metadata_Extend.rds`).
    -   Sample parameters (`per_sample_parameters_df.rds`) file with selected values for QC, dimension reduction, and clustering resolution.
    -   Cell-level metadata (`merged_scRNAseq_metadata.rds`) for the matched scRNAseq samples.
-   [:file_folder: renv](/renv) and `renv.lock` file: Required for the restoration of the R environment set for this analysis.

## How to download the repository and run locally

**Requirements**:\
- [R software](https://cloud.r-project.org/) version 4.1.1 or above installed.\
- [RStudio IDE](https://rstudio.com/products/rstudio/download/).

You can download the repository as a ZIP from this URL:\
[Download master.zip](https://github.com/DSimkin/Spatial.HNSCC/archive/refs/heads/main.zip).

After unzipping:

-   Open the `.Rproj` file in RStudio and follow these configuration steps:
    1.  Open the [`get_data.R` script](/config_scripts/get_data.R) and source it (Command-Shift-S [Mac] or Control-Shift-S [Windows]).\
        This script installs all necessary R packages and downloads the spatial transcriptomic data from GEO.\
        <a name="for-reviewers"></a> **FOR REVIEWERS**: The data is stored on GEO and is currently private prior to publication. The secure token is provided in the manuscript under the "Data Availability" section.\
        To access the data, search for GEO accession **GSE268014** on the [Gene Expression Omnibus website](https://www.ncbi.nlm.nih.gov/geo/). You will be asked to provide the secure token, then redirected to the accession page. Download the file `GSE268014_RAW.tar` and place it in the main directory (where the `.Rproj` file is located). After downloading, open the [`get_data_REVIEWERS.R` script](/config_scripts/get_data_REVIEWERS.R) and source it.

    2.  Open the [`get_supp_files.R` script](/config_scripts/get_supp_files.R) and source it.\
        This script downloads all `.rds`, `.csv`, `.tsv`, and supplementary files generated during the analysis to a newly created folder called `results`.
-   Once the configuration steps are complete, open the [analysis scripts](/analysis_scripts) and run them in the order listed in the **Contents** section to inspect and generate the paper results.
