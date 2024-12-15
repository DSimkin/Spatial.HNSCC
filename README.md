---
editor_options: 
  markdown: 
    wrap: 72
---

# Spatial Head and Neck Squamous Cell Carcinoma Paper

This repository contains the data and code for our paper:

> Dor Simkin, Thomas F Barrett ... Sidharth V Puram, Itay Tirosh,
> (2024). *Spatial analysis of head and neck cancer reveals two
> ecosystems with distinct modes of epithelial-to-mesenchymal
> transition*. NOT PUBLISHED.

## Contents

This **repository** contains:

-   [:file_folder: analysis_scripts](/analysis_scripts):

    -   [:file_folder: functions](/analysis_scripts/functions):
        Functions containing R scripts. These functions are sourced as
        needed in each of the analysis scripts.
    -   Individual analysis scripts. The flow of the analysis is as
        such:
        1.  [MP_Generation.R]{.underline} - run Leiden and NMF
            clustering, pool them and recluster them into metaprograms
            (MPs).
        2.  [Spot_Assignment.R]{.underline} - label spots by their
            dominant MP identity and their composite MP proportions (by
            scoring-based deconvolution), and compute chromosomal
            copy-number aberrations (CNAs).
        3.  [State_Proximity_and_Zonation.R]{.underline} - define
            per-sample malignant zones, compute MPs enrichment over each
            zone and summarize their per-sample and over all spatial
            distribution, and characterize tumor's growth pattern as
            measured by tumor coherence (the degree of spatial
            clustering of malignant spots).
        4.  [Overall_Survival_Analysis.R]{.underline} - calculate the
            5-year overall survival associated with each MP expression
            using Cox proportional hazards regressior.
        5.  [TMA_Analysis.R]{.underline} - assess the agreement of p-EMT
            patterns across the different TMA cores belonging to the
            same tumor (patient).
        6.  [CellLines_RNAseq_Analysis.R]{.underline} - process,
            normalize and run differential gene expression analysis of
            bulk RNAseq from 2 cell lines (JHU006, 93VU-147T), each
            treated with several different conditions (OSM, TGFB, IL1B,
            SPP1).
        7.  [Figures_Generation.R]{.underline} - script that uses the
            generated results files and create the remaining paper
            figures.

-   [:file_folder: config_scripts](/config_scripts): configuration
    scripts:

    1.  [get_data.R]{.underline} - a script that fetchs spatial
        transcriptominc (10X Visium) GEO url, download the data and
        arrange it in a proper folder structure for downstream analysis.
    2.  [get_supp_files.R]{.underline} - a script that gets figshare's
        collection url, download supplementary files that have been
        generated in the analysis steps and order them in their correct
        folder structure.
    3.  [restore_r_env_func.R]{.underline} - helper function to restore
        the R environment provided which is suited for the downstream
        analysis. It is used for the automatic installation of all
        required packages.

-   [:file_folder: metadata](/metadata): contains sample-level
    (samples_metadata.rds) and spot-level metadata
    (merged_metadata_Extend.rds) tables, sample parameters
    (per_sample_parameters_df.rds) file with selected values for QC,
    dimension reduction and clustering resolution, and cell-level
    metadata (merged_scRNAseq_metadata.rds) for the matched scRNAseq
    samples.

-   [:file_folder: renv](/renv) and renv.lock file: required for the
    restoration of the R environment set for this analysis.

## How to download the repository and run locally

[**Requirements**]{.underline}: [R
software](https://cloud.r-project.org/) version 4.1.1 or above
installed. [RStudio
IDE](https://rstudio.com/products/rstudio/download/).

You can download the repository as a zip from from this URL:
[master.zip](https://github.com/DSimkin/Spatial.HNSCC/archive/refs/heads/main.zip).
After unzipping: - open the `.Rproj` file in RStudio and follow these
configuration steps:

1.  Open the [`get_data.R` script](/config_scripts/get_data.R) and
    source run it (type Command-Shift-S [in Mac], or Control-Shift-S [in
    Windows]). This script both installs all the necessary R packages
    and download the spatial transcriptomic data from GEO.\
    [**FOR REVIEWERS**]{.underline}: the data is stored on GEO and is
    currently private, prior to publication. The secure token is
    provided in the manuscript under the section Data Availability. To
    access the data, search this GEO accession - GSE268014, on the [Gene
    Expression Omnibus website](https://www.ncbi.nlm.nih.gov/geo/). You
    will be asked to provide the secure token and will be then
    redirected to the accession display page. On the bottom of that
    page, the file `GSE268014_RAW.tar` should be downloaded. Make sure
    you

### 
