library(tidyverse)
library(scalop)
library(here)
library(survival)
library(survminer)
source(here("analysis_scripts/functions/scRNA_Functions.R"))

# Load metaprograms
MPs <- readRDS(file = here("results/Generated_Data/Final_Metaprograms_Extended.rds"))


# Process TCGA RNAseq and Clinical Data -----------------------------------

# Load clinical data
clinical_data <- openxlsx::read.xlsx(xlsxFile = here("aux_data/Liu_1-s2.0-S0092867418302290-mmc1.xlsx"), sheet = 1)[, -1] %>% 
  dplyr::filter(type == "HNSC")
## optionally add hpv status to the summarized clinical table
clinical_data2 <- read_tsv(here("aux_data/combined_study_clinical_data.tsv")) %>% 
  dplyr::filter(stringr::str_sub(`Sample ID`, 14, 15) == "01", `Study ID` == "hnsc_tcga") %>% 
  dplyr::select(`Patient ID`, `Hpv status p16`, `Hpv status ish`, `Patient Primary Tumor Site`)               # "01" is a code for primary tumor, this removes "06" which are recurrent
tcga_clinical_data <- merge(clinical_data, clinical_data2, by.x = "bcr_patient_barcode", by.y = "Patient ID") %>% 
  dplyr::rename("hpv_status_p16" = "Hpv status p16", "hpv_status_ish" = "Hpv status ish") %>% 
  dplyr::filter(!is.na(OS.time))

# Load RSEM scaled RNAseq counts for TCGA HNSC
tcga_bulk <- read.table(here("aux_data/HNSC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt"), sep = "\t", header = TRUE, row.names = 1)
tcga_bulk <- tcga_bulk[-1, seq(2, ncol(tcga_bulk), 3)] %>% dplyr::mutate(across(where(is.character), as.numeric)) %>% as.matrix()
colnames(tcga_bulk) <- gsub(pattern = "\\.", replacement = "-", colnames(tcga_bulk))

# Correct gene names and replace Entrez id to gene symbol
temp <- data.frame(name = rownames(tcga_bulk))
temp <- data.frame(do.call('rbind', strsplit(as.character(temp$name),'|',fixed = TRUE))) %>% dplyr::mutate(X3 = scalop::x2y("ENTREZID", "SYMBOL")(.$X2))
temp$X1[temp$X1 == "?"] <- temp$X3
table(duplicated(temp$X1[!is.na(temp$X1)]))
temp$X1[!is.na(temp$X1) & duplicated(temp$X1) == TRUE]
grep("SLC35E2", temp$X1)
temp$X1[16302] <- "SLC35E2_2"
temp$X1[is.na(temp$X1)] <- temp$X2[which(is.na(temp$X1))]
rownames(tcga_bulk) <- temp$X1

# Transformation to log2-TPM
log_tpm_tcga <- log2(1 + (tcga_bulk * 1000000))
# Filter lowly expressed genes
filt_tcga <- log_tpm_tcga[rowMeans(log_tpm_tcga) > 2, ]
hist(filt_tcga)

# Remove recurrent tumors and those that have no clinical information
table(stringr::str_sub(colnames(filt_tcga), 14, 15))
recurs <- colnames(filt_tcga)[stringr::str_sub(colnames(filt_tcga), 14, 15) == "06"]
normal_tissue <- colnames(filt_tcga)[stringr::str_sub(colnames(filt_tcga), 14, 15) == "11"]
filt_tcga <- filt_tcga[, !colnames(filt_tcga) %in% c(recurs, normal_tissue)]

# Remove non-exsisting RNAseq samples from the clinical data
tcga_clinical_data <- tcga_clinical_data[tcga_clinical_data$bcr_patient_barcode %in% stringr::str_sub(colnames(filt_tcga), 1, 12), ]

# Score bulk experssion matrix for spatial metaprograms
adj_MPs <- refine_programs(filt_tcga, programs = MPs, var_filter = FALSE, cor_filter = TRUE, cor_thresh = 0.5, conserved.genes = 0.5, return_cor_mat = TRUE, return_diagnostic_plot = TRUE)
tcga_MP_scores <- scalop::sigScores(m = filt_tcga, sigs = adj_MPs$programs, conserved.genes = 0.5)
rownames(tcga_MP_scores) <- stringr::str_sub(rownames(tcga_MP_scores), 1, 12)
# saveRDS(tcga_MP_scores, file = here("results/Generated_Data/survival_analysis_tcga_MP_scores.rds"))

# Control for HPV status by removing samples with HPV positive status either by p16 staining or ISH.
tcga_MP_scores <- readRDS(file = here("results/Generated_Data/survival_analysis_tcga_MP_scores.rds"))
patients2rm <- na.omit(union(tcga_clinical_data$bcr_patient_barcode[tcga_clinical_data$hpv_status_ish == "Positive"], 
                             tcga_clinical_data$bcr_patient_barcode[tcga_clinical_data$hpv_status_p16 == "Positive"]))
hpv_neg_tcga_data <- tcga_clinical_data[!tcga_clinical_data$bcr_patient_barcode %in% patients2rm, ]
tcga_hpv_neg_MP_cat <- tcga_MP_scores[rownames(tcga_MP_scores) %in% hpv_neg_tcga_data$bcr_patient_barcode, ] %>% 
  purrr::map(\(x) ifelse(x >= quantile(x, probs = 0.75), "High", 
                         ifelse(x <= quantile(x, probs = 0.25), "Low", NA))) %>% 
  as.data.frame %>% `rownames<-`(rownames(tcga_MP_scores[rownames(tcga_MP_scores) %in% hpv_neg_tcga_data$bcr_patient_barcode, ]))




# Process Leipzig_HN Microarray and Clinical Data -------------------------
library(illuminaHumanv4.db)

# Load clinical data
leip_clinical_data <- readxl::read_xls(here("aux_data/GSE65858_clinical_data.xls"))

# Load non-normalized micro-array data
leip_bulk <- read.table(here("aux_data/GSE65858_Non-normalized_data.txt"), sep = "\t", header = TRUE, row.names = 1)
detectionpvals <- leip_bulk[, grep('Detection.Pval', colnames(leip_bulk))]
leip_bulk <- data.matrix(leip_bulk[, -grep('Detection.Pval', colnames(leip_bulk))])

# Perform background correction on the fluorescent intensities and quantile normalization + log2 transformation (using limma package)
leip_bulk <- limma::neqc(leip_bulk, detection.p = detectionpvals, offset = 16)

# Convert prob-ID rownames to gene symbol
prob2gene <- data.frame(Gene = unlist(mget(x = rownames(leip_bulk), envir = illuminaHumanv4SYMBOL)))
rownames(leip_bulk) <- prob2gene[match(rownames(leip_bulk), rownames(prob2gene)), ]
leip_bulk_proc <- leip_bulk[!is.na(rownames(leip_bulk)), ] %>% fix_dup_genes(., stat = "mean")

# Score bulk experssion matrix for spatial metaprograms
adj_MPs <- refine_programs(leip_bulk_proc, programs = MPs, var_filter = FALSE, cor_filter = TRUE, cor_thresh = 0.5, conserved.genes = 0.5, return_cor_mat = TRUE, return_diagnostic_plot = TRUE)
leip_MP_scores <- scalop::sigScores(m = leip_bulk_proc, sigs = adj_MPs$programs, conserved.genes = 0.5)
# saveRDS(leip_MP_scores, file = here("results/Generated_Data/survival_analysis_leip_MP_scores.rds"))

# Control for HPV status by removing samples positive for HPV DNA
leip_MP_scores <- readRDS(file = here("results/Generated_Data/survival_analysis_leip_MP_scores.rds"))
hpv_neg_leip_data <- leip_clinical_data %>% dplyr::filter(`characteristics: HPV_DNA` == "Negative")
leip_hpv_neg_MP_cat <- leip_MP_scores[rownames(leip_MP_scores) %in% hpv_neg_leip_data$`Sample name`, ] %>% 
  purrr::map(\(x) ifelse(x >= quantile(x, probs = 0.75), "High", 
                         ifelse(x <= quantile(x, probs = 0.25), "Low", NA))) %>% 
  as.data.frame %>% `rownames<-`(rownames(leip_MP_scores[rownames(leip_MP_scores) %in% hpv_neg_leip_data$`Sample name`, ]))




# Process Fred Hutchinson Cancer Research Center Microarray and Clinical Data -------------------------
library(GEOquery)

# Load data and metadata from GEO server
fhcrc_obj <- getGEO("GSE41613", GSEMatrix = TRUE)
# View(pData(phenoData(fhcrc_obj[[1]])))
show(fhcrc_obj)
# Extract expression data and clinical metadata
expr_mat <- Biobase::exprs(fhcrc_obj[[1]])
fhcrc_clinical_data <- Biobase::pData(fhcrc_obj[[1]]) %>% 
  dplyr::select("geo_accession", "age:ch1", "fu time:ch1", "vital:ch1", "Sex:ch1", "tissue type:ch1", "tissue type:ch1", "treatment:ch1", "tumor stage:ch1") %>% 
  dplyr::rename_all(~scalop::substri(., sep = ":", pos = 1)) %>% 
  dplyr::rename_all(~gsub(" ", "_", .)) %>% 
  dplyr::mutate(ID = geo_accession, OS.time = as.numeric(fu_time) * (365 / 12), OS = as.numeric(ifelse(grepl("Alive", .$vital), 0, 1)), .keep = "unused")
prob2gene <- fhcrc_obj[[1]]@featureData@data %>% 
  dplyr::select(ID, "Gene Symbol", ENTREZ_GENE_ID) %>% 
  dplyr::mutate(across(everything(), ~na_if(., ""))) %>% 
  dplyr::mutate(Select_ENT = scalop::substri(ENTREZ_GENE_ID, sep = " /// ", pos = 1)) %>% 
  dplyr::mutate(Gene_Name = scalop::x2y("ENTREZID", "SYMBOL")(Select_ENT))
convert_gene_names <- prob2gene$Gene_Name[match(prob2gene$ID, rownames(expr_mat))]
expr_mat <- expr_mat[-which(is.na(convert_gene_names)), ]
rownames(expr_mat) <- convert_gene_names[!is.na(convert_gene_names)]

# Filter out lowly expressed genes
filt_expr_mat <- expr_mat[rowMeans(expr_mat) > 2, ] %>% fix_dup_genes(., stat = "mean")

# Score bulk experssion matrix for spatial metaprograms
adj_MPs <- refine_programs(filt_expr_mat, programs = MPs, var_filter = FALSE, cor_filter = TRUE, cor_thresh = 0.5, conserved.genes = 0.5, return_cor_mat = TRUE, return_diagnostic_plot = TRUE)
fhcrc_MP_scores <- scalop::sigScores(m = filt_expr_mat, sigs = adj_MPs$programs, conserved.genes = 0.5)
# saveRDS(fhcrc_MP_scores, file = here("results/Generated_Data/survival_analysis_fhcrc_MP_scores.rds"))

# Categorize scores only to low and high
fhcrc_MP_scores <- readRDS(file = here("results/Generated_Data/survival_analysis_fhcrc_MP_scores.rds"))
fhcrc_MP_cat <- fhcrc_MP_scores %>% 
  purrr::map(\(x) ifelse(x >= quantile(x, probs = 0.75), "High", 
                         ifelse(x <= quantile(x, probs = 0.25), "Low", NA))) %>% 
  as.data.frame %>% `rownames<-`(rownames(fhcrc_MP_scores))





# Merge all sample scores and clinical data -------------------------------

filt_tcga_clinical_data <- tcga_clinical_data %>% 
  dplyr::select(bcr_patient_barcode, OS.time, OS) %>% 
  dplyr::rename(ID = bcr_patient_barcode) %>% 
  `rownames<-`(NULL) %>% tibble::column_to_rownames(var = "ID")
filt_leip_clinical_data <- leip_clinical_data %>% 
  dplyr::select("characteristics: ID", "characteristics: OS", "characteristics: OS_EVENT") %>%
  dplyr::rename_all(~gsub("^characteristics: ", "", .)) %>%
  dplyr::mutate(OS.time = OS, OS = as.numeric(as.logical(OS_EVENT)), .before = 2) %>% 
  dplyr::select(-OS_EVENT) %>% 
  tibble::column_to_rownames(var = "ID")
filt_fhcrc_clinical_data <- fhcrc_clinical_data %>% 
  dplyr::select(ID, OS.time, OS) %>% 
  `rownames<-`(NULL) %>% tibble::column_to_rownames(var = "ID")

tcga_MP_scores <- readRDS(file = here("results/Generated_Data/survival_analysis_tcga_MP_scores.rds"))
leip_MP_scores <- readRDS(file = here("results/Generated_Data/survival_analysis_leip_MP_scores.rds"))
fhcrc_MP_scores <- readRDS(file = here("results/Generated_Data/survival_analysis_fhcrc_MP_scores.rds"))
tcga_hpv_neg_MP_cat <- tcga_MP_scores[rownames(tcga_MP_scores) %in% hpv_neg_tcga_data$bcr_patient_barcode, ] %>% 
  purrr::map(\(x) ifelse(x >= quantile(x, probs = 0.75), "High", 
                         ifelse(x <= quantile(x, probs = 0.25), "Low", NA))) %>% 
  as.data.frame %>% `rownames<-`(rownames(tcga_MP_scores[rownames(tcga_MP_scores) %in% hpv_neg_tcga_data$bcr_patient_barcode, ]))
leip_hpv_neg_MP_cat <- leip_MP_scores[rownames(leip_MP_scores) %in% hpv_neg_leip_data$`Sample name`, ] %>% 
  purrr::map(\(x) ifelse(x >= quantile(x, probs = 0.75), "High", 
                         ifelse(x <= quantile(x, probs = 0.25), "Low", NA))) %>% 
  as.data.frame %>% `rownames<-`(rownames(leip_MP_scores[rownames(leip_MP_scores) %in% hpv_neg_leip_data$`Sample name`, ]))
fhcrc_MP_cat <- fhcrc_MP_scores %>% 
  purrr::map(\(x) ifelse(x >= quantile(x, probs = 0.75), "High", 
                         ifelse(x <= quantile(x, probs = 0.25), "Low", NA))) %>% 
  as.data.frame %>% `rownames<-`(rownames(fhcrc_MP_scores))

merged_hpv_neg_cat_df <- rbind.data.frame(tcga_hpv_neg_MP_cat, leip_hpv_neg_MP_cat, fhcrc_MP_cat) %>% dplyr::rename(INF = Inflammatory)
merged_hpv_neg_clinical_data <- rbind(filt_tcga_clinical_data[rownames(filt_tcga_clinical_data) %in% hpv_neg_tcga_data$bcr_patient_barcode, ], 
                                      filt_leip_clinical_data[rownames(filt_leip_clinical_data) %in% hpv_neg_leip_data$`Sample name`, ], filt_fhcrc_clinical_data)
merged_hpv_neg_clinical_data <- merged_hpv_neg_clinical_data %>% dplyr::mutate(OS.time.days = OS.time, OS.time.months = OS.time / (365 / 12))




# Analyze Survival Associations Limited to 5-year Overall Survival -------------------------------------------

surv_df <- merge(merged_hpv_neg_cat_df, merged_hpv_neg_clinical_data, by = 0) %>% 
  dplyr::mutate(OS_5y = ifelse(OS.time.months >= 60, 0, OS),
                OS_5y_time = ifelse(OS.time.months >= 60, 60, OS.time.months))
surv_df <- surv_df %>% dplyr::select(-c(LowQ, Skeletal_Muscle))

formulae <- lapply(colnames(surv_df[, which(colnames(surv_df) == "Fibroblast"):which(colnames(surv_df) == "Secretory_Norm")]), function(x) as.formula(paste0("Surv(OS_5y_time, OS_5y) ~ ", x)))
fit_mod <- survminer::surv_fit(formulae, data = surv_df)
names(fit_mod) <- stringr::str_replace_all(names(fit_mod), pattern = "surv_df::", replacement = "")
p <- survminer::ggsurvplot_list(fit = fit_mod,
                                data = surv_df,
                                risk.table = FALSE,
                                pval = TRUE,
                                legend.title = "",
                                legend.position = "none",
                                break.time.by = 12,
                                ggtheme = theme_minimal(),
                                # surv.median.line = "hv",
                                surv.median.line = "none",
                                risk.table.y.text.col = FALSE,
                                risk.table.y.text = FALSE,
                                palette = c("#E64B35FF", "#3C5488FF"))
p <- p[c("Secretory_Norm", "T_cell", "Hypoxia", "Fibroblast", "B_cell", "pEMT", "Endothelial",
         "Cell_Cycle", "Secretory_Malig", "Complement", "Epithelial", "Senescence", "Macrophage", "INF")]
dev.off()
surv_plots <- arrange_ggsurvplots(p[1:length(p)], print = TRUE, ncol = 5, nrow = 3, title = "Survival plots")
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_2d.pdf"), device = "pdf", plot = surv_plots, dpi = 300, width = 12, height = 10)

