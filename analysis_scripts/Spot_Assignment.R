library(tidyverse)
library(scalop)
library(Seurat)
library(here)
library(scales)
library(patchwork)
library(egg)
source(here("analysis_scripts/functions/scRNA_Functions.R"))
source(here("analysis_scripts/functions/ST_Functions.R"))


# Load per-sample parameter table
sample_parameters <- readRDS(file = here("metadata/per_sample_parameters_df.rds"))
# Load sample level metadata
samples_metadata <- readRDS(file = here("metadata/samples_metadata.rds"))
# Load metaprograms
MPs <- readRDS(file = here("results/Generated_Data/Final_Metaprograms_Extended.rds"))
# Load Epithelial vs Stromal signatures
epi_stroma_progs <- readRDS(file = here("results/Generated_Data/epi_nonEpi_gene_programs.rds"))
# Load state colors
state_cols <- read.delim(file = here("aux_data/state_col_tab_extend.tsv"), sep = "\t", header = TRUE)
state_cols <- setNames(state_cols$V2, state_cols$V1)



# Create merged matrix with the corrected aliase names shared genes -------
expr_mat_ls <- list()
centered_expr_mat_ls <- list()
key_tab <- readRDS(file =  here("aux_data/gene_aliases_correction_table.rds"))
for(i in seq_along(sample_parameters$Sample)) {
  message(paste0("Processing sample: ", sample_parameters$Sample[[i]]))

  # Load expression matrix
  data_dir <- paste(here("data"), sample_parameters$Site[[i]], sample_parameters$Sample[[i]], sep = "/")
  m <- as.matrix(Matrix::readMM(paste(data_dir, grep("matrix(_[A-Za-z0-9]+)?\\.mtx$", list.files(data_dir), value = TRUE), sep = "/")))
  colnames(m) <- read.table(paste(data_dir, grep("barcodes(_[A-Za-z0-9]+)?\\.tsv$", list.files(data_dir), value = TRUE), sep = "/"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)[, 1]
  colnames(m) <- paste(sample_parameters$Sample[[i]], colnames(m), sep = "_")
  rownames(m) <- read.table(paste(data_dir, grep("features(_[A-Za-z0-9]+)?\\.tsv$", list.files(data_dir), value = TRUE), sep = "/"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)[, 2]

  if(sum(grep("GRCh38", rownames(m))) == 0) {
    gene_alias_idx <- match(rownames(m), key_tab$V2.x)
    rownames(m)[!is.na(gene_alias_idx)] <- key_tab$Symbol[na.omit(gene_alias_idx)]
  } else {
    # Remove HPV genes and GRC38 prefix
    m <- m[-grep("^hpv", rownames(m)), ]
    rownames(m) <- sub("GRCh38_", "", rownames(m))
    gene_alias_idx <- match(rownames(m), key_tab$V2.y)
    rownames(m)[!is.na(gene_alias_idx)] <- key_tab$Symbol[na.omit(gene_alias_idx)]
  }

  # Process raw data (no centering)
  proc_mat <- filter_10x(m, complexity_thresh = sample_parameters$Complexity_cut[[i]], percent_mt_thresh = sample_parameters$MT_cut[[i]],
                         centre = "none", gene_filt_method = "none", merge_method = "sum")

  # Output per-sample centered expression matrix to a list
  expr_mat_ls[[i]] <- proc_mat
  names(expr_mat_ls)[i] <- paste(sample_parameters$Site[[i]], sample_parameters$Sample[[i]], sep = "_")
}

shared_features <- readRDS(file =  here("aux_data/shared_features.rds"))
expr_mat_ls <- lapply(expr_mat_ls, function(x) x[shared_features, ])
merged_expr_mat <- do.call(cbind, expr_mat_ls)
# saveRDS(merged_expr_mat, file = here("results/Generated_Data/merged_expr_mat.rds"))



# Per-Spot Scoring And Scoring-based Deconvolution -----------------------------------------------------

if (!dir.exists(file.path(here("results"), "Deconvolution_mtrxs"))) {
  dir.create(file.path(here("results"), "Deconvolution_mtrxs"))
}

for(i in seq_along(sample_parameters$Sample)) {
  message(paste0("Processing sample: ", sample_parameters$Sample[[i]]))

  # Load expression matrix
  data_dir <- paste(here("data"), sample_parameters$Site[[i]], sample_parameters$Sample[[i]], sep = "/")
  m <- as.matrix(Matrix::readMM(paste(data_dir, grep("matrix(_[A-Za-z0-9]+)?\\.mtx$", list.files(data_dir), value = TRUE), sep = "/")))
  colnames(m) <- read.table(paste(data_dir, grep("barcodes(_[A-Za-z0-9]+)?\\.tsv$", list.files(data_dir), value = TRUE), sep = "/"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)[, 1]
  rownames(m) <- read.table(paste(data_dir, grep("features(_[A-Za-z0-9]+)?\\.tsv$", list.files(data_dir), value = TRUE), sep = "/"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)[, 2]

  if(sum(grep("GRCh38", rownames(m))) == 0) {
    m <- m
  } else {
    # Remove HPV genes and GRC38 prefix
    m <- m[-grep("^hpv", rownames(m)), ]
    rownames(m) <- sub("GRCh38_", "", rownames(m))
  }
  # Load spatial image object
  spatial_image <- load_spatial_image(paste0(data_dir, "/spatial"))

  # Process raw data (no centering)
  proc_mat <- filter_10x(m, complexity_thresh = sample_parameters$Complexity_cut[[i]], percent_mt_thresh = sample_parameters$MT_cut[[i]],
                         centre = "none", gene_filt_method = "exp", merge_method = "sum")

  # Score for metaprograms and plot per-sample spatial heatmap
  spot_MP_score <- metaprog_score(proc_mat, MPs, score_diff = 0, score_cut = 0, bin_size = 100, bin_num = 30, conserved.genes = 0.5, center_rows = TRUE, scale = FALSE)
  plot_scores <- spatial_heatmap(spot_MP_score$mp_score_list, spatial_image, normalize = TRUE, size = 0.5)
  all_heatmaps <- do.call(ggarrange, plot_scores)
  ggplot2::ggsave(filename = paste0(here("results/Per_Sample_Results/"), sample_parameters$Sample[[i]], "/", sample_parameters$Sample[[i]], "_PerSpot_Score_Plot_Extend.png"),
                  plot = all_heatmaps, width = 30, height = 25, units = "cm")

  # Compute scoring-based deconvolution and plot results in spatial scatterpie per sample
  score_vecs <- t(do.call(cbind, spot_MP_score$mp_score_list))
  score_props <- apply(score_vecs, 2, function(x) {
    if(length(x[x > 0]) == 0) {prop_vec <- setNames(rep(0, nrow(score_vecs)), rownames(score_vecs))}
    get_positives <- x[x > 0]
    normalize_sum2one <- get_positives / sum(get_positives)
    prop_vec <- setNames(rep(0, nrow(score_vecs)), rownames(score_vecs))
    prop_vec <- setNames(ifelse(rownames(score_vecs) %in% names(normalize_sum2one), yes = normalize_sum2one[match(names(prop_vec), names(normalize_sum2one))], no = 0), rownames(score_vecs))
    return(prop_vec)
  })
  # Filtering out minimum contributions and updated proportions after
  score_props[score_props < 0.1] <- 0
  filt_score_props <- apply(score_props, 2, function(x) x / sum(x))
  filt_score_props[is.na(filt_score_props)] <- 0

  # Save deconvolution matrix as RDS file
  decon_mtrx <- t(filt_score_props)
  saveRDS(decon_mtrx, file = paste0(here("results/Deconvolution_mtrxs/"), sample_parameters$Sample[[i]], "_score_decon_mtrx_Extend.rds"))

  # Plot and save spatial scatterpie
  plot_scatterpie <- spatial_scatterpie(spatial_image, decon_mtrx, img = FALSE, scatterpie_alpha = 1, pie_scale = 0.4) +
    scale_fill_manual(values = state_cols[names(state_cols) %in% colnames(decon_mtrx)], breaks = colnames(decon_mtrx), name = "Meta-Programs")
  ggplot2::ggsave(filename = paste0(here("results/Per_Sample_Results/"), sample_parameters$Sample[[i]], "/", sample_parameters$Sample[[i]], "_Scoring_Based_Deconv_ScatterPie_Extend.png"),
                  plot = plot_scatterpie, width = 35, height = 35, units = "cm")
}

# Load scoring-based deconvolution matrices and create a list uniting all these
path <- here("results/Deconvolution_mtrxs")
matrices <- list.files(path)[grep("_Extend.rds$", list.files(path))]
all_states <- names(MPs)
decon_mats_ls <- list()
for(i in seq_along(matrices)) {
  decon_mat <- readRDS(paste(path, matrices[i], sep = "/"))
  if(sum(all_states %in% colnames(decon_mat)) == length(all_states)) {decon_mat <- decon_mat}
  else {
    missing_states <- all_states[!all_states %in% colnames(decon_mat)]
    missing_mat <- matrix(data = 0, nrow = nrow(decon_mat), ncol = length(missing_states), dimnames = list(rownames(decon_mat), missing_states))
    decon_mat <- cbind(decon_mat, missing_mat)
  }

  # Output to a list
  decon_mats_ls[[i]] <- decon_mat[, all_states]
  names(decon_mats_ls)[i] <- substri(matrices[i], pos = 1)
}
# saveRDS(decon_mats_ls, file = here("results/Generated_Data/Decon_Mats_List_Extend.rds"))





# Create a Merged Metadata and Compute MP Scores per Sample --------------------------------------------

# Create merged matrix
save_score_vecs <- list()
merged_metadata <- tibble()
for(i in seq_along(sample_parameters$Sample)) {
  message(paste0("Processing sample: ", sample_parameters$Sample[[i]]))

  # Load expression matrix
  data_dir <- paste(here("data"), sample_parameters$Site[[i]], sample_parameters$Sample[[i]], sep = "/")
  m <- as.matrix(Matrix::readMM(paste(data_dir, grep("matrix(_[A-Za-z0-9]+)?\\.mtx$", list.files(data_dir), value = TRUE), sep = "/")))
  colnames(m) <- read.table(paste(data_dir, grep("barcodes(_[A-Za-z0-9]+)?\\.tsv$", list.files(data_dir), value = TRUE), sep = "/"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)[, 1]
  colnames(m) <- paste(sample_parameters$Sample[[i]], colnames(m), sep = "_")
  rownames(m) <- read.table(paste(data_dir, grep("features(_[A-Za-z0-9]+)?\\.tsv$", list.files(data_dir), value = TRUE), sep = "/"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)[, 2]

  if(sum(grep("GRCh38", rownames(m))) == 0) {
    m <- m
  } else {
    # Remove HPV genes and GRC38 prefix
    m <- m[-grep("^hpv", rownames(m)), ]
    rownames(m) <- sub("GRCh38_", "", rownames(m))
  }

  # Process raw data (no centering)
  proc_mat <- filter_10x(m, complexity_thresh = sample_parameters$Complexity_cut[[i]], percent_mt_thresh = sample_parameters$MT_cut[[i]],
                         centre = "none", gene_filt_method = "exp", merge_method = "sum")

  # Score for Epithelial vs Stromal signatures and Metaprogram identity
  epi_stroma_score <- metaprog_score(proc_mat, epi_stroma_progs, score_diff = 0.4, score_cut = 0, bin_size = 100, bin_num = 30, conserved.genes = 0.5, center_rows = TRUE, scale = FALSE)
  spot_MP_score <- metaprog_score(proc_mat, MPs, score_diff = 0, score_cut = 0, bin_size = 100, bin_num = 30, conserved.genes = 0.5, center_rows = TRUE, scale = FALSE)
  save_score_vecs[[i]] <- spot_MP_score$mp_score_list
  names(save_score_vecs)[i] <- sample_parameters$Sample[[i]]

  # Create a metadata tibble to store all collected variables
  cell_stats <- cell_qc(m) %>% dplyr::filter(CellID %in% colnames(proc_mat))
  metadata <- tibble(SpotID = unname(sapply(colnames(proc_mat), function(x) str_split(x, pattern = "_")[[1]][[2]])),
                     Sample = unname(sapply(colnames(proc_mat), function(x) str_split(x, pattern = "_")[[1]][[1]])),
                     Key = colnames(proc_mat),
                     Complexity = unname(cell_stats$Complexity),
                     Percent_Mito = unname(cell_stats$Percent_mt),
                     Percent_HK = unname(cell_stats$HK),
                     MPid = unname(spot_MP_score$meta_program[match(names(spot_MP_score$meta_program), colnames(proc_mat))]),
                     EpiStroma = unname(epi_stroma_score$meta_program[match(names(epi_stroma_score$meta_program), colnames(proc_mat))]))

  # Join metadata to merged metadata tibble
  merged_metadata <- rbind(merged_metadata, metadata)
}
# saveRDS(save_score_vecs, file = here("results/Generated_Data/per_sample_state_score_vectors_Extend.rds"))

# Plot QC parameters for the merged expression data
ggpubr::ggviolin(merged_metadata$Complexity, fill = "lightskyblue", xlab = "Complexity", ylab = FALSE) + ggpubr::ggviolin(merged_metadata$Percent_Mito, fill = "lightskyblue", xlab = "Percent_Mito", ylab = FALSE) +
  ggpubr::ggviolin(merged_metadata$Percent_HK, fill = "lightskyblue", xlab = "Percent_HouseKeeping", ylab = FALSE)

# Save merged metadata
merged_metadata <- merged_metadata %>% tibble::add_column(Percent_Epi = NA, Percent_Stroma = NA, SpotType = NA, binCNAstatus = NA, CNAscore = NA, CNAcorr = NA, CNAsignal = NA, Subclone = NA)
# saveRDS(merged_metadata, file = here("metadata/merged_metadata_Extend.rds"))





# Compute Epithelial Percentage per Spot (based on the deconvolution proportions) --------

# Load merged metadata tibble
metadata <- readRDS(file = here("metadata/merged_metadata_Extend.rds"))
# Load scoring-based deconvolution matrices
decon_mats_ls <- readRDS(file = here("results/Generated_Data/Decon_Mats_List_Extend.rds"))

# Calculate percentage of Epithelial states per spots based on the cell-type proportions of the deconvolution matrix
Epi_states <- c("Senescence", "Epithelial", "Secretory_Malig", "Hypoxia", "pEMT")
Stroma_states <- c("T_cell", "B_cell", "Macrophage", "Fibroblast", "Endothelial", "Complement", "Skeletal_Muscle", "Secretory_Norm")
for(i in seq_along(decon_mats_ls)) {
  spot_epi_prop <- apply(decon_mats_ls[[i]], 1, function(x) sum(x[Epi_states]) / sum(x[c(Epi_states, Stroma_states)]))
  spot_epi_prop[is.na(spot_epi_prop)] <- 0
  metadata$Percent_Epi[metadata$Sample == names(decon_mats_ls[i])] <- spot_epi_prop[match(names(spot_epi_prop), metadata$SpotID[metadata$Sample == names(decon_mats_ls[i])]) &
                                                                                      names(spot_epi_prop) %in% metadata$SpotID[metadata$Sample == names(decon_mats_ls[i])]]
  # In order to avoid calling spots with no deconvolution proportions stromal spots (due to assigning 0 to NA values) repeat the process for stromal proportions (and not spimply 1 - epi_prop)
  spot_stroma_prop <- apply(decon_mats_ls[[i]], 1, function(x) sum(x[Stroma_states]) / sum(x[c(Epi_states, Stroma_states)]))
  spot_stroma_prop[is.na(spot_stroma_prop)] <- 0
  metadata$Percent_Stroma[metadata$Sample == names(decon_mats_ls[i])] <- spot_stroma_prop[match(names(spot_stroma_prop), metadata$SpotID[metadata$Sample == names(decon_mats_ls[i])]) &
                                                                                            names(spot_stroma_prop) %in% metadata$SpotID[metadata$Sample == names(decon_mats_ls[i])]]
}

### OPTIONALLY - add deconvolution proportions to the metadata
# Order decon_mats list to match the order of the samples in the metadata
decon_mats_ls <- decon_mats_ls[unique(metadata$Sample)]
# Add sample name to rownames of each decon_mat
decon_mats_ls <- lapply(names(decon_mats_ls), function(samp) {
  rownames(decon_mats_ls[[samp]]) <- paste0(samp, "_", rownames(decon_mats_ls[[samp]]))
  decon_mats_ls[[samp]]
})
decon_props <- do.call(rbind.data.frame, decon_mats_ls) %>% tibble::rownames_to_column(var = "Key")
metadata <- as_tibble(merge(metadata, decon_props, by = "Key", sort = FALSE))
# saveRDS(metadata, file = here("metadata/merged_metadata_Extend.rds"))





# Gene Expression Heatmap -------------------------------------------------

# Load merged matrix
merged_expr_mat <- readRDS(file = here("results/Generated_Data/merged_expr_mat.rds"))
# Load merged metadata tibble
metadata <- readRDS(file = here("metadata/merged_metadata_Extend.rds"))

# Rearrange metaprograms for expression heatmap plotting
MPs <- MPs[c("Skeletal_Muscle", "Secretory_Norm", "Fibroblast", "Endothelial", "Complement", "Macrophage", "T_cell", "B_cell", "Cell_Cycle", "Epithelial", "Inflammatory", "Hypoxia", "pEMT", "LowQ", "Secretory_Malig", "Senescence")]
metadata$MPid <- factor(metadata$MPid, levels = c(names(MPs)))

# Assign spot to Epithelial, Stromal or Mixed based on the agreement between the scoring of the merged matrix and the per-sample deconvolution proportions
metadata <- metadata %>% mutate(EpiStroma = ifelse(Percent_Epi == 1, "Epithelial",
                                                   ifelse(Percent_Stroma == 1, "Stroma",
                                                          ifelse(Percent_Epi == 0 & Percent_Stroma == 0, "Filtered_Out",
                                                                 ifelse(EpiStroma == "Epithelial" & Percent_Epi >= 0.85, "Epithelial",
                                                                        ifelse(EpiStroma == "Non_Epithelial" & Percent_Stroma >= 0.85, "Stroma", "Mixed"))))))

# Subset the merged centered expression matrix to contain only the genes that appear in the metaprograms
MP_genes <- unique(unname(unlist(MPs)))
filt_merged_expr_mat <-  filter_10x(merged_expr_mat[rownames(merged_expr_mat) %in% MP_genes, ], centre = "mean", log = FALSE, raw_to_CPM = FALSE,
                                    complexity_thresh = NULL, percent_mt_thresh = NULL, gene_filt_method = "none", merge_method = "sum", dbl_remove = FALSE)

# Plot a separate heatmap for each of the groups (Epithelial, Stromal, Mixed)
hm_plots <- map(unique(metadata$EpiStroma[!metadata$EpiStroma %in% "Filtered_Out"]), ~ plot_heatmap(filt_merged_expr_mat, metadata, group.by = "MPid", features = Unlist(MPs), cells = metadata$Key[metadata$EpiStroma == .], disp.max = 4, disp.min = -4, label = FALSE, group.colors = state_cols, sep_line_col = "#000000", anno_leg_name = "Metaprograms"))
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_2b_mixed.pdf"), plot = hm_plots[[1]], device = "pdf", dpi = 300)
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_2b_epithelial.pdf"), plot = hm_plots[[2]], device = "pdf", dpi = 300)
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_2b_non_epithelial.pdf"), plot = hm_plots[[3]], device = "pdf", dpi = 300)

# saveRDS(metadata, file = here("metadata/merged_metadata_Extend.rds"))




# Cell-type Scores Heatmap ------------------------------------------------

# Load spots MP scores
score_vecs <- readRDS(here("results/Generated_Data/per_sample_state_score_vectors_Extend.rds"))
order_MPs <- c("Skeletal_Muscle", "Secretory_Norm", "Fibroblast", "Endothelial", "Complement", "Macrophage", "T_cell", "B_cell",
               "Cell_Cycle", "Epithelial", "Inflammatory", "Hypoxia", "pEMT", "LowQ", "Secretory_Malig", "Senescence")
score_vecs <- lapply(score_vecs, function(x) x[order_MPs])

# Merge score list to matrix that represent all states over all spots
make_equal_lengths <- lapply(score_vecs, function(x) {
  if(sum(is.na(names(x))) != 0) {
    spot_names <- names(x[!is.na(names(x))][[1]])
    na_vec <- setNames(rep(NA, length(spot_names)), spot_names)
    x <- lapply(x, function(i) {
      if (is.null(i)) {
        na_vec
      } else {
        i
      }
    })
  } else {
    x <- x
  }
})
state_score_mat <- sapply(seq_along(names(MPs)), function(i) do.call(c, lapply(make_equal_lengths, function(x) x[[i]]))) %>%
  `row.names<-`(scalop::substri(row.names(.), sep = "\\.", pos = 2)) %>% magrittr::set_colnames(names(MPs))

order_spot_scores_df <- as.data.frame(state_score_mat) %>%
  rowwise() %>%
  dplyr::mutate(MP = names(.)[which.max(c_across(everything()))], max_val = max(c_across(Fibroblast:Secretory_Norm), na.rm = TRUE)) %>%
  `rownames<-`(rownames(state_score_mat)) %>% ungroup

plot_spot_scores_df <- order_spot_scores_df %>%
  tibble::rownames_to_column(var = "Spot") %>%
  dplyr::mutate(MP = factor(MP, levels = order_MPs)) %>%
  dplyr::group_by(MP) %>% arrange(desc(max_val), .by_group = TRUE) %>% dplyr::ungroup() %>%
  dplyr::select(-c("MP", "max_val")) %>%
  tidyr::pivot_longer(!Spot, names_to = "State", values_to = "Score") %>%
  dplyr::mutate(Spot = factor(Spot, levels = unique(Spot)), State = factor(State, levels = order_MPs))
# saveRDS(plot_spot_scores_df, file = here("results/Generated_Data/All_Spots_MP_Scores_for_Heatmap.rds"))

plot_spot_scores_df <- readRDS(here("results/Generated_Data/All_Spots_MP_Scores_for_Heatmap.rds"))
p <- gmap(plot_spot_scores_df, x = Spot, y = State, fill = Score, limits = c(-2, 2), ratio = 0.6, y.name = "", midpoint = 0, geom = "raster",
          axis.rel = 1.5, legend.title.rel = 1.25, legend.rel = 1.0, legend.height = 5, legend.width = 0.6, x.num = F, y.num = F)
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_2a.pdf"), plot = p, dpi = 300, device = "pdf")





# Per-Sample MP proportions -----------------------------------------------

samples_metadata$Ext_Site <- ifelse(samples_metadata$Site == "Oropharynx" & samples_metadata$HPV == "HPVpos", "HPV+ OP",
                                    ifelse(samples_metadata$Site == "Oropharynx" & samples_metadata$HPV == "HPVneg", "HPV- OP", samples_metadata$Site))
per_sample_states <- as.data.frame(table(metadata$MPid[!is.na(metadata$MPid)], metadata$Sample[!is.na(metadata$MPid)])) %>%
  rowwise %>% mutate(Site = as.factor(samples_metadata$Ext_Site[samples_metadata$Sample == Var2]), .before = "Freq") %>% ungroup %>%
  group_by(Var2) %>% mutate(Prop = Freq / sum(Freq), , OrderBy = Prop[Var1 == "Senescence"]) %>% ungroup
p <- per_sample_states %>%
  dplyr::mutate(Group = ifelse(Var1 %in% c("T_cell", "B_cell", "Macrophage", "Fibroblast", "Endothelial", "Complement", "Skeletal_Muscle", "Secretory_Norm"), "Non-Malignant", "Malignant")) %>%
  dplyr::group_by(Var1, Var2) %>%
  summarize(Malignant = sum(Prop[Group == "Malignant"]), `Non-Malignant` = -sum(Prop[Group == "Non-Malignant"])) %>%
  dplyr::group_by(Var2) %>%
  dplyr::mutate(Ord = sum(`Non-Malignant`)) %>% ungroup %>%
  dplyr::mutate(Site = factor(samples_metadata$Ext_Site[match(Var2, samples_metadata$Sample)], levels = c("Oral", "Laryngeal", "HPV+ OP", "HPV- OP"))) %>%
  ggplot(aes(x = fct_reorder(Var2, Ord), y = Malignant, fill = Var1)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  geom_col(aes(y = `Non-Malignant`), position = "stack") +
  geom_hline(aes(yintercept = 0)) +
  scale_fill_manual(name = "States", values = state_cols[levels(per_sample_states$Var1)]) + scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0.01)) +
  theme_classic() + theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.9, size = 16), plot.title = element_text(hjust = 0.5, size = 16),
                          legend.text = element_text(size = 14), legend.title = element_text(size = 16)) +
  facet_grid(.~Site, scales = 'free', space = 'free', drop = F) + guides(fill = guide_legend(override.aes = list(size = 7), reverse = TRUE))
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_2c.pdf"), plot = p, device = "pdf", dpi = 300, width = 10, height = 6)

# Summary of state proportion per site
per_site_states <- metadata %>% dplyr::mutate(Site = as.factor(samples_metadata$Site[match(Sample, samples_metadata$Sample)])) %>% dplyr::select(MPid, Site) %>% table() %>% as.data.frame() %>%
  group_by(Site) %>% mutate(Prop = Freq / sum(Freq)) %>% ungroup
p <- per_site_states %>%
  dplyr::mutate(Group = ifelse(MPid %in% c("T_cell", "B_cell", "Macrophage", "Fibroblast", "Endothelial", "Complement", "Skeletal_Muscle", "Secretory_Norm"), "Non-Malignant", "Malignant")) %>%
  dplyr::group_by(MPid, Site) %>%
  summarize(Malignant = sum(Prop[Group == "Malignant"]), `Non-Malignant` = -sum(Prop[Group == "Non-Malignant"])) %>%
  dplyr::group_by(Site) %>%
  dplyr::mutate(Ord = sum(`Non-Malignant`)) %>%
  ggplot(aes(x = fct_reorder(Site, Ord), y = Malignant, fill = MPid)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  geom_col(aes(y = `Non-Malignant`), position = "stack") +
  geom_hline(aes(yintercept = 0)) +
  scale_fill_manual(name = "States", values = state_cols[levels(per_site_states$MPid)]) + scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  theme_classic() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 16), plot.title = element_text(hjust = 0.5, size = 16),
                          legend.text = element_text(size = 14), legend.title = element_text(size = 16)) +
  guides(fill = guide_legend(override.aes = list(size = 7), reverse = TRUE))
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_2f.pdf"), plot = p, device = "pdf", dpi = 300, width = 7, height = 6)


# Draw statistics
states_per_site_stats <- metadata %>%
  dplyr::mutate(Site = as.factor(samples_metadata$Site[match(Sample, samples_metadata$Sample)]))
prop_stat_tab <- test_state_prop_diff(states_per_site_stats$MPid, states_per_site_stats$Sample, states_per_site_stats$Site, transform = NULL)
write.csv(prop_stat_tab, file = here("results/Generated_Data/State_Proportion_Per_Site_ANOVA_Results_Extend.csv"))




# CNA Inference & Subclones Division --------------------------------

# Load merged metadata tibble
metadata <- readRDS(file = here("metadata/merged_metadata_Extend.rds"))
# Load scoring-based deconvolution matrices
decon_mats_ls <- readRDS(file = here("results/Generated_Data/Decon_Mats_List_Extend.rds"))
sample_parameters <- read.csv(here("metadata/per_sample_parameters.csv")) %>% dplyr::select(-X)

# Iterate through all samples to infere CNAs in each
for(i in seq_along(sample_parameters$Sample)) {
  message(paste0("Processing sample: ", sample_parameters$Sample[[i]]))

  # Load expression matrix
  data_dir <- paste(here("data"), sample_parameters$Site[[i]], sample_parameters$Sample[[i]], sep = "/")
  m <- as.matrix(Matrix::readMM(paste(data_dir, grep("matrix(_[A-Za-z0-9]+)?\\.mtx$", list.files(data_dir), value = TRUE), sep = "/")))
  colnames(m) <- read.table(paste(data_dir, grep("barcodes(_[A-Za-z0-9]+)?\\.tsv$", list.files(data_dir), value = TRUE), sep = "/"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)[, 1]
  rownames(m) <- read.table(paste(data_dir, grep("features(_[A-Za-z0-9]+)?\\.tsv$", list.files(data_dir), value = TRUE), sep = "/"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)[, 2]

  if(sum(grep("GRCh38", rownames(m))) == 0) {
    m <- m
  } else {
    # Remove HPV genes and GRC38 prefix
    m <- m[-grep("^hpv", rownames(m)), ]
    rownames(m) <- sub("GRCh38_", "", rownames(m))
  }
  # Load spatial image object
  spatial_image <- load_spatial_image(paste0(data_dir, "/spatial"))
  # Extract sample specific metadata
  metadata_spots <- metadata[metadata$Sample == sample_parameters$Sample[[i]], ]

  # Process raw data (no centering)
  proc_mat <- filter_10x(m, complexity_thresh = sample_parameters$Complexity_cut[[i]], percent_mt_thresh = sample_parameters$MT_cut[[i]],
                         centre = "none", gene_filt_method = "exp", merge_method = "sum")

  # Subset expression matrix to include only spots present in the merged metadata
  proc_mat <- proc_mat[, metadata_spots$SpotID]

  # Extract Reference vs Epithelial signatures for reference spot selection
  if(sum(decon_mats_ls[[sample_parameters$Sample[[i]]]][, "Skeletal_Muscle"]) != 0) {
    ref_vs_epi_prog <- list(Reference = Unlist(MPs[c("Fibroblast", "Endothelial", "Complement", "Skeletal_Muscle")]),
                            Epithelial = Unlist(MPs[c("Epithelial", "Senescence")]))
  } else {
    ref_vs_epi_prog <- list(Reference = Unlist(MPs[c("Fibroblast", "Endothelial", "Complement")]),
                            Epithelial = Unlist(MPs[c("Epithelial", "Senescence")]))
  }
  epi_ref_scores <- metaprog_score(proc_mat, ref_vs_epi_prog, score_diff = 0.4, score_cut = 0, bin_size = 100, bin_num = 30, conserved.genes = 0.5, center_rows = TRUE, scale = FALSE)
  epi_ref_scores <- do.call(cbind.data.frame, epi_ref_scores$mp_score_list) %>% tibble::rownames_to_column(var = "SpotID")

  # Assign reference and query for CNA inference
  ref <- epi_ref_scores$SpotID[which(epi_ref_scores$Epithelial <= 0 & epi_ref_scores$Reference > quantile(epi_ref_scores$Reference, probs = 0.95))]
  table(metadata_spots$EpiStroma[metadata_spots$SpotID %in% ref])
  ref <- ref[ref %in% metadata_spots$SpotID[metadata_spots$EpiStroma == "Stroma" & metadata_spots$Percent_Stroma == 1]]
  query <- colnames(proc_mat)[!colnames(proc_mat) %in% metadata_spots$SpotID[metadata_spots$EpiStroma == "Stroma" & metadata_spots$Percent_Stroma == 1]]


  # Process non-centered matrix including only the reference and query spots
  epi_ref_m <- filter_10x(m, cells = c(query, ref), centre = "none", merge_method = "sum", gene_filt_method = "exp",
                          complexity_thresh = sample_parameters$Complexity_cut[[i]], percent_mt_thresh = sample_parameters$MT_cut[[i]])

  # Load genome and compute a CNA matrix
  hg38 <- readRDS(here("aux_data/hg38.rds"))
  cna_score_mat <- calc_cna(matrix = epi_ref_m, query = query, ref = ref, genome = hg38, range = c(-3,3), window = 100, noise = 0.15, isLog = TRUE, per_chr = TRUE, scale = NULL, top_genes_num = NULL, verbose = TRUE)
  # Compute and plot CNA-correlation vs CNA-signal graph
  sig_and_cor <- cna_sig_cor(cna_score_mat, epi_cells = query, ref_cells = ref)
  # metadata_spots$SpotType <- ifelse(metadata_spots$SpotID %in% ref, "Reference", metadata_spots$EpiStroma)
  metadata_spots$SpotType <- ifelse(metadata_spots$EpiStroma == "Epithelial" & metadata_spots$Percent_Epi >= 0.85, yes = "Epithelial",
                                    ifelse(metadata_spots$EpiStroma == "Stroma" & metadata_spots$Percent_Stroma == 1 & !(metadata_spots$SpotID %in% ref), yes = "Stroma",
                                           ifelse(metadata_spots$SpotID %in% ref, yes = "Reference", no = "Mixed")))
  metadata_spots$CNAsignal <- sig_and_cor$CNA_Signal[match(metadata_spots$SpotID, names(sig_and_cor$CNA_Signal))]
  metadata_spots$CNAcorr <- sig_and_cor$CNA_Correlation[match(metadata_spots$SpotID, names(sig_and_cor$CNA_Correlation))]
  metadata_spots$CNAscore <- metadata_spots$CNAcorr + (metadata_spots$CNAsignal * max(metadata_spots$CNAcorr[!is.na(metadata_spots$CNAcorr)]) / max(metadata_spots$CNAsignal[!is.na(metadata_spots$CNAsignal)]))
  plot_SigCor <- ggpubr::ggscatter(metadata_spots, x = "CNAcorr", y = "CNAsignal", color = "SpotType")
  ggplot2::ggsave(filename = paste0(here("results/Per_Sample_Results/"), sample_parameters$Sample[[i]], "/", sample_parameters$Sample[[i]], "_CNA_Sig_vs_Cor.png"),
                  plot = plot_SigCor, width = 35, height = 25, units = "cm")
  # Plot CNA matrix ordered by correlation
  order_cells <- metadata_spots %>% dplyr::select(SpotID, CNAcorr) %>% na.omit() %>% .[order(metadata_spots$CNAcorr[!is.na(metadata_spots$CNAcorr)]), ] %>% pull(SpotID)
  ordered_cna_mat <- cna_score_mat[, order_cells]
  cna_plot <- infercna::ggcna(ordered_cna_mat, genome = "hg38", title = "", legend.title = "Inferred CNA\n[log2 ratio]\n", legend.height = 0.6) +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "right") + labs(y = "Spot", x = "Chromosome")
  ggplot2::ggsave(filename = paste0(here("results/Per_Sample_Results/"), sample_parameters$Sample[[i]], "/", sample_parameters$Sample[[i]], "_CNA_Matrix_Plot.png"),
                  plot = cna_plot, width = 35, height = 25, units = "cm")

  # Plot spatial heatmaps for all CNA statistics
  cna_stats <- as.list(metadata_spots[, c("SpotID", "CNAsignal", "CNAcorr", "CNAscore", "SpotType")])
  cna_spatial_hms <- spatial_heatmap(cna_stats, spatial_image, normalize = FALSE, size = 1.8)
  plot_spatial_hm <- egg::ggarrange(cna_spatial_hms$SpotType, cna_spatial_hms$CNAsignal, cna_spatial_hms$CNAcorr, cna_spatial_hms$CNAscore, ncol = 2, nrow = 2)
  ggplot2::ggsave(filename = paste0(here("results/Per_Sample_Results/"), sample_parameters$Sample[[i]], "/", sample_parameters$Sample[[i]], "_CNA_PerSpot_Scoring_Spotmap.png"),
                  plot = plot_spatial_hm, width = 32, height = 22, units = "cm")

  metadata[metadata$Sample == sample_parameters$Sample[[i]], ] <- metadata_spots

  ## CNA Subclone devision
  if(sample_parameters$Subclone_function[[i]] == "spatial_subclones") {
    subclones_clusters <- spatial_subclones(cna_score_mat, query, separate = "arm", genecut = 10, genome = hg38, top_region = 2/3, top_method = "abs", reduction_dims = 30, cluster_k = sample_parameters$Subclones_clust_k[[i]])
  }
  if(sample_parameters$Subclone_function[[i]] == "ST_subclones") {
    subclones_clusters <- ST_subclones(cna_score_mat, query, separate = "arm", genecut = 10, genome = hg38, top_cna_region = 2/3, top_arm_region = 1, amp_del_genes = "all", reduction_dims = 30, cluster_k = sample_parameters$Subclones_clust_k[[i]], resolution = sample_parameters$Subclones_res[[i]],
                                       event_conf_val = 3.5, var_arm_len = 10, expr_diff_cutoff = 0.2, sum_diff_cutoff = 2, cor_cutoff = 0.95, max_iter = 20)
  }

  metadata_spots <- metadata_spots %>% mutate(Subclones = ifelse(metadata_spots$SpotID %in% ref, "Reference",
                                                                 ifelse(metadata_spots$SpotID %in% query, "Query", "Stromal")))
  metadata_spots$Subclones <- ifelse(metadata_spots$Subclones == "Query", yes = subclones_clusters[match(metadata_spots$SpotID, names(subclones_clusters))], no = metadata_spots$Subclones)

  malig_cut <- sample_parameters$CNA_Malig_cut[[i]]
  metadata_spots$test <- ifelse(!is.na(metadata_spots$CNAscore) & metadata_spots$CNAscore > malig_cut, "yes", "no")
  plot_spatial_features(metadata = metadata_spots, image_obj = spatial_image, color_by = "test", pt.size = 1.5, image.alpha = 0, rm_spots = NULL, na.value = "gray", stroke = 0)

  order_cells <- split(metadata_spots$SpotID, metadata_spots$Subclones)
  ref_stroma_spots <- order_cells[c("Reference", "Stromal")]
  order_cells[c("Reference", "Stromal")] <- NULL
  mixed_spots <- lapply(order_cells, function(x) metadata_spots[metadata_spots$SpotID %in% x & !is.na(metadata_spots$CNAscore) & metadata_spots$CNAscore <= malig_cut, ] %>% .[order(.$CNAcorr[!is.na(.$CNAcorr)]), ] %>%
                          pull(SpotID)) %>% unlist(., use.names = FALSE)
  order_cells <- lapply(order_cells, function(x) metadata_spots[metadata_spots$SpotID %in% x & !is.na(metadata_spots$CNAscore) & metadata_spots$CNAscore > malig_cut, ] %>% .[order(.$CNAcorr[!is.na(.$CNAcorr)]), ] %>% pull(SpotID))
  order_cells[lengths(order_cells) < 100] <- NULL
  names(order_cells) <- seq(1, length(order_cells))
  order_cells_full <- c(order_cells, list(Mixed = c(mixed_spots, ref_stroma_spots$Stromal)), list(Non_Malignant = ref_stroma_spots$Reference))
  order_cells_full$Mixed <- metadata_spots %>% dplyr::filter(SpotID %in% order_cells_full$Mixed) %>% dplyr::arrange(desc(CNAscore)) %>% pull(SpotID)
  # saveRDS(order_cells_full, file = paste0(here("results/Per_Sample_Results/"), sample_parameters$Sample[[i]], "/", sample_parameters$Sample[[i]], "_CNA_Subclones.rds"))

  cna_score_mat2 <- calc_cna(matrix = proc_mat, query = c(query, order_cells_full$Mixed), ref = order_cells_full$Non_Malignant, genome = hg38, range = c(-3,3), window = 100, noise = 0.15, isLog = TRUE, per_chr = TRUE, scale = NULL, top_genes_num = NULL, verbose = TRUE)
  sig_and_cor2 <- cna_sig_cor(cna_score_mat2, epi_cells = c(query, order_cells_full$Mixed), ref_cells = order_cells_full$Non_Malignant)
  metadata_spots$CNAsignal <- sig_and_cor2$CNA_Signal[match(metadata_spots$SpotID, names(sig_and_cor2$CNA_Signal))]
  metadata_spots$CNAcorr <- sig_and_cor2$CNA_Correlation[match(metadata_spots$SpotID, names(sig_and_cor2$CNA_Correlation))]
  metadata_spots$CNAscore <- metadata_spots$CNAcorr + (metadata_spots$CNAsignal * max(metadata_spots$CNAcorr[!is.na(metadata_spots$CNAcorr)]) / max(metadata_spots$CNAsignal[!is.na(metadata_spots$CNAsignal)]))
  metadata_spots$Subclone <- ifelse(metadata_spots$SpotID %in% Unlist(order_cells_full), names(Unlist(order_cells_full)[match(metadata_spots$SpotID, Unlist(order_cells_full))]), NA)
  order_cells_full <- lapply(order_cells_full, function(clone) {
    metadata_spots %>%
      dplyr::filter(SpotID %in% clone) %>%
      dplyr::arrange(desc(CNAscore)) %>% pull(SpotID)
  })
  ordered_cna_mat2 <- cna_score_mat2[, unlist(order_cells_full)]
  cna_plot2 <- infercna::ggcna(ordered_cna_mat2, groups = order_cells_full, interorder = FALSE, genome = "hg38", title = "", legend.title = "Inferred CNA\n[log2 ratio]\n", legend.height = 0.6) +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "right") + labs(y = "Spot", x = "Chromosome")
  ggplot2::ggsave(filename = paste0(here("results/Per_Sample_Results/"), sample_parameters$Sample[[i]], "/", sample_parameters$Sample[[i]], "_CNA_Subclones_Plot.png"),
                  plot = cna_plot2, width = 35, height = 25, units = "cm")

  clone_purity <- metadata_spots %>% dplyr::select(SpotID, CNAscore) %>%
    arrange(match(SpotID, Unlist(order_cells_full))) %>% reshape2::melt(.) %>%
    dplyr::mutate(SpotID = factor(SpotID, levels = Unlist(order_cells_full)))
  p <- gmap(clone_purity, x = variable, y = SpotID, fill = value, limits = c(0, round(max(clone_purity$value), digits = 1)), ratio = 15, magma_pal = T, midpoint = malig_cut, geom = "raster",
            axis.rel = 1.5, legend.title.rel = 1.25, legend.rel = 1.0, legend.height = 5, legend.width = 0.6, x.num = F, y.num = F, angle = T) + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank())
  ggplot2::ggsave(filename = paste0(here("results/Per_Sample_Results/"), sample_parameters$Sample[[i]], "/", sample_parameters$Sample[[i]], "_CNA_Subclones_Purity.png"),
                  plot = p, width = 10, height = 25, units = "cm")

  metadata_spots$Subclone <- ifelse(metadata_spots$SpotID %in% Unlist(order_cells_full), names(Unlist(order_cells_full)[match(metadata_spots$SpotID, Unlist(order_cells_full))]), NA)
  metadata_spots %>% dplyr::group_by(Subclone) %>% dplyr::summarise(mean(CNAscore))
  data_dir <- paste(here("data"), sample_parameters$Site[[i]], sample_parameters$Sample[[i]], sep = "/")
  spatial_image <- load_spatial_image(paste0(data_dir, "/spatial"))
  clone_cols <- setNames(c("black", "#0c2a50ff", "#BC80BD", "#88CCEE", "#FCCDE5", "#33A02C", "#FFED6F", "#FE9929", "#E31A1C", "#D9D9D9", "#B15928"),
                         c("Non_Malignant", "Mixed", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
  p <- plot_spatial_features(metadata = metadata_spots, image_obj = spatial_image, color_by = "Subclone", cols = clone_cols, pt.size = 1.5, image.alpha = 0, rm_spots = NULL, na.value = "gray", stroke = 0) +
    guides(fill = guide_legend(override.aes = list(size = 5)))
  ggplot2::ggsave(filename = paste0(here("results/Per_Sample_Results/"), sample_parameters$Sample[[i]], "/", sample_parameters$Sample[[i]], "_Subclones_Spatial_Distribution.png"),
                  plot = p, width = 16, height = 14, units = "cm")

  # Get a threshold dependent malignancy status
  order_cells_full <- readRDS(paste0(here("results/Per_Sample_Results/"), sample_parameters$Sample[[i]], "/", sample_parameters$Sample[[i]], "_CNA_Subclones.rds"))
  malig_spots <- Unlist(order_cells_full[!names(order_cells_full) %in% c("Mixed", "Non_Malignant")])
  subclones_clusters <- setNames(names(Unlist(order_cells_full[names(order_cells_full) != "Non_Malignant"])), unname(Unlist(order_cells_full[names(order_cells_full) != "Non_Malignant"])))
  epi_ref_m <- filter_10x(m, cells = Unlist(order_cells_full), centre = "none", merge_method = "sum", gene_filt_method = "exp",
                          complexity_thresh = sample_parameters$Complexity_cut[[i]], percent_mt_thresh = sample_parameters$MT_cut[[i]])

  malig_assign <- malig_status(epi_ref_m, as.factor(subclones_clusters), ref_cells = order_cells_full$Non_Malignant, top_region = 1/3, top_cells = 1/4, class_cut = 0.95, KNN_correct = FALSE)
  metadata_spots$test <- ifelse(metadata_spots$SpotID %in% malig_spots, "Malignant",
                                ifelse(metadata_spots$SpotID %in% malig_assign[names(malig_assign) %in% c("Normal", "Reference") & !names(malig_assign) %in% malig_spots], "Non_Malignant", "Mixed"))
  plot_spatial_features(metadata = metadata_spots, image_obj = spatial_image, color_by = "test", pt.size = 1.5, image.alpha = 0, rm_spots = NULL, na.value = "gray", stroke = 0)

  metadata_spots$binCNAstatus <- ifelse(metadata_spots$SpotID %in% malig_spots, "Malignant",
                                        ifelse(metadata_spots$SpotID %in% malig_assign[names(malig_assign) %in% c("Normal", "Reference") & !names(malig_assign) %in% malig_spots], "Non_Malignant", "Mixed"))
  metadata[metadata$Sample == sample_parameters$Sample[[i]], ] <- metadata_spots %>% dplyr::select(-test, -Subclones)
}

# saveRDS(metadata, file = here("metadata/merged_metadata_Extend.rds"))


# Plot per sample purity (proportions of malignant, non-malignant and mixed spots) based on CNA inference
samples_metadata$Ext_Site <- ifelse(samples_metadata$Site == "Oropharynx" & samples_metadata$HPV == "HPVpos", "HPV+ OP",
                                    ifelse(samples_metadata$Site == "Oropharynx" & samples_metadata$HPV == "HPVneg", "HPV- OP", samples_metadata$Site))
per_sample_purity <- as.data.frame(table(metadata$binCNAstatus[!is.na(metadata$binCNAstatus)], metadata$Sample[!is.na(metadata$binCNAstatus)])) %>%
  rowwise %>% mutate(Site = as.factor(samples_metadata$Ext_Site[samples_metadata$Sample == Var2]), .before = "Freq") %>% ungroup %>%
  group_by(Var2) %>% mutate(Prop = Freq / sum(Freq), OrderBy = Prop[Var1 == "Malignant"]) %>% ungroup
p <- per_sample_purity %>%
  dplyr::mutate(Group = ifelse(Var1 %in% "Malignant", "Malignant", "Non-Malignant")) %>%
  dplyr::group_by(Var1, Var2) %>%
  summarize(Malignant = sum(Prop[Group == "Malignant"]), `Non-Malignant` = -sum(Prop[Group == "Non-Malignant"])) %>%
  dplyr::group_by(Var2) %>%
  dplyr::mutate(Ord = sum(`Non-Malignant`)) %>% ungroup %>%
  dplyr::mutate(Site = factor(samples_metadata$Ext_Site[match(Var2, samples_metadata$Sample)], levels = c("Oral", "Laryngeal", "HPV+ OP", "HPV- OP"))) %>%
  ggplot(aes(x = fct_reorder(Var2, Ord), y = Malignant, fill = Var1)) +
  geom_col(position = "stack") +
  geom_col(aes(y = `Non-Malignant`), position = position_stack(reverse = TRUE)) +
  geom_hline(aes(yintercept = 0)) +
  scale_fill_manual(name = "Malignancy\nStatus", values = setNames(c("#FCF2B4FF", "#FD9969FF", "#51127CFF"), c("Non_Malignant", "Mixed", "Malignant"))) + scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  theme_classic() + theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.9, size = 16), plot.title = element_text(hjust = 0.5, size = 16),
                          legend.text = element_text(size = 14), legend.title = element_text(size = 16)) +
  facet_grid(.~Site, scales = 'free', space = 'free', drop = F)
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_1c.pdf"), plot = p, device = "pdf", dpi = 300, width = 10, height = 6)

# Percentage of spots at each purity category
as.data.frame(table(metadata$binCNAstatus)) %>% dplyr::mutate(Percentage = round(Freq / sum(Freq), digits = 2) * 100)




# Add Spots Coordinates to Metadata ---------------------------------------

metadata <- readRDS(file = here("metadata/merged_metadata_Extend.rds"))
spat_coord <- list()
for(i in seq_along(sample_parameters$Sample)) {
  data_dir <- paste(here("data"), sample_parameters$Site[[i]], sample_parameters$Sample[[i]], sep = "/")
  # Load sample tissue position table
  spatial_image <- load_spatial_image(paste0(data_dir, "/spatial"))
  spatial_image$Coordinates$SpotID <- paste(sample_parameters$Sample[[i]], spatial_image$Coordinates$SpotID, sep = "_")
  coord_tab <- spatial_image$Coordinates
  # Load sample measurment table (QuPath output)
  spot_measure <- read.delim(file.path(here("results/Spot_measurements"), grep(sample_parameters$Sample[[i]], list.files(here("results/Spot_measurements")), value = TRUE)), header = TRUE)
  spot_measure <- spot_measure[-1, c("array_row", "array_col", "cx", "cy", "Num.Detections")] %>% dplyr::rename("pxl_in_rows" = "cx", "pxl_in_cols" = "cy", "Cell_Count" = "Num.Detections")
  image_data <- left_join(coord_tab, spot_measure, by = c("array_row", "array_col", "pxl_in_rows", "pxl_in_cols"))
  # Output to a list
  spat_coord[[i]] <- image_data
  names(spat_coord)[i] <- sample_parameters$Sample[[i]]
}
spat_coord <- do.call(rbind.data.frame, spat_coord)
metadata <- inner_join(metadata, spat_coord[, -which(colnames(spat_coord) == "in_tissue")], by = c("Key" = "SpotID"))
# saveRDS(metadata, file = here("metadata/merged_metadata_Extend.rds"))





# Calculate data quality statistics -------------------------------------------------------------------

# Spatial Data
cohort_qc_stats <- tibble()
for(i in seq_along(sample_parameters$Sample)) {
  message(paste0("Processing sample: ", sample_parameters$Sample[[i]]))

  # Load expression matrix
  data_dir <- paste(here("data"), sample_parameters$Site[[i]], sample_parameters$Sample[[i]], sep = "/")
  m <- as.matrix(Matrix::readMM(paste(data_dir, grep("matrix(_[A-Za-z0-9]+)?\\.mtx$", list.files(data_dir), value = TRUE), sep = "/")))
  colnames(m) <- read.table(paste(data_dir, grep("barcodes(_[A-Za-z0-9]+)?\\.tsv$", list.files(data_dir), value = TRUE), sep = "/"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)[, 1]
  colnames(m) <- paste(sample_parameters$Sample[[i]], colnames(m), sep = "_")
  rownames(m) <- read.table(paste(data_dir, grep("features(_[A-Za-z0-9]+)?\\.tsv$", list.files(data_dir), value = TRUE), sep = "/"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)[, 2]

  if(sum(grep("GRCh38", rownames(m))) == 0) {
    m <- m
  } else {
    # Remove HPV genes and GRC38 prefix
    m <- m[-grep("^hpv", rownames(m)), ]
    rownames(m) <- sub("GRCh38_", "", rownames(m))
  }

  # Process raw data (no centering)
  proc_mat <- filter_10x(m, complexity_thresh = sample_parameters$Complexity_cut[[i]], percent_mt_thresh = sample_parameters$MT_cut[[i]],
                         centre = "none", gene_filt_method = "exp", merge_method = "sum")

  # Create a metadata tibble to store all collected variables
  cell_stats <- cell_qc(m) %>% dplyr::mutate(Depth = Matrix::colSums(m), .before = Complexity) %>% dplyr::filter(CellID %in% colnames(proc_mat))

  # Join sample QC stats to merged QC stats tibble
  cohort_qc_stats <- rbind(cohort_qc_stats, cell_stats)
}
# readr::write_csv(cohort_qc_stats, file = here("results/Generated_Data/Spatial_Cohort_QC_Data_Extend.csv"))

cohort_qc_stats <- readr::read_csv(file = here("results/Generated_Data/Spatial_Cohort_QC_Data_Extend.csv"))
ranger_report_stats <- readr::read_csv(file = here("results/Generated_Data/Spatial_QC_Stats.csv"))[, -1]
qc_stats_df <- data.frame(Nsamples = length(unique(scalop::substri(cohort_qc_stats$CellID, sep = "_", pos = 1))),
                          Nspots = nrow(cohort_qc_stats),
                          `Median Genes per Spot` = median(cohort_qc_stats$Complexity),
                          `Median UMI Counts per Spot` = median(cohort_qc_stats$Depth),
                          `Mean Reads per Spot` = round(mean(ranger_report_stats$`Mean Reads per Spot`)))
# readr::write_csv(qc_stats_df, file = here("results/Generated_Data/Spatial_Cohort_Summary_QC_Stats.csv"))

