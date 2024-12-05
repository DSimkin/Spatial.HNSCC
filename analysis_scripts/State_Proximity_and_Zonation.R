library(tidyverse)
library(scalop)
library(ggpubr)
library(here)
source(here("analysis_scripts/functions/scRNA_Functions.R"))
source(here("analysis_scripts/functions/ST_Functions.R"))


# Load per-sample parameter table
sample_parameters <- readRDS(file = here("metadata/per_sample_parameters_df.rds"))
# Load sample level metadata
samples_metadata <- readRDS(file = here("metadata/samples_metadata.rds"))
# Load metaprograms
MPs <- readRDS(file = here("results/Generated_Data/Final_Metaprograms_Extended.rds"))
# Load merged metadata tibble
metadata <- readRDS(file = here("metadata/merged_metadata_Extend.rds"))
# Load scoring-based deconvolution matrices
decon_mats_ls <- readRDS(file = here("results/Generated_Data/Decon_Mats_List_Extend.rds"))
# Load state colors
state_cols <- read.delim(file = here("aux_data/state_col_tab_extend.tsv"), sep = "\t", header = TRUE)
state_cols <- setNames(state_cols$V2, state_cols$V1)


# sample_parameters <- sample_parameters %>% dplyr::filter(Sample %in% basename(list.dirs(here("data"), recursive = TRUE)))

# Zone Classification & Enrichment ----------------------------------------

EpiStroma_plots <- list()
Zones_plots <- list()
Enrichment_plots <- list()
abundance_tabs <- list()
state_per_zone_prop <- list()

for(i in seq_along(sample_parameters$Sample)) {
  message(paste0("Processing sample: ", sample_parameters$Sample[[i]]))

  path <- paste(here("data"), sample_parameters$Site[[i]], sample_parameters$Sample[[i]], sep = "/")
  spots_positions <- read.csv(paste0(path, "/spatial/tissue_positions_list.csv"), header = FALSE)
  colnames(spots_positions) <- c("SpotID", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")

  # Add score-based deconvolution results to the metadata
  decon_mat <- decon_mats_ls[[sample_parameters$Sample[[i]]]]
  metadata_spots <- tibble(SpotID = rownames(decon_mat)) %>% left_join(., spots_positions, by = "SpotID")

  # Classify zones
  zones <- classify_zones(metadata_spots, decon_mat, with_plots = TRUE, zones_limits = c(3, 6), stromal_prop_cutoff = 0.7, rm_distant_spots = TRUE)
  EpiStroma_plots[[i]] <- zones$EpiStroma_plot
  names(EpiStroma_plots)[i] <- sample_parameters$Sample[[i]]
  Zones_plots[[i]] <- zones$Zones_plot
  names(Zones_plots)[i] <- sample_parameters$Sample[[i]]
  merged_meta <- zones$merged_meta

  # State enrichment per zone
  zone_enrich <- zone_enrichment(merged_meta, decon_mat, min_spots_per_zone = 70, set_colors = state_cols)
  abundance_tabs[[i]] <- zone_enrich$mean_zone_abund
  names(abundance_tabs)[i] <- sample_parameters$Sample[[i]]
  state_per_zone_prop[[i]] <- zone_enrich$zone_prop
  names(state_per_zone_prop)[i] <- sample_parameters$Sample[[i]]
  Enrichment_plots[[i]] <- zone_enrich$enrichment_plot
  names(Enrichment_plots)[i] <- sample_parameters$Sample[[i]]
}

# saveRDS(EpiStroma_plots, file = here("results/Generated_Data/All_samples_EpiStroma_plots_Extend.rds"))
# saveRDS(Zones_plots, file = here("results/Generated_Data/All_samples_Zones_plots_Extend.rds"))
# saveRDS(Enrichment_plots, file = here("results/Generated_Data/All_samples_StateEnrichment_plots_Extend.rds"))
# saveRDS(abundance_tabs, file = here("results/Generated_Data/All_samples_ZoneAbundance_tables_Extend.rds"))
# saveRDS(state_per_zone_prop, file = here("results/Generated_Data/All_samples_ZoneProportions_tables_Extend.rds"))
EpiStroma_plots <- readRDS(file = here("results/Generated_Data/All_samples_EpiStroma_plots_Extend.rds"))
Zones_plots <- readRDS(file = here("results/Generated_Data/All_samples_Zones_plots_Extend.rds"))
Enrichment_plots <- readRDS(file = here("results/Generated_Data/All_samples_StateEnrichment_plots_Extend.rds"))
abundance_tabs <- readRDS(file = here("results/Generated_Data/All_samples_ZoneAbundance_tables_Extend.rds"))
state_per_zone_prop <- readRDS(file = here("results/Generated_Data/All_samples_ZoneProportions_tables_Extend.rds"))

# Add Zones to metadata
zones_ls <- lapply(seq_along(Zones_plots), function(i) {
  Zones_plots[[i]]$data %>% dplyr::select(SpotID, Zone) %>% dplyr::mutate(Key = paste0(names(Zones_plots)[i], "_", SpotID), .before = 1) %>% dplyr::select(-SpotID)
})
zones_meta <- do.call(rbind, zones_ls)
metadata <- dplyr::left_join(metadata, zones_meta, by = "Key")
# saveRDS(metadata, file = here("metadata/merged_metadata_Extend.rds"))

# Plot per sample states enrichment
p <- Enrichment_plots$P5736 + theme(axis.title.x = element_blank(), text = element_text(size = 14), legend.title = element_text(size = 16),
                                                 panel.border = element_rect(fill = NA, colour = "black", linewidth = 2), aspect.ratio = 1.4, plot.margin = margin(0.2, , 0.2, 0.4, "cm")) +
  ylab("Z-Score") + guides(color = guide_legend(override.aes = list(linewidth = 7))) +
  scale_x_discrete(expand = c(0, 0), labels = function(x) gsub("_", " ", x, fixed = TRUE)) + scale_y_continuous(expand = c(0, 0))  # optionally add + theme_bw() befor theme()
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_4a_right.pdf"), device = "pdf", plot = p, dpi = 300)

# Plot per state zone enrichment
sample_cols <- setNames(samples_metadata$Color, samples_metadata$Sample)
orders_MPs <- c("Skeletal_Muscle", "B_cell", "T_cell", "Macrophage", "Complement", "Fibroblast", "Endothelial",
                "pEMT", "Inflammatory", "LowQ", "Cell_Cycle", "Epithelial", "Secretory_Malig", "Senescence", "Hypoxia")
states_zonation_plots <- lapply(orders_MPs, function(state) {
  state_zone_enrichment(state = state, abundance_tabs, colors = sample_cols, trend_style = "linear", trend_line_size = 2, exlude_samp = "P5666") +
    ggtitle(state) + theme(legend.position = "none", axis.text.x = element_text(size = 8), plot.margin = margin(0.2, 0.6, 0.2, 0.4, "cm"))
})
p <- do.call(ggarrange, states_zonation_plots)
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_4a.pdf"), device = "pdf", plot = p, dpi = 300)



# Plot all states trends
summary_tab <- sapply(names(MPs)[!names(MPs) %in% "Secretory_Norm"], function(state) {
  per_sample_zone_zscores <- lapply(abundance_tabs, function(x) x[state, ]) %>% do.call(plyr::rbind.fill, .) %>%
    dplyr::rename_all(~colnames(abundance_tabs[[1]])) %>% magrittr::set_rownames(names(abundance_tabs)) %>%
    rownames_to_column(var = "Sample") %>% dplyr::mutate(Samples = factor(Sample, levels = unique(Sample)), .keep = "unused", .before = 1) %>% column_to_rownames(var = "Samples") %>% t() %>% scale() %>% t() %>%
    as.data.frame()
  state_per_zone_means <- sapply(per_sample_zone_zscores, function(x) mean(na.omit(x)))
})
summary_plot <- reshape2::melt(summary_tab)
p <- ggplot(summary_plot, aes(x = Var1, y = value, group = Var2)) +
  geom_line(size = 2, aes(color = Var2)) + scale_color_manual(values = state_cols, na.translate = FALSE) +
  theme_bw() + theme(axis.title = element_blank(), text = element_text(size = 14), legend.title = element_text(size = 16), plot.margin = margin(0.2, , 0.2, 0.5, "cm"),
                     panel.border = element_rect(fill = NA, colour = "black", linewidth = 3), aspect.ratio = 1.2) +
  guides(color = guide_legend(override.aes = list(linewidth = 7), title = "State")) + scale_x_discrete(expand = c(0, 0), labels = function(x) gsub("_", " ", x, fixed = TRUE))
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_4b.pdf"), device = "pdf", plot = p, dpi = 300)



# Plot pEMT state zonation
sample_cols <- setNames(samples_metadata$Color, samples_metadata$Sample)
pEMT_zonation <- state_zone_enrichment(state = "pEMT", abundance_tabs, colors = sample_cols, trend_style = "smooth", ratio = 1.2, panel_border_size = 2)
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_4c_top.pdf"), plot = pEMT_zonation, device = "pdf", dpi = 300)

# Different peak patterns of pEMT
pEMT_zonation <- state_zone_enrichment(state = "pEMT", abundance_tabs, colors = sample_cols, with_plot = FALSE) %>% replace(is.na(.), 0)
pattern_1 <- pEMT_zonation %>%  mutate(Class = names(.)[max.col(.)], Sample = sample_parameters$Sample) %>% dplyr::filter(Class == "Zone_1")
pat1_plot <- state_zone_enrichment(state = "pEMT", abundance_tabs, colors = sample_cols, choose_samp = pattern_1$Sample, show_trend = FALSE, ratio = 1.2, panel_border_size = 2, line_size = 2)
pattern_2 <- pEMT_zonation %>%  mutate(Class = names(.)[max.col(.)], Sample = sample_parameters$Sample) %>% dplyr::filter(Class == "Zone_2")
state_zone_enrichment(state = "pEMT", abundance_tabs, colors = sample_cols, choose_samp = pattern_2$Sample, show_trend = FALSE)
pattern_3 <- pEMT_zonation %>%  mutate(Sample = sample_parameters$Sample) %>% dplyr::filter(pEMT_zonation$Zone_3 > pEMT_zonation$Zone_2 & pEMT_zonation$Zone_2 > pEMT_zonation$Zone_1)
pat3_plot <- state_zone_enrichment(state = "pEMT", abundance_tabs, colors = sample_cols, choose_samp = pattern_3$Sample, show_trend = FALSE, ratio = 1.2, panel_border_size = 2, line_size = 2)
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_4c_bottom_left.pdf"), plot = pat1_plot, device = "pdf", dpi = 300)
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_4c_bottom_right.pdf"), plot = pat3_plot, device = "pdf", dpi = 300)

samples_metadata <- samples_metadata %>%
  dplyr::mutate(pEMT_Pattern = dplyr::case_when(Sample %in% c("P5219", "P5624", "P5709", "P5418", "P5659", "P5909") ~ "Core_Enriched",
                                                Sample %in% c("P5448", "P5707", "P5666", "P5684", "P5692", "P5710", "P5764", "P5766", "P5767", "P5756") ~ "Margin_Enriched", .default = "none"))
# saveRDS(samples_metadata, file = here("metadata/samples_metadata.rds"))



# Plot pEMT scores heatmap for margin-enriched sample vs core-enriched sample
# Load spots MP scores
score_vecs <- readRDS(here("results/Generated_Data/per_sample_state_score_vectors_Extend.rds"))

# Merge score list to matrix that represent all states over all spots
make_equal_lengths <- lapply(score_vecs, function(x) {
  if(!all(names(MPs) %in% names(x))) {
    na_vec <- list(setNames(rep(NA, length(x[[1]])), names(x[[1]])))
    for(y in which(!names(MPs) %in% names(x))) {
      x <- append(x, na_vec, after = y - 1)
    }
    x <- x %>% purrr::set_names(names(MPs))
  } else {
    x <- x
  }
})
state_score_mat <- sapply(seq_along(names(MPs)), function(i) do.call(c, lapply(make_equal_lengths, function(x) x[[i]]))) %>%
  `row.names<-`(scalop::substri(row.names(.), sep = "\\.", pos = 2)) %>% magrittr::set_colnames(names(MPs))

order_spot_scores_df <- as.data.frame(state_score_mat) %>%
  rowwise() %>%
  dplyr::mutate(MP = names(.)[which.max(c_across(everything()))], max_val = max(c_across(Fibroblast:Secretory_Norm), na.rm = TRUE)) %>% ungroup() %>%
  dplyr::mutate(SpotID = rownames(state_score_mat), Zone = metadata$Zone[match(SpotID, metadata$Key)], .before = 1)

p5756_hm <- order_spot_scores_df %>% dplyr::mutate(Sample = scalop::substri(SpotID, sep = "_", pos = 1)) %>%
  dplyr::filter(Sample == "P5756", Zone != "Filtered_out", !is.na(Zone)) %>% droplevels() %>%
  dplyr::select(SpotID, Zone, pEMT) %>%
  dplyr::group_by(Zone) %>% arrange(desc(pEMT), .by_group = TRUE) %>% dplyr::ungroup() %>%
  dplyr::mutate(SpotID = factor(SpotID, levels = unique(SpotID)), State = factor("pEMT"))
p <- gmap(p5756_hm, x = SpotID, y = State, fill = pEMT, limits = c(-2, 2), y.name = "", midpoint = 0, tile.size = 0,
          legend.height = 5, legend.width = 0.6, x.num = T, y.num = F, minimal = TRUE, geom = "raster") +
  theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 0.3, axis.ticks = element_blank(), axis.text.x = element_blank(), text = element_text(size = 16)) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  facet_grid(. ~ Zone, scales = "free", labeller = labeller(Zone = function(x) gsub("_", " ", x, fixed = TRUE)))
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_4c_bottom_heatmap_1.pdf"), plot = p, device = "pdf", dpi = 300)

p5659_hm <- order_spot_scores_df %>% dplyr::mutate(Sample = scalop::substri(SpotID, sep = "_", pos = 1)) %>%
  dplyr::filter(Sample == "P5659", Zone != "Filtered_out", !is.na(Zone)) %>% droplevels() %>%
  dplyr::select(SpotID, Zone, pEMT) %>%
  dplyr::group_by(Zone) %>% arrange(desc(pEMT), .by_group = TRUE) %>% dplyr::ungroup() %>%
  dplyr::mutate(SpotID = factor(SpotID, levels = unique(SpotID)), State = factor("pEMT"))
p <- gmap(p5659_hm, x = SpotID, y = State, fill = pEMT, limits = c(-2, 2), y.name = "", midpoint = 0, tile.size = 0,
          legend.height = 5, legend.width = 0.6, x.num = T, y.num = F, minimal = TRUE, geom = "raster") +
  theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 0.3, axis.ticks = element_blank(), axis.text.x = element_blank(), text = element_text(size = 16)) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  facet_grid(. ~ Zone, scales = "free", labeller = labeller(Zone = function(x) gsub("_", " ", x, fixed = TRUE)))
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_4c_bottom_heatmap_2.pdf"), plot = p, device = "pdf", dpi = 300)





# pEMT-Fibro vs pEMT-INF Colocalization -----------------------------------

get_pEMT_grouping <- lapply(sample_parameters$Sample, function(samp) {
  samp_meta <- metadata %>% dplyr::filter(Sample == samp)
  samp_neighb_tab <- .neighbors_table(samp_meta, var = "MPid", spot_id = "Key") %>%
    dplyr::select(-V1) %>% dplyr::rename_all(~ stringr::str_replace_all(., "V", "Neighb_")) %>%
    dplyr::rename_all(~ str_replace_all(., pattern = "\\d", function(num) as.character(as.numeric(num) - 1))) %>%
    tibble::rownames_to_column(var = "Key")
  samp_meta <- left_join(samp_meta, samp_neighb_tab, by = "Key")
  samp_meta <- samp_meta %>% dplyr::filter(pEMT > 0.25) %>%  #dplyr::filter(pEMT > 0 & (Fibroblast > 0 | Inflammatory > 0)) %>%
    dplyr::mutate(pEMT_Coloc = case_when(pEMT > 0.25 & Fibroblast > 0 & Inflammatory == 0 ~ "pEMT-Fibro",
                                         pEMT > 0.25 & Inflammatory > 0 & Fibroblast == 0 ~ "pEMT-INF",
                                         pEMT > 0.25 & Fibroblast > 0 & Inflammatory > 0 ~ "pEMT-Both",
                                         pEMT > 0.25 & Fibroblast == 0 & Inflammatory == 0 ~ "pEMT-None"))

  samp_meta <- samp_meta %>% dplyr::mutate(pEMT_Neighb = ifelse(rowSums(across(starts_with("Neighb_"), ~. %in% "pEMT")) != 0, TRUE, FALSE),
                                           INF_Neighb = ifelse(rowSums(across(starts_with("Neighb_"), ~. %in% "Inflammatory")) != 0, TRUE, FALSE),
                                           Fibro_Neighb = ifelse(rowSums(across(starts_with("Neighb_"), ~. %in% "Fibroblast")) != 0, TRUE, FALSE)) %>%
    dplyr::mutate(pEMT_Neighbors = case_when(pEMT_Neighb & !INF_Neighb & !Fibro_Neighb ~ "pEMT",
                                             INF_Neighb & !pEMT_Neighb & !Fibro_Neighb ~ "INF",
                                             Fibro_Neighb & !pEMT_Neighb & !INF_Neighb ~ "Fibroblast",
                                             pEMT_Neighb & INF_Neighb & !Fibro_Neighb ~ "pEMT-INF",
                                             pEMT_Neighb & Fibro_Neighb & !INF_Neighb ~ "pEMT-Fibro",
                                             INF_Neighb & Fibro_Neighb & !pEMT_Neighb ~ "INF-Fibro",
                                             pEMT_Neighb & INF_Neighb & Fibro_Neighb ~ "pEMT-INF-Fibro",
                                             !pEMT_Neighb & !INF_Neighb & !Fibro_Neighb ~ "None")) %>%
    dplyr::select(-dplyr::ends_with("_Neighb"))
  out <- samp_meta %>% dplyr::select(Key, Sample, MPid, Fibroblast:Secretory_Norm, dplyr::starts_with("Neighb_"), dplyr::starts_with("pEMT_"), Zone)
})
pEMT_coloc <- do.call(rbind.data.frame, get_pEMT_grouping)
pEMT_coloc %>% dplyr::mutate(Adj_Activator = ifelse(pEMT_Coloc == "pEMT-None" & pEMT_Neighbors == "None", FALSE, TRUE)) %>% pull(Adj_Activator) %>% table   # 3709 spots found to have no colocalization or neighboring pEMT//Fibroblasts//INF states
pEMT_coloc %>% dplyr::group_by(Sample, pEMT_Coloc) %>% dplyr::summarise(Freq = n()) %>% tidyr::pivot_wider(names_from = pEMT_Coloc, values_from = Freq)
pEMT_coloc %>% dplyr::group_by(Sample, pEMT_Coloc, pEMT_Neighbors) %>% dplyr::summarise(Freq = n()) %>%
  tidyr::pivot_wider(names_from = c(pEMT_Coloc, pEMT_Neighbors), values_from = Freq)
samp_order <- samples_metadata$Sample[order(samples_metadata$Sample)]
pEMT_coloc %>% dplyr::group_by(Sample, pEMT_Coloc, pEMT_Neighbors) %>% dplyr::summarise(Freq = n()) %>%
  tidyr::pivot_wider(names_from = c(pEMT_Coloc, pEMT_Neighbors), values_from = Freq) %>% tibble::column_to_rownames(var = "Sample") %>%
  dplyr::rowwise() %>% dplyr::mutate(Max = names(.)[which.max(c_across(everything()))]) %>% dplyr::ungroup() %>% dplyr::mutate(Sample = samp_order, .before = 1) %>%
  dplyr::pull(Sample, Max)

pEMT_coloc_df <- pEMT_coloc %>% dplyr::group_by(Sample, pEMT_Coloc) %>% dplyr::summarise(Freq = n()) %>% tidyr::pivot_wider(names_from = pEMT_Coloc, values_from = Freq) %>% ungroup
pEMT_no_coloc_and_neighbs <- pEMT_coloc %>% dplyr::group_by(Sample, pEMT_Coloc, pEMT_Neighbors) %>% dplyr::summarise(Freq = n()) %>%
  tidyr::pivot_wider(names_from = c(pEMT_Coloc, pEMT_Neighbors), values_from = Freq) %>% dplyr::select(Sample, "pEMT-None_None") %>% ungroup
pEMT_coloc_df <- dplyr::left_join(pEMT_coloc_df, pEMT_no_coloc_and_neighbs, by = "Sample") %>%
  dplyr::mutate(pEMT_Pattern = case_when(Sample %in% c("P5219", "P5624", "P5709", "P5418", "P5659", "P5909") ~ "pEMT-Core",
                                         Sample %in% c("P5448", "P5707", "P5666", "P5684", "P5692", "P5710", "P5764", "P5766", "P5767", "P5756") ~ "pEMT-Margin")) %>%
  dplyr::left_join(samples_metadata[, c("Sample", "Site", "HPV", "Color")], by = "Sample")
pEMT_coloc_df <- pEMT_coloc_df %>% dplyr::mutate(`pEMT-None` = `pEMT-None` - `pEMT-None_None`) %>% dplyr::select(-`pEMT-None_None`)
pEMT_coloc_df <- pEMT_coloc_df %>% dplyr::mutate(across(where(is.numeric),
                                                        ~./(rowSums(pEMT_coloc_df %>% dplyr::select(where(is.numeric)), na.rm = TRUE))))

plot_df <- pEMT_coloc_df %>% dplyr::select(-"pEMT-None") %>% dplyr::filter(!is.na(pEMT_Pattern))
plot_df <- reshape2::melt(plot_df)
spec_cols <- setNames(c("lightblue", unname(state_cols[c("Fibroblast", "Inflammatory")])), levels(plot_df$variable))
p <- ggpubr::ggboxplot(plot_df, x = "variable", y = "value", fill = "variable", facet.by = "pEMT_Pattern", ylab = "Number Of Spots Showing Colocalization", xlab = "", legend.title = "", outlier.shape = NA) +
  theme(strip.text = element_text(size = 16), panel.background = element_rect(fill = alpha("antiquewhite", 0.4)), legend.position = "top", axis.text = element_text(size = 14), axis.title.y = element_text(size = 16), legend.text = element_text(size = 16)) +
  scale_fill_manual(values = c("grey30", unname(state_cols[c("Fibroblast", "Inflammatory")])))
to_compare <- list(c("pEMT-Fibro", "pEMT-INF"))
p <- p + stat_compare_means(comparisons = to_compare, method = "t.test")
g <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-t', g$layout$name))
fills <- alpha(c("#B24745FF", "#DF8F44FF"), 0.8)
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(g)
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_5a.pdf"), plot = g, device = "pdf", dpi = 300)



## Per-spot pEMT, INF and Fibroblast coexpression heatmap with margin vs core annotation bar
# Load merged matrix
merged_expr_mat <- readRDS(file = here("results/Generated_Data/merged_expr_mat.rds"))
# Load per-sample spots MP score vectors
score_vecs <- readRDS(here("results/Generated_Data/per_sample_state_score_vectors_Extend.rds"))

score_vecs <- lapply(score_vecs, function(samp) {
  as.data.frame(samp) %>% dplyr::select(pEMT, Fibroblast, Inflammatory)
})
score_vecs <- do.call(rbind.data.frame, score_vecs) %>% `rownames<-`(scalop::substri(rownames(.), sep = "\\.", pos = 2))

feat_select <- Unlist(MPs[c("pEMT", "Fibroblast", "Inflammatory")])
filt_merged_expr_mat <-  filter_10x(merged_expr_mat[rownames(merged_expr_mat) %in% feat_select, ], centre = "mean", log = FALSE, raw_to_CPM = FALSE,
                                    complexity_thresh = NULL, percent_mt_thresh = NULL, gene_filt_method = "none", merge_method = "sum", dbl_remove = FALSE)

spot_select <- unique(c(metadata$Key[metadata$MPid == "pEMT"], pEMT_coloc$Key[pEMT_coloc$MPid %in% c("pEMT", "Fibroblast", "Inflammatory", "Hypoxia")]))
hm_meta <- metadata %>%
  dplyr::mutate(pEMT_Pattern = case_when(Sample %in% c("P5219", "P5624", "P5709", "P5418", "P5659", "P5909") ~ "pEMT-Core",
                                         Sample %in% c("P5448", "P5707", "P5666", "P5684", "P5692", "P5710", "P5764", "P5766", "P5767", "P5756") ~ "pEMT-Margin")) %>%
  dplyr::filter(Key %in% spot_select, !is.na(pEMT_Pattern))
expr_mat <- filt_merged_expr_mat[, hm_meta$Key]
score_df <- score_vecs[hm_meta$Key, ] %>% dplyr::mutate(`CAF-INF` = Fibroblast - Inflammatory)
x_order <- score_df %>% dplyr::arrange(desc(`CAF-INF`)) %>% rownames()
gene_diff_corr <- sapply(rownames(expr_mat), function(gene) cor(x = score_df$`CAF-INF`, y = expr_mat[gene, rownames(score_df)]))

gene_corr <- cor(t(expr_mat))
genes_hc <- hclust(dist(1 - gene_corr), method = "average")
plot_genes_corr <- reshape2::melt(gene_corr[genes_hc$order, genes_hc$order])
corr_plot <- ggplot(data = plot_genes_corr, aes(x = Var1, y = Var2, fill = value)) + geom_raster() +
  scale_fill_gradient2(limits = c(-0.5, 0.5), low = "dodgerblue4", mid = "antiquewhite", high = "red4", midpoint = mean(plot_genes_corr$value), oob = squish, name = "Correlation\nCoefficient") +
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), legend.text.align = 0.5) +
  scale_x_discrete(name = "\nPrograms", labels = NULL, breaks = NULL) + scale_y_discrete(name = "\nPrograms", labels = NULL, breaks = NULL)
genes_ordered <- colnames(gene_corr[genes_hc$order, genes_hc$order])
plot(genes_hc)
rect.hclust(genes_hc , k = 4, border = 2:6)
clusterCut <- stats::cutree(tree = genes_hc, k = 4)
gene_anno_df <- data.frame(row.names = genes_ordered,
                           cluster = as.factor(clusterCut[genes_ordered]),
                           genes = genes_ordered) %>%
  dplyr::mutate(MP = as.factor(case_when(genes %in% MPs$Fibroblast ~ "Fibroblast Genes",
                                         genes %in% MPs$pEMT ~ "pEMT Genes",
                                         genes %in% MPs$Inflammatory ~ "INF Genes")))
table(gene_anno_df$MP, gene_anno_df$cluster)
gene_anno_plot <- ggplot(gene_anno_df, aes(x = factor(genes, levels = genes), y = 1)) +
  geom_raster(aes(fill = cluster)) + theme(legend.position = "bottom", axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), axis.text = element_blank(), axis.line = element_blank(), axis.text.x = element_blank()) +
  labs(x = "", y = "", fill = "")
mp_anno_plot <- ggplot(gene_anno_df, aes(x = factor(genes, levels = genes), y = 1)) +
  geom_raster(aes(fill = MP)) + theme(legend.position = "bottom", axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), axis.text = element_blank(), axis.line = element_blank(), axis.text.x = element_blank()) +
  labs(x = "", y = "", fill = "")
egg::ggarrange(corr_plot, gene_anno_plot, mp_anno_plot, ncol = 1, nrow = 3, heights = c(40, 2, 2))

sort_genes <- names(c(sort(gene_diff_corr[match(gene_anno_df$genes[gene_anno_df$cluster == "1" & gene_anno_df$MP == "Fibroblast Genes"], names(gene_diff_corr))], decreasing = TRUE),
                      sort(gene_diff_corr[match(gene_anno_df$genes[gene_anno_df$cluster == "3" & gene_anno_df$MP == "Fibroblast Genes"], names(gene_diff_corr))]),
                      sort(gene_diff_corr[match(gene_anno_df$genes[gene_anno_df$cluster == "3" & gene_anno_df$MP == "pEMT Genes"], names(gene_diff_corr))]),
                      sort(gene_diff_corr[match(gene_anno_df$genes[gene_anno_df$cluster == "2" & gene_anno_df$MP == "pEMT Genes"], names(gene_diff_corr))]),
                      sort(gene_diff_corr[match(gene_anno_df$genes[gene_anno_df$cluster == "3" & gene_anno_df$MP == "INF Genes"], names(gene_diff_corr))]),
                      sort(gene_diff_corr[match(gene_anno_df$genes[gene_anno_df$cluster == "2" & gene_anno_df$MP == "INF Genes"], names(gene_diff_corr))])))
expr_mat <- expr_mat[rev(sort_genes), x_order]
plotting_df <- reshape2::melt(expr_mat)

# Create annotation bars
plotting_df <- plotting_df %>% dplyr::mutate(Pattern = as.factor(hm_meta$pEMT_Pattern[match(Var2, hm_meta$Key)]),
                                             MP = as.factor(case_when(Var1 %in% MPs$Fibroblast ~ "Fibroblast Genes",
                                                                      Var1 %in% MPs$pEMT ~ "pEMT Genes",
                                                                      Var1 %in% MPs$Inflammatory ~ "INF Genes")))
# saveRDS(plotting_df, file = here("results/Generated_Data/INF_CAF_pEMT_genes_coexprossion.rds"))
annotate_spots <- ggbar(plotting_df$Pattern, dir = "h", legend_title = "pEMT Pattern", cols = c("#B24745FF", "#DF8F44FF")) +
  theme(legend.direction = "horizontal", legend.text = element_text(size = 14, margin = margin(r = 0.8, unit = 'cm')),
        legend.title = element_text(size = 16, margin = margin(r = 0.8, unit = 'cm')))
genes2MP <- rev(ifelse(sort_genes %in% MPs$Fibroblast, "Fibroblast Genes",
                       ifelse(sort_genes %in% MPs$pEMT, "pEMT Genes",
                              ifelse(sort_genes %in% MPs$Inflammatory, "INF Genes", NA))))
annotate_genes <- ggbar(genes2MP, dir = "v", legend_title = "Metaprogram\nGenes", cols = unname(state_cols[c("Fibroblast", "Inflammatory", "pEMT")])) +
  theme(legend.direction = "horizontal", legend.text = element_text(size = 14, margin = margin(r = 0.8, unit = 'cm')),
        legend.title = element_text(size = 16, margin = margin(r = 0.8, unit = 'cm')), plot.margin = unit(c(0, 0, 0, 0), "null"), panel.margin = unit(c(0, 0, 0, 0), "null"))

lej1 = cowplot::get_plot_component(annotate_spots, 'guide-box-top', return_all = TRUE)
pdf(file = here("results/Paper_Figures/Main_figs/Fig_5b_pattern_legend.pdf"), width = 12, height = 3)
cowplot::ggdraw(lej1)
dev.off()
lej2 = cowplot::get_plot_component(annotate_genes, 'guide-box-top', return_all = TRUE)
pdf(file = here("results/Paper_Figures/Main_figs/Fig_5b_genes_legend.pdf"), width = 12, height = 3)
cowplot::ggdraw(lej2)
dev.off()

## Merge heatmap with annotation bars
p_hm <- gmap(plotting_df, x = Var2, y = Var1, fill = value, limits = c(-4, 4), ratio = NULL, y.name = "", x.labels = NULL, geom = "raster", y.labels = NULL,
     axis.rel = 1.5, legend.title.rel = 1.25, legend.rel = 1.0, legend.height = 5, legend.width = 0.6, x.num = F, y.num = F, angle = T) +
  theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), panel.border = element_rect(linewidth = 1),
        panel.spacing = unit(0, "lines"), plot.margin = unit(c(0, 0, 0, 0), "null"), panel.margin = unit(c(0, 0, 0, 0), "null"))
p_gene_anno <- annotate_genes + cowplot::theme_nothing()
p_spot_anno <- annotate_spots + cowplot::theme_nothing()

p_comb1 <- p_gene_anno + plot_spacer() + p_hm + patchwork::plot_layout(nrow = 1, ncol = 3, widths = c(0.1, -0.02, 1), guides = "collect")
p_comb2 <- plot_spacer() + p_spot_anno + patchwork::plot_layout(nrow = 1, ncol = 2, widths = c(0.092, 1))
p_final <- p_comb2 / plot_spacer() + p_comb1 + patchwork::plot_layout(nrow = 3, ncol = 1, heights = c(0.1, -0.04, 1), guides = "collect")
# ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_5b.pdf"), plot = p_final, device = "pdf", dpi = 300)



## Plot the differential MP abundace between pEMT core and margin samples
# Load abundances (as generated for the consensus interaction network)
pEMT_patterns_abund <- read.csv(file = here("results/Generated_Data/Node_Abundance_Neighborhood_Network.csv")) %>%
  dplyr::select(State, pEMTcore, pEMTmargin) %>% reshape2::melt(.) %>%
  dplyr::mutate(State = factor(stringr::str_replace_all(State, pattern = "TNF_Signaling", replacement = "Inflammatory"), levels = levels(metadata$MPid)))

per_sample_states <- as.data.frame(table(metadata$MPid[!is.na(metadata$MPid)], metadata$Sample[!is.na(metadata$MPid)])) %>%
  rowwise %>% mutate(pEMT_Pattern = as.factor(samples_metadata$pEMT_Pattern[samples_metadata$Sample == Var2]), .before = "Freq") %>% ungroup %>%
  group_by(Var2) %>% mutate(Prop = Freq / sum(Freq), , OrderBy = Prop[Var1 == "Senescence"]) %>% ungroup %>%
  dplyr::filter(pEMT_Pattern != "none")
p <- per_sample_states %>%
  dplyr::mutate(Group = ifelse(Var1 %in% c("T_cell", "B_cell", "Macrophage", "Fibroblast", "Endothelial", "Complement", "Skeletal_Muscle", "Secretory_Norm"), "Non-Malignant", "Malignant")) %>%
  dplyr::group_by(Var1, Var2) %>%
  summarize(Malignant = sum(Prop[Group == "Malignant"]), `Non-Malignant` = -sum(Prop[Group == "Non-Malignant"])) %>%
  dplyr::group_by(Var2) %>%
  dplyr::mutate(Ord = sum(`Non-Malignant`)) %>% ungroup %>%
  dplyr::mutate(pEMT_Pattern = factor(samples_metadata$pEMT_Pattern[match(Var2, samples_metadata$Sample)], levels = c("Margin_Enriched", "Core_Enriched"))) %>%
  ggplot(aes(x = fct_reorder(Var2, Ord), y = Malignant, fill = Var1)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  geom_col(aes(y = `Non-Malignant`), position = "stack") +
  geom_hline(aes(yintercept = 0)) +
  scale_fill_manual(name = "States", values = state_cols[levels(per_sample_states$Var1)]) + scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0.01)) +
  theme_classic() + theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.9, size = 16), plot.title = element_text(hjust = 0.5, size = 16),
                          legend.text = element_text(size = 14), legend.title = element_text(size = 16)) +
  facet_grid(.~pEMT_Pattern, scales = 'free', space = 'free', drop = F) + guides(fill = guide_legend(override.aes = list(size = 7), reverse = TRUE))
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_5e.pdf"), plot = p, device = "pdf", dpi = 300, width = 12, height = 8)


p <- pEMT_patterns_abund %>%
  dplyr::mutate(Group = ifelse(State %in% c("T_cell", "B_cell", "Macrophage", "Fibroblast", "Endothelial", "Complement"), "Non-Malignant", "Malignant")) %>%
  dplyr::group_by(State, variable) %>%
  summarize(Malignant = sum(value[Group == "Malignant"]), `Non-Malignant` = -sum(value[Group == "Non-Malignant"])) %>%
  dplyr::group_by(variable) %>%
  dplyr::mutate(Ord = sum(`Non-Malignant`)) %>%
  ggplot(aes(x = fct_reorder(variable, Ord), y = Malignant, fill = State)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  geom_col(aes(y = `Non-Malignant`), position = "stack") +
  geom_hline(aes(yintercept = 0)) +
  scale_fill_manual(name = "States", values = state_cols[levels(pEMT_patterns_abund$State)]) + scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.3, 0.2)) +
  theme_classic() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 16), plot.title = element_text(hjust = 0.5, size = 16),
                          legend.text = element_text(size = 14), legend.title = element_text(size = 16), aspect.ratio = 1.4) +
  guides(fill = guide_legend(override.aes = list(size = 7), reverse = TRUE))
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_5c_left.pdf"), plot = p, device = "pdf", dpi = 300)

# Draw statistics
MPs_per_pEMT_pattern_stats <- metadata %>% dplyr::mutate(pEMT_Pattern = as.factor(samples_metadata$pEMT_Pattern[match(Sample, samples_metadata$Sample)])) %>%
  dplyr::filter(!is.na(MPid) & !MPid %in% c("Skeletal_Muscle", "Secretory_Norm", "Secretory_Malig", "LowQ"), pEMT_Pattern %in% c("Core_Enriched", "Margin_Enriched")) %>% droplevels
prop_stat_tab <- test_state_prop_diff(MPs_per_pEMT_pattern_stats$MPid, MPs_per_pEMT_pattern_stats$Sample, MPs_per_pEMT_pattern_stats$pEMT_Pattern, transform = NULL)
# write.csv(prop_stat_tab, file = here("results/Generated_Data/State_Proportion_pEMTcore_vs_pEMTmargin_Extend.csv"))


pEMT_patterns_abund <- read.csv(file = here("results/Generated_Data/Node_Abundance_Neighborhood_Network.csv")) %>%
  dplyr::select(State, pEMTcore, HPVpos_pEMTmargin, HPVneg_pEMTmargin) %>% reshape2::melt(.) %>%
  dplyr::mutate(State = factor(stringr::str_replace_all(State, pattern = "TNF_Signaling", replacement = "Inflammatory"), levels = levels(metadata$MPid)))

p <- pEMT_patterns_abund %>%
  dplyr::mutate(Group = ifelse(State %in% c("T_cell", "B_cell", "Macrophage", "Fibroblast", "Endothelial", "Complement"), "Non-Malignant", "Malignant")) %>%
  dplyr::group_by(State, variable) %>%
  summarize(Malignant = sum(value[Group == "Malignant"]), `Non-Malignant` = -sum(value[Group == "Non-Malignant"])) %>%
  dplyr::group_by(variable) %>%
  dplyr::mutate(Ord = sum(`Non-Malignant`)) %>%
  ggplot(aes(x = fct_reorder(variable, Ord), y = Malignant, fill = State)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  geom_col(aes(y = `Non-Malignant`), position = "stack") +
  geom_hline(aes(yintercept = 0)) +
  scale_fill_manual(name = "States", values = state_cols[levels(pEMT_patterns_abund$State)]) + scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.15, 0.2), labels = c("pEMT margin\nHPV negative", "pEMT margin\nHPV positive", "pEMT core")) +
  theme_classic() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 16), plot.title = element_text(hjust = 0.5, size = 16),
                          legend.text = element_text(size = 14), legend.title = element_text(size = 16), aspect.ratio = 1.4) +
  guides(fill = guide_legend(override.aes = list(size = 7), reverse = TRUE))
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_4f_summary.pdf"), plot = p, device = "pdf", dpi = 300, width = 12, height = 8)

per_sample_states <- as.data.frame(table(metadata$MPid[!is.na(metadata$MPid)], metadata$Sample[!is.na(metadata$MPid)])) %>%
  rowwise %>% mutate(pEMT_Pattern = as.factor(dplyr::case_when(samples_metadata$pEMT_Pattern[samples_metadata$Sample == Var2] == "Core_Enriched" ~ "Core_Enriched",
                                                               samples_metadata$pEMT_Pattern[samples_metadata$Sample == Var2] == "Margin_Enriched" ~ paste0(samples_metadata$pEMT_Pattern[samples_metadata$Sample == Var2], ".", samples_metadata$HPV[samples_metadata$Sample == Var2]),
                                                               samples_metadata$pEMT_Pattern[samples_metadata$Sample == Var2] == "none" ~ NA)), .before = "Freq") %>% ungroup %>%
  dplyr::filter(!is.na(pEMT_Pattern)) %>%
  group_by(Var2) %>% mutate(Prop = Freq / sum(Freq), , OrderBy = Prop[Var1 == "Senescence"]) %>% ungroup
p <- per_sample_states %>%
  dplyr::mutate(Group = ifelse(Var1 %in% c("T_cell", "B_cell", "Macrophage", "Fibroblast", "Endothelial", "Complement", "Skeletal_Muscle", "Secretory_Norm"), "Non-Malignant", "Malignant")) %>%
  dplyr::group_by(Var1, Var2) %>%
  summarize(Malignant = sum(Prop[Group == "Malignant"]), `Non-Malignant` = -sum(Prop[Group == "Non-Malignant"])) %>%
  dplyr::group_by(Var2) %>%
  dplyr::mutate(Ord = sum(`Non-Malignant`)) %>% ungroup %>% droplevels %>% rowwise %>%
  dplyr::mutate(pEMT_Pattern = factor(dplyr::case_when(samples_metadata$pEMT_Pattern[samples_metadata$Sample == Var2] == "Core_Enriched" ~ "Core_Enriched",
                                                       samples_metadata$pEMT_Pattern[samples_metadata$Sample == Var2] == "Margin_Enriched" ~ paste0(samples_metadata$pEMT_Pattern[samples_metadata$Sample == Var2], ".", samples_metadata$HPV[samples_metadata$Sample == Var2]),
                                                       samples_metadata$pEMT_Pattern[samples_metadata$Sample == Var2] == "none" ~ NA),
                                             levels = c("Margin_Enriched.HPVneg", "Margin_Enriched.HPVpos", "Core_Enriched"))) %>%
  ggplot(aes(x = fct_reorder(Var2, Ord), y = Malignant, fill = Var1)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  geom_col(aes(y = `Non-Malignant`), position = "stack") +
  geom_hline(aes(yintercept = 0)) +
  scale_fill_manual(name = "States", values = state_cols[levels(per_sample_states$Var1)]) + scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0.01)) +
  theme_classic() + theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.9, size = 16), plot.title = element_text(hjust = 0.5, size = 16),
                          legend.text = element_text(size = 14), legend.title = element_text(size = 16)) +
  facet_grid(.~pEMT_Pattern, scales = 'free', space = 'free', drop = F, ) + guides(fill = guide_legend(override.aes = list(size = 7), reverse = TRUE))
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_4f_full.pdf"), plot = p, device = "pdf", dpi = 300, width = 12, height = 8)





# Coherence Analysis ------------------------------------------------------

source(here("analysis_scripts/functions/Connectivity_Functions.R"))
# Load merged metadata tibble
metadata <- readRDS(file = here("metadata/merged_metadata_Extend.rds"))


# Classify sample histological growth pattern using coherence at the level of malignancy status
malig_coherence_scores <- calc_spatial_coherence(metadata = metadata, spot_class = "binCNAstatus", zone = "All", site = "All", samples = "All",
                                                 n_cores = 20, n_perm = 100)
# saveRDS(malig_coherence_scores, file = here("results/Generated_Data/all_samples_malignancy_coherence_scores.rds"))

# Plot all sample's malignancy status
malig_coherence_scores <- readRDS(file = here("results/Generated_Data/all_samples_malignancy_coherence_scores.rds"))
samp_levels <- names(sort(do.call(cbind.data.frame, malig_coherence_scores) %>% colMeans(.[!is.na(.)]), decreasing = TRUE))
all_CNAstat_plots <- lapply(seq_along(sample_parameters$Sample), function(i) {
  # print the name of the sample been currently processed
  message(paste0("Processing sample: ", sample_parameters$Sample[[i]]))

  # Load spatial image object
  data_dir <- paste(here("data"), sample_parameters$Site[[i]], sample_parameters$Sample[[i]], sep = "/")
  spatial_image <- load_spatial_image(paste0(data_dir, "/spatial"))

  # Subset metadata to include only the sample you need
  samp_metadata <- metadata[metadata$Sample == sample_parameters$Sample[[i]], ]
  p <- plot_spatial_features(metadata = samp_metadata, image_obj = spatial_image, color_by = "binCNAstatus", cols = setNames(c("#FCF2B4FF", "#FD9969FF", "#51127CFF"), c("Non_Malignant", "Mixed", "Malignant")),
                             pt.size = 1.7, image.alpha = 0, rm_spots = NULL, na.value = "gray", stroke = 0) + theme(legend.position = "none")
})
names(all_CNAstat_plots) <- sample_parameters$Sample
all_CNAstat_plots <- all_CNAstat_plots[samp_levels]
merged_plot <- do.call(plot_grid, all_CNAstat_plots)
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_6b.pdf"), plot = merged_plot, device = "pdf", dpi = 300)



# Run coherence quantification for malignant spots only (so that the genetic and non-genetic components can be compared)
malig_spots_coherence <- calc_malig_spatial_coherence(metadata = metadata, samples = "All", n_cores = 20, n_perm = 100)
# saveRDS(malig_spots_coherence, file = here("results/Generated_Data/all_samples_malignant_spots_MPS_vs_Subclones_coherence_scores_Extend.rds"))
malig_spots_coherence <- readRDS(file = here("results/Generated_Data/all_samples_malignant_spots_MPS_vs_Subclones_coherence_scores_Extend.rds"))

malig_coher_df <- do.call(rbind, malig_spots_coherence) %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::mutate(Sample = factor(Sample, levels = do.call(rbind, malig_spots_coherence) %>% dplyr::arrange(desc(Tumor)) %>% rownames()),
                Subclones_Mean = ifelse(is.na(Subclones_Sd), Tumor, Subclones_Mean))

# Plot boxplot of Tumor bulkyness, Subclone coherence and MP coherence
malig_coher_df <- do.call(rbind, malig_spots_coherence) %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::mutate(Subclones_Mean = ifelse(is.na(Subclones_Sd), Tumor, Subclones_Mean)) %>%
  dplyr::select(Sample, Tumor, MPs_Mean, Subclones_Mean) %>% dplyr::rename(MPs = MPs_Mean, Subclones = Subclones_Mean) %>%
  dplyr::mutate(Color = samples_metadata$Color[match(Sample, samples_metadata$Sample)]) %>% reshape2::melt(.) %>%
  dplyr::mutate(Sample = factor(Sample, levels = unique(Sample)), Color = factor(Color, levels = unique(Color)), variable = factor(variable, levels = c("MPs", "Subclones", "Tumor")))
sample_cols <- setNames(samples_metadata$Color, samples_metadata$Sample)
p <- ggplot(data = malig_coher_df, aes(x = variable, y = value)) +
  geom_boxplot(size = 1) +
  geom_dotplot(aes(fill = Sample), binaxis = "y", stackdir = "center", method="dotdensity", stackgroups = T, binpositions="all", stroke = 2, dotsize = 1) +
  scale_fill_manual(values = sample_cols) + labs(x = "", y = "Mean Coherence Score") +
  theme_classic2() + theme(text = element_text(size = 16), legend.title = element_text(size = 18), axis.title = element_text(size = 18)) +
  guides(fill = guide_legend(override.aes = list(size = 8)))
p$layers[[1]]$aes_params$fill <- "#E7E9EB"
ggplot2::ggsave(filename = here("results/Diagnostic_Plots/Malignant_Spots_Mean_Samples_Subclones_vs_MPs_Coherence_Scores_Boxplot.png"), plot = p)


malig_coher_df <- malig_coher_df %>% dplyr::mutate(pEMT_Pattern = samples_metadata$pEMT_Pattern[match(Sample, samples_metadata$Sample)])
pattern_cols <- setNames(c("#B24745FF", "#DF8F44FF", "#374E55ff"), c("Core_Enriched", "Margin_Enriched", "none"))
p <- malig_coher_df %>% dplyr::filter(variable == "Tumor") %>% dplyr::arrange(desc(value)) %>%
  ggplot(aes(x = factor(Sample, levels = Sample), y = value, group = 1)) +
  geom_line(color = "#0c2a50ff", size = 3) +
  geom_point(aes(fill = pEMT_Pattern), shape = 21, size = 6, stroke = 1.5) +
  theme_classic() + theme(axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust= 1), axis.text.y = element_text(size = 14),
                          axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_text(size = 16)) +
  ylab("Mean Cohernece Score") + geom_hline(yintercept = c(0.65, 0.4), linetype = "dotted") +
  scale_fill_manual(name = "pEMT Pattern", values = pattern_cols)
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_5b.pdf"), plot = p, device = "pdf", dpi = 300)

samples_metadata <- samples_metadata %>% dplyr::mutate(Growth_Pattern = dplyr::case_when(Sample %in% c("P5764", "P5217", "P5456", "P5418", "P5766", "P5684") ~ "Compartmentalized",
                                                                                         Sample %in% c("P5692", "P5219", "P5659", "P5903", "P5624", "P5735") ~ "Intermixed", .default = "Average"))
# saveRDS(samples_metadata, file = here("metadata/samples_metadata.rds"))

# Check if there is a significant association between pEMT pattern and Growth pattern
contig_tab <- table(samples_metadata$Growth_Pattern, samples_metadata$pEMT_Pattern)
chisq.test(contig_tab)
chisq.test(contig_tab, simulate.p.value = TRUE)
chi.test <- function(a, b) {
  return(chisq.test(cbind(a, b)))
}
plot_chi <- as.data.frame(contig_tab)
p <- ggplot(plot_chi, aes(x = Var1, y = Freq, fill = Var2)) +
  geom_col() + scale_fill_manual(name = "pEMT Pattern", values = pattern_cols) + theme_bw() + ylab("Number of Samples") +
  scale_y_continuous(expand = expansion(add = c(0, 0.2)), breaks = seq(2, 12, 2)) + scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +
  theme(axis.text = element_text(size = 14), axis.title.x = element_blank(), axis.title = element_text(size = 16), legend.text = element_text(size = 14), legend.title = element_text(size = 16)) +
  annotate("text", x = 3.24, y = 11.5, label = paste0("Chi-squared\np-value: ", round(chisq.test(contig_tab, simulate.p.value = TRUE)$p.value, digits = 3))) +
  geom_signif(comparisons = list(c("Compartmentalized","Intermixed"), c("Compartmentalized","Intermixed", "Average")), test = "chi.test", y_position = 6)
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_5d.pdf"), plot = p, device = "pdf", dpi = 300)



# Table of metaprograms proportions at each sample
prop_df <- metadata %>% dplyr::filter(!is.na(MPid)) %>%
  dplyr::group_by(Sample, MPid) %>% dplyr::summarise(n = length(MPid)) %>%
  dplyr::group_by(Sample) %>% dplyr::mutate(Prop = n / sum(n)) %>%
  dplyr::select(-n) %>%
  tidyr::pivot_wider(names_from = MPid, values_from = Prop) %>% ungroup %>%
  dplyr::mutate(Sample = factor(Sample, levels = samp_levels))
prop_mat <- prop_df %>% tibble::column_to_rownames("Sample") %>% as.matrix
malig_coher_vec <- malig_coher_df %>% dplyr::filter(variable == "Tumor") %>% dplyr::select(Sample, value) %>%
  dplyr::left_join(., samples_metadata[, c("Sample", "Growth_Pattern")], by = "Sample")

per_sample_states <- as.data.frame(table(metadata$MPid[!is.na(metadata$MPid)], metadata$Sample[!is.na(metadata$MPid)])) %>%
  rowwise %>% mutate(GP = as.factor(samples_metadata$Growth_Pattern[samples_metadata$Sample == Var2]), .before = "Freq") %>% ungroup %>%
  group_by(Var2) %>% mutate(Prop = Freq / sum(Freq), , OrderBy = Prop[Var1 == "Senescence"]) %>% ungroup %>% dplyr::filter(GP != "Average") %>% droplevels()
p <- per_sample_states %>%
  dplyr::mutate(Group = ifelse(Var1 %in% c("T_cell", "B_cell", "Macrophage", "Fibroblast", "Endothelial", "Complement", "Skeletal_Muscle", "Scretory_Norm"), "Non-Malignant", "Malignant")) %>%
  dplyr::group_by(Var1, Var2) %>%
  summarize(Malignant = sum(Prop[Group == "Malignant"]), `Non-Malignant` = -sum(Prop[Group == "Non-Malignant"])) %>%
  dplyr::group_by(Var2) %>%
  dplyr::mutate(Ord = sum(`Non-Malignant`)) %>% ungroup %>%
  dplyr::mutate(GP = factor(samples_metadata$Growth_Pattern[match(Var2, samples_metadata$Sample)], levels = c("Compartmentalized", "Intermixed"))) %>%
  ggplot(aes(x = fct_reorder(Var2, Ord), y = Malignant, fill = Var1)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  geom_col(aes(y = `Non-Malignant`), position = "stack") +
  geom_hline(aes(yintercept = 0)) +
  scale_fill_manual(name = "States", values = state_cols[levels(per_sample_states$Var1)]) + scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  theme_classic() + theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.9, size = 16), plot.title = element_text(hjust = 0.5, size = 16),
                          legend.text = element_text(size = 14), legend.title = element_text(size = 16)) +
  facet_grid(.~GP, scales = 'free', space = 'free', drop = F) + guides(fill = guide_legend(override.aes = list(size = 7), reverse = TRUE))
# ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_5e.pdf"), plot = p, device = "pdf", dpi = 300, width = 10, height = 6)

per_GP_states <- metadata %>% dplyr::mutate(GP = as.factor(samples_metadata$Growth_Pattern[match(Sample, samples_metadata$Sample)])) %>%
  dplyr::select(MPid, GP) %>% table() %>% as.data.frame() %>%
  group_by(GP) %>% mutate(Prop = Freq / sum(Freq)) %>% ungroup %>% dplyr::filter(GP != "Average") %>% droplevels()
p <- per_GP_states %>%
  dplyr::mutate(Group = ifelse(MPid %in% c("T_cell", "B_cell", "Macrophage", "Fibroblast", "Endothelial", "Complement", "Skeletal_Muscle", "Scretory_Norm"), "Non-Malignant", "Malignant")) %>%
  dplyr::group_by(MPid, GP) %>%
  summarize(Malignant = sum(Prop[Group == "Malignant"]), `Non-Malignant` = -sum(Prop[Group == "Non-Malignant"])) %>%
  dplyr::group_by(GP) %>%
  dplyr::mutate(Ord = sum(`Non-Malignant`)) %>%
  ggplot(aes(x = fct_reorder(GP, Ord), y = Malignant, fill = MPid)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  geom_col(aes(y = `Non-Malignant`), position = "stack") +
  geom_hline(aes(yintercept = 0)) +
  scale_fill_manual(name = "States", values = state_cols[levels(per_GP_states$MPid)]) + scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  theme_classic() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 16), plot.title = element_text(hjust = 0.5, size = 16),
                          legend.text = element_text(size = 14), legend.title = element_text(size = 16)) +
  guides(fill = guide_legend(override.aes = list(size = 7), reverse = TRUE))
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_5c_right.pdf"), plot = p, device = "pdf", dpi = 300, width = 7, height = 6)

# Draw statistics
states_per_GP_stats <- metadata %>%
  dplyr::mutate(GP = as.factor(samples_metadata$Growth_Pattern[match(Sample, samples_metadata$Sample)])) %>%
  dplyr::filter(GP != "Average") %>% droplevels()
prop_stat_tab <- test_state_prop_diff(states_per_GP_stats$MPid, states_per_GP_stats$Sample, states_per_GP_stats$GP, transform = NULL)
write.csv(prop_stat_tab, file = here("results/Generated_Data/State_Proportion_Per_Growth_Pattern_Ttest_Results.csv"))

