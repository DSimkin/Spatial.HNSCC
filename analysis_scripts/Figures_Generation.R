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
# Load single-cell cell-type colors
sc_cellType_cols <- readRDS(file = here("aux_data/sc_cellType_cols.rds"))
sc_cellType_cols <- setNames(sc_cellType_cols$V2, sc_cellType_cols$V1)
# Load state colors
state_cols <- read.delim(file = here("aux_data/state_col_tab_extend.tsv"), sep = "\t", header = TRUE)
state_cols <- setNames(state_cols$V2, state_cols$V1)





# Plot CNA Scores per MP as Boxplot ---------------------------------------
cna_score <- data.frame(Sample = factor(metadata$Sample),
                        CNAscore = metadata$CNAscore,
                        Percent_Epi = metadata$Percent_Epi,
                        MPid = metadata$MPid,
                        Zone = metadata$Zone)

# Scale CNA scores from all samples in order to be able to compare all
cna_score <- cna_score %>% dplyr::group_by(Sample) %>% dplyr::mutate(CNAscore = as.numeric(scale(CNAscore, center = TRUE, scale = FALSE))) %>% dplyr::ungroup() %>%
  dplyr::filter(!is.na(MPid), !is.na(CNAscore))

p0 <- ggplot(cna_score[!is.na(cna_score$MPid), ], aes(x = reorder(MPid, CNAscore, FUN = median), y = CNAscore, fill = MPid)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE, show.legend = FALSE, coef = 0) + scale_fill_manual(name = "Metaprograms", values = state_cols) +
  scalop::theme_scalop() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust= 0.9), axis.title.x = element_blank()) + ylab("Normalized CNA Score")
ylim1 <- c(quantile(cna_score$CNAscore[!is.na(cna_score$MPid) & cna_score$MPid == "Complement"], probs = 0.25),
           quantile(cna_score$CNAscore[!is.na(cna_score$MPid) & cna_score$MPid == "Epithelial"], probs = 0.75))
p <- p0 + coord_cartesian(ylim = ylim1 * 1.05)

# same plot colored by malignant non-malignant and macrophage/inflammatory
state_cols2 <- setNames(c(rep("#FCF2B4FF", 7), "pink2", "tan3", rep("#51127CFF", 7)), levels(reorder(cna_score$MPid, cna_score$CNAscore, FUN = median)))
cna_score <- cna_score %>% dplyr::filter(MPid != "LowQ")
p0 <- ggplot(cna_score[!is.na(cna_score$MPid), ], aes(x = reorder(MPid, CNAscore, FUN = median), y = CNAscore, fill = MPid)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE, show.legend = FALSE, coef = 0) + scale_fill_manual(name = "Metaprograms", values = state_cols2) +
  scalop::theme_scalop() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust= 0.9), axis.title.x = element_blank()) + ylab("Normalized CNA Score")
p <- p0 + coord_cartesian(ylim = ylim1 * 1.05)
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_2c_top.pdf"), device = "pdf", plot = p, dpi = 300, width = 9, height = 6)


# Plot CNA-scores versus Zone
p <- ggpubr::ggboxplot(data = cna_score[!is.na(cna_score$MPid) & !is.na(cna_score$Zone) & cna_score$Zone != "Filtered_out", ],
                       x = "Zone", y = "CNAscore", fill = "Zone", palette = "jama", bxp.errorbar = TRUE, bxp.errorbar.width = 0.2, outlier.shape = NA) +
  ylab("CNA\nScore") + theme(text = element_text(size = 14), axis.text = element_text(size = 14), axis.title = element_text(size = 16),
                             axis.title.x = element_blank(), axis.title.y = element_text(angle = 0, vjust = 0.5), legend.position = "none") +
  scale_x_discrete(labels = function(x) gsub("_", " ", x, fixed = TRUE)) + scale_y_continuous(limits = c(-1, 0.6))
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_4a_middle_1.pdf"), device = "pdf", plot = p, dpi = 300, width = 6, height = 4)





# INF Spots vs Macrophage Spots Colocalization and Malignant Frequencies ----------------------
prop_tab <- as.data.frame(table(metadata$MPid[!is.na(metadata$MPid)])) %>% dplyr::mutate(Sum = sum(Freq), Prop = Freq/Sum)
malig_states <- c("Senescence", "Cell_Cycle", "Epithelial", "Hypoxia", "pEMT", "LowQ", "Secretory_Malig")
non_malig_types <- c("T_cell", "B_cell", "Macrophage", "Fibroblast", "Endothelial", "Complement", "Skeletal_Muscle", "Secretory_Norm")
col_pal <- c(state_cols, setNames(c("#FCF2B4FF", "#51127CFF", scales::alpha("grey", alpha = 0.4)), c("Non-Malignant", "Malignant", "Unknown")))

coloc_states_ls <- list()
malig_fracs_ls <- list()
pie_charts_ls <- list()
for(state in names(MPs)) {

  # Calculate percentage of cell types/states that colocalize with a specific cell type/state
  state_coloc <- metadata %>% dplyr::filter(MPid == state) %>% dplyr::select(Fibroblast:Secretory_Norm) %>% colSums %>% as.data.frame %>%
    tibble::rownames_to_column() %>% `colnames<-`(c("Var1", "Coloc")) %>%
    dplyr::left_join(., prop_tab[, c("Var1", "Prop")]) %>%
    dplyr::mutate(Group = dplyr::case_when(Var1 %in% malig_states ~ "Malignant",
                                           Var1 %in% non_malig_types ~ "Non-Malignant",
                                           Var1 == "Inflammatory" ~ "Unknown"), .before = "Coloc") %>%
    dplyr::mutate(Group = factor(Group, levels = c(Group[Var1 == state], sort(unique(Group[!Group %in% Group[Var1 == state]])))), Coloc = Coloc / (1 / -log2(Prop))) %>%
    dplyr::mutate(frac = (Coloc / sum(Coloc)) * 100) %>%
    dplyr::group_by(Group) %>% dplyr::arrange(desc(frac), .by_group = TRUE) %>% ungroup() %>%
    dplyr::mutate(ymax = cumsum(frac), lab_ypos = cumsum(frac) - 0.5 * frac) %>%
    dplyr::mutate(ymin = c(0, head(ymax, n = -1)), text_size = dplyr::case_when(frac < 1 ~ 0,
                                                                                frac >= 10 ~ 10,
                                                                                frac > 1 & frac < 10 ~ round(frac)))
  coloc_states_ls[[which(names(MPs) == state)]] <- state_coloc
  names(coloc_states_ls)[which(names(MPs) == state)] <- state

  # Calculate fraction of malignant cell states vs non-malognant cell types for cell that colocalize with a specific cell type/state
  state_malig_frac <- state_coloc %>% dplyr::group_by(Group) %>% dplyr::reframe(Malig_frac = sum(frac)) %>%
    dplyr::mutate(ymax = cumsum(Malig_frac), lab_ypos = cumsum(Malig_frac) - 0.5 * Malig_frac) %>% dplyr::mutate(ymin = c(0, head(ymax, n = -1)))
  malig_fracs_ls[[which(names(MPs) == state)]] <- state_malig_frac
  names(malig_fracs_ls)[which(names(MPs) == state)] <- state

  # Plot Pie-chart
  p <- ggplot() +
    geom_rect(data = state_coloc, aes(fill = factor(Var1, levels = rev(Var1)), ymax = ymax, ymin = ymin, xmax = 6, xmin = 0), color = "black") +
    geom_rect(data = state_malig_frac, aes(fill = Group, ymax = ymax, ymin = ymin, xmax = 7, xmin = 6), color = "black") +
    coord_polar(theta = "y", start = 0) +
    geom_text(data = state_coloc, aes(x = 3.5, y = lab_ypos, label = paste0(round(frac), "%")), color = "black", size = state_coloc$text_size) +
    geomtextpath::geom_textpath(data = state_malig_frac, aes(x = 6.5, y = lab_ypos, label = paste0(round(Malig_frac), "%")), angle = 90, size = 5, color = "#7D7D7D", fontface = "bold") +
    scale_fill_manual(values = col_pal, name = "State") + theme_void() +
    theme(legend.text = element_text(size = 16), legend.title = element_text(size = 18), plot.title = element_text(size = 8, hjust = 0.5, face = "bold", margin = margin(5, 0, -10, 0))) +
    guides(fill = guide_legend(override.aes = list(size = 10), reverse = TRUE)) + ggtitle(state)
  pie_charts_ls[[which(names(MPs) == state)]] <- p
  names(pie_charts_ls)[which(names(MPs) == state)] <- state
}

orders_MPs <- c("Skeletal_Muscle", "Secretory_Norm", "Endothelial", "Fibroblast", "Complement", "B_cell", "T_cell", "Macrophage",
                "Inflammatory", "pEMT", "Hypoxia", "Cell_Cycle", "Epithelial", "Secretory_Malig", "Senescence", "LowQ")
pie_charts_ls_plot <- lapply(pie_charts_ls, function(x) {
  plt <- x + theme(legend.position = "none") + theme(plot.margin = margin(t = 0, r = 0, b = -10, l = 0))
  plt$layers[[3]]$aes_params$size <- plt$layers[[3]]$aes_params$size / 2.8
  plt$layers[[4]]$aes_params$size <- 3
  return(plt)
})
pie_charts_ls_plot <- pie_charts_ls_plot[orders_MPs]
p <- do.call(ggarrange, pie_charts_ls_plot)
ggplot2::ggsave(filename = here("results/Diagnostic_Plots/PerState_Colocalization_PieChart.png"), plot = p, width = 14)


# Same plot with seperation to Malignant, Non-malignant Epithelial, Stromal and Immune compartments
prop_tab <- as.data.frame(table(metadata$MPid[!is.na(metadata$MPid)])) %>% dplyr::mutate(Sum = sum(Freq), Prop = Freq/Sum)
malig_states <- c("Senescence", "Cell_Cycle", "Epithelial", "Hypoxia", "pEMT", "LowQ", "Secretory_Malig")
norm_epi <- "Secretory_Norm"
stromal_types <- c("Fibroblast", "Endothelial", "Complement", "Skeletal_Muscle")
immune_types <- c("T_cell", "B_cell", "Macrophage", "Inflammatory")
col_pal <- c(state_cols, setNames(c("#FCF2B4FF", "#51127CFF", "#FC8E64FF", "#C73D73FF"), c("Non-Malignant_Epithelial", "Malignant", "Immune", "Stroma")))

coloc_states_ls <- list()
malig_fracs_ls <- list()
pie_charts_ls <- list()
for(state in names(MPs)) {

  # Calculate percentage of cell types/states that colocalize with a specific cell type/state
  state_coloc <- metadata %>% dplyr::filter(MPid == state) %>% dplyr::select(Fibroblast:Secretory_Norm) %>% colSums %>% as.data.frame %>%
    tibble::rownames_to_column() %>% `colnames<-`(c("Var1", "Coloc")) %>%
    dplyr::left_join(., prop_tab[, c("Var1", "Prop")]) %>%
    dplyr::mutate(Group = dplyr::case_when(Var1 %in% malig_states ~ "Malignant",
                                           Var1 == norm_epi ~ "Non-Malignant_Epithelial",
                                           Var1 %in% stromal_types ~ "Stroma",
                                           Var1 %in% immune_types ~ "Immune"), .before = "Coloc") %>%
    dplyr::mutate(Coloc = Coloc / (1 / -log2(Prop)), frac = (Coloc / sum(Coloc)) * 100) %>%
    dplyr::filter(frac >= 1) %>% dplyr::mutate(frac = frac * (100 / sum(frac)))
  order_groups <- state_coloc %>% dplyr::group_by(Group) %>% dplyr::reframe(order_fracs = sum(frac)) %>% dplyr::arrange(desc(order_fracs)) %>% dplyr::pull(Group)
  state_coloc <- state_coloc %>% dplyr::mutate(Group = factor(Group, levels = order_groups)) %>%
    dplyr::group_by(Group) %>% dplyr::arrange(desc(frac), .by_group = TRUE) %>% ungroup() %>%
    dplyr::mutate(ymax = cumsum(frac), lab_ypos = cumsum(frac) - 0.5 * frac) %>%
    dplyr::mutate(ymin = c(0, head(ymax, n = -1)), text_size = dplyr::case_when(frac < 1 ~ 0,
                                                                                frac >= 10 ~ 10,
                                                                                frac > 1 & frac < 10 ~ round(frac)))

  coloc_states_ls[[which(names(MPs) == state)]] <- state_coloc
  names(coloc_states_ls)[which(names(MPs) == state)] <- state

  # Calculate fraction of malignant cell states vs non-malognant cell types for cell that colocalize with a specific cell type/state
  state_malig_frac <- state_coloc %>% dplyr::group_by(Group) %>% dplyr::reframe(Malig_frac = sum(frac)) %>%
    dplyr::arrange(desc(Malig_frac)) %>% dplyr::mutate(Group = factor(Group, levels = Group), MPid = state)
  if(!all(c("Malignant", "Immune", "Stroma", "Non-Malignant_Epithelial") %in% levels(state_malig_frac$Group))) {
    state_malig_frac <- state_malig_frac %>% dplyr::add_row(Group = setdiff(c("Malignant", "Immune", "Stroma", "Non-Malignant_Epithelial"), levels(state_malig_frac$Group)),
                                                            Malig_frac = 0, MPid = state)
  }
  malig_fracs_ls[[which(names(MPs) == state)]] <- state_malig_frac
  names(malig_fracs_ls)[which(names(MPs) == state)] <- state
  state_malig_frac <- state_malig_frac %>% dplyr::mutate(ymax = cumsum(Malig_frac)) %>% dplyr::mutate(ymin = c(0, head(ymax, n = -1))) %>%
    dplyr::mutate(lab_ypos = ifelse(Malig_frac < 3, NA, cumsum(Malig_frac) - 0.5 * Malig_frac))

  # Plot Pie-chart
  p <- ggplot() +
    geom_rect(data = state_coloc, aes(fill = factor(Var1, levels = rev(Var1)), ymax = ymax, ymin = ymin, xmax = 6, xmin = 0), color = "black") +
    geom_rect(data = state_malig_frac, aes(fill = Group, ymax = ymax, ymin = ymin, xmax = 7, xmin = 6), color = "black") +
    coord_polar(theta = "y", start = 0) +
    geom_text(data = state_coloc, aes(x = 3.5, y = lab_ypos, label = paste0(round(frac), "%")), color = "black", size = state_coloc$text_size) +
    geomtextpath::geom_textpath(data = state_malig_frac, aes(x = 6.5, y = lab_ypos, label = paste0(round(Malig_frac), "%")), angle = 90, size = 5, color = "#7D7D7D", fontface = "bold") +
    scale_fill_manual(values = col_pal, name = "State") + theme_void() +
    theme(legend.text = element_text(size = 16), legend.title = element_text(size = 18), plot.title = element_text(size = 8, hjust = 0.5, face = "bold", margin = margin(5, 0, -10, 0))) +
    guides(fill = guide_legend(override.aes = list(size = 10), reverse = TRUE)) + ggtitle(state)
  pie_charts_ls[[which(names(MPs) == state)]] <- p
  names(pie_charts_ls)[which(names(MPs) == state)] <- state
}
p <- pie_charts_ls$Inflammatory + theme(legend.position = "none") + ggtitle("")
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_2e_right.pdf"), device = "pdf", plot = p, dpi = 300)
p <- pie_charts_ls$Macrophage + theme(legend.position = "none") + ggtitle("")
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_2e_left.pdf"), device = "pdf", plot = p, dpi = 300)
lej1 <- cowplot::get_legend(pie_charts_ls$Epithelial)
lej2 <- cowplot::get_legend(pie_charts_ls$Skeletal_Muscle)
p <- plot_grid(lej1, lej2, ncol = 2)
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_2e_legends.pdf"), device = "pdf", plot = p, dpi = 300)


orders_MPs <- c("Secretory_Norm", "Skeletal_Muscle", "Endothelial", "Fibroblast", "Complement", "T_cell", "B_cell", "Macrophage",
                "Inflammatory", "pEMT", "Hypoxia", "Cell_Cycle", "Epithelial", "Secretory_Malig", "Senescence", "LowQ")
pie_charts_ls_plot <- lapply(pie_charts_ls, function(x) {
  plt <- x + theme(legend.position = "none") + theme(plot.margin = margin(t = 0, r = 0, b = -10, l = 0))
  plt$layers[[3]]$aes_params$size <- plt$layers[[3]]$aes_params$size / 2.8
  plt$layers[[4]]$aes_params$size <- 3
  return(plt)
})
pie_charts_ls_plot <- pie_charts_ls_plot[orders_MPs]
p <- do.call(ggarrange, pie_charts_ls_plot)
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_3a.pdf"), device = "pdf", plot = p, dpi = 300, width = 14)



# Plot heatmap of the different fractions of colocalizing stromal/immune/malignant of each cell-type/state
orders_MPs <- metadata %>% dplyr::filter(!is.na(MPid), !is.na(CNAscore)) %>%
  dplyr::group_by(Sample) %>% dplyr::mutate(CNAscore = as.numeric(scale(CNAscore, center = TRUE, scale = FALSE))) %>% dplyr::ungroup() %>%
  dplyr::reframe(MP_order = reorder(MPid, CNAscore, FUN = median)) %>% pull(MP_order) %>% levels() %>% .[-grep("LowQ", .)]
hm_df <- do.call(rbind.data.frame, malig_fracs_ls) %>%
  dplyr::mutate(MPid = factor(MPid, levels = orders_MPs), Group = str_replace_all(Group, pattern = "Non-Malignant_Epithelial", replacement = "Non-Malignant\nEpithelial")) %>%
  dplyr::mutate(Group = factor(Group, levels = rev(c("Malignant", "Immune", "Stroma", "Non-Malignant\nEpithelial")))) %>% dplyr::filter(MPid != "LowQ")
p <- gmap(hm_df, x = MPid, y = Group, fill = Malig_frac, limits = c(0, 100), midpoint = 50, ratio = 0.2, y.name = "",
     axis.rel = 1.5, legend.title.rel = 1.25, legend.rel = 1.0, legend.height = 5, legend.width = 0.6, legend.breaks = c(0, 50, 100), x.num = F, y.num = F, magma_pal = TRUE) +
  theme(axis.text = element_text(size = 16), axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1))
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_3c_bottom.pdf"), device = "pdf", plot = p, dpi = 300, width = 12, height = 6)




# Matched Single-cell OSM, OSMR, SPP1 and IL1B Expression in Macrophage, Neutrophil and pEMT ------------------

# Load matched merged scRNA-seq expression matrix
sc_merged_mat <- readRDS(file = here("aux_data/merged_matched_scRNAseq_expr_mat.rds"))
# Load single cell HTAN metadata and process it
sc_metadata <- readRDS(file = here("metadata/merged_scRNAseq_metadata.rds"))

# Extract 3 groups of cells - pEMT, IS-Macrophage and Neutrophil
group_meta <- sc_metadata %>% dplyr::filter(grepl("pEMT", Label) | Label %in% c("Neutrophil", "Immunosuppr.-Macrophage")) %>%
  dplyr::mutate(Label = stringr::str_replace_all(Label, pattern = "Cancer_pEMT_Hypoxic", replacement = "Cancer_pEMT"))
filt_expr_mat <- as.matrix(sc_merged_mat[c("OSM", "OSMR", "IL1B", "SPP1"), group_meta$Key])
plot_df <- reshape2::melt(filt_expr_mat) %>% dplyr::mutate(MPid = factor(group_meta$Label[match(Var2, group_meta$Key)], levels = c("Cancer_pEMT", "Immunosuppr.-Macrophage", "Neutrophil")))
p <- ggplot(data = plot_df, aes(MPid, value, fill = MPid)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE, linewidth = 0.8) +
  scale_y_continuous(expand = c(0, 0), position = "right") +
  facet_grid(rows = vars(Var1), scales = "free", switch = "y") +
  scalop::theme_scalop() +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"), panel.background = element_rect(fill = NA, color = "black", linewidth = 1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), strip.background = element_blank(), strip.text = element_text(face = "bold"), strip.text.y.left = element_text(angle = 0),
        axis.ticks.x = element_blank(), aspect.ratio = 0.3) +
  scale_fill_manual(values = setNames(c("#E4CAB5", "#FC9272", "#DDCC77"), c("Cancer_pEMT", "Immunosuppr.-Macrophage", "Neutrophil"))) +
  scale_x_discrete(expand = c(0.2, 0.1), labels = function(x) gsub("_", " ", x, fixed = TRUE)) + xlab("") + ylab(expression("Log"[2]*" Expression Level"))
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_5e.pdf"), plot = p, device = "pdf", dpi = 300, width = 8, height = 10)

# Extract 4 groups of cells - pEMT, IS-Macrophage, Neutrophil and Fibroblast (to also check TGFB)
group_meta <- sc_metadata %>% dplyr::filter(grepl("pEMT", Label) | Label %in% c("Neutrophil", "Immunosuppr.-Macrophage", "Fibroblast")) %>%
  dplyr::mutate(Label = stringr::str_replace_all(Label, pattern = "Cancer_pEMT_Hypoxic", replacement = "Cancer_pEMT"))
filt_expr_mat <- as.matrix(sc_merged_mat[c("OSM", "OSMR", "IL1B", "SPP1", "TGFB1", "TGFBR1", "TGFB2", "TGFBR2", "TGFB3", "TGFBR3"), group_meta$Key])
plot_df <- reshape2::melt(filt_expr_mat) %>% dplyr::mutate(MPid = factor(group_meta$Label[match(Var2, group_meta$Key)], levels = c("Cancer_pEMT", "Immunosuppr.-Macrophage", "Neutrophil", "Fibroblast")))
p <- ggplot(data = plot_df, aes(MPid, value, fill = MPid)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE, linewidth = 0.8) +
  scale_y_continuous(expand = c(0, 0), position = "right") +
  facet_grid(rows = vars(Var1), scales = "free", switch = "y") +
  scalop::theme_scalop() +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"), panel.background = element_rect(fill = NA, color = "black", linewidth = 1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), strip.background = element_blank(), strip.text = element_text(face = "bold"), strip.text.y.left = element_text(angle = 0),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values = setNames(c("#E4CAB5", "#FC9272", "#DDCC77", "#88CCEE"), c("Cancer_pEMT", "Immunosuppr.-Macrophage", "Neutrophil", "Fibroblast"))) +
  scale_x_discrete(expand = c(0.15, 0.05)) + xlab("") + ylab(expression("Log"[2]*" Expression Level"))
ggplot2::ggsave(filename = here("results/Diagnostic_Plots/scRNAseq_Ligand_Receptor_Expression_ViolinPlot_with_TGFB.png"), plot = p, width = 9.5, height = 14)

# Test LR expression in all cell types
filt_expr_mat <- as.matrix(sc_merged_mat[c("OSM", "OSMR", "IL1B", "SPP1", "TGFB1", "TGFBR1", "TGFB2", "TGFBR2", "TGFB3", "TGFBR3"), colnames(sc_merged_mat) %in% sc_metadata$Key])
plot_df <- reshape2::melt(filt_expr_mat) %>% dplyr::mutate(MPid = sc_metadata$Label[match(Var2, sc_metadata$Key)])
p <- ggplot(data = plot_df, aes(MPid, value, fill = MPid)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE, linewidth = 0.8) +
  scale_y_continuous(expand = c(0, 0), position = "right") +
  facet_grid(rows = vars(Var1), scales = "free", switch = "y") +
  scalop::theme_scalop() +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"), panel.background = element_rect(fill = NA, color = "black", linewidth = 1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), strip.background = element_blank(), strip.text = element_text(face = "bold"), strip.text.y.left = element_text(angle = 0),
        axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.9)) +
  scale_x_discrete(expand = c(0., 0.0)) + xlab("") + ylab(expression("Log"[2]*" Expression Level"))
sapply(unique(sc_metadata$Label[grepl("Cancer_", sc_metadata$Label)]), function(x) rowMeans(filt_expr_mat[, sc_metadata$Key[sc_metadata$Label == x]]))






# HPV Cellular Density  ---------------------------------------------------

metadata <- metadata %>% mutate(Site = samples_metadata$Site[match(metadata$Sample, samples_metadata$Sample)],
                                HPV = samples_metadata$HPV[match(metadata$Sample, samples_metadata$Sample)])

# Plot the difference in cellular density in HPV-associated vs HPV-negative HNSCC tumors
p <- ggplot(data = metadata) +
  ggridges::geom_density_ridges(mapping = aes(x = Cell_Count,
                                              y = HPV,
                                              fill = HPV),
                                color = "black", stat = "binline", bins = max(metadata$Cell_Count) + 1, alpha = 0.8) +
  ggridges::geom_density_ridges(mapping = aes(x = Cell_Count,
                                              y = HPV,
                                              fill = HPV),
                                stat = "density_ridges", quantile_lines = TRUE, quantile_fun = median, alpha = 0, size = 0.75, linetype = "dashed", vline_color = "grey50", color = NA) +
  scale_fill_manual(values = setNames(c("#FD8D3C", "#74C476"), c("HPVneg", "HPVpos"))) +
  ggridges::theme_ridges() + xlab("Numer of Cells per Spot") + ylab("") + scale_y_discrete(expand = expand_scale(mult = c(0.01, 0.4))) +
  theme(text = element_text(size = 14), axis.text = element_text(size = 14), axis.title = element_text(size = 16), axis.title.x = element_text(hjust = 0.5), legend.position = "none")
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_3e.pdf"), plot = p, device = "pdf", dpi = 300)


# Difference between HPVneg, HPVneg-OPSCC and HPVpos-OPSCC
metadata <- metadata %>%
  dplyr::mutate(Label = ifelse(Site != "Oropharynx", "HPVneg",
                               ifelse(Site == "Oropharynx" & HPV == "HPVpos", "HPVpos OPSCC", "HPVneg OPSCC")))
p <- ggplot(data = metadata %>% dplyr::rename("Cell Count" = "Cell_Count"), aes(x = Label, y = `Cell Count`)) +
  geom_violin(aes(fill = Label), draw_quantiles = c(0.25, 0.75), linetype = "dashed") + geom_violin(fill = "transparent", draw_quantiles = 0.5, linewidth = 1) +
  theme_classic() + theme(axis.ticks.x = element_blank(), axis.text = element_text(size = 14), axis.title.x = element_blank(), axis.title.y = element_text(size = 16), legend.position = "none") +
  scale_y_continuous(expand = c(0.005, 0.005)) + scale_fill_manual(values = setNames(c("#6BAED6", "#FD8D3C", "#74C476"), c("HPVneg OPSCC", "HPVneg", "HPVpos OPSCC")))
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_3f.pdf"), plot = p, device = "pdf", dpi = 300)

# Difference between subsite of HPVneg tumors
metadata <- metadata %>%
  dplyr::mutate(Label = dplyr::case_when(Site == "Oral" ~ "Oral Cavity",
                                         Site == "Laryngeal" ~ "Laryngeal",
                                         Site == "Oropharynx" & HPV == "HPVneg" ~ "HPVneg OPSCC", .default = NA))
p <- ggplot(data = metadata %>% dplyr::rename("Cell Count" = "Cell_Count") %>% dplyr::filter(!is.na(Label)), aes(x = Label, y = `Cell Count`)) +
  geom_violin(aes(fill = Label), draw_quantiles = c(0.25, 0.75), linetype = "dashed") + geom_violin(fill = "transparent", draw_quantiles = 0.5, linewidth = 1) +
  theme_classic() + theme(axis.ticks.x = element_blank(), axis.text = element_text(size = 14), axis.title.x = element_blank(), axis.title.y = element_text(size = 16), legend.position = "none") +
  scale_y_continuous(expand = c(0.005, 0.005)) + scale_fill_manual(values = setNames(c("#756BB1", "#E31A1C", "#238B45"), c("HPVneg OPSCC", "Oral Cavity", "Laryngeal")))
ggplot2::ggsave(filename = here("results/Diagnostic_Plots/HPVneg_Subsite_Cellular_Density_ViolinPlot.png"), plot = p)


# By HPV status splitted to each component (Epithelial, Stromal & Mixed)
p_list <- map(unique(metadata$EpiStroma[!metadata$EpiStroma %in% "Filtered_Out"]),
              ~ ggplot(data = metadata[metadata$EpiStroma == ., ], aes(x = HPV, y = Cell_Count)) +
                geom_violin(aes(fill = HPV), draw_quantiles = c(0.25, 0.75), linetype = "dashed") + geom_violin(fill = "transparent", draw_quantiles = 0.5, linewidth = 1) +
                theme_classic() + theme(axis.ticks.x = element_blank(), axis.text = element_text(size = 14), axis.title.x = element_blank(), axis.title.y = element_text(size = 16), legend.position = "none") +
                scale_fill_manual(values = setNames(c("#FD8D3C", "#74C476"), c("HPVneg", "HPVpos"))) + scale_y_continuous(expand = c(0.005, 0.005), breaks = c(20, 40, 60), limits = c(0, 60)) +
                ylab("Numer of Cells\nper Spot") + ggtitle(paste0(., " Spots")))
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_3g_stroma.pdf"), plot = p_list[[1]], device = "pdf", dpi = 300, width = 6, height = 4.5)
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_3g_epithelial.pdf"), plot = p_list[[2]], device = "pdf", dpi = 300, width = 6, height = 4.5)
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_3g_mixed.pdf"), plot = p_list[[3]], device = "pdf", dpi = 300, width = 6, height = 4.5)



# HPV Status Affect MP Abundance ------------------------------------------

hpv_status_states_prop <- as.data.frame(table(droplevels(metadata$MPid[!is.na(metadata$MPid) & metadata$MPid != "Skeletal_Muscle"]),
                                              metadata$Sample[!is.na(metadata$MPid) & metadata$MPid != "Skeletal_Muscle"])) %>%
  dplyr::mutate(HPVstat = as.factor(samples_metadata$HPV[match(Var2, samples_metadata$Sample)]), .before = "Freq") %>%
  group_by(Var2) %>% mutate(Prop = Freq / sum(Freq), , OrderBy = Prop[Var1 == "Senescence"]) %>% ungroup

p <- hpv_status_states_prop %>%
  dplyr::mutate(Group = ifelse(Var1 %in% c("T_cell", "B_cell", "Macrophage", "Fibroblast", "Endothelial", "Complement", "Secretory_Norm"), "Non-Malignant", "Malignant")) %>%
  dplyr::group_by(Var1, Var2) %>%
  summarize(Malignant = sum(Prop[Group == "Malignant"]), `Non-Malignant` = -sum(Prop[Group == "Non-Malignant"]), HPVstat = HPVstat) %>%
  dplyr::group_by(Var2) %>%
  dplyr::mutate(Ord = sum(`Non-Malignant`)) %>% ungroup %>%
  ggplot(aes(x = fct_reorder(Var2, Ord), y = Malignant, fill = Var1)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  geom_col(aes(y = `Non-Malignant`), position = "stack") +
  geom_hline(aes(yintercept = 0)) +
  scale_fill_manual(name = "States", values = state_cols[levels(hpv_status_states_prop$Var1)]) + scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0.01)) +
  theme_classic() + theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.9, size = 16), plot.title = element_text(hjust = 0.5, size = 16),
                          legend.text = element_text(size = 14), legend.title = element_text(size = 16)) +
  facet_grid(.~HPVstat, scales = 'free', space = 'free', drop = F) + guides(fill = guide_legend(override.aes = list(size = 7), reverse = TRUE))
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_4b.pdf"), plot = p, device = "pdf", dpi = 300)


sum_hpv_status_MPs_prop <- metadata %>% dplyr::mutate(HPVstat = as.factor(samples_metadata$HPV[match(Sample, samples_metadata$Sample)])) %>%
  dplyr::filter(!is.na(MPid) & MPid != "Skeletal_Muscle") %>% droplevels %>%
  dplyr::select(MPid, HPVstat) %>% table() %>% as.data.frame() %>%
  group_by(HPVstat) %>% mutate(Prop = Freq / sum(Freq)) %>% ungroup
p <- sum_hpv_status_MPs_prop %>%
  dplyr::mutate(Group = ifelse(MPid %in% c("T_cell", "B_cell", "Macrophage", "Fibroblast", "Endothelial", "Complement", "Secretory_Norm"), "Non-Malignant", "Malignant")) %>%
  dplyr::group_by(MPid, HPVstat) %>%
  summarize(Malignant = sum(Prop[Group == "Malignant"]), `Non-Malignant` = -sum(Prop[Group == "Non-Malignant"])) %>%
  dplyr::group_by(HPVstat) %>%
  dplyr::mutate(Ord = sum(`Non-Malignant`)) %>%
  ggplot(aes(x = fct_reorder(HPVstat, Ord), y = Malignant, fill = MPid)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  geom_col(aes(y = `Non-Malignant`), position = "stack") +
  geom_hline(aes(yintercept = 0)) +
  scale_fill_manual(name = "States", values = state_cols[levels(sum_hpv_status_MPs_prop$MPid)]) + scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.3, 0.2)) +
  theme_classic() + theme(axis.title = element_blank(), axis.text.x = element_text(size = 16), plot.title = element_text(hjust = 0.5, size = 16),
                          legend.text = element_text(size = 14), legend.title = element_text(size = 16), aspect.ratio = 1.4) +
  guides(fill = guide_legend(override.aes = list(size = 7), reverse = TRUE))
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_3c.pdf"), plot = p, device = "pdf", dpi = 300)

# Draw statistics
MPs_per_hpv_stats <- metadata %>% dplyr::mutate(HPVstat = as.factor(samples_metadata$HPV[match(Sample, samples_metadata$Sample)])) %>%
  dplyr::filter(!is.na(MPid) & MPid != "Skeletal_Muscle") %>% droplevels
prop_stat_tab <- test_state_prop_diff(MPs_per_hpv_stats$MPid, MPs_per_hpv_stats$Sample, MPs_per_hpv_stats$HPVstat, transform = NULL)
write.csv(prop_stat_tab, file = here("results/Generated_Data/State_Proportion_Per_HPV_Status_Ttest_Results_Extend.csv"))


# Zoom in on pEMT, Hypoxia and Inflammatory
samples_metadata$Ext_Site <- ifelse(samples_metadata$Site == "Oropharynx" & samples_metadata$HPV == "HPVpos", "HPV+ OP",
                                    ifelse(samples_metadata$Site == "Oropharynx" & samples_metadata$HPV == "HPVneg", "HPV- OP", samples_metadata$Site))
samples_metadata <- samples_metadata %>% dplyr::mutate(Ext_Site = stringr::str_replace_all(Ext_Site, pattern = "Oral", replacement = "Oral Cavity"))
per_sample_states <- as.data.frame(table(metadata$MPid[grep("pEMT|Hypoxia|Inflammatory", metadata$MPid), drop = TRUE], metadata$Sample[grep("pEMT|Hypoxia|Inflammatory", metadata$MPid)])) %>%
  dplyr::mutate(Site = factor(samples_metadata$Ext_Site[match(Var2, samples_metadata$Sample)], levels = c("Laryngeal", "HPV+ OP", "HPV- OP", "Oral Cavity")), .before = "Freq") %>%
  group_by(Var2) %>% mutate(Prop = Freq / sum(Freq), OrderBy = Prop[Var1 == "pEMT"])
p <- ggplot(data = per_sample_states, aes(x = reorder(Var2, OrderBy), y = Prop, fill = Var1)) + geom_col() +
  scale_fill_manual(name = "States", values = state_cols[levels(per_sample_states$Var1)]) + scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0), limits = c(0, NA)) +
  theme_classic() + theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.9, size = 16), plot.title = element_text(hjust = 0.5, size = 16),
                          legend.text = element_text(size = 14), legend.title = element_text(size = 16),
                          strip.background = element_blank(), strip.text = element_text(size = 14)) +
  facet_grid(.~Site, scales = 'free', space = 'free', drop = F) + guides(fill = guide_legend(override.aes = list(size = 7), reverse = TRUE))
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_3d.pdf"), plot = p, device = "pdf", dpi = 300, width = 12.2, height = 8)




# Plot Cellular Density Histogram and Complexity over Cell Count ----------------

summary(metadata$Cell_Count[metadata$Cell_Count <= 47 & metadata$Cell_Count >= 1])
summary(metadata$Cell_Count[metadata$Cell_Count <= 47 & metadata$Cell_Count >= 1 & metadata$EpiStroma == "Epithelial"])
summary(metadata$Cell_Count[metadata$Cell_Count <= 47 & metadata$Cell_Count >= 1 & metadata$EpiStroma == "Stroma"])
add_med_x <- median(metadata$Cell_Count[metadata$Cell_Count <= 47 & metadata$Cell_Count >= 1])
add_med_y <- nrow(metadata[metadata$Cell_Count == add_med_x, ])
p <- ggpubr::gghistogram(data = metadata[metadata$Cell_Count <= 47 & metadata$Cell_Count >= 1 & metadata$EpiStroma != "Filtered_Out", ],
                         x = "Cell_Count", fill = "EpiStroma", bins = 47, position = "stack") +
  scale_fill_manual(name = "Spot Type", values = setNames(c("#94b594", "#df7e66", "#224b5e"), c("Mixed", "Epithelial", "Stroma"))) +
  scalop::theme_scalop(legend.position = "right") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA, colour = "black", linewidth = 2),
                                                          aspect.ratio = 1.2, panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                                          axis.ticks.length.y = grid::unit(0.15, "cm"), legend.position = "bottom", legend.justification = "left",
                                                          axis.text.y = element_text(margin = margin(r = 0.1, unit = "cm"))) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 3100)) + scale_x_continuous(expand = c(0, 0), limits = c(0, 46)) + scale_x_continuous(expand = c(0, 0)) +
  ylab("Number of Spots") + xlab("Number of Cells per Spot") +
  geom_segment(aes(x = add_med_x, y = 0, xend = add_med_x, yend = add_med_y), linetype = 2)
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_1c.pdf"), plot = p, device = "pdf", dpi = 300)


# Plot Complexity as a function of cellular density:
p <- ggplot(metadata, aes(x = Cell_Count, y = Complexity)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", colour = "white") +
  scalop::theme_scalop(legend.position = "right") + scale_x_continuous(expand = c(0.01, 0.01))
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_1a.pdf"), device = "pdf", plot = p, dpi = 300, width = 11, height = 10)




# pEMT Patterns Differential Abundance and Coherence Heatmap --------------

# Load Node attributes tables
node_attr_tab <- read.csv(file = here("results/Generated_Data/Node_Abundance_Neighborhood_Network.csv"))
node_attr_tab <- node_attr_tab %>% dplyr::select(State, pEMTcore, pEMTmargin, Coherence_Core, Coherence_Margin) %>%
  dplyr::mutate(State = stringr::str_replace_all(State, "TNF_Signaling", "Inflammatory"))

# Compute differential abundace and coherence and plot
plot_hm <- node_attr_tab %>% dplyr::reframe(State = State,
                                            Abundance = pEMTcore - pEMTmargin,
                                            Coherence = Coherence_Core - Coherence_Margin) %>%
  dplyr::mutate(Abundance = Abundance * 100) %>%
  reshape2::melt(.) %>% dplyr::mutate(State = factor(State, levels = c("T_cell", "B_cell", "Macrophage", "Complement", "Fibroblast", "Endothelial", "Cell_Cycle", "pEMT", "Epithelial", "Hypoxia", "Inflammatory", "Senescence")))
p1 <- gmap(plot_hm[plot_hm$variable == "Abundance", ], x = State, y = variable, fill = value, limits = c(-7, 7), midpoint = c(-2, 2), geom = "tile", tile.size = 0,
          y.name = "", legend.height = 0.6, legend.width = 8, x.num = F, y.num = F, minimal = TRUE, legend.title = "Percent Difference") +
  theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 0.5, axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 14), title = element_text(size = 14), legend.direction = "horizontal") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  facet_grid(. ~ State, scales = "free") + theme(strip.text.x = element_text(size = 14, face = "bold"))
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_4e_abundance.pdf"), plot = p1, device = "pdf", dpi = 300, width = 20, height = 8)

p2 <- gmap(plot_hm[plot_hm$variable == "Coherence", ], x = State, y = variable, fill = value, limits = c(-0.3, 0.3), midpoint = c(-0.085, 0.085), geom = "tile", tile.size = 0,
           y.name = "", legend.height = 0.6, legend.width = 8, x.num = F, y.num = F, minimal = TRUE, legend.title = "Score Difference") +
  theme_bw() + theme(panel.grid = element_blank(), aspect.ratio = 0.5, axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 14), title = element_text(size = 14), legend.direction = "horizontal") +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  facet_grid(. ~ State, scales = "free") + theme(strip.text.x = element_text(size = 14, face = "bold"))
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_4e_coherence.pdf"), plot = p2, device = "pdf", dpi = 300, width = 20, height = 8)





# pEMT DEGs between pEMT-core & pEMT-margin  ------------------------------

# Load merged expression matrix
merged_expr_mat <- readRDS(file = here("results/Generated_Data/merged_expr_mat.rds"))
malig_states <- c("Senescence", "Cell_Cycle", "Epithelial", "Hypoxia", "pEMT", "LowQ", "Secretory_Malig")
norm_epi <- "Secretory_Norm"
stromal_types <- c("Fibroblast", "Endothelial", "Complement", "Skeletal_Muscle")
immune_types <- c("T_cell", "B_cell", "Macrophage", "Inflammatory")
col_pal <- c(state_cols, setNames(c("#FCF2B4FF", "#51127CFF", "#FC8E64FF", "#C73D73FF"), c("Non-Malignant_Epithelial", "Malignant", "Immune", "Stroma")))

# Select spots to test DEG on
spot_compos <- metadata %>% dplyr::select(Key, Fibroblast:Secretory_Norm, MPid, Zone)
INF_pEMT_spots <- spot_compos %>% dplyr::filter(pEMT > 0 & Inflammatory > 0 & Fibroblast == 0 & Zone == "Zone_3") %>%
  dplyr::rowwise(.) %>% dplyr::mutate(pEMT_INF = sum(pEMT, Inflammatory)) %>%
  dplyr::filter(pEMT_INF >= 0.5) %>% dplyr::select(Key, pEMT, Inflammatory, Fibroblast, MPid) %>% ungroup %>%
  dplyr::filter(MPid == "pEMT")
CAF_pEMT_spots <- spot_compos %>% dplyr::filter(pEMT > 0 & Fibroblast > 0 & Inflammatory == 0 & Zone == "Zone_1") %>%
  dplyr::rowwise(.) %>% dplyr::mutate(pEMT_CAF = sum(pEMT, Fibroblast)) %>%
  dplyr::filter(pEMT_CAF >= 0.5) %>% dplyr::select(Key, pEMT, Inflammatory, Fibroblast, MPid) %>% ungroup %>%
  dplyr::filter(MPid == "pEMT")


# Check for colocalizing cell types/states in the choosen spots:
coloc_states_ls <- list()
malig_fracs_ls <- list()
pie_charts_ls <- list()
for(pattern in c("Core", "Margin")) {
  if(pattern == "Core") {
    select_spots <- INF_pEMT_spots$Key
    i <- 1
  } else if (pattern == "Margin") {
    select_spots <- CAF_pEMT_spots$Key
    i <- 2
  }

  # Calculate percentage of cell types/states that colocalize with a specific cell type/state
  state_coloc <- metadata %>% dplyr::filter(Key %in% select_spots) %>% dplyr::select(Fibroblast:Secretory_Norm) %>% colMeans %>% as.data.frame %>%
    tibble::rownames_to_column() %>% `colnames<-`(c("Var1", "Coloc")) %>%
    dplyr::mutate(Group = dplyr::case_when(Var1 %in% malig_states ~ "Malignant",
                                           Var1 == norm_epi ~ "Non-Malignant_Epithelial",
                                           Var1 %in% stromal_types ~ "Stroma",
                                           Var1 %in% immune_types ~ "Immune"), .before = "Coloc") %>%
    dplyr::mutate(frac = (Coloc / sum(Coloc)) * 100) %>%
    dplyr::filter(frac >= 1) %>% dplyr::mutate(frac = frac * (100 / sum(frac)))
  order_groups <- state_coloc %>% dplyr::group_by(Group) %>% dplyr::reframe(order_fracs = sum(frac)) %>% dplyr::arrange(desc(order_fracs)) %>% dplyr::pull(Group)
  state_coloc <- state_coloc %>% dplyr::mutate(Group = factor(Group, levels = order_groups)) %>%
    dplyr::group_by(Group) %>% dplyr::arrange(desc(frac), .by_group = TRUE) %>% ungroup() %>%
    dplyr::mutate(ymax = cumsum(frac), lab_ypos = cumsum(frac) - 0.5 * frac) %>%
    dplyr::mutate(ymin = c(0, head(ymax, n = -1)), text_size = dplyr::case_when(frac < 1 ~ 0,
                                                                                frac >= 10 ~ 10,
                                                                                frac > 1 & frac < 10 ~ round(frac)))

  coloc_states_ls[[i]] <- state_coloc
  names(coloc_states_ls)[i] <- pattern

  # Calculate fraction of malignant cell states vs non-malognant cell types for cell that colocalize with a specific cell type/state
  state_malig_frac <- state_coloc %>% dplyr::group_by(Group) %>% dplyr::reframe(Malig_frac = sum(frac)) %>%
    dplyr::arrange(desc(Malig_frac)) %>% dplyr::mutate(Group = factor(Group, levels = Group))
  if(!all(c("Malignant", "Immune", "Stroma", "Non-Malignant_Epithelial") %in% levels(state_malig_frac$Group))) {
    state_malig_frac <- state_malig_frac %>% dplyr::add_row(Group = setdiff(c("Malignant", "Immune", "Stroma", "Non-Malignant_Epithelial"), levels(state_malig_frac$Group)),
                                                            Malig_frac = 0)
  }
  malig_fracs_ls[[i]] <- state_malig_frac
  names(malig_fracs_ls)[i] <- pattern
  state_malig_frac <- state_malig_frac %>% dplyr::mutate(ymax = cumsum(Malig_frac)) %>% dplyr::mutate(ymin = c(0, head(ymax, n = -1))) %>%
    dplyr::mutate(lab_ypos = ifelse(Malig_frac < 3, NA, cumsum(Malig_frac) - 0.5 * Malig_frac))

  # Plot Pie-chart
  p <- ggplot() +
    geom_rect(data = state_coloc, aes(fill = factor(Var1, levels = rev(Var1)), ymax = ymax, ymin = ymin, xmax = 6, xmin = 0), color = "black") +
    geom_rect(data = state_malig_frac, aes(fill = Group, ymax = ymax, ymin = ymin, xmax = 7, xmin = 6), color = "black") +
    coord_polar(theta = "y", start = 0) +
    geom_text(data = state_coloc, aes(x = 3.5, y = lab_ypos, label = paste0(round(frac), "%")), color = "black", size = state_coloc$text_size) +
    geomtextpath::geom_textpath(data = state_malig_frac, aes(x = 6.5, y = lab_ypos, label = paste0(round(Malig_frac), "%")), angle = 90, size = 5, color = "#7D7D7D", fontface = "bold") +
    scale_fill_manual(values = col_pal, name = "State") + theme_void() +
    theme(legend.text = element_text(size = 16), legend.title = element_text(size = 18), plot.title = element_text(size = 8, hjust = 0.5, face = "bold", margin = margin(5, 0, -10, 0)), legend.position = "none") +
    guides(fill = guide_legend(override.aes = list(size = 10), reverse = TRUE)) + ggtitle(pattern)
  pie_charts_ls[[i]] <- p
  names(pie_charts_ls)[i] <- pattern
}
p <- do.call(ggarrange, pie_charts_ls)
ggplot2::ggsave(filename = here("results/Diagnostic_Plots/CAF-pEMT_vs_INF-pEMT_Spots_Colocalization.png"), plot = p)


# Remove genes with zero variance across all cells
var_filter <- apply(merged_expr_mat, 1, var)
filt_merged_expr_mat <- merged_expr_mat[var_filter != 0, ]
# Filter out lowly expressed genes and subset column to contain only comparison groups
filt_merged_expr_mat <- filt_merged_expr_mat[rowMeans(filt_merged_expr_mat) >= 0.1, ]
cent_merged_expr_mat <- filter_10x(filt_merged_expr_mat, centre = "mean", log = FALSE, raw_to_CPM = FALSE, complexity_thresh = NULL, percent_mt_thresh = NULL, gene_filt_method = "none")

# Select spots to test DEG on
spot_compos <- metadata %>% dplyr::select(Key, Fibroblast:Secretory_Norm, MPid, Zone)
INF_pEMT_spots <- spot_compos %>% dplyr::filter(pEMT > 0 & Inflammatory > 0 & Fibroblast == 0 & Zone == "Zone_3") %>%
  dplyr::rowwise(.) %>% dplyr::mutate(pEMT_INF = sum(pEMT, Inflammatory)) %>%
  dplyr::filter(pEMT_INF >= 0.75) %>% dplyr::select(Key, pEMT, Inflammatory, Fibroblast, MPid) %>% ungroup %>%
  dplyr::filter(MPid == "pEMT")
CAF_pEMT_spots <- spot_compos %>% dplyr::filter(pEMT > 0 & Fibroblast > 0 & Inflammatory == 0 & Zone == "Zone_1") %>%
  dplyr::rowwise(.) %>% dplyr::mutate(pEMT_CAF = sum(pEMT, Fibroblast)) %>%
  dplyr::filter(pEMT_CAF >= 0.75) %>% dplyr::select(Key, pEMT, Inflammatory, Fibroblast, MPid) %>% ungroup %>%
  dplyr::filter(MPid == "pEMT")

# INF_spots <- spot_compos %>% dplyr::filter(Inflammatory >= 0.6) %>% dplyr::pull(Key)
INF_spots <- spot_compos %>% dplyr::filter(Inflammatory >= 0.75) %>% dplyr::pull(Key)
CAF_spots <- spot_compos %>% dplyr::filter(Fibroblast >= 0.75) %>% dplyr::pull(Key)

# cent_CAF_INF_expr_mat <- cent_merged_expr_mat[, c(INF_spots, CAF_spots)]
cent_pemt_expr_mat <- cent_merged_expr_mat[, c(INF_pEMT_spots$Key, CAF_pEMT_spots$Key)]
cent_CAF_INF_expr_mat <- cent_merged_expr_mat[, c(INF_spots, CAF_spots)]      # Excluding HPV-positive samples
plot_df <- data.frame(INF_CAF = rowMeans(cent_CAF_INF_expr_mat[, is.element(colnames(cent_CAF_INF_expr_mat), INF_spots)]) - rowMeans(cent_CAF_INF_expr_mat[,is.element(colnames(cent_CAF_INF_expr_mat), CAF_spots)]),
                      INF_pEMT_CAF_pEMT = rowMeans(cent_pemt_expr_mat[,is.element(colnames(cent_pemt_expr_mat), INF_pEMT_spots$Key)]) - rowMeans(cent_pemt_expr_mat[,is.element(colnames(cent_pemt_expr_mat), CAF_pEMT_spots$Key)]))
plot_df <- plot_df %>% dplyr::mutate(Label = ifelse(abs(INF_pEMT_CAF_pEMT) >= 2 & abs(INF_CAF) <= 0.5, rownames(.), NA)) %>%
  dplyr::mutate(Group = dplyr::case_when(rownames(.) %in% MPs$Fibroblast ~ "CAF Gene",
                                         rownames(.) %in% MPs$Inflammatory ~ "INF Gene",
                                         !is.na(Label) ~ "DEG", .default = "Not Significant"))

lr_network <- readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
lr_network <- lr_network %>% dplyr::filter(database == "ramilowski" | database == "kegg")
ligands <- lr_network %>% dplyr::pull(from) %>% unique()
receptors <- lr_network %>% dplyr::pull(to) %>% unique()

plot_df <- plot_df %>% dplyr::mutate(Label = ifelse(Group != "Not Significant" & rownames(.) %in% c(ligands, receptors), rownames(.), NA),
                                      LR = dplyr::case_when(Group != "Not Significant" & rownames(.) %in% ligands ~ "Ligands",
                                                            Group != "Not Significant" & rownames(.) %in% receptors ~ "Receptors", .default = "None")) %>%
  dplyr::mutate(Group = dplyr::case_when(rownames(.) %in% MPs$Fibroblast ~ "CAF Gene",
                                         rownames(.) %in% MPs$Inflammatory ~ "INF Gene",
                                         .default = "Other"),
                LR = factor(LR, levels = c("Ligands", "Receptors", "None")))

plot_df <- plot_df %>% dplyr::mutate(Lab_Type = ifelse(!is.na(Label), "Label", NA),
                                       Label = ifelse(rownames(.) %in% c("OSM", "OSMR", "LIFR", "IL6ST", "IL1A", "IL1B", "IL1RN", "IL1R1", "IL1R2", "SPP1", "CD44", "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3"), rownames(.), Label),
                                       LR = dplyr::case_when(!is.na(Label) & rownames(.) %in% ligands ~ "Ligands",
                                                             !is.na(Label) & rownames(.) %in% receptors ~ "Receptors", .default = "None"))
p <- ggplot(plot_df, aes(x = INF_CAF, y = INF_pEMT_CAF_pEMT, label = Label, fill = Group, color = LR)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_point(shape = 21,
             stroke = dplyr::case_when(plot_df$Group == "Other" & is.na(plot_df$Label) ~ 0.8,
                                       plot_df$Group %in% c("CAF Gene", "INF Gene") & is.na(plot_df$Label) ~ 1, .default = 1.5),
             size = dplyr::case_when(plot_df$Group == "Other" & is.na(plot_df$Label) ~ 1,
                                     plot_df$Group %in% c("INF Gene", "CAF Gene") & is.na(plot_df$Label) ~ 2,
                                     !is.na(plot_df$Label) ~ 3),
             alpha = ifelse(plot_df$Group == "Other" & is.na(plot_df$Label), 0.2, 1)) +
  scale_fill_manual(name = "Group", values = setNames(c(alpha("grey30", 0.2), alpha(unname(state_cols[c("Fibroblast", "Inflammatory")]), 0.6)), c("Other", "CAF Gene", "INF Gene"))) +
  scale_color_manual(name = "Ligand-Receptor", values = setNames(c("grey30", "#053061", "#b2182b"), c("None", "Ligands", "Receptors"))) +
  theme_classic() + labs(x = expression("INF - CAF (log"[2]*"FC)"), y = expression("INF-pEMT - CAF-pEMT (log"[2]*"FC)")) +
  ggrepel::geom_text_repel(data = plot_df[is.na(plot_df$Lab_Type) & !is.na(plot_df$Label), ], max.overlaps = Inf, color = "black", show.legend = F, fontface = "bold") +
  ggrepel::geom_label_repel(data = plot_df[plot_df$Lab_Type == "Label" & !is.na(plot_df$Label), ], max.overlaps = Inf, color = "black", show.legend = FALSE, size = 5) +
  scalop::theme_scalop(legend.position = "right") +
  guides(fill = guide_legend(override.aes = list(size = 5)),
         color = guide_legend(override.aes = list(stroke = 3, size = 3)))
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_5d.pdf"), plot = p, device = "pdf", dpi = 300, width = 14, height = 12)

