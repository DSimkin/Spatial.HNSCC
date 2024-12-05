library(Seurat)
library(tidyverse)
library(scalop)
library(patchwork)
library(reshape2)
library(ggpubr)
library(scales)
library(egg)
library(MetBrewer)
library(here)
source(here("analysis_scripts/functions/scRNA_Functions.R"))
source(here("analysis_scripts/functions/ST_Functions.R"))




# Step 1 - Extract Leiden Clusters Gene Programs --------------------------

# Load per-sample parameter table
sample_parameters <- readRDS(file = here("metadata/per_sample_parameters_df.rds"))

if (!dir.exists(file.path(here("results"), "Leiden_programs"))) {
  dir.create(file.path(here("results"), "Leiden_programs"))
}

# Iterate through all samples to get Leiden clusters and their programs
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
  m <- fix_dup_genes(m, stat = "sum")
  HN_obj <- CreateSeuratSpatial(dgCMatrix = m, data.dir = data_dir)

  # Filter spots
  HN_obj[["percent.mt"]] <- PercentageFeatureSet(HN_obj, pattern = "^MT-")
  HN_obj <- subset(HN_obj, subset = percent.mt < sample_parameters$MT_cut[[i]] & nFeature_Spatial > sample_parameters$Complexity_cut[[i]])

  # Process object
  all_genes <- rownames(HN_obj)
  HN_obj <- NormalizeData(HN_obj, normalization.method = "LogNormalize", scale.factor = 10000) %>%
    FindVariableFeatures(assay = "Spatial", selection.method = "vst", nfeatures = 8000) %>%
    ScaleData(features = all_genes) %>%
    RunPCA(assay = "Spatial", verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:sample_parameters$Dim_Red[[i]], k.param = 20, prune.SNN = 1/15, n.trees = 100) %>%
    FindClusters(algorithm = 4, resolution = sample_parameters$Clust_Res[[i]], n.iter = 100, verbose = FALSE) %>%
    FindAllMarkers(min.pct = 0.25, logfc.threshold = 0.35, return.thresh = .001, only.pos = TRUE)

  # Extract top 50 genes per cluster and save these programs
  clusters_top_genes <- as.data.frame(matrix(data = 0, nrow = 50, ncol = length(unique(HN_obj$cluster)))) %>% rename_all(~ stringr::str_replace_all(., "V", "Cluster_"))
  for (j in seq_along(unique(HN_obj$cluster))){
    if(nrow(HN_obj %>% filter(p_val_adj < 0.001 & cluster == j)) >= 50){
      clusters_top_genes[, j] <- HN_obj %>% filter(p_val_adj < 0.001 & cluster == j) %>% .[order(.$avg_log2FC, decreasing = TRUE), ] %>% head(50) %>% dplyr::select(gene)
    } else{
      de_genes <- HN_obj %>% filter(p_val_adj < 0.001 & cluster == j) %>% .[order(.$avg_log2FC, decreasing = TRUE), ] %>% dplyr::select(gene)
      fill_na <- rep(NA, 50 - nrow(de_genes))
      clusters_top_genes[, j] <- c(de_genes$gene, fill_na)
    }
  }
  write_tsv(clusters_top_genes, file = paste0(here(), "/results/Leiden_programs/", sample_parameters$Sample[[i]], "_Clusters_Program.tsv"))
}




# Step 2 - Construct Primary Epithelial vs Non-Epithelial Metaprograms From Leiden Programs Only -------

# Load all samples gene programs
all_progs <- list()
toLoad <- list.files(here("results/Leiden_programs"))
for(i in seq_along(toLoad)){
  prog <- read.table(file = paste(here("results/Leiden_programs"), toLoad[[i]], sep = "/"), header = TRUE)
  prog <- prog %>% rename_all(~ stringr::str_replace_all(., "Cluster", substri(toLoad[[i]], pos = 1, sep = "_")))
  all_progs <- do.call(c, list(all_progs, prog))
}
rm_na <- lapply(all_progs, function(x) na.omit(x))
all_progs <- all_progs[lengths(rm_na) > 40]

# Hierarchical clustering to split into 2 groups - Epithelial and Non-Epithelial programs
jac_mat <- scalop::jaccard(all_progs)
hc <- hclust(dist(1 - jac_mat), method = "average")
jac_plot <- melt(jac_mat[hc$order, hc$order])

p1 <- ggplot(data = jac_plot, aes(x = Var1, y = Var2, fill = value)) + geom_raster() +
  scale_fill_gradient2(limits = c(0, 0.2), low = "dodgerblue4", mid = "antiquewhite", high = "red4", midpoint = mean(jac_plot$value), oob = squish, name = "Jaccard\nIndex") +
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), legend.text.align = 0.5) +
  scale_x_discrete(name = "\nPrograms", labels = NULL, breaks = NULL) + scale_y_discrete(name = "\nPrograms", labels = NULL, breaks = NULL)

prog_ordered <- colnames(jac_mat[hc$order, hc$order])
plot(hc)
clusterCut <- stats::cutree(tree = hc, k = 6)

metaprog_df <- data.frame(row.names = prog_ordered,
                          cluster = as.factor(clusterCut[prog_ordered]),
                          cells = prog_ordered)
metaprog_annotation_plot <- ggplot(metaprog_df, aes(x = factor(cells, levels = cells), y = 1)) +
  geom_raster(aes(fill = cluster)) + theme(legend.position = "bottom", axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), axis.text = element_blank(), axis.line = element_blank(), axis.text.x = element_blank()) +
  labs(x = "", y = "", fill = "")
MP_plot <- egg::ggarrange(p1, metaprog_annotation_plot, ncol = 1, nrow = 2, heights = c(40, 2))
if (!dir.exists(file.path(here("results"), "Diagnostic_Plots"))) {
  dir.create(file.path(here("results"), "Diagnostic_Plots"))
}
ggplot2::ggsave(filename = here("results/Diagnostic_Plots/Initial_Leiden_Program_Clusters.png"), plot = MP_plot)

# Check cluster identities
clust_list <- split(names(clusterCut), clusterCut)
mp_freqs <- sapply(clust_list, function(k) sort(table(unlist(all_progs[k])), decreasing = T), simplify = F)
metaprograms <- sapply(mp_freqs, function(tab) head(names(tab)[tab >= (max(tab)/2)], 100), simplify = F)

## Epithelial programs
epi_prog_names <- names(clusterCut[which(clusterCut == 1 | clusterCut == 5)])
epi_progs_HR <- all_progs[epi_prog_names]
epi_mp_freqs <- mp_freqs[c(1, 5)]
all_epi_gene_freqs <- sort(table(unlist(sapply(epi_mp_freqs, untable))), decreasing = TRUE)
all_epi_gene_freqs <- all_epi_gene_freqs[-grep("^MT-", names(all_epi_gene_freqs))]
epi_metaprograms <- head(names(all_epi_gene_freqs), 150)

## Non-Epithelial programs (Stromal + Immune)
non_epi_prog_names <- names(clusterCut[which(clusterCut != 1 | clusterCut != 5)])
non_epi_progs_HR <- all_progs[non_epi_prog_names]
non_epi_mp_freqs <- mp_freqs[-c(1, 5)]
all_non_epi_gene_freqs <- sort(table(unlist(sapply(non_epi_mp_freqs, untable))), decreasing = TRUE)
all_non_epi_gene_freqs <- all_non_epi_gene_freqs[-grep("^MT-", names(all_non_epi_gene_freqs))]
non_epi_metaprograms <- head(names(all_non_epi_gene_freqs), 150)

# Combine into a list
rm_shared_genes <- intersect(epi_metaprograms, non_epi_metaprograms)
epi_nonEpi_metaprograms <- list(Epithelial = epi_metaprograms[!epi_metaprograms %in% rm_shared_genes],
                                Non_Epithelial = non_epi_metaprograms[!non_epi_metaprograms %in% rm_shared_genes])
if (!dir.exists(file.path(here("results"), "Generated_Data"))) {
  dir.create(file.path(here("results"), "Generated_Data"))
}
saveRDS(epi_nonEpi_metaprograms, file = here("results/Generated_Data/epi_nonEpi_gene_programs.rds"))




# Step 3 - Split Data to Epithelial, Non-Epithelial & Mixed Based on Scoring: get processed expression matrix for NMF --------

# Load per-sample parameter table
sample_parameters <- readRDS(file = here("metadata/per_sample_parameters_df.rds"))
# Load Epithelial vs Non-Epithelial gene programs for scoring
epi_stroma_progs <- readRDS(file = here("results/Generated_Data/epi_nonEpi_gene_programs.rds"))

for(dirname in c("Per_Sample_Results", "NMF_programs")) {
  if (!dir.exists(file.path(here("results"), dirname))) {
    dir.create(file.path(here("results"), dirname))
  }
}
for(group in c("Epithelial", "Non_Epithelial", "Mixed")) {
  if (!dir.exists(file.path(here("results/NMF_programs"), paste0(group, "_mats")))) {
    dir.create(file.path(here("results/NMF_programs"), paste0(group, "_mats")))
  }
}

# Iterate through all samples classify spots (Epi vs Non-Epi vs Mixed) and output matrices for NMF
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

  # Filter and preprocess - first without centering (for scoring)
  m_proc <- filter_10x(m, centre = "none", complexity_thresh = sample_parameters$Complexity_cut[[i]],
                       percent_mt_thresh = sample_parameters$MT_cut[[i]], merge_method = "sum", gene_filt_method = "exp")

  # Score data-set for Epi vs Non-Epi gene programs
  epi_stroma_sigs <- metaprog_score(m_proc, epi_stroma_progs, score_diff = 0.4, score_cut = 0, bin_size = 100, bin_num = 30, conserved.genes = 0.5, center_rows = TRUE, scale = FALSE)

  # Plot programs score distribution on a spatial heatmap
  epi_stroma_plots <- spatial_heatmap(sig_scores = epi_stroma_sigs$mp_score_list, image_obj = spatial_image, normalize = FALSE)
  save_p <- egg::ggarrange(epi_stroma_plots$Epithelial, epi_stroma_plots$Non_Epithelial, ncol = 2, nrow = 1)
  if (!dir.exists(file.path(here("results/Per_Sample_Results"), sample_parameters$Sample[[i]]))) {
    dir.create(file.path(here("results/Per_Sample_Results"), sample_parameters$Sample[[i]]))
  }
  ggplot2::ggsave(filename = paste0(here("results/Per_Sample_Results/"), sample_parameters$Sample[[i]], "/", sample_parameters$Sample[[i]], "_Epi_vs_Non-Epi_plot.png"),
                  plot = save_p, width = 30, height = 15, units = "cm")

  # Split spots to 3 groups and for each create & save a processed centered matrix with non-negative expression values - for NMF
  spot_groups <- split(names(epi_stroma_sigs$meta_program), epi_stroma_sigs$meta_program)
  names(spot_groups)[names(spot_groups) == "Unresolved"] <- "Mixed"
  save_mats <- lapply(names(spot_groups), function(group) {
    m_centered <- filter_10x(m, cells = spot_groups[[group]], complexity_thresh = sample_parameters$Complexity_cut[[i]],
                             percent_mt_thresh = sample_parameters$MT_cut[[i]], centre = "mean", merge_dup_genes = TRUE, gene_filt_method = "exp", merge_method = "sum")
    m_centered[m_centered < 0] <- 0     # Assign zero to negative values
    if (length(which(rowSums(m_centered) == 0)) > 0) {m_centered <- m_centered[-which(rowSums(m_centered) == 0), ]}   # NMF doesn't work if all the values in one row are zero
    # Save as sparse matrix
    cent_mat <- Matrix(m_centered)
    saveRDS(cent_mat, file = paste0(here("results/NMF_programs/"), group, "_mats/", sample_parameters$Sample[[i]], ".RDS"))
  })
}

# Write all available sample names to a file (for NMF engine script)
write_lines(sample_parameters$Sample, here("results/NMF_programs/samples.txt"), sep = "\n")
# Create output directories to store NMF basis and coefficients
for(group in c("Epithelial", "Non_Epithelial", "Mixed")) {
  if (!dir.exists(file.path(here("results/NMF_programs"), paste0(group, "_nmf_output")))) {
    dir.create(file.path(here("results/NMF_programs"), paste0(group, "_nmf_output")))
  }
}



# Step 4 - Process NMF programs -------------------------------------------
library(NMF)

# Iterate over spot groups (Epi, Non_Epi, Mix) to load all nmf gene (W matrix / Basis) ranks
spot_groups <- c("Epithelial", "Non_Epithelial", "Mixed")
prog_genes_ls <- lapply(spot_groups, function(group) {

  # Load all group's NMF objects for each sample
  path <- paste0(here("results/NMF_programs/"), group, "_nmf_output")
  sample_ls <- list.files(path)

  # Create list of NMF matrices where each sample is an entry
  nmf_mats_ls <- list()
  group_mark <- switch(group,
                       Epithelial = "epithel_",
                       Non_Epithelial = "stroma_",
                       Mixed = "mixed_")

  for(i in seq_along(sample_ls)) {
    nmf_obj <- readRDS(paste(path, sample_ls[[i]], sep = "/"))
    samp_name <- stringr::str_split(sample_ls[[i]], pattern = "_")[[1]][[1]]
    nmf_mats <- c()
    for(j in names(nmf_obj$fit)) {
      get_basis <- basis(nmf_obj$fit[[j]])
      colnames(get_basis) <- paste0(group_mark, samp_name, ".", j, ".", 1:j)
      nmf_mats <- cbind(nmf_mats, get_basis)
    }
    nmf_mats_ls[[i]] <- nmf_mats
    names(nmf_mats_ls)[i] <- paste0(group_mark, samp_name)
    rm(nmf_obj, nmf_mats, get_basis)
  }
  return(nmf_mats_ls)
})

prog_genes_ls <- unlist(prog_genes_ls, recursive = FALSE)
saveRDS(prog_genes_ls, file = here("results/Generated_Data/all_nmf_wBasis_progs.rds"))




# Step 5 - Generate Metaprograms: Consensus Clustering of Both NMF and Leiden programs -------------------------
library(NMF)
samples_metadata <- readRDS(file = here("metadata/samples_metadata.rds"))

# Load Louvain programs
all_progs <- list()
toLoad <- list.files(here("results/Leiden_programs"))
for(i in seq_along(toLoad)){
  prog <- read.table(file = paste(here("results/Leiden_programs"), toLoad[[i]], sep = "/"), header = TRUE)
  new_names <- paste0("louvain_", scalop::substri(toLoad[[i]], pos = 1, sep = "_"), ".1.")
  prog <- prog %>% rename_all(~ stringr::str_replace_all(., "Cluster_", new_names))
  all_progs <- do.call(c, list(all_progs, prog))
}
rm_na <- lapply(all_progs, function(x) na.omit(x))
all_progs <- all_progs[lengths(rm_na) > 40]
all_progs <- do.call(cbind, all_progs)

# Load NMF programs (only those that passed robust_NMF filter)
nmf_programs_sig <- readRDS(file = here("results/Generated_Data/all_nmf_wBasis_progs.rds"))
nmf_programs_sig <- lapply(nmf_programs_sig, function(x) apply(x, 2, function(y) names(sort(y, decreasing = TRUE))[1:50]))
nmf_filter <- robust_nmf_programs(nmf_programs_sig, intra_min = 35, intra_max = 10, inter_filter = TRUE, inter_min = 10)
nmf_programs_sig <- lapply(nmf_programs_sig, function(x) x[, is.element(colnames(x), nmf_filter), drop = FALSE])
nmf_programs_sig <- do.call(cbind, nmf_programs_sig)

# Merge programs
merged_progs <- cbind(nmf_programs_sig, all_progs)

# Calculate similarity between programs
progs_intersect <- apply(merged_progs , 2, function(x) apply(merged_progs, 2, function(y) length(intersect(x, y))))

# Hierarchical clustering of the similarity matrix
progs_intersect_hc <- hclust(as.dist(50 - progs_intersect), method = "average")
progs_intersect_hc <- reorder(as.dendrogram(progs_intersect_hc), colMeans(progs_intersect))
progs_intersect <- progs_intersect[order.dendrogram(progs_intersect_hc), order.dendrogram(progs_intersect_hc)]

# Plot similarity for sanity check
jac_plot <- reshape2::melt(progs_intersect)
ggplot(data = jac_plot, aes(x = Var1, y = Var2, fill = value)) + geom_raster() +
  scale_fill_gradient2(limits = c(0, 30), low = "dodgerblue4", mid = "antiquewhite", high = "red4", midpoint = 0 , oob = squish, name = "Intersecting\nGenes") +
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), legend.text.align = 0.5) +
  scale_x_discrete(name = "\nPrograms", labels = NULL, breaks = NULL) + scale_y_discrete(name = "\nPrograms", labels = NULL, breaks = NULL)

# Save merged-filtered Leiden and NMF gene programs and programs intersect values
# saveRDS(progs_intersect, file = here("results/Generated_Data/merged_programs_intersect.rds"))
# saveRDS(merged_progs, file = here("results/Generated_Data/merged_programs_signatures.rds"))

# Load merged-filtered Leiden and NMF gene programs and programs intersect values
progs_intersect <- readRDS(file = here("results/Generated_Data/merged_programs_intersect.rds"))
progs_signatures <- readRDS(file = here("results/Generated_Data/merged_programs_signatures.rds"))

# Replace ambient genes with a random sequence in order to remove those later (these sequences will not repeat enough times to be included in the final metaprograms)
rm_ambient_igs <- str_c(c("^IGLC3$", "^IGKC$", "^IGLC2$", "^IGHA1$", "^IGHG1$", "^IGHG3$", "^IGHG4$", "^IGHM$", "^IGLC1$", "^IGLC7$", "^IGHG2$"), collapse = "|")
temp_progs <- matrix(nrow = nrow(progs_signatures), ncol = ncol(progs_signatures), dimnames = list(c(1:50), colnames(progs_signatures)))
for(i in 1:ncol(progs_signatures)) {
  if(length(grep(rm_ambient_igs, progs_signatures[, i])) > 0) {
    curr_prog <- progs_signatures[, i]
    rand_str <- stringi::stri_rand_strings(length(grep(rm_ambient_igs, curr_prog)), 5, pattern = "[a-z]")
    curr_prog[grep(rm_ambient_igs, curr_prog)] <- rand_str
    temp_progs[, i] <- curr_prog
  } else {
    temp_progs[, i] <- progs_signatures[, i]
  }
}
progs_signatures <- temp_progs

# Cluster merged programs based on a clustering approach that updates metaprograms in each iteration
nmf_wBasis <- readRDS(file = here("results/Generated_Data/all_nmf_wBasis_progs.rds"))      # NMF programs scoring help in resolving genes at the border of a cluster
MP_clust <- cluster_metaprograms(progs_signatures, progs_intersect, nmf_wBasis, min_init_intersect = 13, min_clust_intersect = 9, min_group_size = 5)   # 13 MPs
# saveRDS(MP_clust, file = here("results/Generated_Data/MP_clust.rds"))

# Sort Jaccard similarity plot according to new clusters:
idxs_sorted <- c()
for(i in 1:length(MP_clust$Cluster_list)) {
  idxs_sorted <- c(idxs_sorted, match(MP_clust$Cluster_list[[i]], colnames(progs_intersect)))
}
idxs_new <- c(idxs_sorted, which(is.na(match(1:ncol(progs_intersect), idxs_sorted))))     # combine indexes in clusters with indexes without clusters

# Plot re-ordered similarity matrix heatmap
MP_plot <- reshape2::melt(progs_intersect[idxs_new, idxs_new])
p1 <- ggplot(data = MP_plot, aes(x = Var1, y = Var2, fill = value)) + geom_raster() +
  scale_fill_gradient2(limits = c(0, 25), low = "dodgerblue4", mid = "antiquewhite", high = "red4", midpoint = 0 , oob = squish, name = "Jaccard\nIndex") +
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), legend.text.align = 0.5, axis.text = element_blank()) +
  scale_x_discrete(name = "\nPrograms") +
  scale_y_discrete(name = "\nPrograms")

## Annotate metaprograms heatmap
metaprog_df <- data.frame(Programs = colnames(progs_intersect)[idxs_new]) %>%
  mutate(Cluster = as.factor(ifelse(Programs %in% Unlist(MP_clust$Cluster_list), yes = names(Unlist(MP_clust$Cluster_list))[match(Programs, Unlist(MP_clust$Cluster_list))], no = "No_Cluster")),
         Sample = scalop::substri(Programs, pos = 2),
         Site = as.factor(samples_metadata$Site[match(Sample, samples_metadata$Sample)]),
         Class = scalop::substri(Programs, pos = 1))
levels(metaprog_df$Cluster) <- list(Fibroblast = "Cluster_1", Senescence = "Cluster_2", Macrophage = "Cluster_3", T_cell = "Cluster_4", Skeletal_Muscle = "Cluster_5", B_cell = "Cluster_6",
                                    Hypoxia_pEMT = "Cluster_7", LowQ = "Cluster_8", Endothelial = "Cluster_9", Cell_Cycle = "Cluster_10", Complement = "Cluster_11", Epithelial = "Cluster_12",
                                    Secretory = "Cluster_13", Unclustered = "No_Cluster")

annotate_clusters <- ggbar(metaprog_df$Cluster, dir = "h", cols = MetBrewer::met.brewer("Derain", nlevels(metaprog_df$Cluster)), legend_title = "Cluster") + theme(legend.direction = "horizontal")
annotate_samples <- ggbar(metaprog_df$Sample, dir = "h", cols = scalop::discrete_colours, legend_title = "Sample") + theme(legend.direction = "horizontal")
annotate_location <- ggbar(metaprog_df$Site, dir = "h", cols = MetBrewer::met.brewer("Degas", 3), legend_title = "Tumor\nLocation") + theme(legend.direction = "horizontal")
annotate_class <- ggbar(metaprog_df$Class, dir = "h", cols = met.brewer("Hokusai1", 5), legend_title = "Class") + theme(legend.direction = "horizontal")

MP_plot <- annotate_location + theme(legend.position = "none") +
  annotate_class + theme(legend.position = "none") +
  annotate_samples + theme(legend.position = "none") +
  annotate_clusters + theme(legend.position = "none") +
  p1 + plot_layout(nrow = 5, heights = c(0.05, 0.05, 0.05, 0.05, 1), guides = "collect")
# ggplot2::ggsave(filename = here("results/Diagnostic_Plots/MetaProgram_Clusters.png"), plot = MP_plot)

lej1 = cowplot::get_plot_component(annotate_clusters, 'guide-box-top', return_all = TRUE)
cowplot::ggdraw(lej1)
lej2 = cowplot::get_plot_component(annotate_samples, 'guide-box-top', return_all = TRUE)
cowplot::ggdraw(lej2)
lej3 = cowplot::get_plot_component(annotate_location, 'guide-box-top', return_all = TRUE)
cowplot::ggdraw(lej3)
lej4 = cowplot::get_plot_component(annotate_class, 'guide-box-top', return_all = TRUE)
cowplot::ggdraw(lej4)

# Save generated metaprograms as csv file
MP <-  do.call(cbind.data.frame, MP_clust$MP_list)
colnames(MP) <- c("Fibroblast", "Senescence", "Macrophage", "T_cell", "Skeletal_Muscle", "B_cell", "Hypoxia_pEMT", "LowQ", "Endothelial", "Cell_Cycle", "Complement", "Epithelial", "Secretory")
MP <- sapply(MP, function(x) c(x[!is.na(x)], x[is.na(x)]))
write.csv(MP, file = here("results/Generated_Data/Merged_MetaProgram_Clusters.csv"), row.names = FALSE, na = '')

# Enrichment for metaprograms annotation:
hg38 <- readRDS(here("aux_data/hg38.rds"))
universe <- hg38$symbol
MP_clusters <- MP_clust$MP_list
# msigdb modules to include
msigdb_sigs <- msigdb(category = c("H","C2","C5","C8"))
marker_genes_list <- readRDS(file = here("aux_data/marker_genes_list.rds"))
# flatten
msigdb_sigs <- unlist(msigdb_sigs, recursive = FALSE)
# run hypergeometric test
enrichments <- enricher(test_gene_sets = MP_clusters,
                        universe = universe,
                        ref_gene_sets = c(msigdb_sigs, marker_genes_list))








# Further Resolve Mixed Custers ---------------------------------------
MP_clust <- readRDS(file = here("results/Generated_Data/MP_clust.rds"))
progs_signatures <- readRDS(file = here("results/Generated_Data/merged_programs_signatures.rds"))

# Separate Hypoxia-pEMT Metaprogram
# Hypoxia vs pEMT
str(MP_clust)
hypox_pEMT_clust <- MP_clust$Cluster_list[[7]]
hypox_pEMT_mat <- progs_signatures[, colnames(progs_signatures) %in% hypox_pEMT_clust]
hypox_pEMT_mat <- setNames(split(hypox_pEMT_mat, rep(1:ncol(hypox_pEMT_mat), each = nrow(hypox_pEMT_mat))), colnames(hypox_pEMT_mat))

jac_mat <- scalop::jaccard(hypox_pEMT_mat)
hc <- hclust(dist(1 - jac_mat), method = "ward.D")
jac_plot <- melt(jac_mat[hc$order, hc$order])

p1 <- ggplot(data = jac_plot, aes(x = Var1, y = Var2, fill = value)) + geom_raster() +
  scale_fill_gradient2(limits = c(0, 0.4), low = "dodgerblue4", mid = "antiquewhite", high = "red4", midpoint = 0 , oob = squish, name = "Jaccard\nIndex") +
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), legend.text.align = 0.5) +
  scale_x_discrete(name = "\nPrograms", labels = NULL, breaks = NULL) + scale_y_discrete(name = "\nPrograms", labels = NULL, breaks = NULL)

prog_ordered <- colnames(jac_mat[hc$order, hc$order])
plot(hc)
rect.hclust(hc , k = 4, border = 2:6)
clusterCut <- stats::cutree(tree = hc, k = 4)

metaprog_df <- data.frame(row.names = prog_ordered,
                          cluster = as.factor(clusterCut[prog_ordered]),
                          cells = prog_ordered)

metaprog_annotation_plot <- ggplot(metaprog_df, aes(x = factor(cells, levels = cells), y = 1)) +
  geom_raster(aes(fill = cluster)) + theme(legend.position = "bottom", axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), axis.text = element_blank(), axis.line = element_blank(), axis.text.x = element_blank()) +
  labs(x = "", y = "", fill = "")
hypox_pEMT_plot <- egg::ggarrange(p1, metaprog_annotation_plot, ncol = 1, nrow = 2, heights = c(40, 2))
# ggplot2::ggsave(filename = here("results/Diagnostic_Plots/Hypoxia_pEMT_Cluster_ZoomIn.png"), plot = hypox_pEMT_plot)

clust_list <- split(names(clusterCut), clusterCut)
mp_freqs <- sapply(clust_list, function(k) sort(table(unlist(hypox_pEMT_mat[k])), decreasing = T), simplify = F)
metaprograms <- sapply(mp_freqs, function(tab) head(names(tab)[tab >= 3], 50), simplify = F)

enrichments <- enricher(test_gene_sets = metaprograms,
                        universe = universe,
                        ref_gene_sets = c(msigdb_sigs, marker_genes_list))

hypox_vs_pEMT_progs <- do.call(cbind.data.frame, metaprograms)
colnames(hypox_vs_pEMT_progs) <- c("Shared", "Hypoxia", "Inflammatory", "pEMT")

# Assign intersecting genes to the the program with highest gene rank
intersect(hypox_vs_pEMT_progs[, 2], hypox_vs_pEMT_progs[, 4])
intersect(hypox_vs_pEMT_progs[, 2], hypox_vs_pEMT_progs[, 3])
intersect(hypox_vs_pEMT_progs[, 3], hypox_vs_pEMT_progs[, 4])
to_pEMT <- c("SERPINE1", "PTHLH", "LAMC2", "INHBA", "IGFBP3")
to_Hypoxia <- c("ADM", "SLC2A1", "ERRFI1", "DDIT4")
to_INF <- c("MMP9", "ANGPTL4", "CXCL8", "PPP1R15A")
hypox_vs_pEMT_progs_final <- list(Hypoxia = hypox_vs_pEMT_progs$Hypoxia[!hypox_vs_pEMT_progs$Hypoxia %in% c(to_pEMT, to_INF)],
                                  pEMT = hypox_vs_pEMT_progs$pEMT[!hypox_vs_pEMT_progs$pEMT %in% c(to_Hypoxia, to_INF)],
                                  Inflammatory = hypox_vs_pEMT_progs$Inflammatory[!hypox_vs_pEMT_progs$Inflammatory %in% c(to_Hypoxia, to_pEMT)])


# Separate Secretory Metaprogram to Malignant vs Non-Malignant Programs
secret_progs <- progs_signatures[, colnames(progs_signatures) %in% MP_clust$Cluster_list$Cluster_13]
secret_progs <- setNames(split(secret_progs, rep(1:ncol(secret_progs), each = nrow(secret_progs))), colnames(secret_progs))

jac_mat <- scalop::jaccard(secret_progs)
hc <- hclust(dist(1 - jac_mat), method = "ward.D")
jac_plot <- melt(jac_mat[hc$order, hc$order])

p1 <- gmap(jac_plot, x = Var1, y = Var2, fill = value, limits = c(0, 1), midpoint = 0.25, ratio = 1, y.name = "Programs", x.name = "Programs",
           axis.rel = 1.5, legend.title.rel = 1.25, legend.rel = 1.0, legend.height = 5, legend.width = 0.6, x.num = F, y.num = F, angle = T)
# ggplot2::ggsave(filename = here("results/Diagnostic_Plots/Resolve_Secretory_MP_to_Malignant_vs_NonMalignant.png"), plot = p1)

prog_ordered <- colnames(jac_mat[hc$order, hc$order])
plot(hc)
rect.hclust(hc , k = 2, border = 2:6)
clusterCut <- stats::cutree(tree = hc, k = 2)

metaprog_df <- data.frame(row.names = prog_ordered,
                          cluster = as.factor(clusterCut[prog_ordered]),
                          cells = prog_ordered,
                          sample = scalop::substri(prog_ordered, pos = 2))

anno_clusts <- ggplot(metaprog_df, aes(x = factor(cells, levels = cells), y = 1)) +
  geom_raster(aes(fill = cluster)) + theme(legend.position = "none", axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), axis.text = element_blank(), axis.line = element_blank(), axis.text.x = element_blank()) +
  labs(x = "", y = "", fill = "")
anno_samps <- ggplot(metaprog_df, aes(x = factor(cells, levels = cells), y = 1)) +
  geom_raster(aes(fill = sample)) + theme(legend.position = "bottom", axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), axis.text = element_blank(), axis.line = element_blank(), axis.text.x = element_blank()) +
  labs(x = "", y = "", fill = "")
p1 <- gmap(jac_plot, x = Var1, y = Var2, fill = value, limits = c(0, 1), midpoint = 0.25, ratio = 1, y.name = "", x.name = "",
           axis.rel = 1.5, legend.title.rel = 1.25, legend.rel = 1.0, legend.height = 5, legend.width = 0.6, x.num = T, y.num = T, angle = T, x.labels = NULL, y.labels = NULL) + theme(axis.ticks = element_blank())
p2 <- egg::ggarrange(p1, anno_clusts, anno_samps, ncol = 1, nrow = 3, heights = c(1, 0.05, 0.05))
# ggplot2::ggsave(filename = here("results/Diagnostic_Plots/Resolve_Secretory_MP_to_Malignant_vs_NonMalignant_Annotated.png"), plot = p2)

clust_list2 <- split(names(clusterCut), clusterCut)
mp_freqs <- sapply(clust_list2, function(k) sort(table(unlist(secret_progs[k])), decreasing = T), simplify = F)
malig_secret <- mp_freqs$`1` %>% .[. >= 3] %>% names()
epithel_secret <- mp_freqs$`2` %>% .[. > 2] %>% names()

intersect(malig_secret, epithel_secret)
setdiff(malig_secret, epithel_secret)
setdiff(malig_secret, names(mp_freqs$`2`))
setdiff(epithel_secret, malig_secret)
setdiff(epithel_secret, names(mp_freqs$`1`))

secretory_MPs <- list(Malignant_Secretory = malig_secret,
                      Epithelial_Secretory = epithel_secret)


# Save the complete consensus metaprograms to a RDS file for downstream use
MPs <- read.csv(file = here("results/Generated_Data/Merged_MetaProgram_Clusters.csv"))
MPs <- as.list(MPs)
MPs <- MPs[!names(MPs) %in% c("Hypoxia_pEMT", "Secretory")]
MPs$Hypoxia <- hypox_vs_pEMT_progs_final$Hypoxia
MPs$pEMT <- hypox_vs_pEMT_progs_final$pEMT
MPs$Inflammatory <- hypox_vs_pEMT_progs_final$Inflammatory
MPs$Secretory_Malig <- secretory_MPs$Malignant_Secretory
MPs$Secretory_Norm <- secretory_MPs$Epithelial_Secretory
saveRDS(MPs, file = here("results/Generated_Data/Final_Metaprograms_Extended.rds"))

MPs_df <- as.data.frame(matrix(nrow = 50, ncol = length(MPs), dimnames = list(1:50, names(MPs))))
for(i in seq_along(MPs)) {
  if(length(MPs[[i]]) == 50) {MPs_df[, i] <- MPs[[i]]}
  else {MPs_df[, i] <- c(MPs[[i]], rep(NA, 50 - length(MPs[[i]])))}
}
write.csv(MPs_df, file = here("results/Generated_Data/Final_Metaprograms_Extended.csv"), row.names = FALSE, na = '')



# Update MP_clust$Cluster_list and MP_clust$MP_list to include the splited mixed (#7th) program
names(MP_clust$MP_list) <- c("Fibroblast", "Senescence", "Macrophage", "T_cell", "Skeletal_Muscle", "B_cell", "Mixed", "LowQ", "Endothelial", "Cell_Cycle", "Complement", "Epithelial", "Secretory")
MP_clust$MP_list <- MP_clust$MP_list[!names(MP_clust$MP_list) %in% c("Mixed", "Secretory")]
MP_clust$MP_list$Hypoxia <- hypox_vs_pEMT_progs_final$Hypoxia
MP_clust$MP_list$pEMT <- hypox_vs_pEMT_progs_final$pEMT
MP_clust$MP_list$Inflammatory <- hypox_vs_pEMT_progs_final$Inflammatory
MP_clust$MP_list$Malignant_Secretory <- secretory_MPs$Malignant_Secretory
MP_clust$MP_list$Epithelial_Secretory <- secretory_MPs$Epithelial_Secretory

new_clust_list <- list(Hypoxia = clust_list$`2`, pEMT = clust_list$`4`, Inflammatory = clust_list$`3`)
new_clust_list2 <- list(Secretory_Malig = clust_list2$`1`, Secretory_Norm = clust_list2$`2`)
MP_clust$Cluster_list <- c(MP_clust$Cluster_list[1:6], new_clust_list, MP_clust$Cluster_list[8:12], new_clust_list2)
names(MP_clust$Cluster_list) <- c("Fibroblast", "Senescence", "Macrophage", "T_cell", "Skeletal_Muscle", "B_cell", "Hypoxia", "pEMT", "Inflammatory", "LowQ", "Endothelial", "Cell_Cycle", "Complement_Fibro", "Epithelial", "Secretory_Malig", "Secretory_Norm")
saveRDS(MP_clust, file = here("results/Generated_Data/MP_clust_all_programs_Extended.rds"))





# Final MP plot -----------------------------------------------------------

MP_clust <- readRDS(file = here("results/Generated_Data/MP_clust_all_programs_Extended.rds"))
progs_intersect <- readRDS(file = here("results/Generated_Data/merged_programs_intersect.rds"))
progs_signatures <- readRDS(file = here("results/Generated_Data/merged_programs_signatures.rds"))

# Sort Jaccard similarity plot according to new clusters:
idxs_sorted <- c()
for(i in 1:length(MP_clust$Cluster_list)) {
  idxs_sorted <- c(idxs_sorted, match(MP_clust$Cluster_list[[i]], colnames(progs_intersect)))
}
idxs_new <- c(idxs_sorted, which(is.na(match(1:ncol(progs_intersect), idxs_sorted))))     # combine indexes in clusters with indexes without clusters

# Plot re-ordered similarity matrix heatmap
prog_jacc_mat <- progs_intersect[idxs_new, idxs_new]
prog_plot_df <- reshape2::melt(prog_jacc_mat)
p1 <- ggplot(data = prog_plot_df, aes(x = Var1, y = Var2, fill = value)) + geom_raster() +
  scale_fill_gradient2(limits = c(0, 25), low = "dodgerblue4", mid = "antiquewhite", high = "red4", midpoint = 0 , oob = squish, name = "Jaccard\nIndex") +
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), legend.text.align = 0.5,
        axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank())

## Annotate metaprograms heatmap
metaprog_df <- data.frame(Programs = colnames(progs_intersect)[idxs_new]) %>%
  mutate(Cluster = factor(ifelse(Programs %in% Unlist(MP_clust$Cluster_list), yes = names(Unlist(MP_clust$Cluster_list))[match(Programs, Unlist(MP_clust$Cluster_list))], no = "Unclustered"),
                          levels = c(names(MP_clust$Cluster_list), "Unclustered")),
         Sample = factor(scalop::substri(Programs, pos = 2), levels = samples_metadata$Sample),
         Site = as.factor(samples_metadata$Site[match(Sample, samples_metadata$Sample)]),
         Class = factor(scalop::substri(Programs, pos = 1)))
levels(metaprog_df$Class) <- list("NMF\nEpithelial" = "epithel", "Leiden" = "louvain", "NMF\nMixed" = "mixed", "NMF\nStroma" = "stroma")

# Load state colors
state_cols <- read.delim(file = here("aux_data/state_col_tab_extend.tsv"), sep = "\t", header = TRUE)
state_cols <- state_cols %>%
  tibble::add_row(V1 = "Unclustered", V2 = "#0c2a50ff")
state_cols <- setNames(state_cols$V2, state_cols$V1)
# Get color palette for samples annotation
sample_cols <- setNames(samples_metadata$Color, samples_metadata$Sample)
# Create annotation bars
annotate_clusters <- ggbar(metaprog_df$Cluster, dir = "h", cols = state_cols, legend_title = "Cluster") +
  theme(legend.direction = "horizontal", legend.text = element_text(size = 14, margin = margin(r = 0.8, unit = 'cm')),
        legend.title = element_text(size = 16, margin = margin(r = 0.8, unit = 'cm')))
annotate_samples <- ggbar(metaprog_df$Sample, dir = "h", cols = sample_cols, legend_title = "Sample") +
  theme(legend.direction = "horizontal", legend.text = element_text(size = 14, margin = margin(r = 0.5, unit = 'cm')),
        legend.title = element_text(size = 16, margin = margin(r = 0.5, unit = 'cm')))
annotate_location <- ggbar(metaprog_df$Site, dir = "h", cols = setNames(c("#94b594", "#df7e66", "#224b5e"), c("Laryngeal", "Oral", "Oropharynx")), legend_title = "Tumor\nLocation") +   #  "#41AB5D", "#FEB24C", "#2171B5"   "#005A32", "#FD8D3C", "#6BAED6"
  theme(legend.direction = "horizontal", legend.text = element_text(size = 14, margin = margin(r = 0.5, unit = 'cm')),
        legend.title = element_text(size = 16, margin = margin(r = 0.5, unit = 'cm')))
annotate_class <- ggbar(metaprog_df$Class, dir = "h", cols = setNames(c("#9E4A11", "#E19F33", "#408777", "#2A604D"), c("Leiden", "NMF\nEpithelial", "NMF\nMixed", "NMF\nStroma")), legend_title = "Class") +
  theme(legend.direction = "horizontal", legend.text = element_text(size = 14, margin = margin(r = 0.5, unit = 'cm')),
        legend.title = element_text(size = 16, margin = margin(r = 0.5, unit = 'cm')))

MP_plot <- annotate_class + theme(legend.position = "none") +
  annotate_samples + theme(legend.position = "none") +
  annotate_location + theme(legend.position = "none") +
  annotate_clusters + theme(legend.position = "none") +
  p1 + plot_layout(nrow = 5, heights = c(0.05, 0.05, 0.05, 0.05, 1), guides = "collect")
# ggplot2::ggsave(filename = here("results/Diagnostic_Plots/MetaProgram_Clusters_V2.png"), plot = MP_plot)

lej1 = cowplot::get_plot_component(annotate_clusters, 'guide-box-top', return_all = TRUE)
cowplot::ggdraw(lej1)
lej2 = cowplot::get_plot_component(annotate_samples, 'guide-box-top', return_all = TRUE)
cowplot::ggdraw(lej2)
lej3 = cowplot::get_plot_component(annotate_location, 'guide-box-top', return_all = TRUE)
cowplot::ggdraw(lej3)
lej4 = cowplot::get_plot_component(annotate_class, 'guide-box-top', return_all = TRUE)
cowplot::ggdraw(lej4)


# MP plot without the non-clustered part:
progs2rm <- metaprog_df %>% dplyr::filter(Cluster == "Unclustered") %>% dplyr::pull(Programs)
filt_prog_jacc_mat <- prog_jacc_mat[!rownames(prog_jacc_mat) %in% progs2rm, !colnames(prog_jacc_mat) %in% progs2rm]
filt_metaprog_df <- metaprog_df %>% dplyr::filter(Cluster != "Unclustered")

prog_plot_df <- reshape2::melt(filt_prog_jacc_mat)
p1 <- ggplot(data = prog_plot_df, aes(x = Var1, y = Var2, fill = value)) + geom_raster() +
  scale_fill_gradient2(limits = c(0, 25), low = "dodgerblue4", mid = "antiquewhite", high = "red4", midpoint = 0 , oob = squish, name = "Jaccard\nIndex") +
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), legend.text.align = 0.5,
        axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank())

# Create annotation bars
annotate_clusters <- ggbar(filt_metaprog_df$Cluster, dir = "h", cols = state_cols, legend_title = "Cluster") +
  theme(legend.direction = "horizontal", legend.text = element_text(size = 14, margin = margin(r = 0.8, unit = 'cm')),
        legend.title = element_text(size = 16, margin = margin(r = 0.8, unit = 'cm')))
annotate_samples <- ggbar(filt_metaprog_df$Sample, dir = "h", cols = colorspace::darken(sample_cols, 0.2), legend_title = "Sample") +
  theme(legend.direction = "horizontal", legend.text = element_text(size = 14, margin = margin(r = 0.5, unit = 'cm')),
        legend.title = element_text(size = 16, margin = margin(r = 0.5, unit = 'cm')))
annotate_location <- ggbar(filt_metaprog_df$Site, dir = "h", cols = setNames(c("#94b594", "#df7e66", "#224b5e"), c("Laryngeal", "Oral", "Oropharynx")), legend_title = "Tumor\nLocation") +   #  "#41AB5D", "#FEB24C", "#2171B5"   "#005A32", "#FD8D3C", "#6BAED6"
  theme(legend.direction = "horizontal", legend.text = element_text(size = 14, margin = margin(r = 0.5, unit = 'cm')),
        legend.title = element_text(size = 16, margin = margin(r = 0.5, unit = 'cm')))
annotate_class <- ggbar(filt_metaprog_df$Class, dir = "h", cols = setNames(c("#9E4A11", "#E19F33", "#408777", "#2A604D"), c("Leiden", "NMF\nEpithelial", "NMF\nMixed", "NMF\nStroma")), legend_title = "Class") +
  theme(legend.direction = "horizontal", legend.text = element_text(size = 14, margin = margin(r = 0.5, unit = 'cm')),
        legend.title = element_text(size = 16, margin = margin(r = 0.5, unit = 'cm')))

MP_plot <- annotate_class + theme(legend.position = "none") +
  annotate_samples + theme(legend.position = "none") +
  annotate_location + theme(legend.position = "none") +
  annotate_clusters + theme(legend.position = "none") +
  p1 + plot_layout(nrow = 5, heights = c(0.05, 0.05, 0.05, 0.05, 1), guides = "collect")
# ggplot2::ggsave(filename = here("results/Diagnostic_Plots/MetaProgram_Clusters_V3.png"), plot = MP_plot)


# Custom color palette
library(RColorBrewer)
library(viridis)
custom_pal <- c(grDevices::colorRampPalette(c(colorspace::lighten("antiquewhite", .75), "antiquewhite", rev(magma(30, begin = .25, end = .9))[1]))(5), rev(magma(30, begin = .25, end = .9)))
lims = c(2, 30)
grad_steps <- (max(lims) - min(lims)) / (length(custom_pal) - 1)
col_vals <- seq(from = min(lims), to = max(lims), by = grad_steps)

p2 <- ggplot(data = prog_plot_df, aes(x = Var1, y = Var2, fill = value)) + geom_raster() +
  scale_fill_gradientn(limits = lims, colours = custom_pal, values = scales::rescale(col_vals), oob = squish, name = "Jaccard\nIndex") +
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), legend.text.align = 0.5,
        axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank())
MP_plot <- annotate_class + theme(legend.position = "none") +
  annotate_samples + theme(legend.position = "none") +
  annotate_location + theme(legend.position = "none") +
  annotate_clusters + theme(legend.position = "none") +
  p2 + plot_layout(nrow = 5, heights = c(0.05, 0.05, 0.05, 0.05, 1), guides = "collect")
if (!dir.exists(file.path(here("results"), "Paper_Figures"))) {
  dir.create(file.path(here("results"), "Paper_Figures"))
}
for(dirtype in c("Main_figs", "Supp_figs")) {
  if (!dir.exists(file.path(here("results/Paper_Figures"), dirtype))) {
    dir.create(file.path(here("results/Paper_Figures"), dirtype))
  }
}
# ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_2a.pdf"), plot = MP_plot, device = "pdf", dpi = 300)

lej1 = cowplot::get_plot_component(annotate_clusters, 'guide-box-top', return_all = TRUE)
pdf(file = here("results/Paper_Figures/Supp_figs/SFig_2a_clust_legend.pdf"), width = 12, height = 3)
cowplot::ggdraw(lej1)
dev.off()
lej2 = cowplot::get_plot_component(annotate_samples, 'guide-box-top', return_all = TRUE)
pdf(file = here("results/Paper_Figures/Supp_figs/SFig_2a_samp_legend.pdf"), width = 12, height = 3)
cowplot::ggdraw(lej2)
dev.off()
lej3 = cowplot::get_plot_component(annotate_location, 'guide-box-top', return_all = TRUE)
pdf(file = here("results/Paper_Figures/Supp_figs/SFig_2a_location_legend.pdf"), width = 12, height = 3)
cowplot::ggdraw(lej3)
dev.off()
lej4 = cowplot::get_plot_component(annotate_class, 'guide-box-top', return_all = TRUE)
pdf(file = here("results/Paper_Figures/Supp_figs/SFig_2a_class_legend.pdf"), width = 12, height = 3)
cowplot::ggdraw(lej4)
dev.off()





# Overlap Between Spatial and Single-Cell Metaprograms --------------------

MP_clust <- readRDS(file = here("results/Generated_Data/MP_clust_all_programs_Extended.rds"))
# Load single-cell HNSCC epithelial metaprograms
scMPs <- read.csv(here("aux_data//HTAN_Cancer_metaprograms.csv"))
scMPs <- split(scMPs$Gene, scMPs$Program)
scMPs <- scMPs[-grep("LQ", names(scMPs))]
# Load single-cell HNSCC TME signatures
sc_sigs <- read.csv(here("aux_data/HTAN_full_signatures.csv"))
sc_sigs <- split(sc_sigs$Gene, sc_sigs$Celltype)
sc_sigs <- sc_sigs[-grep("T-NK-cell", names(sc_sigs))]
nk_vs_t <- read.csv(here("aux_data/nksplit2.csv"))
nk_vs_t <- split(nk_vs_t$Gene, nk_vs_t$Cluster)
sc_sigs <- c(sc_sigs, nk_vs_t, scMPs)


MP_jacc <- scalop::jaccard(MP_clust$MP_list[c("Fibroblast", "Endothelial", "Macrophage", "T_cell", "B_cell", "Epithelial", "Senescence", "Cell_Cycle", "Secretory", "Hypoxia", "pEMT", "Inflammatory")],
                           sc_sigs[c("Fibroblast", "Endothelial", "Macrophage", "T-cell", "B-cell", "Plasma_cell", "Epithelial", "EpiSen_A", "EpiSen_B", "G1S", "G2M", "Mucin", "Hypoxia", "pEMT", "Stress")])
plot_jacc <- reshape2::melt(MP_jacc)

lims <- c(0, max(MP_jacc))
custom_pal <- c(grDevices::colorRampPalette(c(colorspace::lighten("antiquewhite", .75), "antiquewhite", rev(viridis::magma(30, begin = .25, end = .9))[1]))(15), rev(viridis::magma(15, begin = .25, end = .9)))
midpoint <- mean(lims)
zero_col_val <- round(length(custom_pal) / 2)
low_grad_steps <- (midpoint - min(lims)) / (zero_col_val - 1)
high_grad_steps <- (max(lims) - midpoint) / (zero_col_val - 1)
col_vals <- c(seq(from = min(lims), to = midpoint, by = low_grad_steps), seq(from = midpoint, to = max(lims), by = high_grad_steps)[-1])
order_x_axis <- c("Fibroblast", "Endothelial", "T_cell", "B_cell", "Macrophage", "Inflammatory", "pEMT", "Hypoxia", "Cell_Cycle", "Secretory", "Senescence", "Epithelial")
order_y_axis <- c("Fibroblast", "Endothelial", "T-cell", "B-cell", "Plasma_cell", "Macrophage", "Stress", "pEMT", "Hypoxia", "G1S", "G2M", "Mucin", "EpiSen_A", "EpiSen_B", "Epithelial")

p <- ggplot(data = plot_jacc, aes(x = factor(Var2, levels = order_x_axis), y = factor(Var1, levels = order_y_axis), fill = value)) + geom_tile() +
  scale_fill_gradientn(limits = lims, colours = custom_pal, values = scales::rescale(col_vals), oob = squish, name = "Jaccard\nIndex") +
  ylab("Single-cell Metaprograms") + xlab("Spatial Metaprograms") +
  theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.line = element_blank(), aspect.ratio = 1.2,
        legend.text.align = 0.5, text = element_text(size = 14), axis.title.x = element_text(vjust = -0.6), axis.text.x = element_text(angle = 90, hjust = 1), legend.title = element_text(size = 16), axis.title = element_text(size = 16)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))
# ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_2b.pdf"), device = "pdf", plot = p, dpi = 300)


