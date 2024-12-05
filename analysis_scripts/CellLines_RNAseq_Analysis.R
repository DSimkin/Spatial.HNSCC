library(tidyverse)
library(scalop)
library(here)
source(here("analysis_scripts/functions/scRNA_Functions.R"))

# Load RNAseq expression matrix
raw_mat <- read.table(here("aux_data/spatial_transcriptomics_rename_featurecounts.txt"), comment.char = "#", header=TRUE, sep = "\t")
# Remove chromosomal information and move Geneid to rownames
raw_mat <- raw_mat %>% dplyr::select(-c("Chr", "Start", "End", "Strand", "Length")) %>% tibble::column_to_rownames("Geneid")
# Rename samples
colnames(raw_mat) <- scalop::substri(colnames(raw_mat), sep = "_", pos = 3) %>% stringr::str_remove_all(., ".RNA")

# Create metadata + run cell QC and add its parameters to the metadata
cell_stats <- cell_qc(raw_mat)

metadata <- tibble::tibble(Sample = colnames(raw_mat),
                           Cell_Line = scalop::substri(colnames(raw_mat), sep = "\\.", pos = 1),
                           Condition = stringr::str_remove_all(colnames(raw_mat), "93VU.|JHU6.|.Rep1|.Rep2"),
                           Tech_Rep = str_sub(colnames(raw_mat), start = -4),
                           Complexity = unname(cell_stats$Complexity),
                           Percent_Mito = unname(cell_stats$Percent_mt),
                           Percent_HK = unname(cell_stats$HK),)

# Plot QC parameters
ggpubr::ggviolin(metadata$Complexity, fill = "lightskyblue", xlab = "Complexity", ylab = FALSE) + ggpubr::ggviolin(metadata$Percent_Mito, fill = "lightskyblue", xlab = "Percent_Mito", ylab = FALSE) +
  ggpubr::ggviolin(metadata$Percent_HK, fill = "lightskyblue", xlab = "Percent_HouseKeeping", ylab = FALSE)



# Preprocessing and QC ----------------------------------------------------

# Convert raw counts to TPM   ##NOTE: CHECK IF GENE LENGTH NORMALIZATION IS NEEDED
scaling_factor <- 1000000/colSums(raw_mat)
m_CPM <- sweep(raw_mat, MARGIN = 2, STATS = scaling_factor, FUN = "*")

# Filter genes (choose between 2 options: filter based soley on highly expressed genes, or filter based on variability ensuring that the gene is expressed highly in at least 2 samples)
m_filt <- m_CPM[-grep("MT-|MTND1P23|^RPS|^RPL|^RNA|^AC0|^AL3|^AL5|^FP\\d", rownames(m_CPM)), ]
top_var_genes <- apply(log2(m_filt + 1) > 5, 1, function(gene) sum(gene) > 2)
m_filt <- m_filt[top_var_genes, ]

# Log2 transform
m_proc <- log2(m_filt + 1) 
hist(rowMeans(m_proc), breaks = 100)

# Plot post-processing QC parameters
post_process_qc <- cell_qc(m_proc)
ggpubr::ggviolin(post_process_qc$Complexity, fill = "lightskyblue", xlab = "Complexity", ylab = FALSE) + ggpubr::ggviolin(post_process_qc$Percent_mt, fill = "lightskyblue", xlab = "Percent_Mito", ylab = FALSE) +
  ggpubr::ggviolin(post_process_qc$HK, fill = "lightskyblue", xlab = "Percent_HouseKeeping", ylab = FALSE)

# Center the data
get_avg <- rowMeans(m_proc)
m_cent <- sweep(m_proc, 1, get_avg)
hist(as.matrix(m_cent), breaks = 100)



# Samples Correlation and Clustering --------------------------------------

# Create correlation matrix - per sample data
corr_mat <-  cor(as.matrix(m_cent))

# Cluster the correlation matrix
corr_hclust  <- hclust(dist(1 - corr_mat, method = "euclidean"), method = "ward.D2")
corr_clusts <- corr_mat[corr_hclust$order, corr_hclust$order]
corr_clusts_melt  <- reshape2::melt(corr_clusts)

# Get a named (samples) vector of clusters and add it to the metadata
plot(corr_hclust)
rect.hclust(corr_hclust , k = 5, border = 2:6)
clusters_vector <- as.factor(dendextend:::cutree(corr_hclust, k = 5, order_clusters_as_data = FALSE))
clusters_vector <- clusters_vector[corr_hclust$order]
metadata <- clusters_vector %>% 
  enframe(name = "Sample", value = "hclust_cluster") %>% left_join(metadata, .) 

# Get a tibble with the order of the samples from the hclust clustering in order to add annotation bar for the plot
sample_id_orders <- data.frame(Sample = colnames(corr_clusts)) %>% tibble::rowid_to_column(var = "Sample_order")
ordered_metadata <- metadata %>% left_join(., sample_id_orders)

# Generate dendrogram
heatmap_dendogram <- corr_hclust %>%
  as.dendrogram %>% 
  dendextend::set("branches_k_color", k = 5) %>% 
  dendextend::set("branches_lwd", 1) %>% 
  dendextend::set("labels", "none") %>% 
  ggplot() +
  ggtitle("Corr matrix") + #add title
  theme(plot.title = element_text(hjust = 0.5)) + #put the title in the middle
  ggplot2::coord_cartesian(expand = F) # remove empty background around the plot

# Generate heatmap
corr_HM <- ggplot(corr_clusts_melt, aes(Var2, Var1, fill = value)) + 
  geom_tile(color = "black", linewidth = 0.001) + 
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "darkred", oob = scales::squish, midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name = "Correlation\ncoefficient", guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  theme(axis.text.y = element_blank(), axis.text.x =  element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        axis.ticks = element_blank(), axis.line = element_blank(), legend.title = element_text(size = 12), legend.text = element_text(size = 8), 
        plot.margin = margin(t = -10, r = 0, b = -10, l = 0, unit = "pt"), legend.key.size = unit(25, 'pt')) +
  labs(x = "", y = "")

# Make the annotation bar
anno_p <- ggplot(ordered_metadata, aes(x = factor(Sample, levels = Sample), y = 1)) +
  geom_tile(aes(fill = as.factor(Cell_Line)), color = "black", size = 0.02) +
  theme(legend.position = "right", axis.ticks = element_blank(), panel.background = element_rect(fill = "white"), axis.line = element_blank(),
    axis.text.x = element_text(angle = 90, size = 8), axis.text.y = element_blank(), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), legend.key.size = unit(25, 'pt')) +
  scale_x_discrete(labels = ordered_metadata$Sample) + 
  labs(x = "", y = "", fill = "Cell Line")

# Bind the legends
HM_legend <- ggpubr::get_legend(corr_HM)
anno_legend <- ggpubr::get_legend(anno_p)
legends <- ggpubr::ggarrange(HM_legend, anno_legend, nrow = 2)

# Bind all together and save the plot
p <- ggpubr::ggarrange(heatmap_dendogram, corr_HM, anno_p, ncol = 1, heights = c(30, 165, 60), align = "v", legend = "none") %>%
  ggpubr::ggarrange(., legends, widths = c(0.9, 0.1)) 
# ggplot2::ggsave(filename = here("results/Diagnostic_Plots/Bulk_RNAseq_Samples_Correlation_Heatmap.png"), plot = p)




# KMeans on PCA -----------------------------------------------------------

# For the kneeplot function (from pcatools package) the rownames in the metadata needs to be the same as the colnames of the expression matrix
pca_metadata <- metadata %>% as.data.frame %>% tibble::column_to_rownames("Sample") %>% .[colnames(m_cent), ]
# Run PCA
pca_obj <- stats::prcomp(t(m_cent))
# Add PCA component to metadata
pca_metadata <- pca_metadata %>% dplyr::mutate(PC1 = pca_obj$x[, 1],
                                               PC2 = pca_obj$x[, 2],
                                               PC3 = pca_obj$x[, 3], 
                                               PC4 = pca_obj$x[, 4])

p <- ggplot(data = pca_metadata, aes(x = PC1, y = PC2, color = Condition, shape = Cell_Line)) +
  geom_point(size = 5) + xlab(paste0("PC1 (", as.data.frame(summary(pca_obj)$importance)$PC1[[2]] * 100, "%)")) + ylab(paste0("PC2 (", as.data.frame(summary(pca_obj)$importance)$PC2[[2]] * 100, "%)")) +
  scalop::theme_scalop(legend.position = "right") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_8a_left.pdf"), plot = p, device = "pdf", dpi = 300, width = 10, height = 8)

# RUN `KMeans on PCA` PART AGAIN AFTER CENTERING WITHIN EACH CELL LINE SEPERATLEY
VU_mat <- m_proc[, metadata$Sample[metadata$Cell_Line == "93VU"]]
JHU_mat <- m_proc[, metadata$Sample[metadata$Cell_Line == "JHU6"]]
m_cent <- cbind.data.frame(sweep(VU_mat, 1, rowMeans(VU_mat)), sweep(JHU_mat, 1, rowMeans(JHU_mat)))

pca_metadata <- metadata %>% as.data.frame %>% tibble::column_to_rownames("Sample") %>% .[colnames(m_cent), ]
# Run PCA
pca_obj <- stats::prcomp(t(m_cent))
# Add PCA component to metadata
pca_metadata <- pca_metadata %>% dplyr::mutate(PC1 = pca_obj$x[, 1],
                                               PC2 = pca_obj$x[, 2],
                                               PC3 = pca_obj$x[, 3], 
                                               PC4 = pca_obj$x[, 4])

p <- ggplot(data = pca_metadata, aes(x = PC1, y = PC2, color = Condition, shape = Cell_Line)) +
  geom_point(size = 5) + xlab(paste0("PC1 (", as.data.frame(summary(pca_obj)$importance)$PC1[[2]] * 100, "%)")) + ylab(paste0("PC2 (", as.data.frame(summary(pca_obj)$importance)$PC2[[2]] * 100, "%)")) +
  scalop::theme_scalop(legend.position = "right") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_8a_right.pdf"), plot = p, device = "pdf", dpi = 300, width = 10, height = 8)


# Calculate KMeans
kmeans_obj <- stats::kmeans(pca_obj$x[, 1:2], centers = 3)
# Extract the kmeans cluster to a df
kmeans_df <- kmeans_obj$cluster %>% 
  enframe(name = "Sample", value = "kmeans_cluster")

# Combine kmeans to my metadata
metadata <- left_join(metadata, kmeans_df, by = "Sample")

top_PC1_genes <- pca_obj$rotation %>% .[, 1, drop = F] %>% 
  rowMeans() %>% 
  sort(decreasing = T) %>% 
  head(50) %>% 
  names 
bottom_PC1_genes <- pca_obj$rotation %>% .[, 1, drop = F] %>% 
  rowMeans() %>% 
  sort(decreasing = F) %>% 
  head(50) %>% 
  names 
top_PC2_genes <- pca_obj$rotation %>% .[, 2, drop = F] %>% 
  rowMeans() %>% 
  sort(decreasing = T) %>% 
  head(50) %>% 
  names 
bottom_PC2_genes <- pca_obj$rotation %>% .[, 2, drop = F] %>% 
  rowMeans() %>% 
  sort(decreasing = F) %>% 
  head(50) %>% 
  names
PC_genes <- data.frame(Top_PC1 = top_PC1_genes, Bottom_PC1 = bottom_PC1_genes,
                       Top_PC2 = top_PC2_genes, Bottom_PC2 = bottom_PC2_genes,
                       Shared_Top = c(intersect(top_PC1_genes, top_PC2_genes), rep(NA, 50 - length(intersect(top_PC1_genes, top_PC2_genes)))), 
                       Shared_Bottom = c(intersect(bottom_PC1_genes, bottom_PC2_genes), rep(NA, 50 - length(intersect(bottom_PC1_genes, bottom_PC2_genes)))))
write.csv(PC_genes, file = here("results/Generated_Data/Bulk_RNAseq_Top_and_Bottom_50_Genes_for_PC1_and_PC2.csv"))




# Differently Expressed Genes ---------------------------------------------

# Split the uncentered matrix by cell-lines and center within each cell-line separately
VU_mat <- m_proc[, metadata$Sample[metadata$Cell_Line == "93VU"]]
JHU_mat <- m_proc[, metadata$Sample[metadata$Cell_Line == "JHU6"]]
m_cent <- cbind.data.frame(sweep(VU_mat, 1, rowMeans(VU_mat)), sweep(JHU_mat, 1, rowMeans(JHU_mat)))


### OSM vs CONTROL
upregulated_group <- c("93VU.OSM.Rep1", "93VU.OSM.Rep2", "JHU6.OSM.Rep1", "JHU6.OSM.Rep2")  
downregulated_group <- c("93VU.CTR.Rep1", "93VU.CTR.Rep2", "JHU6.CTR.Rep1", "JHU6.CTR.Rep2")
  
m.k1 <- as.matrix(m_cent[, upregulated_group])
m.other <- as.matrix(m_cent[, downregulated_group])
genes <- rownames(m_cent)
lfc <- rowMeans(m.k1) - rowMeans(m.other)
pval <- sapply(rownames(m_cent), function(gene) {
  t.test(m.k1[gene, ], m.other[gene, ], alternative = 'two.sided', paired = T, var.equal = T)$p.value
})
padj <- p.adjust(pval, method = "BH")
ord <- order(lfc, decreasing = T)
de_meta <- data.frame(Lfc = lfc,
                      Pval = pval,
                      Padj = padj,
                      Gene = names(lfc))

de_meta <- de_meta %>% dplyr::mutate(Signif_DE = dplyr::case_when(Pval < 0.01 & Lfc >= 1 ~ "Upregulated by OSM", 
                                                                  Pval < 0.01 & Lfc <= -1 ~ "Upregulated in the control", .default = "Not Significant"),
                                     Label = ifelse(Pval < 0.01 & abs(Lfc) >= 1, Gene, NA))

hg38 <- readRDS(here("aux_data/hg38.rds"))
universe <- rownames(m_cent)
osm_gene_set <- de_meta$Gene[de_meta$Signif_DE == "Upregulated by OSM"]
# msigdb modules and 3CA MPs to include
msigdb_sigs <- msigdb(category = c("H","C2","C5"))
mps_3ca <- as.list(readxl::read_xlsx(path = here("aux_data/3CA_Metaprograms.xlsx"), sheet = 1))
# flatten
msigdb_sigs <- unlist(msigdb_sigs, recursive = FALSE)
# run hypergeometric test
enrichments <- enricher(test_gene_sets = osm_gene_set,
                        universe = universe,
                        ref_gene_sets = c(msigdb_sigs, mps_3ca)) %>% .[[1]]


plot_enrichment_df <- enrichments %>% dplyr::arrange(qvalue) %>% dplyr::slice_head(n = 20) %>% dplyr::select("Description", "qvalue") %>% 
  dplyr::mutate(Collection = dplyr::case_when(scalop::substri(Description, sep = "\\.", pos = 1) == "C2" ~ "Curated",
                                              scalop::substri(Description, sep = "\\.", pos = 1) == "C5" ~ "Ontology",
                                              scalop::substri(Description, sep = "\\.", pos = 1) == "H" ~ "Hallmark", .default = "3CA")) %>% 
  dplyr::mutate(Source = scalop::substri(scalop::substri(Description, sep = "\\.", pos = 2), sep = "_", pos = 1),
                Description = ifelse(Collection == "3CA", Description, 
                                     str_replace_all(scalop::substri(scalop::substri(Description, sep = "\\.", pos = 2), sep = "_", pos = -1), pattern = "_", replacement = " ")),
                qvalue = -log10(qvalue)) %>%
  dplyr::rename("Enricher program" = "Description")
plot_enrichment_df$`Enricher program` <- ave(plot_enrichment_df$`Enricher program`, plot_enrichment_df$`Enricher program`, FUN = function(i) 
  if(length(i) > 1) {
    i[-1] <- paste0(i[-1], " (", letters[seq(length(i[-1])) + 1], ")"); i} else i)
plot_enrichment_df <- plot_enrichment_df %>% dplyr::arrange(qvalue)
# Load state colors 
state_cols <- read.delim(file = here("aux_data/state_col_tab_extend.tsv"), sep = "\t", header = TRUE)
state_cols <- setNames(state_cols$V2, state_cols$V1)
plot_enrichment <- ggpubr::ggbarplot(plot_enrichment_df, x = "Enricher program", y = "qvalue", orientation = "horizontal", fill = "Collection") +
  ggtitle("Enriched programs in OSM induction") + xlab("") + scale_y_continuous(expand = c(0, NA)) +
  scale_fill_manual(values = setNames(c("#9dbac5", "#44b7c2", "#024B7ACC", "#0c2a50ff"), c("Curated", "Hallmark", "Ontology", "3CA"))) +
  theme(title = element_text(size = 16), text = element_text(size = 14)) + ylab("-log10(qvalue)")
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_7b.pdf"), plot = plot_enrichment, device = "pdf", dpi = 300, width = 14, height = 12)

emt_path <- head(enrichments$Description, 20) %>% .[grep("EMT|EPITHELIAL_MESENCHYMAL_TRANSITION", ., ignore.case = TRUE)]
emt_genes <- unique(unlist(str_split(enrichments$geneID[enrichments$Description %in% emt_path], pattern = "/")))
emt_genes_3ca <- c(mps_3ca$`EMT-I`, mps_3ca$`EMT-II`)
emt_genes_3ca <- emt_genes_3ca[emt_genes_3ca %in% de_meta$Gene[de_meta$Signif_DE == "Upregulated by OSM"]]
emt_genes <- unique(c(emt_genes, emt_genes_3ca))
inf_path <- head(enrichments$Description, 20) %>% .[grep("INFLAMMATORY_RESPONSE", ., ignore.case = TRUE)]
inf_genes <- unique(unlist(str_split(enrichments$geneID[enrichments$Description %in% inf_path], pattern = "/")))

de_meta <- de_meta %>% dplyr::mutate(Path = dplyr::case_when(Gene %in% intersect(emt_genes, inf_genes) ~ "Both",
                                                             Gene %in% emt_genes ~ "EMT",
                                                             Gene %in% inf_genes ~ "Inflammatory", .default = NA)) %>% 
  dplyr::mutate(nudge = dplyr::case_when(is.na(Path) & !is.na(Label) & Signif_DE == "Upregulated by OSM" ~ 5.5 - Lfc,
                                         is.na(Path) & !is.na(Label) & Signif_DE == "Upregulated in the control" ~ -4.5 - Lfc, .default = NA))
                
p <- ggplot(data = de_meta, aes(x = Lfc, y = -log10(Pval), label = Label)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') +
  geom_point(shape = 21,
             fill = dplyr::case_when(de_meta$Signif_DE == "Not Significant" ~ alpha("grey30", 0.2),
                                     de_meta$Signif_DE == "Upregulated by OSM" ~ "#b2182b",
                                     de_meta$Signif_DE == "Upregulated in the control" ~ "#2166ac"),
             color = dplyr::case_when(is.na(de_meta$Path) ~ alpha("grey30", 0.5),
                                     de_meta$Path == "Both" ~ "#117733",
                                     de_meta$Path == "EMT" ~ state_cols["pEMT"],
                                     de_meta$Path == "Inflammatory" ~ state_cols["Inflammatory"]),
             stroke = dplyr::case_when(is.na(de_meta$Label) ~ 0.8,
                                       is.na(de_meta$Path) & !is.na(de_meta$Label) ~ 1,
                                       !is.na(de_meta$Path) & !is.na(de_meta$Label) ~ 1.5),
             size = dplyr::case_when(is.na(de_meta$Label) ~ 1,
                                     is.na(de_meta$Path) & !is.na(de_meta$Label) ~ 2,
                                     !is.na(de_meta$Path) & !is.na(de_meta$Label) ~ 3),
             alpha = ifelse(is.na(de_meta$Label), 0.2, 1)) + 
  ggrepel::geom_text_repel(data = de_meta[is.na(de_meta$Path) & !is.na(de_meta$Label), ], 
                           max.overlaps = Inf, color = "black", show.legend = F, force = 0.5, nudge_x = de_meta[is.na(de_meta$Path) & !is.na(de_meta$Label), "nudge"], direction = "y", 
                           hjust = 1, segment.size = 0.3, size = ifelse(de_meta[is.na(de_meta$Path) & !is.na(de_meta$Label), "Lfc"] > 0, 4, 3)) +
  ggrepel::geom_label_repel(data = de_meta[!is.na(de_meta$Path) & !is.na(de_meta$Label), ], aes(fill = Path), max.overlaps = Inf, color = "black", show.legend = FALSE, size = 5) +
  scale_fill_manual(name = "Pathway", values = setNames(c(alpha("#117733", 0.6), alpha(unname(state_cols[c("pEMT", "Inflammatory")]), 0.8)), c("Both", "EMT", "Inflammatory")), na.value = alpha("grey30", 0.5)) +
  theme_classic() + labs(x = expression("INF - CAF (log"[2]*"FC)"), y = expression("INF-pEMT - CAF-pEMT (log"[2]*"FC)")) +
  scalop::theme_scalop(legend.position = "right")
ggplot2::ggsave(filename = here("results/Diagnostic_Plots/Bulk_RNAseq_OSM_vs_Control_DEG_VolcanoPlot.png"), plot = p, width = 14, height = 12)

## Simplified version without gene names, keeping only selected labels labels
de_meta <- de_meta %>% dplyr::mutate(Select_genes = ifelse(Gene %in% c("VEGFC", "CD44", "IL1B", "PDPN", "COL8A2", "LAMC2",  "AREG", "OSMR", "ABCA1", "HIF1A", "CD274", "FAP", "TGFA"), Gene, NA),
                                     nudge = ifelse(!is.na(Select_genes), 4.5 - Lfc, NA))
p <- ggplot(data = de_meta, aes(x = Lfc, y = -log10(Pval), label = Select_genes)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') +
  geom_point(shape = 21,
             fill = dplyr::case_when(de_meta$Signif_DE == "Not Significant" ~ alpha("grey30", 0.2),
                                     de_meta$Signif_DE != "Not Significant" & is.na(de_meta$Path) ~ "#0c2a50ff",
                                     de_meta$Signif_DE != "Not Significant" & de_meta$Path == "Both" ~ "#117733",
                                     de_meta$Signif_DE != "Not Significant" & de_meta$Path == "EMT" ~ state_cols["pEMT"],
                                     de_meta$Signif_DE != "Not Significant" & de_meta$Path == "Inflammatory" ~ state_cols["Inflammatory"]),
             color = dplyr::case_when(de_meta$Signif_DE == "Not Significant" ~ alpha("grey30", 0.5),
                                      de_meta$Signif_DE != "Not Significant" ~ "black"),
             stroke = dplyr::case_when(is.na(de_meta$Label) ~ 0.8,
                                       !is.na(de_meta$Label) ~ 1),
             size = dplyr::case_when(is.na(de_meta$Label) ~ 1,
                                     is.na(de_meta$Path) & !is.na(de_meta$Label) ~ 2,
                                     !is.na(de_meta$Path) & !is.na(de_meta$Label) ~ 3),
             alpha = ifelse(is.na(de_meta$Label), 0.2, 1)) + 
  ggrepel::geom_label_repel(data = de_meta[!is.na(de_meta$Select_genes), ], aes(fill = Path), max.overlaps = Inf, color = "black", show.legend = FALSE, size = 5,
                            nudge_x = de_meta[!is.na(de_meta$Select_genes), "nudge"], direction = "y", hjust = 1) +
  scale_fill_manual(name = "Pathway", values = setNames(c(alpha("#117733", 0.6), alpha(unname(state_cols[c("pEMT", "Inflammatory")]), 0.8)), c("Both", "EMT", "Inflammatory")), na.value = alpha("grey30", 0.5)) +
  theme_classic() + labs(x = expression("Response to OSM (log"[2]*"[Fold Change] average of two cell lines - JHU6 & 93VU)"), y = expression("-log"[10]*"[P-value]")) +
  scalop::theme_scalop(legend.position = "right") + scale_y_continuous(expand = c(0, 0)) + expand_limits(y = max(-log10(de_meta$Pval)) + 0.15) + 
  scale_x_continuous(expand = expansion(add = 0.15), limits = c(min(de_meta$Lfc[!is.na(de_meta$Label)]), max(de_meta$Lfc[!is.na(de_meta$Label)]))) + theme(aspect.ratio = 0.5)
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_5f.pdf"), plot = p, device = "pdf", dpi = 300, width = 14, height = 10)


osm_ctr_cent <- m_cent[, c(upregulated_group, downregulated_group)]
corr_mat <- cor(as.matrix(osm_ctr_cent[rownames(osm_ctr_cent) %in% c(emt_genes, inf_genes), ]))
corr_hclust  <- hclust(dist(1 - corr_mat, method = "euclidean"), method = "ward.D2")
gene_corr_mat <- cor(t(as.matrix(osm_ctr_cent[rownames(osm_ctr_cent) %in% c(emt_genes, inf_genes), ])))
gene_hclust  <- hclust(dist(1 - gene_corr_mat, method = "euclidean"), method = "ward.D2")
plot_df <- osm_ctr_cent[rownames(osm_ctr_cent) %in% c(emt_genes, inf_genes), ] %>% tibble::rownames_to_column("Gene") %>% reshape2::melt() %>% 
  dplyr::mutate(variable = factor(variable, levels = colnames(osm_ctr_cent)[corr_hclust$order]),
                Gene = factor(Gene, levels = rownames(osm_ctr_cent[rownames(osm_ctr_cent) %in% c(emt_genes, inf_genes), ])[gene_hclust$order]))

p <- gmap(plot_df, x = variable, y = Gene, fill = value, limits = c(-2, 2), midpoint = 0, ratio = 1, y.name = "",
          axis.rel = 1.5, legend.title.rel = 1.25, legend.rel = 1.0, legend.height = 5, legend.width = 0.6, x.num = F, y.num = F, magma_pal = FALSE) +
  theme(axis.text = element_text(size = 16), axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1))
ggplot2::ggsave(filename = here("results/Diagnostic_Plots/Bulk_RNAseq_OSM_vs_Control_DEG_Heatmap_short.png"), plot = p, width = 14, height = 12)

order_genes <- de_meta[de_meta$Signif_DE == "Upregulated by OSM", ] %>% dplyr::arrange(desc(Lfc)) %>% dplyr::pull(Gene)
plot_df <- osm_ctr_cent[order_genes, ] %>% tibble::rownames_to_column("Gene") %>% reshape2::melt() %>% 
  dplyr::mutate(Gene = factor(Gene, levels = order_genes))
p <- gmap(plot_df, x = variable, y = Gene, fill = value, limits = c(-2, 2), midpoint = 0, ratio = 1, y.name = "",
          axis.rel = 1.5, legend.title.rel = 1.25, legend.rel = 1.0, legend.height = 5, legend.width = 0.6, x.num = F, y.num = F, magma_pal = FALSE) +
  theme(axis.text = element_text(size = 16), axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1))
ggplot2::ggsave(filename = here("results/Diagnostic_Plots/Bulk_RNAseq_OSM_vs_Control_DEG_Heatmap_long.png"), plot = p, width = 14, height = 12)

# Save OSM vs CTR DEGs to csv file
de_df <- de_meta %>% dplyr::filter(Signif_DE != "Not Significant") %>% 
  dplyr::arrange(desc(Lfc)) %>% dplyr::select(Gene, Lfc, Pval, Signif_DE) %>% dplyr::rename("Group" = "Signif_DE")
write.csv(de_df, file = here("results/Generated_Data/Bulk_RNAseq_OSM_vs_Control_DEGs.csv"), row.names = FALSE)




### TGFB vs CONTROL
upregulated_group <- c("93VU.TGFB.Rep1", "93VU.TGFB.Rep2", "JHU6.TGFB.Rep1", "JHU6.TGFB.Rep2")
downregulated_group <- c("93VU.CTR.Rep1", "93VU.CTR.Rep2", "JHU6.CTR.Rep1", "JHU6.CTR.Rep2")

m.k1 <- as.matrix(m_cent[, upregulated_group])
m.other <- as.matrix(m_cent[, downregulated_group])
genes <- rownames(m_cent)
lfc <- rowMeans(m.k1) - rowMeans(m.other)
pval <- sapply(rownames(m_cent), function(gene) {
  t.test(m.k1[gene, ], m.other[gene, ], alternative = 'two.sided', paired = F, var.equal = T)$p.value
})
padj <- p.adjust(pval, method = "BH")
ord <- order(lfc, decreasing = T)
de_meta <- data.frame(Lfc = lfc,
                      Pval = pval,
                      Padj = padj,
                      Gene = names(lfc))

de_meta <- de_meta %>% dplyr::mutate(Signif_DE = dplyr::case_when(Pval < 0.01 & Lfc >= 1 ~ "Upregulated by TGFB", 
                                                                  Pval < 0.01 & Lfc <= -1 ~ "Upregulated in the control", .default = "Not Significant"),
                                     Label = ifelse(Pval < 0.01 & abs(Lfc) >= 1, Gene, NA))

tgfb_gene_set <- de_meta$Gene[de_meta$Signif_DE == "Upregulated by TGFB"]
# run hypergeometric test
enrichments <- enricher(test_gene_sets = tgfb_gene_set,
                        universe = universe,
                        ref_gene_sets = c(msigdb_sigs, mps_3ca)) %>% .[[1]]
plot_enrichment_df <- enrichments %>% dplyr::arrange(qvalue) %>% dplyr::slice_head(n = 20) %>% dplyr::select("Description", "qvalue") %>% 
  dplyr::mutate(Collection = dplyr::case_when(scalop::substri(Description, sep = "\\.", pos = 1) == "C2" ~ "Curated",
                                              scalop::substri(Description, sep = "\\.", pos = 1) == "C5" ~ "Ontology",
                                              scalop::substri(Description, sep = "\\.", pos = 1) == "H" ~ "Hallmark", .default = "3CA")) %>% 
  dplyr::mutate(Source = scalop::substri(scalop::substri(Description, sep = "\\.", pos = 2), sep = "_", pos = 1),
                Description = ifelse(Collection == "3CA", Description, 
                                     str_replace_all(scalop::substri(scalop::substri(Description, sep = "\\.", pos = 2), sep = "_", pos = -1), pattern = "_", replacement = " ")),
                qvalue = -log10(qvalue)) %>%
  dplyr::rename("Enricher program" = "Description")
plot_enrichment_df$`Enricher program` <- ave(plot_enrichment_df$`Enricher program`, plot_enrichment_df$`Enricher program`, FUN = function(i) 
  if(length(i) > 1) {
    i[-1] <- paste0(i[-1], " (", letters[seq(length(i[-1])) + 1], ")"); i} else i)
plot_enrichment_df <- plot_enrichment_df %>% dplyr::arrange(qvalue)
plot_enrichment <- ggpubr::ggbarplot(plot_enrichment_df, x = "Enricher program", y = "qvalue", orientation = "horizontal", fill = "Collection") +
  ggtitle("Enriched programs in TGFB induction") + xlab("") + scale_y_continuous(expand = c(0, NA)) +
  scale_fill_manual(values = setNames(c("#9dbac5", "#44b7c2", "#024B7ACC", "#0c2a50ff"), c("Curated", "Hallmark", "Ontology", "3CA"))) +
  theme(title = element_text(size = 16), text = element_text(size = 14)) + ylab("-log10(qvalue)")
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_7c.pdf"), plot = plot_enrichment, device = "pdf", dpi = 300, width = 14, height = 12)

tgfb_emt_path <- head(enrichments$Description, 20) %>% .[grep("EMT|EPITHELIAL_MESENCHYMAL_TRANSITION", ., ignore.case = TRUE)]
tgfb_emt_genes <- unique(unlist(str_split(enrichments$geneID[enrichments$Description %in% tgfb_emt_path], pattern = "/")))

col_ord <- colnames(m_cent)[c(grep("CTR", colnames(m_cent)), grep("OSM.Rep", colnames(m_cent)), grep("TGFB", colnames(m_cent)))]
osm_vs_tgfb_emt <- sort(rowMeans(m_cent[unique(c(emt_genes, tgfb_emt_genes)), colnames(m_cent)[grep("OSM.Rep", colnames(m_cent))]]) - rowMeans(m_cent[unique(c(emt_genes, tgfb_emt_genes)), colnames(m_cent)[grep("TGFB", colnames(m_cent))]]), decreasing = TRUE)
test_groups <- split(osm_vs_tgfb_emt, cut(osm_vs_tgfb_emt, c(Inf, 0.7, -0.7, -Inf)))
names(test_groups) <- c("TGFB", "Shared", "OSM.Rep")
test <- sapply(names(test_groups), function(x) {
  if(x == "Shared") {
    names(sort(rowMeans(m_cent[names(test_groups[[which(names(test_groups) == x)]]), colnames(m_cent)[grep("OSM.Rep|TGFB", colnames(m_cent))]]), decreasing = TRUE))
  } else {
    names(sort(rowMeans(m_cent[names(test_groups[[which(names(test_groups) == x)]]), colnames(m_cent)[grep(x, colnames(m_cent))]]), decreasing = TRUE))
  }
})
row_ord <- rev(unlist(test[c("Shared", "OSM.Rep", "TGFB")], use.names = FALSE))

plot_df <- m_cent[rownames(m_cent) %in% row_ord, col_ord] %>% 
  tibble::rownames_to_column("Gene") %>% reshape2::melt() %>% 
  dplyr::mutate(variable = factor(variable, levels = col_ord),
                Gene = factor(Gene, levels = row_ord))
anno_df <- metadata %>% dplyr::filter(Sample %in% plot_df$variable) %>% dplyr::arrange(match(Sample, levels(plot_df$variable))) %>% 
  dplyr::select(Sample:Tech_Rep) %>% dplyr::rename(Replicate = Tech_Rep) %>% 
  tidyr::pivot_longer(cols = c("Cell_Line", "Condition", "Replicate")) %>% dplyr::mutate(Sample = factor(Sample, levels = levels(plot_df$variable)),
                                                                                      name = factor(name, levels = c("Replicate", "Cell_Line", "Condition")))
set_cols <- setNames(c("lightgrey", "#0c2a50ff", "#D2aF81FF", "#5A9599FF", "#FAFD7CFF", "#B7E4F9FF", "#FB6467FF"), c("Rep1", "Rep2", "93VU", "JHU6", "CTR", "OSM", "TGFB"))
p <- gmap(plot_df, x = variable, y = Gene, fill = value, limits = c(-2, 2), midpoint = 0, ratio = NULL, y.name = "",
          axis.rel = 1.5, legend.title.rel = 1.25, legend.rel = 1.0, legend.height = 5, legend.width = 0.6, x.num = F, y.num = F, magma_pal = FALSE) +
  geom_hline(yintercept = c(which(levels(plot_df$Gene) == "ABCA1") - 0.5, which(levels(plot_df$Gene) == "SERPINE1") + 0.5), size = 1) +
  geom_vline(xintercept = c(4.5, 8.5), size = 0.2) +
  theme(axis.text = element_text(size = 16), axis.ticks.x = element_blank(), axis.text.x = element_blank())
p_anno <- ggplot(anno_df, aes(x = Sample, y = name, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = set_cols) +
  labs(x = "", y = "") +
  geom_vline(xintercept = c(4.5, 8.5), size = 0.2) +
  theme_bw(base_size = 16) +
  theme(aspect.ratio = NULL, panel.grid = element_blank(), legend.position = "right", axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.text = element_text(size = 16, face = "bold")) +
  scale_y_discrete(expand = c(0, 0), labels = function(x) gsub("_", " ", x, fixed = TRUE)) + scale_x_discrete(expand = c(0, 0))
p1 <- p_anno + plot_spacer() + p + patchwork::plot_layout(nrow = 3, heights = c(0.1, -0.05, 1), guides = "collect")
ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_5g.pdf"), plot = p1, device = "pdf", dpi = 300, width = 8, height = 14)

de_meta <- de_meta %>% dplyr::mutate(Path = ifelse(Gene %in% tgfb_emt_genes, "EMT", NA))
p <- ggplot(data = de_meta, aes(x = Lfc, y = -log10(Pval), label = Label)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') +
  geom_point(shape = 21,
             fill = dplyr::case_when(de_meta$Signif_DE == "Not Significant" ~ alpha("grey30", 0.2),
                                     de_meta$Signif_DE == "Upregulated by TGFB" ~ "#b2182b",
                                     de_meta$Signif_DE == "Upregulated in the control" ~ "#2166ac"),
             color = dplyr::case_when(is.na(de_meta$Path) ~ alpha("grey30", 0.5),
                                      de_meta$Path == "EMT" ~ state_cols["pEMT"]),
             stroke = dplyr::case_when(is.na(de_meta$Label) ~ 0.8,
                                       is.na(de_meta$Path) & !is.na(de_meta$Label) ~ 1,
                                       !is.na(de_meta$Path) & !is.na(de_meta$Label) ~ 1.5),
             size = dplyr::case_when(is.na(de_meta$Label) ~ 1,
                                     is.na(de_meta$Path) & !is.na(de_meta$Label) ~ 2,
                                     !is.na(de_meta$Path) & !is.na(de_meta$Label) ~ 3),
             alpha = ifelse(is.na(de_meta$Label), 0.2, 1)) + 
  ggrepel::geom_label_repel(data = de_meta[!is.na(de_meta$Path) & !is.na(de_meta$Label), ], aes(fill = Path), max.overlaps = Inf, color = "black", show.legend = FALSE, size = 5) +
  scale_fill_manual(name = "Pathway", values = setNames(c(alpha("#117733", 0.6), alpha(unname(state_cols[c("pEMT", "Inflammatory")]), 0.8)), c("Both", "EMT", "Inflammatory")), na.value = alpha("grey30", 0.5)) +
  theme_classic() + labs(x = expression("log"[2]*"[Fold Change]"), y = expression("-log"[10]*"[P-value]")) +
  scalop::theme_scalop(legend.position = "right") + scale_y_continuous(expand = c(0, 0)) + expand_limits(y = max(-log10(de_meta$Pval)) + 0.15) + scale_x_continuous(expand = expansion(add = 0.15))
ggplot2::ggsave(filename = here("results/Diagnostic_Plots/Bulk_RNAseq_TGFB_vs_Control_DEG_VolcanoPlot.pdf"), plot = p, device = "pdf", dpi = 300, width = 12, height = 12)

# Another version of the plot
de_meta <- de_meta %>% dplyr::mutate(Select_genes = ifelse(Gene %in% c("VIM", "TGFBI", "SNAI2", "LAMC2", "SERPINE1", "MMP2", "COL1A1", "INHBA", "TAGLN", "FSTL3", "PMEPA1"), Gene, NA),
                                     nudge = ifelse(!is.na(Select_genes), 3.75 - Lfc, NA))
p <- ggplot(data = de_meta, aes(x = Lfc, y = -log10(Pval), label = Label)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') +
  geom_point(shape = 21,
             fill = dplyr::case_when(de_meta$Signif_DE == "Not Significant" ~ alpha("grey30", 0.2),
                                     de_meta$Signif_DE != "Not Significant" & is.na(de_meta$Path) ~ "#0c2a50ff",
                                     de_meta$Signif_DE != "Not Significant" & de_meta$Path == "EMT" ~ state_cols["pEMT"]),
             color = dplyr::case_when(de_meta$Signif_DE == "Not Significant" ~ alpha("grey30", 0.5),
                                      de_meta$Signif_DE != "Not Significant" ~ "black"),
             stroke = dplyr::case_when(is.na(de_meta$Label) ~ 0.8,
                                       !is.na(de_meta$Label) ~ 1),
             size = dplyr::case_when(is.na(de_meta$Label) ~ 1,
                                     is.na(de_meta$Path) & !is.na(de_meta$Label) ~ 2,
                                     !is.na(de_meta$Path) & !is.na(de_meta$Label) ~ 3),
             alpha = ifelse(is.na(de_meta$Label), 0.2, 1)) + 
  ggrepel::geom_label_repel(data = de_meta[!is.na(de_meta$Path), ], aes(fill = Path), max.overlaps = Inf, color = "black", show.legend = FALSE, size = 5,
                            nudge_x = de_meta[!is.na(de_meta$Path), "nudge"], direction = "y", hjust = 1) +
  scale_fill_manual(name = "Pathway", values = setNames(c(alpha("#117733", 0.6), alpha(unname(state_cols[c("pEMT", "Inflammatory")]), 0.8)), c("Both", "EMT", "Inflammatory")), na.value = alpha("grey30", 0.5)) +
  theme_classic() + labs(x = expression("Response to TGFB (log"[2]*"[Fold Change] average of two cell lines - JHU6 & 93VU)"), y = expression("-log"[10]*"[P-value]")) +
  scalop::theme_scalop(legend.position = "right") + scale_y_continuous(expand = c(0, 0)) + expand_limits(y = max(-log10(de_meta$Pval)) + 0.15) + 
  scale_x_continuous(expand = expansion(add = 0.15), limits = c(min(de_meta$Lfc[!is.na(de_meta$Label)]), max(de_meta$Lfc[!is.na(de_meta$Label)]))) + theme(aspect.ratio = 0.6)
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_7a_top.pdf"), plot = p, device = "pdf", dpi = 300, width = 14, height = 12)

# Save TGFB vs CTR DEGs to csv file
de_df <- de_meta %>% dplyr::filter(Signif_DE != "Not Significant") %>% 
  dplyr::arrange(desc(Lfc)) %>% dplyr::select(Gene, Lfc, Pval, Signif_DE) %>% dplyr::rename("Group" = "Signif_DE")
write.csv(de_df, file = here("results/Generated_Data/Bulk_RNAseq_TGFB_vs_Control_DEGs.csv"), row.names = FALSE)




### IL1B vs CONTROL
upregulated_group <- c("93VU.IL1B.Rep1", "93VU.IL1B.Rep2", "JHU6.IL1B.Rep1", "JHU6.IL1B.Rep2")
downregulated_group <- c("93VU.CTR.Rep1", "93VU.CTR.Rep2", "JHU6.CTR.Rep1", "JHU6.CTR.Rep2")

m.k1 <- as.matrix(m_cent[, upregulated_group])
m.other <- as.matrix(m_cent[, downregulated_group])
genes <- rownames(m_cent)
lfc <- rowMeans(m.k1) - rowMeans(m.other)
pval <- sapply(rownames(m_cent), function(gene) {
  t.test(m.k1[gene, ], m.other[gene, ], alternative = 'two.sided', paired = F, var.equal = T)$p.value
})
padj <- p.adjust(pval, method = "BH")
ord <- order(lfc, decreasing = T)
de_meta <- data.frame(Lfc = lfc,
                      Pval = pval,
                      Padj = padj,
                      Gene = names(lfc))

de_meta <- de_meta %>% dplyr::mutate(Signif_DE = dplyr::case_when(Pval < 0.01 & Lfc >= 1 ~ "Upregulated by IL1B", 
                                                                  Pval < 0.01 & Lfc <= -1 ~ "Upregulated in the control", .default = "Not Significant"),
                                     Label = ifelse(Pval < 0.01 & abs(Lfc) >= 1, Gene, NA))
p <- ggplot(data = de_meta, aes(x = Lfc, y = -log10(Pval), col = Signif_DE, label = Label)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(name = "Group", values = setNames(c(alpha("grey30", 0.2), "#2166ac", "#b2182b"), c("Not Significant", "Upregulated in the control", "Upregulated by IL1B"))) +
  theme_classic() + labs(x = expression("log"[2]*"[Fold Change]"), y = expression("-log"[10]*"[P-value]")) +
  scalop::theme_scalop(legend.position = "right") + scale_y_continuous(expand = c(0, 0)) + expand_limits(y = max(-log10(de_meta$Pval)) + 0.15) + scale_x_continuous(expand = expansion(add = 0.15))
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_7a_bottom_right.pdf"), plot = p, device = "pdf", dpi = 300, width = 12, height = 12)




### SPP1 vs CONTROL
upregulated_group <- c("93VU.SPP1.Rep1", "93VU.SPP1.Rep2", "JHU6.SPP1.Rep1", "JHU6.SPP1.Rep2")
downregulated_group <- c("93VU.CTR.Rep1", "93VU.CTR.Rep2", "JHU6.CTR.Rep1", "JHU6.CTR.Rep2")

m.k1 <- as.matrix(m_cent[, upregulated_group])
m.other <- as.matrix(m_cent[, downregulated_group])
genes <- rownames(m_cent)
lfc <- rowMeans(m.k1) - rowMeans(m.other)
pval <- sapply(rownames(m_cent), function(gene) {
  t.test(m.k1[gene, ], m.other[gene, ], alternative = 'two.sided', paired = F, var.equal = T)$p.value
})
padj <- p.adjust(pval, method = "BH")
ord <- order(lfc, decreasing = T)
de_meta <- data.frame(Lfc = lfc,
                      Pval = pval,
                      Padj = padj,
                      Gene = names(lfc))

de_meta <- de_meta %>% dplyr::mutate(Signif_DE = dplyr::case_when(Pval < 0.01 & Lfc >= 1 ~ "Upregulated by SPP1", 
                                                                  Pval < 0.01 & Lfc <= -1 ~ "Upregulated in the control", .default = "Not Significant"),
                                     Label = ifelse(Pval < 0.01 & abs(Lfc) >= 1, Gene, NA))
p <- ggplot(data = de_meta, aes(x = Lfc, y = -log10(Pval), col = Signif_DE, label = Label)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(name = "Group", values = setNames(c(alpha("grey30", 0.2), "#2166ac", "#b2182b"), c("Not Significant", "Upregulated in the control", "Upregulated by IL1B"))) +
  theme_classic() + labs(x = expression("log"[2]*"[Fold Change]"), y = expression("-log"[10]*"[P-value]")) +
  scalop::theme_scalop(legend.position = "right") + scale_y_continuous(expand = c(0, 0)) + expand_limits(y = max(-log10(de_meta$Pval)) + 0.15) + scale_x_continuous(expand = expansion(add = 0.15))
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_7a_bottom_left.pdf"), plot = p, device = "pdf", dpi = 300, width = 12, height = 12)

