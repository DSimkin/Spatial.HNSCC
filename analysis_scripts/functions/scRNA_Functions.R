library(tidyverse)
library(reshape2)
library(scales)
library(patchwork)
library(here)
library(Matrix)
library(biomaRt)
library(S4Vectors)
library(ComplexHeatmap)

source(here("analysis_scripts/functions/gmap_discrete_binary.R"))
# Sys.setenv(RETICULATE_PYTHON = "/Users/dorsi/.conda/envs/leiden/bin/python")    ### DON'T FORGET TO FIRST ACTIVATE `leiden` CONDA ENV: on terminal - `conda activate leiden`
library(reticulate)


# QC tables and functions --------------------------------------------------

compute_housekeeping <- function(m, return_sorted = FALSE, cell_subset = NULL, custom_hk = NULL, log = TRUE, custom_gene_align = FALSE, verbose = FALSE){
  stopifnot(is.logical(return_sorted), (is.null(cell_subset) | is.vector(cell_subset)), (is.null(custom_hk) | is.vector(custom_hk)), is.logical(log))

  if(!is.null(cell_subset))
    m <- m[, cell_subset]

  if(isTRUE(custom_gene_align)){
    genes <- unname(unlist(sapply(rownames(m), function(x) stringr::str_split(x, pattern = "_", simplify = TRUE)[-1])))
    prefix <- unname(unlist(sapply(rownames(m), function(x) stringr::str_split(x, pattern = "_", simplify = TRUE)[1])))
  } else{
    genes <- rownames(m)
  }

  HOUSEKEEPING_GENES_LIST <- c("ACTB", "B2M", "HNRPLL", "HPRT", "PSMB2", "PSMB4", "PPIA", "PRPS1", "PRPS1L1", "PRPS1L3", "PRPS2", "PRPSAP1", "PRPSAP2",
                               "RPL10", "RPL10A", "RPL10L", "RPL11", "RPL12", "RPL13", "RPL14", "RPL15", "RPL17", "RPL18", "RPL19", "RPL21", "RPL22",
                               "RPL22L1", "RPL23", "RPL24", "RPL26", "RPL27", "RPL28", "RPL29", "RPL3", "RPL30", "RPL32", "RPL34", "RPL35", "RPL36", "RPL37",
                               "RPL38", "RPL39", "RPL39L", "RPL3L", "RPL4", "RPL41", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL7L1", "RPL8", "RPL9", "RPLP0",
                               "RPLP1", "RPLP2", "RPS10", "RPS11", "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", "RPS17", "RPS18", "RPS19", "RPS20",
                               "RPS21", "RPS24", "RPS25", "RPS26", "RPS27", "RPS27A", "RPS27L", "RPS28", "RPS29", "RPS3", "RPS3A", "RPS4X", "RPS5", "RPS6",
                               "RPS6KA1", "RPS6KA2", "RPS6KA3", "RPS6KA4", "RPS6KA5", "RPS6KA6", "RPS6KB1", "RPS6KB2", "RPS6KC1", "RPS6KL1", "RPS7", "RPS8",
                               "RPS9", "RPSA", "TRPS1", "UBB")
  if(is.null(custom_hk))
    g <- genes[genes %in% HOUSEKEEPING_GENES_LIST]
  else
    g <- genes[genes %in% custom_hk]

  if(isTRUE(custom_gene_align)){
    g <- paste(prefix, g, sep = "_")
    h <- base::apply(m[rownames(m) %in% g, ], 2, mean)
  } else {
    h <- base::apply(m[g, ], 2, mean)
  }

  if(isTRUE(return_sorted)){
    h <- base::sort(h)
  }

  if(isTRUE(log))
    h <- base::log2(h + 1)

  return (h)
}


percent_features <- function(m, pattern = NULL, features = NULL, custom_gene_align = FALSE){
  if(!is.null(features) && !is.null(pattern)){
    warning("Both pattern and features provided. Pattern is being ignored.")
  }

  if(isTRUE(custom_gene_align)){
    genes <- unname(unlist(sapply(rownames(m), function(x) stringr::str_split(x, pattern = "_", simplify = TRUE)[-1])))
    prefix <- unname(unlist(sapply(rownames(m), function(x) stringr::str_split(x, pattern = "_", simplify = TRUE)[1])))
  } else{
    genes <- rownames(m)
  }

  features <- features %||% grep(pattern = pattern, x = genes, value = TRUE)
  if(isTRUE(custom_gene_align)){
    features <- paste(prefix, features, sep = "_")
    percent.featureset <- colSums(m[rownames(m) %in% features, , drop = FALSE]) / colSums(m) * 100
  } else {
    percent.featureset <- colSums(m[features, , drop = FALSE]) / colSums(m) * 100
  }

  return(percent.featureset)
}


cell_qc <- function(m, log = TRUE, ...){
  stopifnot(is.logical(log))

  res <- tibble::tibble(CellID = colnames(m))
  res$Complexity <- colSums(m != 0)
  res$HK <- compute_housekeeping(m, log = log, ...)
  res$Percent_mt <- percent_features(m, pattern = "^MT-", ...)

  return (res)
}


gene_qc <- function(m, scaling_factor = 10, pseudo_count = 1, log_base = 2, norm = c("umi", "tpm", "cpm")){
  if(length(norm) > 1){
    warning("Please choose only one data normalization method to start with.")
  }

  if(norm == "cpm" | norm == "tpm"){m <- as.matrix(m)}
  if(norm == "umi"){
    # UMI to CPM
    scaling_factor <- 1000000/colSums(m)
    m <- sweep(m, MARGIN = 2, STATS = scaling_factor, FUN = "*")
  }

  res <- tibble::tibble(Gene = rownames(m),
                        logMean = base::log2(base::rowMeans(m) + 1))
  m <- base::log(x = m / scaling_factor + pseudo_count, base = log_base)
  res$Mean <- base::rowMeans(m)
  res$Variance <- matrixStats::rowVars(m)
  res$SD <- base::apply(m, 1, stats::sd)

  return (res)
}

.sparsity <- function(x) length(which(x == 0)) / (dim(x)[1] * dim(x)[2]) * 100

qc_stats <- function(m_preQC, m_postQC){

  cells_pre_qc <- ncol(m_preQC)
  cells_post_qc <- ncol(m_postQC)
  genes_pre_qc <- nrow(m_preQC)
  genes_post_qc <- nrow(m_postQC)

  DataFrame("Cells pre-QC" = cells_pre_qc,
            "Cells post-QC" = cells_post_qc,
            "Cells dropped" = sprintf("%.2f%%", (1 - (cells_post_qc / cells_pre_qc)) * 100),
            "Genes pre-QC" = genes_pre_qc,
            "Genes post-QC" = genes_post_qc,
            "Genes dropped" = sprintf("%.2f%%", (1 - (genes_post_qc / genes_pre_qc)) * 100),
            "Sparsity pre QC" = sprintf("%.2f%%",.sparsity(m_preQC)),
            "Sparsity post QC" = sprintf("%.2f%%",.sparsity(m_postQC)),
            row.names = unique(cell_QC$Sample))
}



# Utility functions -------------------------------------------------------

untable <- function(x) {
  stopifnot(is.table(x))
  obs <- as.data.frame(x)[rep(1:prod(dim(x)),c(x)),-length(dim(x))-1]
  rownames(obs) <- NULL; obs
}

Melt <- function(x) {
  if (!is.data.frame(x = x)) {
    x <- as.data.frame(x = x)
  }
  return(data.frame(
    rows = rep.int(x = rownames(x = x), times = ncol(x = x)),
    cols = unlist(x = lapply(X = colnames(x = x), FUN = rep.int, times = nrow(x = x))),
    vals = unlist(x = x, use.names = FALSE)
  ))
}


Unlist = function(L, nested.names = FALSE) {
  if (nested.names) {
    Names = unlist(sapply(L, names), use.names = F)
  } else {
    Names = rep(names(L), lengths(L))
  }
  stats::setNames(unlist(L), Names)
}

mat2list <- function(x) {
  stopifnot(is.matrix(x))
  lapply(seq_len(ncol(x)), function(i) x[, i])
}

cbind.fill <- function(...){
  nm <- list(...)
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow))
  do.call(cbind, lapply(nm, function (x)
    rbind(x, matrix(, n-nrow(x), ncol(x)))))
}

grouped_rowcenter = function(m, groups, by = 'mean') {
  res <- sapply(groups, function(cells) {
    scalop::rowcenter(m[, cells], by = by)},
    simplify = F)
  Names <- Unlist(sapply(res, colnames, simplify = F))
  res = do.call(cbind.data.frame, res)
  colnames(res) <- Names
  res <- as.matrix(res)
  rownames(res) <- rownames(m)
  return(res)
}


get_top_genes <- function(expr_mat,
                          cells = colnames(expr_mat),
                          n_top = 50,
                          rm_ribo_mito = TRUE) {

  ordered_mat <- expr_mat[order(rowMeans(expr_mat[, cells]), decreasing = TRUE), cells]
  if(isTRUE(rm_ribo_mito)) {
    ordered_mat <- ordered_mat[-grep("^MT-|^RPL|^RPS", rownames(ordered_mat)), ]
  }
  head(setNames(round(rowMeans(ordered_mat), 2), rownames(ordered_mat)), n_top)
}


ident_by_feats <- function(matrix, features, cells = colnames(matrix), umi_thresh = 1, filter_thresh = 0.5, return_freqs = FALSE) {
  feat_expr_cell_ids <- lapply(features, function(feat) matrix[feat, cells] %>% .[. >= umi_thresh] %>% names)
  collect_ids <- do.call(c, feat_expr_cell_ids)
  freq_of_intersecting_feats <- table(collect_ids) %>% as.data.frame(stringsAsFactors = FALSE) %>%
    dplyr::mutate(Freq = Freq / length(features)) %>%
    dplyr::filter(Freq >= filter_thresh)

  if(return_freqs) {
    return(freq_of_intersecting_feats)
  } else {
    ids <- freq_of_intersecting_feats %>% dplyr::pull(collect_ids)
    return(ids)
  }
}


# Use UMI or Log2-CPM expression matrix and vector of marker genes to classify cells/spots to those containing each combination of these markers
class_by_markers <- function(expr_mat,
                             markers,
                             filter_thresh = 1,
                             ...) {

  # Spots expressing all markers
  out_list <- list(ident_by_feats(expr_mat, features = markers, filter_thresh = filter_thresh, ...))
  names(out_list) <- str_c(paste0(markers, "+"), collapse = ", ")

  # Spots expressing all but one marker
  chk_n_m1 <- lapply(seq_along(markers), function(i) {
    n_m1 <- ident_by_feats(expr_mat, features = markers[-i], filter_thresh = filter_thresh, ...)
    n_m1_uniq <- n_m1[!n_m1 %in% unlist(out_list, use.names = FALSE)]
  })
  names(chk_n_m1) <- sapply(seq_along(markers), function(i) str_c(paste0(markers[-i], "+"), collapse = ", "))
  chk_n_m1 <- chk_n_m1[lapply(chk_n_m1, length) > 0]
  out_list <- c(out_list, chk_n_m1)

  for(iter in seq_len(length(markers) - 2)) {
    group_markers <- as.list(as.data.frame(combn(markers, 1 + iter)))
    chk_n_miter <- lapply(seq_along(group_markers), function(i) {
      n_miter <- ident_by_feats(expr_mat, features = markers[!markers %in% group_markers[[i]]], filter_thresh = filter_thresh, ...)
      n_miter_uniq <- n_miter[!n_miter %in% unlist(out_list, use.names = FALSE)]
    })
    names(chk_n_miter) <- sapply(seq_along(group_markers), function(i) str_c(paste0(markers[!markers %in% group_markers[[i]]], "+"), collapse = ", "))
    chk_n_miter <- chk_n_miter[lapply(chk_n_miter, length) > 0]
    out_list <- c(out_list, chk_n_miter)
  }

  return(out_list)
}



# 10X Preprocessing -------------------------------------------------------


filter_10x <- function(matrix,
                       cells = colnames(matrix),
                       complexity_thresh = 1000,
                       percent_mt_thresh = 20,
                       dbl_remove = FALSE,
                       umi_thresh = 5,
                       cells_thresh = 20,
                       aver_log_thresh = 4,
                       var_thresh = 0.2,
                       centre = "mean",
                       gene_filt_method = c("exp", "var", "none"),
                       merge_dup_genes = FALSE,
                       merge_method = c("sum", "mean", "max", "median"),
                       ref = NULL,
                       log = TRUE,
                       raw_to_CPM = TRUE){
  m <- matrix[, colnames(matrix) %in% cells]

  # Deal with duplicated gene names in the expression matrix
  if(isTRUE(merge_dup_genes) & sum(duplicated(rownames(m))) > 0){
    m <- fix_dup_genes(m, stat = merge_method)
  }

  # Filter cells by complexity cutoff
  if(!is.null(complexity_thresh)){
    complexity <- colSums(m != 0)
    m <- m[, complexity >= complexity_thresh]
  } else if(is.null(complexity_thresh)){m <- m}

  # Filter out doublets if specified in the arrguments
  if(isFALSE(dbl_remove)){m <- m}
  if(isTRUE(dbl_remove)){
    set.seed(5045410)
    d1 <- detect_doublets(m, method = "scdblfinder", scale = 0.006)
    d2 <- detect_doublets(m, method = "scds", scale = 0.006)
    d3 <- detect_doublets(m, method = "scran", scale = 0.006)

    dbl <- cbind(d1, d2, d3)[, c(2, 4, 6)]
    dbl <- ifelse(dbl == "doublet", 1, 0)
    dbl_sum <- apply(dbl, 1, sum)       # table(dbl_sum)
    dbl <- ifelse(dbl_sum > 1, "doublet", "singlet")
    m <- m[, dbl == "singlet"]
  }

  # Filter cells with >20% mitochondrial genes
  if(!is.null(percent_mt_thresh)){
    per_mt <- percent_features(m, pattern = "^MT-")
    m <- m[, per_mt < percent_mt_thresh]
  } else if(is.null(percent_mt_thresh)){m <- m}

  # Return Matrix after removing doublets and filtering low complexity cells
  if(gene_filt_method == "none" & centre == "none" & isFALSE(log)){return(m)}

  # UMI to CPM
  if(isTRUE(raw_to_CPM)) {
    scaling_factor <- 1000000/colSums(m)
    m_CPM <- sweep(m, MARGIN = 2, STATS = scaling_factor, FUN = "*")
  } else {m_CPM <- m}

  # Filter genes with mean log expression > 4 OR passing UMI threshold
  if(gene_filt_method == "none"){m <- m_CPM}
  if(gene_filt_method == "exp"){
    sum_gene_exp <- rowSums(m >= umi_thresh)
    m <- m_CPM[log2(rowMeans(m_CPM) + 1) > aver_log_thresh | sum_gene_exp >= cells_thresh, ]
  }
  if(gene_filt_method == "var"){
    m_tmp <- base::log(x = m_CPM / 10 + 1, base = 2)
    vari <- matrixStats::rowVars(m_tmp)
    m <- m_CPM[log2(rowMeans(m_CPM) + 1) > aver_log_thresh & vari > var_thresh, ]
  }

  # Log2 transformation
  if(isFALSE(log)){m <- m}
  if(isTRUE(log)){m <- log2(1 + (m/10))}

  # Centering to get relative expression
  if(centre == "mean"){avg <- rowMeans(m)}
  if(centre == "median"){avg <- apply(m, 1, median)}
  if(centre == "median" | centre == "mean"){m_cent <- sweep(m, 1, avg)}
  if(centre == "ref"){centmat <- sweep(m, 1, ref[names(ref) %in% rownames(m)])}
  if(centre != "none"){return(m_cent)}
  if(centre == "none"){return(m)}
}


##### log and un_log functions:
log_norm <- function(x, centre = FALSE){
  logged <- log2(1 + (x/10))
  if(isTRUE(centre)){
    avg <- rowMeans(logged)
    logged <- sweep(logged, 1, avg)
  }
  return(logged)
}

un_log <- function(x){
  tpmed <- 10 * (2^x - 1)
  return(tpmed)
}

#Perform hierarchical clustering with average distance and reorder correlation matrix
order_Hclust <- function(x){
  d <- as.dist((1 - x) / 2)
  hc <- hclust(d, method = "average")
  x <- x[hc$order, hc$order]
}

minMax_normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}

### summary statistics for dealing with duplicated gene names in exprssion matrix
fix_dup_genes <- function(m, stat = c("sum", "mean", "max", "median")){

  fix_dup_gene  <- function(gene){
    gene_mat <- m[rownames(m) == gene, ]
    statFUN(gene_mat)
  }

  statFUN <- switch(match.arg(stat),
                   sum = matrixStats::colSums2,
                   mean = matrixStats::colMeans2,
                   max = matrixStats::colMaxs,
                   median = matrixStats::colMedians)

  # gene order, keeping only first instance of duplicated genes
  gene_ord <- rownames(m)[!duplicated(rownames(m))]
  # find duplicated genes
  dup_genes <- unique(rownames(m)[duplicated(rownames(m))])
  # for each duplicated gene, apply summary statistic and return numeric vector
  # rbind duplicate gene vectors
  new_mat <- do.call(rbind, sapply(dup_genes, fix_dup_gene, simplify = F))
  m <- m[!rownames(m) %in% dup_genes, ]
  m <- rbind(m, new_mat)
  m[gene_ord, ]
}



# Doublets detection wrapper function -------------------------------------


detect_doublets <- function(matrix,
                            method = "scdblfinder",
                            scransubset = rownames(matrix),
                            doubleprop = NULL,
                            scale = 0.004,
                            sctrans = F, scpc = 10){

  if(is.null(doubleprop)){doubleprop <- scale * ncol(matrix) / 500}

  if(method != "seurat"){
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = matrix))
  }

  if(method == "scdblfinder"){
    sce <- scDblFinder::scDblFinder(sce, dbr = doubleprop)
    out <- cbind.data.frame(sce$scDblFinder.score, sce$scDblFinder.class)
  }

  if(method == "scds"){
    sce <- scds::cxds_bcds_hybrid(sce)
    sce$bin <- ifelse(sce$hybrid_score > quantile(sce$hybrid_score, probs = 1-doubleprop), "doublet", "singlet")
    out <- cbind.data.frame(sce$hybrid_score, sce$bin)
  }

  if(method == "scran"){
    scd <- scDblFinder::computeDoubletDensity(sce, subset.row = intersect(scransubset, rownames(matrix)))
    scd <- log2(scd+1)
    sce$scran_dblscore <- scd
    sce$bin <- ifelse(sce$scran_dblscore > quantile(sce$scran_dblscore, probs = 1-doubleprop), "doublet", "singlet")
    out <- cbind.data.frame(sce$scran_dblscore, sce$bin)
  }

  if(method == "seurat"){
    library("Seurat")
    library("DoubletFinder")
    library("sctransform")
    seu <- CreateSeuratObject(counts = matrix, project = "x", min.cells = 20, min.features = 1000)
    if(isFALSE(sctrans)){
      seu <- NormalizeData(seu)
      seu <- ScaleData(seu)
      seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
    }
    if(isTRUE(sctrans)){
      seu <- SCTransform(seu)
    }
    seu <- RunPCA(seu)
    seu <- FindNeighbors(seu, reduction = "pca", dims = 1:scpc)
    seu <- FindClusters(seu, resolution = 0.5)
    sweeplist <- paramSweep_v3(seu, PCs = 1:scpc, sct = sctrans)
    sweepstats <- summarizeSweep(sweeplist, GT = FALSE)
    bcmvn <- find.pK(sweepstats)
    pk <- as.character(bcmvn[which.max(bcmvn$BCmetric),"pK"])
    pk <- as.numeric(pk)
    homotypic.prop <- modelHomotypic(seu@meta.data$seurat_clusters)
    nExp_poi <- round(doubleprop * ncol(seu))
    nExp_poi.adj <- round(nExp_poi * (1-homotypic.prop))
    seu <- doubletFinder_v3(seu, PCs = 1:scpc, pN = 0.25, pK = pk, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = sctrans)
    n1 <- ncol(seu@meta.data)
    out <- seu@meta.data[, (n1-1):n1]
  }

  colnames(out) <- c("DoubletScore", "DoubletClass")
  rownames(out) <- colnames(matrix)
  return(out)
}


### Remove doublets
rm_doublets <- function(m) {
  set.seed(123)
  d1 <- detect_doublets(m, method = "scdblfinder", scale = 0.006)
  d2 <- detect_doublets(m, method = "scds", scale = 0.006)
  d3 <- detect_doublets(m, method = "scran", scale = 0.006)

  dbl <- cbind(d1, d2, d3)[, c(2, 4, 6)]
  dbl <- ifelse(dbl == "doublet", 1, 0)
  dbl_sum <- apply(dbl, 1, sum)       # table(dbl_sum)
  dbl <- ifelse(dbl_sum > 1, "doublet", "singlet")
  m <- m[, dbl == "singlet"]
}



# Dimensionality Reduction and Clustering ------------------------------------------------

standard_umap<-function(matrix,metric="correlation",
                        n_neighbors=50,spread=10,min_dist=0.01,var_genes=NULL,...){
  if(!is.null(var_genes)){
    s1<-apply(matrix,1,sd)
    n1<-names(tail(sort(s1),n=var_genes))
    matrix<-matrix[n1,]
  }


  tm2<-t(as.matrix(matrix))
  um2<-uwot::umap(tm2,metric=metric,n_neighbors=n_neighbors,
                  spread=spread,min_dist=min_dist,
                  pca_center = FALSE, fast_sgd = TRUE,n_threads = 100,...)
  colnames(um2)<-c("V1","V2")
  rownames(um2)<-rownames(tm2)
  dtsne22<-as.data.frame(um2)
  return(dtsne22)
}

cluster_louvain <- function(data, k) {
  knn <- FNN::get.knn(as.matrix(data), k = k)
  knn <- data.frame(from = rep(1:nrow(knn$nn.index), k), to = as.vector(knn$nn.index), weight = 1 / (1 + as.vector(knn$nn.dist)))
  simple_igraph <- igraph::simplify(igraph::graph_from_data_frame(knn, directed = FALSE))
  louvain_clusts <- igraph::cluster_louvain(simple_igraph)
  clusters <- igraph::membership(louvain_clusts)
  return (clusters)
}

### Extract clusters by DBScan, finding the knee point of knn distance graph
cluster_dbscan <- function(matrix, min_pts = log(nrow(matrix)), probs = 3/4, diff_cut = 10){
  # Calculate KNN distance
  dist <- dbscan::kNNdist(matrix, min_pts)
  # Find knee point through derivative of KNN graph
  dist_cutoff <- quantile(sort(dist), probs = probs)
  dist <- dist[dist > dist_cutoff]
  # Order and scale result
  dist <- dist[order(dist)]
  dist <- dist / max(dist)
  # Derivative
  der_dist <- diff(dist) / (1 / length(dist))
  # Get first point where derivative is higher than 1
  knee <- dist[length(der_dist) - length(der_dist[der_dist > diff_cut])] + dist_cutoff
  # Use knee point in plot as `eps` for dbscan algorithm
  db_run <- dbscan::dbscan(matrix, eps = knee, minPts = min_pts)
  # Extract clusters
  clusters <- as.factor(db_run$cluster)
  return(clusters)
}


### Input should be either embeddings (such as PCA matrix) or centered (!!) expression matrix on which PCA will be done
cluster_leiden <- function(matrix,
                           red_dim = TRUE,
                           dims = 1:30,
                           resolution = 0.4,
                           clust_k = 20,
                           iter_n = 10,
                           verbose = FALSE,
                           ...) {
  if (is.null(dim(matrix))) {
    stop("Input is not a matrix!!")
  }
  if(ncol(matrix) > 150 & isTRUE(red_dim)) {
    # Perform dimensionality reduction first
    pca <- irlba::prcomp_irlba(t(matrix), center = FALSE, scale = TRUE, n = dims, retx = TRUE)
    attr(pca$rotation, "dimnames")[[1]] <- rownames(matrix)
    attr(pca$x, "dimnames")[[1]] <- colnames(matrix)
    matrix <- pca$x
  }
  if(!is.null(dims)) {mat <- matrix[, dims]}
  if(is.null(dims) & isFALSE(red_dim)) {mat <- matrix}
  graphs <- Seurat::FindNeighbors(mat, k.param = clust_k, verbose = verbose, ...)
  if(isTRUE(verbose)){message("Clustering with Leiden")}
  clusters <- Seurat::FindClusters(graphs$snn, algorithm = 4, resolution = resolution, n.iter = iter_n, verbose = verbose, ...)
  out <- clusters[, grep("res.", colnames(clusters))]
}

# Cluster coordinate matrix with cells as rows
cluster_coord <- function(matrix, method = c("dbscan", "louvain"), louvain = 25, ...){
  if(method == "dbscan"){
    clusters <- cluster_dbscan(matrix, ...)
  }
  if(method == "louvain"){
    clusters <- cluster_louvain(matrix, k = louvain)
  }

  rename_clusts <- function(x){paste("cluster", x, sep = "")}
  clusters <- rename_clusts(clusters)
  names(clusters) <- rownames(matrix)
  return(clusters)
}

is_hpv <- function(m){
  if(length(rownames(m)[grep("^hpv", rownames(m))]) == 0){
    print("No HPV genes are found")
  } else{
    table(substri(rownames(m), pos = 1)) %>% .[grep("^hpv", names(.))]
  }
}



# Differential Expression -------------------------------------------------


dea <- function(x, g1, g2, name) {
  de_res <- tibble(Gene = rownames(x),
                   log2FC = rowMeans(x[, g1]) - rowMeans(x[, g2]),
                   p.val = sapply(1:nrow(x), function(i) wilcox.test(x = x[i, g1], y = x[i, g2])$p.value),
                   p.adj = p.adjust(p.val, "fdr"),
                   Name = name)
  de_res
}

library(BiocParallel)
plapply <- function(X, FUN, ..., pack = F, profile = F, n_workers = 15) {

  start_time <- Sys.time()
  res <- bplapply(X = X, FUN = FUN, ... = ..., BPPARAM = MulticoreParam(workers = n_workers))
  end_time <- Sys.time()
  end_time - start_time
  # message(end_time - start_time)
  print(end_time - start_time)

  if(isTRUE(pack))
    res <- unlist(lapply(res, function(x) x), recursive = F)

  res
}

btw_samples_dea <- function(samp1, samp2, name) {
  # TODO: should add a condition to narrow samples gene names (rows) to shared genes
  if(!all(rownames(samp1) == rownames(samp2))) {
    samp2 <- samp2[rownames(samp1), ]
  }

  de_res <- tibble(Gene = rownames(samp1),
                   log2FC = sapply(rownames(samp1), function(gene) mean(samp1[gene, ]) - mean(samp2[gene, ])),
                   p.val = sapply(rownames(samp1), function(gene) {
                     sapply(1:ncol(samp1), function(i) {
                       sapply(1:ncol(samp2), function(j) wilcox.test(x = samp1[gene, i], y = samp2[gene, j])$p.value)
                     })
                   }),
                   p.adj = p.adjust(p.val, "fdr"),
                   Name = name)
  de_res
}


dea2 <- function(x, g1, g2, name) {
  de_res <- tibble(Gene = rownames(x),
                   log2FC = rowMeans(x[, g1]) - rowMeans(x[, g2]),
                   p.val = sapply(1:nrow(x), function(i) wilcox.test(x = x[i, g1], y = x[i, g2])$p.value),
                   p.adj = p.adjust(p.val, "fdr"),
                   avg.expr <- rowMeans(x[, g1]),
                   Name = name)
  de_res
}


DEGenes <- function(expr_mat, cl1, cl2, name, min_pct = 0.25, min_diff_pct = -Inf) {
  if(!is.null(min_pct)) {
    # Calculate per-gene percent expressed
    pct1 <- round(rowSums(expr_mat[, cl1, drop = FALSE] > 0) / length(cl1), digits = 3)
    pct2 <- round(rowSums(expr_mat[, cl2, drop = FALSE] > 0) / length(cl2), digits = 3)
    pct_res <- as.data.frame(cbind(pct1, pct2))
    colnames(pct_res) <- c("pct_1", "pct_2")

    # feature selection (based on percentages)
    alpha_min <- pmax(pct_res$pct_1, pct_res$pct_2)
    names(alpha_min) <- rownames(pct_res)
    features <- names(which(alpha_min >= min_pct))
    if (length(features) == 0) {
      warning("No features pass min_pct threshold; returning empty data.frame")
      return(pct_res[features, ])
    }
    alpha_diff <- alpha_min - pmin(pct_res$pct_1, pct_res$pct_2)
    features <- names(which(alpha_min >= min_pct & alpha_diff >= min_diff_pct))
    if (length(features) == 0) {
      warning("No features pass min_diff_pct threshold; returning empty data.frame")
      return(pct_res[features, ])
    }
  }

  de_res <- tibble(Gene = rownames(expr_mat[features, ]),
                   log2FC = rowMeans(expr_mat[features, cl1]) - rowMeans(expr_mat[features, cl2]),
                   p.val = sapply(1:nrow(expr_mat[features, ]), function(i) wilcox.test(x = expr_mat[i, cl1], y = expr_mat[i, cl2])$p.value),
                   p.adj = p.adjust(p.val, "fdr"),
                   Name = name)

  return(de_res)
}



# Scoring -----------------------------------------------------------------


####Add scores
score <- function(matrix, genes, bin_size = 100, bin_num = 30, center_rows = TRUE, conserved.genes = 0.7){

  filt_mat <- matrix[rownames(matrix) %in% genes, ]
  if(is.null(bin_size) || is.null(bin_num)){
    score <- colMeans(filt_mat)
    } else {
    score <- scalop::sigScores(as.matrix(matrix), genes, center.rows = center_rows, expr.binsize = bin_size, expr.nbin = bin_num, conserved.genes = conserved.genes)[, 1]
  }
  names(score) <- colnames(matrix)
  return(score)
}


metaprog_score <- function(matrix, metaprog_genes, score_diff = 0, score_cut = NULL, bin_size = 100, bin_num = 30, center_rows = TRUE, conserved.genes = 0.7, scale = TRUE) {

  # Prepare input as dataframe
  if(!is.list(metaprog_genes)) {stop("Input should be list of metaprograms containing genes")}
  metaprog_df <- reshape2::melt(unlist(metaprog_genes), value.name = "Gene")
  metaprog_df$Program <- unlist(lapply(seq_along(metaprog_genes), function(x) {
    prog_name <- names(metaprog_genes[x])
    leng <- length(metaprog_genes[[x]])
    out <- c(rep(prog_name, leng))
  }))

  # Subset metaprogram data frame to include only genes that appear in the matrix being scored
  # metaprog_df <- metaprog_df[metaprog_df$Gene %in% rownames(matrix),]

  #Score cells
  programs <- unique(metaprog_df$Program)
  score_ls <- parallel::mclapply(programs, function(x){
    meta_score <- score(matrix, metaprog_df[metaprog_df$Program == x, "Gene"], bin_size = bin_size, bin_num = bin_num, center_rows = center_rows, conserved.genes = conserved.genes)
  }, mc.cores = length(programs))
  names(score_ls) <- programs
  score_ls <- score_ls[sapply(score_ls, is.numeric)]
  scores <- do.call(cbind, score_ls)

  # Classify cells by the highest ranking metaprogram
  cell_membership <- character(length = nrow(scores))
  names(cell_membership) <- rownames(scores)
  for(i in 1:nrow(scores)){
    first_scoring_mp <- which.max(scores[i, ])
    second_scoring_mp <- max(scores[i, -first_scoring_mp])
    abs_maximum <- max(scores[i, ])
    if(score_diff != 0){
      scaled_range <- ifelse(isTRUE(scale), score_diff * abs_maximum, score_diff)
      if(!is.null(score_cut)){
        cell_membership[i] <- ifelse(abs_maximum - second_scoring_mp > scaled_range & abs_maximum > score_cut, names(first_scoring_mp), "Unresolved")
      } else {
        cell_membership[i] <- ifelse(abs_maximum - second_scoring_mp > scaled_range, names(first_scoring_mp), "Unresolved")
      }
    } else {
      if(!is.null(score_cut)){
        cell_membership[i] <- ifelse(abs_maximum > score_cut, names(first_scoring_mp), "Unresolved")
      } else {
        cell_membership[i] <- names(first_scoring_mp)
      }
    }
  }
  out <- list(meta_program = cell_membership, mp_score_list = score_ls)
  return(out)
}


refine_programs <- function(mat, programs, subset_cells = colnames(mat),
                            var_filter = TRUE, cor_filter = TRUE, var_quant_thresh = 0.5, cor_thresh = 0.4, cor_diff = 0,
                            return_cor_mat = FALSE, return_diagnostic_plot = FALSE, ...) {
  # Check input
  if(is.character(programs)) {programs <- list(programs)}
  if(isFALSE(var_filter) & isFALSE(cor_filter)) {stop("No filters for refinment are turned on.")}

  if(isTRUE(var_filter)) {
    subset_mat <- mat[rownames(mat) %in% unique(unlist(programs)), subset_cells]
    gene_var <- apply(subset_mat, 1, var)
    ordered_gene_var <- sort(gene_var, decreasing = TRUE)
    gene_filt <- ordered_gene_var[ordered_gene_var > quantile(ordered_gene_var, probs = var_quant_thresh)]
    set_class <- unlist(sapply(names(gene_filt), function(x) names(Unlist(programs)[Unlist(programs) == x]), simplify = TRUE))
    programs <- split(names(set_class), unname(set_class))
  }

  if(isTRUE(cor_filter)) {
    # Initial scoring - create score vectors to correlate
    score_df <- as.data.frame(scalop::sigScores(mat[, subset_cells], programs, ...))
    programs <- lapply(programs, function(x) x[x %in% rownames(mat)])    # check that all genes in the signature exist in the matrix

    # For each gene -- run correlation between its expression values (gene expression in all selected cells/spots) and the score vector for the programs it belongs to
    prog_cor <- sapply(seq_len(ncol(score_df)), function(idx) {
      run_cor <- sapply(programs[[idx]], function(gene) {
        cor(as.numeric(mat[gene, ]), score_df[[idx]])
      })
    })
    names(prog_cor) <- colnames(score_df)

    if(return_cor_mat || cor_diff != 0) {
      cor_mat <- sapply(seq_len(length(prog_cor)), function(idx) {
        test_cor <- sapply(unique(unname(unlist(programs))), function(x) {
          cor(as.numeric(mat[x, ]), score_df[[idx]])
        })
      })
      colnames(cor_mat) <- colnames(score_df)
    }

    # Filter programs based on their correlation coefficient in regrads to selected threshold
    programs <- sapply(prog_cor, function(x) names(x)[x >= cor_thresh], simplify = F)
  }

  if(isTRUE(cor_filter) && cor_diff != 0) {
    # Filter based on the difference between the top and second best correlating gene
    filtered_cor_mat <- cor_mat[unique(unlist(programs)), ]
    gene_membership <- character(length = nrow(filtered_cor_mat))
    names(gene_membership) <- rownames(filtered_cor_mat)
    for(i in 1:nrow(filtered_cor_mat)) {
      first_cor_prog <- which.max(filtered_cor_mat[i, ])
      second_cor_prog <- max(filtered_cor_mat[i, -first_cor_prog])
      abs_maximum <- max(filtered_cor_mat[i, ])

      gene_membership[i] <- ifelse(abs_maximum - second_cor_prog > cor_diff, names(first_cor_prog), "Unresolved")
    }

    distinct_genes <- gene_membership[gene_membership != "Unresolved"]
    programs <- split(names(distinct_genes), distinct_genes)
  }

  if(return_diagnostic_plot) {
    filtered_cor_mat <- cor_mat[unique(unlist(programs)), ]
    plot_cor_mat <- reshape2::melt(filtered_cor_mat)
    cor_plot <- gmap(dat = plot_cor_mat, x = Var2, y = Var1, fill = value, type = "continuous", limits = c(-1, 1),
                     legend.height = 5, legend.width = 1, ratio = 1) + theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1))
  }

  if(isTRUE(return_cor_mat) & isFALSE(return_diagnostic_plot)) {return(list(programs = programs, cor_mat = cor_mat))}
  else if(return_cor_mat && return_diagnostic_plot) {return(list(programs = programs, cor_mat = cor_mat, diagnostic_plot = cor_plot))}
  else if(isTRUE(return_diagnostic_plot) & isFALSE(return_cor_mat)) {return(list(programs = programs, diagnostic_plot = cor_plot))}
  else {return(programs)}
}



# Enrichment Analysis Functions -------------------------------------------


msigdb = function(category = c('H', 'C2', 'C5'),
                  subcategory = NULL,
                  exclude.subcategory = 'HPO', # Human Phenotype Ontology (in GO)
                  split.by.subcategory = F,
                  species = 'Homo sapiens',
                  annotation = c('gene_symbol','entrez_gene'),
                  return.dataframe = F) {

  annotation = match.arg(annotation)

  result = list()

  for (cat in category) {
    d = msigdbr::msigdbr(category = cat, species = species)
    d = d %>% dplyr::filter(!gs_subcat %in% exclude.subcategory)

    if (!is.null(subcategory) && any(subcategory %in% d$gs_subcat)) {
      d = d %>% dplyr::filter(gs_subcat %in% subcategory)
    }

    if (split.by.subcategory) {
      dl = split(d, d$gs_subcat)
    } else {
      dl = list(d)
    }

    if (return.dataframe) sigs = dl

    else {
      sigs = sapply(dl, function(di) {
        split(di[[annotation]], di$gs_name)},
        simplify = F)

      # just in case.
      sigs = sapply(sigs, function(s) s[!duplicated(s)], simplify = F)
    }

    result = c(result, setNames(sigs, cat))
  }

  result
}



enricher = function(test_gene_sets,
                    ref_gene_sets,
                    universe = NULL,
                    minGSSize = 5,
                    maxGSSize = 500,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2) {

  if (is.character(test_gene_sets)) {
    test_gene_sets = list(test_gene_sets)
  }

  test_gene_sets = sapply(test_gene_sets, function(s) limma::alias2Symbol(s), simplify = F)
  universe = limma::alias2Symbol(universe)
  ref_universe = unique(unlist(ref_gene_sets))
  universe = intersect(universe, ref_universe)
  test_gene_sets = sapply(test_gene_sets, function(x) x[x %in% universe], simplify = F)
  ref_gene_sets = sapply(ref_gene_sets, function(x) x[x %in% universe], simplify = F)

  term2gene = data.frame(term = names(Unlist(ref_gene_sets)),
                         gene = Unlist(ref_gene_sets),
                         stringsAsFactors = F)

  # enrichment test
  result = sapply(test_gene_sets, function(gene_set) {
    clusterProfiler::enricher(gene = gene_set,
                              universe = universe,
                              TERM2GENE = term2gene,
                              minGSSize = minGSSize,
                              maxGSSize = maxGSSize,
                              pAdjustMethod = pAdjustMethod,
                              pvalueCutoff = pvalueCutoff,
                              qvalueCutoff = pvalueCutoff)},
    simplify = F)

  sapply(1:length(result), function(i) {
    as.data.frame(result[[i]]) %>% dplyr::mutate(name=rep(names(result)[i],nrow(result[[i]])))},
    simplify = F)
}




# CNA Inference -----------------------------------------------------------

#### genome_sort - with gene list as input, reorder genes by chromosome  position and get breaks for chromosomes and arms
genome_sort <- function(genes, genome = NULL, attributes = c("hgnc_symbol", "start_position", "chromosome_name", "band"), downURL = "uswest.ensembl.org"){

  # Choose which species to use and server to download from
  if(is.null(genome)){
    mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = downURL, port = 80)
    #Get attributes for gene list
    cna_attr <- biomaRt::getBM(attributes = attributes, filters = "hgnc_symbol", values = genes, mart = mart, uniqueRows = TRUE)
  }
  if(!is.null(genome)){
    cna_attr <- as.data.frame(genome)
    rownames(cna_attr) <- cna_attr[, 1]
    cna_attr <- cna_attr[rownames(cna_attr) %in% genes, ]
    colnames(cna_attr)[1] <- "hgnc_symbol"
  }

  # Remove genes not mapped to proper chromosomes
  chr_names <- c(seq(1:22), "X", "Y")
  genes2rm <- setdiff(unique(cna_attr$chromosome_name), chr_names)
  rmgenes <- lapply(genes2rm, function(x){
    gene_drop <- cna_attr[grep(x, cna_attr$chromosome_name), ]
  })
  gene2rm <- do.call(rbind, rmgenes)
  cna_attr <- cna_attr[!is.element(rownames(cna_attr), rownames(gene2rm)), ]
  # Remove doubles
  cna_attr <- cna_attr[!duplicated(cna_attr$hgnc_symbol), ]
  # Remove NA's
  cna_attr <- na.omit(cna_attr)

  # Change X and Y chromosomes to numbers for ordering
  cna_attr$chromosome_name <- gsub("X", 23, cna_attr$chromosome_name)
  cna_attr$chromosome_name <- gsub("Y", 24, cna_attr$chromosome_name)
  cna_attr$chromosome_name <- as.numeric(cna_attr$chromosome_name)
  # Order
  cna_attr <- cna_attr[order(cna_attr$chromosome_name, cna_attr$start_position, decreasing = FALSE), ]
  # Chromosome number as sorting vector
  chr_names <- as.character(unique(cna_attr$chromosome_name))
  # Chromosome length
  chr_length <- lapply(chr_names, function(x){
    l <- nrow(subset(cna_attr, cna_attr$chromosome_name == x))
  })
  chr_length <- do.call(rbind, chr_length)
  # Break vector at each chromosome end
  chr_breaks <- cumsum(chr_length)
  # Set all genes as p or q
  cna_attr$band <- gsub("[[:digit:]]", "", cna_attr$band)
  cna_attr$band <- gsub("p.", "p", cna_attr$band)
  cna_attr$band <- gsub("q.", "q", cna_attr$band)
  # Chromosome arm lenghts
  arm_breaks <- integer(length = length(chr_breaks))
  for(i in seq_along(chr_breaks)){
    starting <- chr_breaks[i]
    ending <- nrow(cna_attr[cna_attr$chromosome_name == i & cna_attr$band == "q", ])
    gap_length <- starting - ending
    arm_breaks[i] <- round(gap_length, digits = 0)
  }

  # Breaks and labels for arms and full chromosomes
  full_breaks <- sort(c(1, chr_breaks, arm_breaks))
  # Add p and q labels
  q2q <- paste(seq(1, length(chr_breaks)), "q", sep = "")
  p2p <- paste(seq(1, length(chr_breaks)), "p", sep = "")
  full_labels <- sort(c(seq(1, length(chr_breaks)), seq(1, length(chr_breaks))))
  full_labels[seq(2, length(full_labels), 2)] <- q2q
  full_labels[seq(1, length(full_labels), 2)] <- p2p
  # Empty label at end of Y chromosome
  full_labels <- c(full_labels, " ")
  # Name X and Y chromosomes
  full_labels[45:48] <- c("Xp", "Xq", "Yp", "Yq")
  chr_breaks <- c(1, chr_breaks)
  chr_names[23:25] <- c("X", "Y", " ")

  output <- list(cna_genes = cna_attr$hgnc_symbol,
                 chr_breaks = chr_breaks,
                 chr_names = chr_names,
                 arm_breaks = arm_breaks,
                 full_breaks = full_breaks,
                 full_labels = full_labels,
                 all = cna_attr)
  return(output)
}


####

genome_break <- function(genes, genecut = 10, separate = "chr", genome = NULL){

  chr_sort <- genome_sort(genes, genome = genome)
  genes <- chr_sort$cna_genes

  if(separate == "chr") {breaks = chr_sort$chr_breaks
  labels = chr_sort$chr_names}
  if(separate == "arm") {breaks = chr_sort$full_breaks
  labels = chr_sort$full_labels}
  if(length(breaks) != length(labels)) {labels = labels[1:length(breaks)]}

  # Breaks as factor
  break_idx <- character(length = length(genes))
  for(i in 1:(length(breaks) - 1)){
    if(i == (length(breaks) - 1)) {end = length(genes)}
    if(i != (length(breaks) - 1)) {end = breaks[i + 1]}
    chr_leng <- end - breaks[i]
    if(chr_leng < genecut & chr_leng > 0 & i != 1) {break_idx[breaks[i] + 1:end] = NA}
    if(chr_leng < genecut & chr_leng > 0 & i == 1) {break_idx[breaks[i]:end] = NA}
    if(chr_leng == 0 & i != 1) {break_idx[breaks[i] + 1] = NA}
    if(chr_leng == 0 & i == 1) {break_idx[breaks[i]] = NA}
    if(chr_leng > genecut & i != 1) {break_idx[breaks[i] + 1:end] = i}
    if(chr_leng > genecut & i == 1) {break_idx[breaks[i]:end] = i}
  }

  # Remove genes on chromosomes with too few genes
  dispose <- which(is.na(break_idx))
  break_idx <- na.omit(break_idx)
  break_idx <- sort(as.numeric(break_idx))
  genes <- chr_sort$cna_genes[-dispose]
  names(break_idx) <- genes

  # Keep only labels with chromosomes that have enough genes
  breaks_df <- cbind.data.frame(breaks, labels)
  rownames(breaks_df) <- seq(1, nrow(breaks_df))
  breaks_df <- breaks_df[unique(break_idx), ]

  out <- list(break_idx = break_idx, breaks_df = breaks_df, dispose = dispose)
  return(out)
}


#### calc_cna - calculate CNA scores per cell/spot given a gene x cell matrix, query cells and reference cell list/vector
# NOTE: matrix should NOT be row(gene)-centered
calc_cna <- function(matrix,
                     query,
                     ref,
                     top_genes_num = NULL,
                     top_genes_vec = NULL,
                     window = 100,
                     range = c(-3, 3),
                     per_chr = TRUE,
                     scale = 0.05,
                     noise = NULL,
                     isLog = FALSE,
                     genome = NULL,
                     min_ref_leng = 10,
                     verbose = FALSE){

  if(all(round(range(rowMeans(matrix)), 3) == 0)) {
    stop(print("Matrix is row-centered. Please provide non-centered data."))
  }
  if(is.list(ref)){
    # Remove small references
    ref_leng <- lapply(ref, length)
    ref <- ref[ref_leng > min_ref_leng]
    # Prepare CNA matrix to work on
    cna_mat <- matrix[, c(query, unlist(ref))]
  }
  else if(!is.list(ref)){
    if(length(ref) < min_ref_leng){
      stop(print("There are not enough reference cells!"))
    }
    cna_mat <- matrix[, c(query, ref)]
    }

  # Order chromosomes
  if (verbose) message("Ordering the genes by their genomic position.")
  order_genes <- genome_sort(rownames(cna_mat), genome = genome)
  genes <- order_genes$cna_genes
  chr_breaks <- order_genes$chr_breaks

  # Optional list of top expressed genes
  if(!is.null(top_genes_num)){
    if(verbose) message("Filtering the expression matrix to include only top ", length(top_genes_num), "genes...")
    # if(isTRUE(isLog)){
    #   cna_mat <- un_log(cna_mat)}
    cna_mat <- cna_mat[genes, ]
    mean_exp <- apply(cna_mat, 1, mean)
    top_exp <- names(head(sort(mean_exp, decreasing = TRUE), n = top_genes_num))
    order_genes <- genome_sort(top_exp, genome = genome)
    genes <- order_genes$cna_genes
    chr_breaks <- order_genes$chr_breaks
  }
  if(!is.null(top_genes_vec)){
    if(verbose) message("Filtering the expression matrix to include only top ", length(top_genes_vec), " selected genes...")
    cna_mat <- cna_mat[top_genes_vec, ]
    order_genes <- genome_sort(rownames(cna_mat), genome = hg38)
    genes <- order_genes$cna_genes
    chr_breaks <- order_genes$chr_breaks
  }

  # Reorder
  ordered_mat <- cna_mat[genes, ]
  # Log before first centering
  if(isFALSE(isLog)){
    ordered_mat <- log_norm(ordered_mat)}
  # First row centering step
  if (verbose) message("Performing mean-centering of the genes.")
  avg <- apply(ordered_mat, 1, mean)
  ordered_mat <- sweep(ordered_mat, 1, avg)
  # Set 3 and -3 as extreme values (as set by the argument "range")
  if (verbose) message("Restricting expression matrix values to between ", range[[1]], " and ", range[[2]], ".")
  ordered_mat <- apply(ordered_mat, 2, function(x) pmax(x, range[1]))
  ordered_mat <- apply(ordered_mat, 2, function(x) pmin(x, range[2]))
  # Unlog to CPM/TPM again for moving average
  ordered_mat <- un_log(ordered_mat)

  # Calculate moving average per chromosome by window of 100 (as set by the argument "window")
  if(isTRUE(per_chr)){
    if (verbose) message("Calculating rolling means with a window size of ", window, " genes, on each chromosome in turn.")
    num <- seq(1:(length(chr_breaks) -1))
    perchr <- lapply(num, function(y){
      if(y == length(num)){
        end <- nrow(ordered_mat)
        }
      if(y != length(num)){
        end <- chr_breaks[y + 1] - 1
        }
      chr <- ordered_mat[chr_breaks[y]:end, ]
      chr_mat <- apply(chr, 2, function(x) caTools::runmean(x, k = window, endrule = "mean"))
    })
    calc_mat <- do.call(rbind, perchr)
  }
  # Calculate moving average for all genes as one chromosome
  if(isFALSE(per_chr)){
    if (verbose) message("Calculating rolling means with a window size of ", window, " genes, on all genes as one chromosome.")
    calc_mat <- apply(ordered_mat, 2, function(x) caTools::runmean(x, k = window, endrule = "mean"))
    }

  # Log before second centering
  if (verbose) message("Converting CNA score values to log(2) space.")
  calc_mat <- log_norm(calc_mat)
  # Substract median per cell
  if (verbose) message("Performing median-centering of the cells.")
  cell_med <- apply(calc_mat, 2, median)
  calc_mat <- sweep(calc_mat, 2, cell_med)
  # Unlog to CPM/TPM again for reference removal
  calc_mat <- un_log(calc_mat)

  # Create max/min values per gene from reference cells
  if(is.list(ref)){
    ref_leng <- seq(1, length(ref))
    mean_mat <- lapply(ref_leng, function(x){
      idx_ref <- ref[[x]]
      m1 <- apply(calc_mat[, idx_ref], 1, mean)
      m1
      })
    ref_mat <- do.call(cbind, mean_mat)
    # Log references
    ref_mat <- log_norm(ref_mat)
    ref_max <- apply(ref_mat, 1, max)
    ref_min <- apply(ref_mat, 1, min)
  }
  else if(!is.list(ref)){
    ref_mat <- apply(calc_mat[, ref], 1, mean)
    ref_mat <- log_norm(ref_mat)
    ref_max <- ref_mat
    ref_min <- ref_mat
  }

  # Expand reference boundaries by scaling percentage
  if(!is.null(scale)){
    rmax <- ref_max + scale * abs(ref_max)
    rmin <- ref_min - scale * abs(ref_min)
  }
  # Or expand by fixed noise factor
  if(!is.null(noise)){
    rmax <- ref_max + noise
    rmin <- ref_min - noise
  }

  # Log CNA matrix
  calc_mat <- log_norm(calc_mat)
  # Centre by reference
  if (verbose) message("Correcting CNA profiles using CNA values from reference cells.")
  score_mat <- ifelse(calc_mat > rmax, calc_mat - rmax,
                     ifelse(calc_mat < rmin, calc_mat - rmin, 0))
  rownames(score_mat) <- rownames(ordered_mat)

  if (verbose) message("Done!")
  return(score_mat)
}




#### Find the malignant and non-malignant subsets of cells from a gene-by-cell matrix of CNA values.
#######Define CNA from score matrix
assign_cna_status <- function(cna_matrix,
                              epi_cells,
                              ref_cells,
                              cna_sig = "abs",
                              top_region = 1/3,
                              top_cells = 1/3,
                              clone_name = NULL,
                              thresh_method = c("DFM", "percentile"),
                              sd_cut = 2.5,
                              class_cut = 0.99,
                              upper_cut = 1,
                              lower_cut = 1,
                              knn_prob = 0.9,
                              global_thr = NULL) {

  if(is.list(ref_cells)){ref_cells <- unlist(ref_cells)}
  if(!is.null(intersect(ref_cells, epi_cells))) {
    dispose <- intersect(ref_cells, epi_cells)
    ref_cells <- setdiff(ref_cells, dispose)
  }
  cna_matrix <- cna_matrix[, c(epi_cells, ref_cells)]
  all_cells <- colnames(cna_matrix)
  epi_matrix <- cna_matrix[, epi_cells]

  # Absolute CNA value per gene for all epithelial cells
  abs_vals <- apply(epi_matrix, 1, function(x) mean(abs(x)))

  # Select top X-percentage genes with highest absolute value
  top_abs <- ifelse(abs_vals > quantile(abs_vals, probs = 1 - top_region), abs_vals, 0)
  # Define matrices where to calculate cell averages - only the regions with most CNA
  score_cna <- cna_matrix[top_abs > 0, ]

  # CNA score for each tumor cell across relevant region
  if(cna_sig == "abs") {cna_score <- apply(score_cna, 2, function(x) mean(abs(x)))}
  if(cna_sig == "square") {cna_score <- apply(score_cna, 2, function(x) mean(x^2))}
  if(cna_sig == "sqrt") {cna_score <- apply(score_cna, 2, function(x) mean(sqrt(x^2)))}

  # Calculate correlation only with epithelial cells with strongest signal
  def_cutoff <- quantile(cna_score, probs = 1 - top_cells)
  epi_cna <- epi_matrix[top_abs > 0, cna_score[epi_cells] > def_cutoff]
  # Correlation vector
  cor_vec <- cor(score_cna, epi_cna)
  cna_cor <- apply(cor_vec, 1, mean)

  # Format in data frame
  df <- as.data.frame(t(rbind.data.frame(cna_cor, cna_score)))
  rownames(df) <- all_cells
  colnames(df) <- c("cna_cor", "cna_score")
  df$Type <- with(df, ifelse(rownames(df) %in% epi_cells, "Epithelial", "Reference"))

  ## Define thresholds and create cell classification score
  if(thresh_method == "percentile"){
    score_thr <- quantile(df[df$Type == "Reference", "cna_score"], probs = class_cut)
    cor_thr <- quantile(df[df$Type == "Reference", "cna_cor"], probs = class_cut)
  }
  else if(thresh_method == "DFM"){
    score_thr <- mean(df[df$Type == "Reference", "cna_score"]) + sd_cut * sd(df[df$Type == "Reference", "cna_score"])
    cor_thr <- mean(df[df$Type == "Reference", "cna_cor"]) + sd_cut * sd(df[df$Type == "Reference", "cna_cor"])
  }
  else if(!is.null(global_thr)){
    score_thr <- max(score_thr, global_thr$score)
    cor_thr <- max(cor_thr, glob_thr$cor)
  }

  score_bin <- ifelse(df[epi_cells, "cna_score"] > score_thr, 1, 0)
  cor_bin <- ifelse(df[epi_cells, "cna_cor"] > cor_thr, 1, 0)
  cna_bin <- score_bin + cor_bin
  df[epi_cells, "base_cna_status"] <- ifelse(cna_bin > upper_cut, "Malignant", ifelse(cna_bin < lower_cut, "Normal", "Unresolved"))
  df$tot_score <- with(df, ifelse(rownames(df) %in% epi_cells, cna_bin, -1))
  df$cna_status <- with(df, ifelse(rownames(df) %in% epi_cells, df$base_cna_status, "Reference"))

  # Plot definition with cutoffs
  p2 <- ggplot(df, aes(x = cna_cor, y = cna_score, colour = cna_status)) + theme_classic() +
    geom_point(aes(size = Type)) +
    scale_size_manual(name = "Type", values = c(2, 0.8)) +
    geom_vline(xintercept = cor_thr) +
    geom_hline(yintercept = score_thr)  +
    ggtitle(paste(clone_name, "CNA Correlation vs Signature Initial Cutoffs", sep = " ")) + labs(x = "CNA Correlation", y = "CNA Score")

  ## Add KNN as parameter
  unres_cells <- rownames(df[df$cna_status == "Unresolved", ])
  if(length(unres_cells) > 2){
    knn_cna <- class::knn(train = df[df$cna_status != "Unresolved" & df$cna_status != "Reference", 1:2],
                          test = df[df$cna_status == "Unresolved", 1:2], k = max(1, round(log(length(epi_cells)))),
                          cl = df$cna_status[df$cna_status != "Unresolved" & df$cna_status != "Reference"], prob = TRUE)
    df[unres_cells, c("knn", "prob")] <- c(as.character(knn_cna), attr(knn_cna, "prob"))
    df$final_cna_status <- ifelse(!rownames(df) %in% unres_cells, df$cna_status,
                                  ifelse(rownames(df) %in% unres_cells & df$knn == "Malignant" & df$prob >= knn_prob, "Malignant",
                                         ifelse(rownames(df) %in% unres_cells & df$knn == "Normal" & df$prob >= knn_prob, "Normal", "Unresolved")))
  }
  if(length(unres_cells) < 3) {df$final_cna_status <- df$cna_status}

  # Plot final CNA definition
  p3 <- ggplot(df, aes(x = cna_cor, y = cna_score, colour = final_cna_status)) + theme_classic() +
    geom_point(aes(size = Type)) +
    scale_size_manual(name = "Type", values = c(2, 0.8)) +
    geom_vline(xintercept = cor_thr) +
    geom_hline(yintercept = score_thr)  +
    ggtitle(paste(clone_name, "CNA Correlation vs Signature Final Definition", sep = " ")) + labs(x = "CNA Correlation", y = "CNA Score")


  plots <- egg::ggarrange(p2, p3, ncol = 2, nrow = 1)
  def_mat <- df[epi_cells, ]
  output <- list(plots = plots, def_mat = def_mat, thresholds = c(cor_thr, score_thr))
  return(output)
}



#### Create amplification/deletion matrix
admat.fun<-function(matrix,breaks,vector,lab,cut=0.15,genes="top",region=1/2,refmat=NULL){
  bn<-unique(breaks)
  tn<-sort(unique(vector))
  #Split matrix per subclone
  cut.mat<-lapply(tn,function(x){
    n1<-names(vector[vector==x])
    m1<-matrix[,n1]
    # if(is.null(ncol(m1))){
    #   m1 <- NULL} else{m1}
  })
  # cut.mat <- cut.mat[-which(sapply(cut.mat, is.null))]
  # tn <- unique(vector)[1:length(cut.mat)]
  #Mean arm expression per subclone
  epmat<-matrix(nrow = length(bn),ncol = length(tn))
  for(i in 1:length(tn)){
    for(j in 1:length(bn)){
      m1<-cut.mat[[i]]
      m2<-m1[breaks==bn[j],]
      absvals<-apply(m2,1,function(x)mean(abs(x)))
      m3<-m2[absvals>quantile(absvals, prob = 1-region),]
      if(!is.null(dim(m3))){c<-mean(colMeans(m3))}else{c=0}
      d<-mean(colMeans(m2))
      if(genes=="top"){epmat[j,i]<-c}
      if(genes=="all"){epmat[j,i]<-d}
      if(genes=="ref" && !is.null(refmat)){
        mr<-refmat[breaks==bn[j],!colnames(refmat) %in% cells]
        mc<-mean(colMeans(mr))
        epmat[j,i]<-c-mc}

    }
  }
  rownames(epmat)<-lab
  colnames(epmat)<-tn
  #Amplification/deletion matrix
  ad<-matrix(ncol = length(tn),nrow = length(bn))

  for(i in 1:length(bn)){
    for(j in 1:length(tn)){
      ad[i,j]<-ifelse(epmat[i,j]>cut,"amp",ifelse(epmat[i,j]< -1*cut,"del","none"))
    }}
  rownames(ad)<-lab
  colnames(ad)<-tn

  ad<-as.data.frame(ad)
  ad[is.na(ad)]<-"none"
  epmat<-as.data.frame(epmat)
  epmat[is.na(epmat)]<-0
  out<-list(ad=ad,epmat=epmat)
  return(out)
}


#### Filter arms by variance cutoff or top variant arms
varfilt<-function(epmat,ad,cut=0.01,n=10){
  armvar<-apply(epmat,1,var)
  if(is.null(cut)){
    z<-names(tail(sort(armvar),n))
    epmat<-epmat[z,]
    ad<-ad[z,]}
  if(!is.null(cut)){
    epmat<-epmat[armvar>cut,]
    ad<-ad[armvar>cut,]
  }
  out<-list(ad=ad,epmat=epmat)
  return(out)
}

### Arms variance:
var_filt <- function(chr_exp_mat, cutoff = 0.01, n = 15) {
  arm_var <- apply(chr_exp_mat, 1, var)
  if(is.null(cutoff)) {
    top_var_arms <- chr_exp_mat[names(tail(sort(arm_var), n)), ]
  }
  if(!is.null(cutoff)) {
    top_var_arms <- chr_exp_mat[arm_var > cutoff, ]
  }
  return(top_var_arms)
}

#### Calculate CNAsignal & CNAcorrelation scores
cna_sig_cor <- function(cna_matrix,
                              epi_cells,
                              ref_cells,
                              top_region = 1/3,
                              top_cells = 1/3,
                              cna_sig = "abs"){

  if(is.list(ref_cells)){ref_cells <- unlist(ref_cells)}
  if(!is.null(intersect(ref_cells, epi_cells))){
    rmcells <- intersect(ref_cells, epi_cells)
    ref_cells <- setdiff(ref_cells, rmcells)
  }
  cna_matrix <- cna_matrix[, c(epi_cells, ref_cells)]
  all_cells <- colnames(cna_matrix)
  epi_matrix <- cna_matrix[, epi_cells]

  # Absolute CNA value per gene for all epithelial cells
  absvals <- apply(epi_matrix, 1, function(x) mean(abs(x)))
  # Top one-third of genes with highest absolute value (as set by the argument "top_region")
  top_abs <- ifelse(absvals > quantile(absvals, probs = 1 - top_region), absvals, 0)
  # Define matrices where to calculate cell averages - only the regions with most CNA
  score_cna <- cna_matrix[top_abs > 0, ]
  # CNA score for each tumor cell across relevant region
  if(cna_sig == "abs"){cna_score <- apply(score_cna, 2, function(x) mean(abs(x)))}
  if(cna_sig == "sqrt"){cna_score <- apply(score_cna, 2, function(x) mean(sqrt(x^2)))}
  if(cna_sig == "square"){cna_score <- apply(score_cna, 2, function(x) mean(x^2))}
  # Calculate correlation only with epithelial cells with strongest signal (defined by a cutoff set by the argument "top_cells")
  cellcut <- quantile(cna_score, probs = 1 - top_cells)
  epi_cna <- epi_matrix[top_abs > 0, cna_score[epi_cells] > cellcut]
  # Correlation vector
  cor_vec <- cor(score_cna, epi_cna)
  cna_cor <- apply(cor_vec, 1, mean)

  out <- list(CNA_Signal = cna_score,
              CNA_Correlation = cna_cor)
  return(out)
}



##########Separate groups by permutation testing
permutation.score<-function(matrix,signature,trials=100,correction="fdr",name="x",method="max",
                            mode="positive",cutoff=0.01,bin=100,cna=NULL,cna_cells=NULL){


  g1<-signature
  if(is.null(cna)){s1 <- score(matrix,g1,bin=bin)}else{
    if(cna=="abs"){s1 <- apply(matrix[g1,],2,function(x)mean(abs(x)))}
    if(cna=="sqrt"){s1 <- apply(matrix[g1,],2,function(x)mean(sqrt(x^2)))}
    if(cna=="cor"){s1<-rowMeans(cor(matrix[g1,],matrix[g1,cna_cells]))}
  }
  N <- trials

  perm_test <- mclapply(1:N, function(i) {
    permuted <- t(apply(matrix, 1, gtools::permute))

    colnames(permuted) <- colnames(matrix); rownames(permuted) <- rownames(matrix)

    if(is.null(cna)){scx <- score(permuted,g1,bin=bin)}else{
      if(cna=="abs"){scx <- apply(permuted[g1,],2,function(x)mean(abs(x)))}
      if(cna=="sqrt"){scx <- apply(permuted[g1,],2,function(x)mean(sqrt(x^2)))}
      if(cna=="cor"){scx<-rowMeans(cor(permuted[g1,],permuted[g1,cna_cells]))}
    }

    if(mode=="positive" & method=="max"){
      res = s1 > max(scx)}
    if(mode=="negative" &method=="max"){
      res = s1 < min(scx)
    }
    if(mode=="positive" & method=="sd"){
      res = s1 > (mean(scx) + 2*sd(scx))
    }
    if(mode=="negative" & method=="sd"){
      res = s1 < (mean(scx) - 2*sd(scx))
    }
    res

  })
  perm_test<-do.call(cbind,perm_test)
  # Count the number of times each cell was an outlier
  perm_test_res <- apply(perm_test, 1, function(x) length(which(x == TRUE)))
  prob_success <- sum(perm_test_res) / (N * length(perm_test_res))
  pv1 <- sapply(perm_test_res,
                function(x) binom.test(x = x, n = N, p = prob_success,
                                       alternative = "greater")$p.value)
  pv2<-p.adjust(p = pv1, method = correction)
  out<-ifelse(pv2<cutoff,name,"rest")
  return(out)
}



cna_subclones <- function(matrix,
                          epi_cells,
                          separate = "arm",
                          genome = NULL,
                          genecut = 10,
                          top_region = 1/3,
                          top_method = "abs",
                          cluster_by_cell = FALSE,
                          reduction_method = "uwot",
                          reduction_dims = 10,
                          clone_size = 10,
                          umap_pca = FALSE,
                          cluster_method = "dbscan",
                          louvain_k = 10,
                          dbs_minpts = NULL,
                          dbs_ptscale = NULL,
                          diffcut = 10,
                          hclust_k = 15,
                          cell_matrix = "numeric",
                          armcut = 0,
                          knn = TRUE,
                          merge = TRUE,
                          merge_method = "ad",
                          corcut = 0.95,
                          adcut = 0.2,
                          adgenes = "all",
                          adregion = 1,
                          armfilt = FALSE,
                          arm_n = 10,
                          dcut = 0.2,
                          name=NULL,
                          perm_bins=50,
                          perm_cut=0.01,
                          perm_trials=100,...){



  ## Subset by chromosome
  breaks <- genome_break(genes = rownames(matrix), separate = separate,
                         genecut = genecut, genome = genome)
  break_idx <- breaks$break_idx
  labels <- breaks$breaks_df$labels
  dispose <- breaks$dispose
  if(length(dispose) != 0) {matrix <- matrix[-dispose, ]}
  break_num <- unique(break_idx)

  ## Set working and reference matrices
  ref_mat <- matrix
  matrix <- matrix[, epi_cells]
  names(break_idx) <- rownames(matrix)

  ## Select the most CNA-rich regions
  if(top_method == "abs"){
    abs_vals <- apply(matrix, 1, function(x) mean(abs(x)))
    mat <- matrix[abs_vals > quantile(abs_vals, probs = 1 - top_region), ]
  }
  ## or alternatively - Select the most variant gene regions
  if(top_method == "sd"){
    stand_dev <- apply(matrix, 1, sd)
    mat <- matrix[stand_dev > quantile(stand_dev, probs = 1 - top_region), ]
  }

  trans_mat1 <- as.matrix(t(mat))

  #### Option 1 - Dimensionality reduction & Clustering
  if(isFALSE(cluster_by_cell)){
    ## Dimensional reduction methods:
    if(reduction_method == "none"){
      trans_mat3 <- trans_mat1
    }
    if(reduction_method == "pca"){
      pca <- prcomp(trans_mat1, center = FALSE, scale. = FALSE)
      trans_mat3 <- pca$x[, 1:reduction_dims]
    }
    if(reduction_method == "svd"){
      svd <- svd::propack.svd(trans_mat1, neig = reduction_dims)
      trans_mat3 <- as.data.frame(svd$u)
    }
    if(reduction_method == "ica"){
      ica <- fastICA::fastICA(trans_mat1, n.comp = reduction_dims, method = "C")
      trans_mat3 <- ica$S
    }
    if(reduction_method == "tsne"){
      trans_mat2 <- Rtsne::Rtsne(trans_mat1, check_duplicates = FALSE, pca = TRUE, perplexity = 30, theta = 0.5, dims = 2, pca_center = FALSE, pca_scale = FALSE, ...)
      trans_mat3 <- as.data.frame(trans_mat2$Y)
    }
    if(reduction_method == "umap"){
      if(isTRUE(umap_pca)){
        pca <- prcomp(trans_mat1, center = FALSE, scale. = FALSE)
        trans_mat1 <- pca$x[, 1:reduction_dims]
        }
      umap <- umap::umap(trans_mat1, metric = "pearson", spread = 5, min_dist = 0.01, n_neighbors = clone_size)
      trans_mat3 <- cbind.data.frame(umap$layout[, 1], umap$layout[, 2])
      colnames(trans_mat3) <- c("V1", "V2")
    }
    if(reduction_method == "uwot"){
      trans_mat3 <- standard_umap(t(trans_mat1), n_neighbors = clone_size, spread = 5, min_dist = 0.01, pca = reduction_dims)
    }

    rownames(trans_mat3) <- colnames(matrix)

    ## Clustering methods:
    if(cluster_method == "louvain"){
      clusters <- cluster_coord(trans_mat3, method = "louvain", louvain = louvain_k)
    }
    if(cluster_method == "dbscan"){
      # Select minpts for knn graph
      if(is.null(dbs_minpts) & is.null(dbs_ptscale)) {minpts <- log(nrow(trans_mat3))}
      if(!is.null(dbs_ptscale)) {minpts <- dbs_ptscale * log(nrow(trans_mat3))}
      if(!is.null(dbs_minpts)) {minpts <- dbs_minpts}
      # Run dbscan with set minpts
      clusters <- cluster_dbscan(trans_mat3, min_pts = minpts, diff_cut = diffcut)
    }
    if(cluster_method == "hclust"){
      hc_cor <- cor(t(trans_mat3))
      hc_dist <- as.dist((1 - hc_cor) / 2)
      hc <- hclust(hc_dist, method = "average")
      clusters <- cutree(hc, k = hclust_k)
    }
  }

  #### Option 2 - Clustering by-cell (not relying on dimensional reduction):
  if(isTRUE(cluster_by_cell)){
    # Subset break vector by filtered matrix
    filt_break_idx <- break_idx[rownames(mat)]
    filt_break_num <- unique(filt_break_idx)
    gene_leng_perArm <- lapply(filt_break_num, function(x) length(filt_break_num[filt_break_num == x]))
    if(length(filt_break_num[gene_leng_perArm > genecut]) != length(filt_break_num)){
      filt_break_num <- filt_break_num[gene_leng_perArm > genecut]
      filt_break_idx <- filt_break_idx[filt_break_idx %in% filt_break_num]
      mat <- mat[names(filt_break_idx), ]
    }
    # Chromosome arms label data-frame
    arm_label_df <- cbind.data.frame(break_num, labels)
    arm_label_df <- arm_label_df[arm_label_df$break_num %in% filt_break_num, ]

    # Select arms with highest expression and standard deviation      TO FIX SEPARATION!!!
    if(armcut != 0){
      chr_arm_stats <- lapply(filt_break_num, function(x){
        chr_arm_genes <- mat[filt_break_idx == x, ]
        avg_arm_exp <- mean(abs(colMeans(chr_arm_genes)))
        avg_arm_sd <- mean(apply(chr_arm_genes, 2, sd))
        c(avg_arm_exp, avg_arm_sd)
      })
      chr_arm_stats <- do.call(rbind.data.frame, chr_arm_stats)
      colnames(chr_arm_stats) <- c("mean_abs_exp","mean_sd")
      rownames(chr_arm_stats) <- filt_break_num
      chr_arm_stats$arm <- arm_label_df$labels
      most_diff_arms <- chr_arm_stats[chr_arm_stats$mean_abs_exp %in% head(sort(chr_arm_stats$mean_abs_exp, decreasing = TRUE), n = armcut) |
                                      chr_arm_stats$mean_sd %in% head(sort(chr_arm_stats$mean_sd, decreasing = TRUE), n = armcut), ]

      filt_break_idx <- filt_break_idx[filt_break_idx %in% rownames(most_diff_arms)]
      for(i in seq_along(most_diff_arms)){
        filt_break_idx[filt_break_idx == rownames(most_diff_arms)[i]] <- most_diff_arms$arm[i]
      }
      filt_break_num <- unique(filt_break_idx)
      mat <- mat[names(filt_break_idx), ]
    }

    ## Create chromosome expression matrics
    # Binary (amp / del / none) Cells - on - Arms matrix
    per_arm_exp <- lapply(filt_break_num, function(x){unname(colMeans(mat[filt_break_idx == x, ]))})
    per_arm_exp_vec <- as.numeric(do.call(cbind, per_arm_exp))
    exp_up_cut <- unname(quantile(per_arm_exp_vec, probs = 0.9))
    exp_low_cut <- unname(quantile(per_arm_exp_vec, probs = 0.1))
    binary_mat <- lapply(filt_break_num, function(x){
      avg_exp <- colMeans(mat[filt_break_idx == x, ])
      classify <- ifelse(avg_exp > exp_up_cut, "amp",
                         ifelse(avg_exp < exp_low_cut, "del", "none"))
      })
    binary_mat <- do.call(cbind.data.frame, binary_mat)
    colnames(binary_mat) <- filt_break_num

    # Numeric Cells - on - Arms matrix
    numeric_mat <- do.call(cbind.data.frame, per_arm_exp)
    colnames(numeric_mat) <- filt_break_num

    if(cluster_method == "hclust"){
      if(cell_matrix == "binary"){
        binary_mat[binary_mat == "amp"] <- 1
        binary_mat[binary_mat == "del"] <- -1
        binary_mat[binary_mat == "none"] <- 0
        binary_mat <- apply(binary_mat, 1, as.numeric)
        hc <- hclust(dist(t(binary_mat)), method = "average")
      }
      if(cell_matrix == "numeric"){
        hc_cor <- cor(t(numeric_mat))
        hc_dist <- as.dist((1 - hc_cor) / 2)
        hc <- hclust(hc_dist, method = "average")
      }
      clusters <- cutree(hc, k = hclust_k)
    }

    trans_mat3 <- numeric_mat
  }

  ## Rename clusters
  ren.fun2<-function(x){paste("subclone",x,sep = "")}
  clusters <- ren.fun2(clusters)
  names(clusters) <- colnames(matrix)

  small_clusts <- table(clusters) <= 2
  small_clusts <- small_clusts[small_clusts == TRUE]
  ## Assign cells in small clusters by knn method if no clusters are larger than min size
  if(isTRUE(knn) & max(table(clusters)) < clone_size){
    trans_mat4 <- trans_mat3
    trans_mat4$subclone <- clusters

    unres_cells <- rownames(trans_mat4[trans_mat4$subclone %in% names(small_clusts), ])
    trans_mat4 <- trans_mat4[, -ncol(trans_mat4)]

    knn_classify <- class::knn(train = trans_mat4[!rownames(trans_mat4) %in% unres_cells, ],
                               test = trans_mat4[rownames(trans_mat4) %in% unres_cells, ], k = max(1, round(log(nrow(trans_mat4)))),
                               cl = clusters[!names(clusters) %in% unres_cells], prob = TRUE)

    trans_mat4[unres_cells, c("knn", "prob")] <- c(as.character(knn_classify), attr(knn_classify, "prob"))

    clusters <- ifelse(names(clusters) %in% unres_cells, trans_mat4$knn, clusters)
    names(clusters) <- rownames(trans_mat4)
  }
  ## Assign cells in small clusters by knn method
  if(isTRUE(knn) && length(small_clusts) > 0){
    trans_mat3$subclone <- clusters
    unres_cells <- rownames(trans_mat3[trans_mat3$subclone %in% names(small_clusts), ])
    trans_mat3 <- trans_mat3[, -ncol(trans_mat3)]

    knn_classify <- class::knn(train = trans_mat3[!rownames(trans_mat3) %in% unres_cells, ],
                               test = trans_mat3[rownames(trans_mat3) %in% unres_cells, ], k = max(1, round(log(nrow(trans_mat3)))),
                               cl = clusters[!names(clusters) %in% unres_cells], prob = TRUE)

    trans_mat3[unres_cells, c("knn","prob")] <- c(as.character(knn_classify), attr(knn_classify, "prob"))

    clusters <- ifelse(names(clusters) %in% unres_cells, trans_mat3$knn, clusters)
    names(clusters) <- rownames(trans_mat3)
  }

  clust_sort <- sort(unique(clusters))
  ltrs <- letters[1:length(clust_sort)]
  for(i in seq_along(clust_sort)){
    clusters[clusters == clust_sort[i]] <- paste(ltrs[i], clust_sort[i], sep = "_")
  }

  sorted_clusts <- clusters

  #### Merge Clusters
  if(isTRUE(merge) && length(unique(clusters)) > 1){
    ## Merge clusters based on correlation of mean expression across chromosome arms:
    if(merge_method == "cor"){
      cormax <- 1
      while(corcut < cormax && length(unique(clusters)) > 1){
        amp_del_mats <- admat.fun(matrix, break_idx, clusters, labels, cut = adcut, genes = adgenes, region = adregion)
        arms_exp_mat <- amp_del_mats$epmat
        amp_del_mat <- amp_del_mats$ad
        clust_combns <- data.frame(t(combn(colnames(arms_exp_mat), 2)))
        if(isTRUE(armfilt)){
          amp_del_mat <- varfilt(arms_exp_mat, amp_del_mat, cut = NULL, n = arm_n)$ad
          arms_exp_mat <- varfilt(arms_exp_mat, amp_del_mat, cut = NULL, n = arm_n)$epmat
        }
        for(i in 1:nrow(clust_combns)){
          clust_a <- as.character(clust_combns[i, 1])
          clust_b <- as.character(clust_combns[i, 2])
          clust_combns[i, "cor"] <- cor(arms_exp_mat[, clust_a], arms_exp_mat[, clust_b])
          cormax <- max(clust_combns$cor)
          clusters[clusters == clust_a] <- ifelse(clust_combns[i, "cor"] > corcut, clust_b, clust_a)
        }
      }
    }
    ## Merge subclones where the maximal distance between two chromosome arms is lower than cut
    if(merge_method == "dist" | merge_method == "both"){
      dist_min <- 0
      while(dist_min < dcut && length(unique(clusters)) > 1){
        amp_del_mats <- admat.fun(matrix, break_idx, clusters, labels, cut = adcut, genes = adgenes, region = adregion)
        arms_exp_mat <- amp_del_mats$epmat
        amp_del_mat <- amp_del_mats$ad
        clust_combns <- data.frame(t(combn(colnames(arms_exp_mat), 2)))
        if(isTRUE(armfilt)){
          amp_del_mat <- varfilt(arms_exp_mat, amp_del_mat, cut = NULL, n = arm_n)$ad
          arms_exp_mat <- varfilt(arms_exp_mat, amp_del_mat, cut = NULL, n = arm_n)$epmat
        }
        for(i in 1:nrow(clust_combns)){
          clust_a <- as.character(clust_combns[i, 1])
          clust_b <- as.character(clust_combns[i, 2])
          clust_combns[i, "dist"] <- max(abs(arms_exp_mat[, clust_a] - arms_exp_mat[, clust_b]))
          dist_min <- min(clust_combns$dist)
          clusters[clusters == clust_a] <- ifelse(clust_combns[i, "dist"] < dcut, clust_b, clust_a)
        }
      }
    }
    ## Matrix of difference to reference per subclone and arm
    if(merge_method == "ad"){clusters <- sorted_clusts}
    if(merge_method == "ad" | merge_method == "both"){
      if(length(unique(clusters)) > 1){
        string_min <- 0
        while(string_min < 1 && length(unique(clusters)) > 1){
          amp_del_mats <- admat.fun(matrix, break_idx, clusters, labels, cut = adcut, genes = adgenes, region = adregion)
          arms_exp_mat <- amp_del_mats$epmat
          amp_del_mat <- amp_del_mats$ad
          clust_combns <- data.frame(t(combn(colnames(arms_exp_mat), 2)))
          if(isTRUE(armfilt)){
            amp_del_mat <- varfilt(arms_exp_mat, amp_del_mat, cut = NULL, n = arm_n)$ad
            arms_exp_mat <- varfilt(arms_exp_mat, amp_del_mat, cut = NULL, n = arm_n)$epmat
          }
          for(i in 1:nrow(clust_combns)){
            clust_a <- as.character(clust_combns[i, 1])
            clust_b <- as.character(clust_combns[i, 2])
            clust_combns[i, "string_dist"] <- sum(stringdist::stringdist(amp_del_mat[, clust_a], amp_del_mat[, clust_b]))
            clusters[clusters == clust_a] <- ifelse(clust_combns[i, "string_dist"] == 0, clust_b, clust_a)
            string_min <- min(clust_combns$string_dist)
          }
        }
      }
    }
  }

  ## Get back to full amp_del matrix if using only top arms
  if(isTRUE(armfilt) || isFALSE(merge) & length(unique(clusters)) > 1){
    amp_del_mats <- admat.fun(matrix, break_idx, clusters, labels, cut = adcut, genes = adgenes, region = adregion)
    arms_exp_mat <- amp_del_mats$epmat
    amp_del_mat <- amp_del_mats$ad
  }

  ## Rename merged clusters
  if(length(unique(clusters)) > 1){
    clust_sort <- sort(unique(clusters))
    amp_del_mat <- amp_del_mat[, clust_sort]
    arms_exp_mat <- arms_exp_mat[, clust_sort]
    ltrs <- LETTERS[1:length(clust_sort)]
    for(i in 1:length(clust_sort)){
      clusters <- ifelse(clusters == clust_sort[i], paste("subclone", ltrs[i], sep = "_"), clusters)
    }
    if(!is.null(name)){clusters <- paste(name, clusters, sep = "_")}
    names(clusters) <- rownames(trans_mat1)
    colnames(amp_del_mat) <- sort(unique(clusters))
    colnames(arms_exp_mat) <- sort(unique(clusters))

    out <- list(amp_del_mat = amp_del_mat, clusters = clusters, arms_exp_mat = arms_exp_mat, dist_mat = trans_mat3)
    }

  if(length(unique(clusters)) < 2){
    clusters <- paste(name, "subclone_1", sep = "_")
    out <- clusters
  }
  return(out)
}





#### Dived the CNA matrix to subclones by reducing feature dimension and clustering with Louvain:
spatial_subclones <- function(cna_matrix,
                              epi_cells,
                              separate = c("arm", "chr"),
                              genome = NULL,
                              genecut = 10,
                              top_region = 1/3,
                              top_method = c("abs", "sd"),
                              reduction_dims = 30,
                              cluster_k = 10){

  ## Subset by chromosome
  breaks <- genome_break(genes = rownames(cna_matrix), separate = separate,
                         genecut = genecut, genome = genome)
  break_idx <- breaks$break_idx
  labels <- breaks$breaks_df$labels
  dispose <- breaks$dispose
  if(length(dispose) != 0) {
    matrix <- cna_matrix[-dispose, ]
  } else {matrix <- cna_matrix}
  break_num <- unique(break_idx)

  ## Set working matrix
  matrix <- matrix[, epi_cells]
  names(break_idx) <- rownames(matrix)

  ## Select the most CNA-rich regions
  if(top_method == "abs"){
    abs_vals <- apply(matrix, 1, function(x) mean(abs(x)))
    mat <- matrix[abs_vals > quantile(abs_vals, probs = 1 - top_region), ]
  }
  ## or alternatively - Select the most variant gene regions
  if(top_method == "sd"){
    stand_dev <- apply(matrix, 1, sd)
    mat <- matrix[stand_dev > quantile(stand_dev, probs = 1 - top_region), ]
  }
  transpose_mat <- as.matrix(t(mat))

  # Dimred
  pca <- prcomp(transpose_mat, center = FALSE, scale. = FALSE)
  dim_red_mat <- pca$x[, 1:reduction_dims]
  rownames(dim_red_mat) <- colnames(matrix)

  ## Clustering methods:
  clusters <- cluster_coord(dim_red_mat, method = "louvain", louvain = cluster_k)
  names(clusters) <- colnames(matrix)
  return(clusters)
}


#### Create amplification/deletion matrix
make_ampdel_mat <- function(mat, break_idx, clusters, labels, genes = "top", region = 1/2, refmat = NULL) {
  break_num <- unique(break_idx)
  clust <- sort(unique(clusters))

  # Split matrix per subclone
  cut_mat <- lapply(clust, function(x) {
    spot_names <- names(clusters[clusters == x])
    clust_mat <- mat[, spot_names]})

  # Mean arm expression per subclone
  arm_exp_mat <- matrix(nrow = length(break_num), ncol = length(clust))
  for(i in 1:length(clust)) {
    for(j in 1:length(break_num)) {
      m1 <- cut_mat[[i]]
      m2 <- m1[which(break_idx == break_num[j]), ]
      abs_vals <- apply(m2, 1, function(x) mean(abs(x)))
      m3 <- m2[abs_vals > quantile(abs_vals, probs = 1 - region), ]
      if(!is.null(dim(m3))) {c <- mean(colMeans(m3))} else {c = 0}
      d <- mean(colMeans(m2))
      if(genes == "top") {arm_exp_mat[j, i] <- c}
      if(genes == "all") {arm_exp_mat[j, i] <- d}
      if(genes == "ref" && !is.null(refmat)) {
        mr <- refmat[break_idx == break_num[j], !colnames(refmat) %in% cells]
        mc <- mean(colMeans(mr))
        arm_exp_mat[j, i] <- c - mc}
    }
  }
  rownames(arm_exp_mat) <- labels
  colnames(arm_exp_mat) <- clust

  arm_exp_mat <- as.data.frame(arm_exp_mat)
  arm_exp_mat[is.na(arm_exp_mat)] <- 0
  return(arm_exp_mat)
}


#### Detect unique chromosomal amplification / deletion events
detect_unique_event <- function(amp_del_mat, conf_val = 4, var_arm_len = 10) {
  arm_var <- apply(amp_del_mat, 1, var)
  top_var_arms <- names(tail(sort(arm_var), var_arm_len))
  trans_arm_mat <- t(amp_del_mat)
  unique_clones <- list()
  for(i in seq_along(top_var_arms)) {
    chr_arm <- top_var_arms[[i]]
    low_OL_bound <- median(trans_arm_mat[, chr_arm]) - conf_val * mad(trans_arm_mat[, chr_arm], constant = 1)
    high_OL_bound <- median(trans_arm_mat[, chr_arm]) + conf_val * mad(trans_arm_mat[, chr_arm], constant = 1)
    unique_clones[[i]] <- list(deletion = rownames(trans_arm_mat)[which(trans_arm_mat[, chr_arm] < low_OL_bound)],
                               amplification = rownames(trans_arm_mat)[which(trans_arm_mat[, chr_arm] > high_OL_bound)])
    names(unique_clones)[i] <- chr_arm
    if(sum(S4Vectors::isEmpty(unique_clones[[i]])) == 2) {unique_clones[i] <- NULL}
  }
  unique_clones <- plyr::compact(unique_clones)
  return(unique_clones)
}



#### Optimized subclone detection function
ST_subclones <- function(cna_matrix,
                         epi_cells,
                         separate = "arm",
                         genome = hg38,
                         genecut = 10,
                         top_cna_region = 1/3,
                         top_arm_region = 1/2,
                         amp_del_genes = "top",
                         reduction_dims = 30,
                         resolution = 2,
                         cluster_k = 10,
                         event_conf_val = 4,
                         var_arm_len = 10,
                         expr_diff_cutoff = 0.15,
                         sum_diff_cutoff = 1.5,
                         cor_cutoff = 0.95,
                         max_iter = 20) {

  ## Subset by chromosome
  breaks <- genome_break(genes = rownames(cna_matrix), separate = separate, genecut = genecut, genome = genome)
  break_idx <- breaks$break_idx
  labels <- breaks$breaks_df$labels
  dispose <- breaks$dispose
  if(length(dispose) != 0) {
    cna_matrix <- cna_matrix[-dispose, ]
  } else {cna_matrix <- cna_matrix}
  break_num <- unique(break_idx)
  cna_matrix <- cna_matrix[, epi_cells]
  # names(break_idx) <- rownames(cna_matrix)

  ## Select the most CNA-rich regions
  abs_vals <- apply(cna_matrix, 1, function(x) mean(abs(x)))
  mat <- cna_matrix[abs_vals > quantile(abs_vals, probs = 1 - top_cna_region), ]

  transpose_mat <- as.matrix(t(mat))

  # Dimensionality reduction
  pca <- prcomp(transpose_mat, center = FALSE, scale. = FALSE)
  dim_red_mat <- pca$x[, 1:reduction_dims]
  rownames(dim_red_mat) <- colnames(mat)

  ## Over-clustering:
  clusters <- cluster_leiden(dim_red_mat, red_dim = FALSE, dims = NULL, resolution = resolution, clust_k = cluster_k, iter_n = 10)
  names(clusters) <- colnames(mat)
  levels(clusters) <- paste("Subclone", sort(as.integer(unique(clusters)), decreasing = FALSE), sep = "_")

  ## Create per over-clustered subclones amp/del matrix
  chr_arm_stats <- make_ampdel_mat(cna_matrix, break_idx, clusters, labels, genes = amp_del_genes, region = top_arm_region)

  ## Detect unique chromosomal events
  uniq_events <- detect_unique_event(chr_arm_stats, conf_val = event_conf_val, var_arm_len = var_arm_len)

  ## Check if clones sharing same unique event are similar and can be merged
  uniq_clusters <- list()
  merged_uniq_clusters <- list()
  for(i in seq_along(uniq_events)) {
    if(lengths(uniq_events[[i]]["deletion"]) >= 2 || lengths(uniq_events[[i]]["amplification"]) >= 2) {
      if(sum(lengths(uniq_events[[i]]) >= 2) == 1) {
        clones2check <- unname(unlist(uniq_events[[i]][which(lengths(uniq_events[[i]]) >= 2)]))
        if(!S4Vectors::isEmpty(uniq_events[[i]][which(lengths(uniq_events[[i]]) < 2)])) {
          uniq_clone <- unname(unlist(uniq_events[[i]][which(lengths(uniq_events[[i]]) < 2)]))
          if(!uniq_clone %in% levels(uniq_clusters) && !uniq_clone %in% attr(merged_uniq_clusters, "merged_subclones")){
            uniq_clusters <- unlist(list(uniq_clusters, droplevels(clusters[which(clusters %in% uniq_clone)])))
          }
        }
      }
      if(sum(lengths(uniq_events[[i]]) > 2) == 2) {
        clones2check <- list(del = unname(unlist(uniq_events[[i]][["deletion"]])),
                             amp = unname(unlist(uniq_events[[i]][["amplification"]])))
      }
      if(is.character(clones2check)) {
        clust_combns <- data.frame(t(combn(clones2check, 2)))
        for(j in seq_len(nrow(clust_combns))) {
          clust_a <- as.character(clust_combns[j, 1])
          clust_b <- as.character(clust_combns[j, 2])
          clust_combns[j, "dist"] <- max(abs(chr_arm_stats[, clust_a] - chr_arm_stats[, clust_b]))
          clust_combns[j, "sum"] <- sum(abs(chr_arm_stats[, clust_a] - chr_arm_stats[, clust_b]))
          clust_combns[j, "cor"] <- cor(x = chr_arm_stats[, clust_a], y = chr_arm_stats[, clust_b], method = "pearson")
        }
        avg_stats <- c(mean(clust_combns$dist), mean(clust_combns$sum), mean(clust_combns$cor))
        if(avg_stats[[1]] < expr_diff_cutoff && avg_stats[[2]] < sum_diff_cutoff && avg_stats[[3]] > cor_cutoff) {
          merge_clones <- factor(names(clusters[which(clusters %in% clones2check)]))
          levels(merge_clones) <- rep(paste0("Merged_", i), length(merge_clones))
          names(merge_clones) <- names(clusters[which(clusters %in% clones2check)])
          merged_uniq_clusters <- unlist(list(merged_uniq_clusters, merge_clones))
          attr(merged_uniq_clusters, "merged_subclones") <- clones2check
        }
        ## add conditions to check also pairs similarity (not only all clusters combined)
        if(avg_stats[[1]] > expr_diff_cutoff || avg_stats[[2]] > sum_diff_cutoff || avg_stats[[3]] < cor_cutoff) {
          uniq_clone <- clones2check[!clones2check %in% attr(merged_uniq_clusters, "merged_subclones") & !clones2check %in% levels(uniq_clusters)]
          uniq_clusters <- unlist(list(uniq_clusters, droplevels(clusters[which(clusters %in% uniq_clone)])))
        }
      }
      #if(is.list(clones2check)) {} ## should be filled
    }
    if(lengths(uniq_events[[i]]["deletion"]) < 2 && lengths(uniq_events[[i]]["amplification"]) < 2) {
      if(!S4Vectors::isEmpty(uniq_events[[i]]["deletion"])) {
        uniq_clone <- unname(unlist(uniq_events[[i]]["deletion"]))
        if(!uniq_clone %in% levels(uniq_clusters) && !uniq_clone %in% attr(merged_uniq_clusters, "merged_subclones")){
          uniq_clusters <- unlist(list(uniq_clusters, droplevels(clusters[which(clusters %in% uniq_clone)])))
        }
      }
      if(!S4Vectors::isEmpty(uniq_events[[i]]["amplification"])) {
        uniq_clone <- unname(unlist(uniq_events[[i]]["amplification"]))
        if(!uniq_clone %in% levels(uniq_clusters) && !uniq_clone %in% attr(merged_uniq_clusters, "merged_subclones")){
          uniq_clusters <- unlist(list(uniq_clusters, droplevels(clusters[which(clusters %in% uniq_clone)])))
        }
      }
    }
  }

  used_clusts <- unique(unlist(uniq_events))
  remaining_clusts <- droplevels(clusters[!clusters %in% used_clusts])

  iter <- 1
  repeat {
    chr_arm_stats <- make_ampdel_mat(cna_matrix, break_idx, remaining_clusts, labels, genes = amp_del_genes, region = top_arm_region)
    clust_combns <- data.frame(t(combn(colnames(chr_arm_stats), 2)))

    for(i in 1:nrow(clust_combns)) {
      clust_a <- as.character(clust_combns[i, 1])
      clust_b <- as.character(clust_combns[i, 2])
      clust_combns[i, "dist"] <- max(abs(chr_arm_stats[, clust_a] - chr_arm_stats[, clust_b]))
      clust_combns[i, "sum"] <- sum(abs(chr_arm_stats[, clust_a] - chr_arm_stats[, clust_b]))
      clust_combns[i, "cor"] <- cor(x = chr_arm_stats[, clust_a], y = chr_arm_stats[, clust_b], method = "pearson")
    }
    filt_clust_combns <- clust_combns[which(clust_combns$dist < expr_diff_cutoff & clust_combns$sum < sum_diff_cutoff & clust_combns$cor > cor_cutoff), ]
    max_cor <- filt_clust_combns[which.max(filt_clust_combns$cor), ]
    if(nrow(filt_clust_combns) > 1 && sum(is.na(max_cor)) == 0) {
      remaining_clusts <- setNames(stringr::str_replace_all(remaining_clusts, pattern = paste0("^", max_cor$X2, "$"), replacement = max_cor$X1), names(remaining_clusts))
      remaining_clusts <- factor(remaining_clusts, levels = paste("Subclone", sort(as.integer(gsub("Subclone_", "", unique(remaining_clusts))), decreasing = FALSE), sep = "_"))
    }
    else if(nrow(filt_clust_combns) == 1 && sum(is.na(max_cor)) != 0) {
      remaining_clusts <- setNames(stringr::str_replace_all(remaining_clusts, pattern = paste0("^", filt_clust_combns$X2, "$"), replacement = filt_clust_combns$X1), names(remaining_clusts))
      remaining_clusts <- factor(remaining_clusts, levels = paste("Subclone", sort(as.integer(gsub("Subclone_", "", unique(remaining_clusts))), decreasing = FALSE), sep = "_"))
      break
    }
    else if(sum(is.na(filt_clust_combns)) != 0 && sum(is.na(max_cor)) != 0) {break}
    iter <- iter + 1
    if(iter > max_iter) {break}
  }

  clusters <- unlist(list(uniq_clusters, merged_uniq_clusters, remaining_clusts))
  levels(clusters) <- paste("Subclone", seq(unique(clusters)), sep = "_")
  return(clusters)
}


### Assign malignant and non-malignant status to all cna sub-clones:
malig_status <- function(expr_mat,
                         subclones,
                         ref_cells,
                         KNN_correct = FALSE,
                         top_region = 1/3,
                         top_cells = 1/4,
                         thresh_method = "percentile",
                         sd_cut = 2.5,
                         class_cut = 0.99, ...) {

  if(!is.factor(subclones)) {stop("Subclones must be provided in as named factor, where the names are the spot IDs.")}
  def_ls <- list()
  plots_ls <- list()
  for(i in levels(subclones)) {
    query <- names(subclones[subclones == i])
    cna_score <- calc_cna(matrix = expr_mat, query = query, ref = ref_cells, genome = hg38, range = c(-3,3), window = 100, noise = 0.15, isLog = TRUE, per_chr = TRUE, scale = NULL, top_genes_num = NULL, verbose = FALSE)
    #cna_status <- define.cna(cna_score, query, ref, tumour_score = NULL, top_cells = 1/4, gs = NULL, no_ts = TRUE, name = NULL, print_ref = TRUE, extra_score = NULL)
    cna_status <- assign_cna_status(cna_matrix = cna_score, epi_cells = query, ref_cells = ref_cells, top_region = top_region, top_cells = top_cells,
                                    thresh_method = thresh_method, sd_cut = sd_cut, class_cut = class_cut, clone_name = i, ...)
    if(isTRUE(KNN_correct)) {
      def_ls[[i]] <- setNames(rownames(cna_status$def_mat), cna_status$def_mat$final_cna_status)
    }
    else if(isFALSE(KNN_correct)) {
      def_ls[[i]] <- setNames(rownames(cna_status$def_mat), cna_status$def_mat$cna_status)
    }
    plots_ls[[i]] <- cna_status$plots
  }
  assignments <- unlist(unname(def_ls))
  reference <- setNames(ref_cells, rep("Reference", length(ref_cells)))
  definitions <- c(assignments[!duplicated(assignments)], reference)
  # out <- list(definitions = definitions, plots = plots_ls)
  return(definitions)
}




# Cell to Cell Communication inference ------------------------------------


#### input: a list of metaprograms, setting argument "detect" to either `known` or `prediction` and the argument "compare" to either `within` or `between`
#### output: list of lists - if compare = within: a list of 4 entries (ligands, receptors, expected_interactions & existing_interactions) for each metaprograms. if compare = between: 2 entries Description and Interaction
detect_interactions <- function(metaprograms_ls,
                                detect = c("known", "prediction"),
                                compare = c("within", "between")){
  # Load ligand-receptor database
  lr_network <- readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
  lr_network <- lr_network %>% mutate(interaction = paste0(lr_network$from, sep = "-", lr_network$to))
  if(detect != "known" & detect != "prediction"){stop("The argument detect must be either 'known' (for high comfidence interactions) or 'prediction' (to includ also predicted interactions).")}
  if(detect == "known"){
    lr_network <- lr_network %>% filter(database == "ramilowski" | database == "kegg")
    ligands <- lr_network %>% pull(from) %>% unique()
    receptors <- lr_network %>% pull(to) %>% unique()
  }
  if(detect == "prediction"){
    ligands <- lr_network %>% pull(from) %>% unique()
    receptors <- lr_network %>% pull(to) %>% unique()
  }

  # detect ligands, receptors, expected interactions and existing interactions for each metaprogram
  if(compare != "within" & compare != "between"){stop("The argument compare must be either 'within' (to check for interaction within a metaprogram) or 'between' (to look for interactions between metaprograms).")}
  if(compare == "within"){
    metaprograms_int_ls <- lapply(metaprograms_ls, function(prog){
      prog_ligs <- ligands[ligands %in% prog]
      prog_receps <- receptors[receptors %in% prog]

      ligands_vec <- ifelse(length(prog_ligs) > 0, yes = str_c(prog_ligs, collapse = "|"), no = "NULL")
      a <- lr_network$interaction[grep(ligands_vec, scalop::substri(lr_network$interaction, pos = 1))]
      receptors_vec <- ifelse(length(prog_receps) > 0, yes = str_c(prog_receps, collapse = "|"), no = "NULL")
      b <- lr_network$interaction[grep(receptors_vec, scalop::substri(lr_network$interaction, pos = 2))]
      expected_ints <- unique(c(a, b))
      existing_ints <- intersect(a,b)

      prog_lr_int <- list(ligands = prog_ligs, receptors = prog_receps, expected_interactions = expected_ints, existing_interactions = existing_ints)
    })
  }
  if(compare == "between"){
    # detect interactions between metaprograms:
    mp_combns <- combn(names(metaprograms_ls), m = 2)
    inter_mp_ls <- list()
    for(i in 1:ncol(mp_combns)){
      mp_a <- mp_combns[1, i]
      mp_b <- mp_combns[2, i]
      mp_a_ligs <- ligands[ligands %in% unlist(metaprograms_ls[mp_a])]
      mp_a_receps <- receptors[receptors %in% unlist(metaprograms_ls[mp_a])]
      mp_b_ligs <- ligands[ligands %in% unlist(metaprograms_ls[mp_b])]
      mp_b_receps <- receptors[receptors %in% unlist(metaprograms_ls[mp_b])]

      ligs_a_vec <- ifelse(length(mp_a_ligs) > 0, yes = str_c(mp_a_ligs, collapse = "|"), no = "NULL")
      lig_a_ints <- lr_network$interaction[grep(ligs_a_vec, scalop::substri(lr_network$interaction, pos = 1))]
      receps_b_vec <- ifelse(length(mp_b_receps) > 0, yes = str_c(mp_b_receps, collapse = "|"), no = "NULL")
      recep_b_ints <- lr_network$interaction[grep(receps_b_vec, scalop::substri(lr_network$interaction, pos = 2))]
      ligs_a_receps_b <- intersect(lig_a_ints, recep_b_ints)

      ligs_b_vec <- ifelse(length(mp_b_ligs) > 0, yes = str_c(mp_b_ligs, collapse = "|"), no = "NULL")
      lig_b_ints <- lr_network$interaction[grep(ligs_b_vec, scalop::substri(lr_network$interaction, pos = 1))]
      receps_a_vec <- ifelse(length(mp_a_receps) > 0, yes = str_c(mp_a_receps, collapse = "|"), no = "NULL")
      recep_a_ints <- lr_network$interaction[grep(receps_a_vec, scalop::substri(lr_network$interaction, pos = 2))]
      ligs_b_receps_a <- intersect(lig_b_ints, recep_a_ints)

      if(length(ligs_a_receps_b) == 0 & length(ligs_b_receps_a) == 0){
        inter_mp_ls[[i]] <- NULL
        next
      } else if(length(ligs_a_receps_b) > 0 & length(ligs_b_receps_a) == 0){
        inter_mp_ls[[i]] <- list(Description = paste("Ligands from", mp_a, "- Receptors from", mp_b, sep = " "),
                                 Interactions = ligs_a_receps_b)
      } else if(length(ligs_a_receps_b) == 0 & length(ligs_b_receps_a) > 0){
        inter_mp_ls[[i]] <- list(Description = paste("Ligands from", mp_b, "- Receptors from", mp_a, sep = " "),
                                 Interactions = ligs_b_receps_a)
      } else {
        message_a <- paste("Interaction_A = Ligands from", mp_a, "- Receptors from", mp_b, sep = " ")
        message_b <- paste("Interaction_B = Ligands from", mp_b, "- Receptors from", mp_a, sep = " ")
        inter_mp_ls[[i]] <- list(Description = paste(message_a, message_b, sep = "      "),
                                 Interaction_A = ligs_a_receps_b,
                                 Interaction_B = ligs_b_receps_a)
      }
      names(inter_mp_ls)[i] <- paste(mp_a, mp_b, sep = "-")
    }
    inter_mp_ls <- inter_mp_ls[-which(sapply(inter_mp_ls, is.null))]
  }
}





# Differential Cell-Type Proportion Tests ---------------------------------

## Calculates and transforms cell type proportions
.trans_props <- function(states = states, sample = sample, transform = NULL) {
  if(is.null(transform)) transform <- "logit"

  tab <- table(sample, states)
  props <- tab/rowSums(tab)
  if(transform=="asin"){
    message("Performing arcsin square root transformation of proportions")
    prop.trans <- asin(sqrt(props))
  }
  else if(transform=="logit"){
    message("Performing logit transformation of proportions")
    props.pseudo <- (tab+0.5)/rowSums(tab+0.5)
    prop.trans <- log(props.pseudo/(1-props.pseudo))
  }

  return(list(Counts = t(tab), TransformedProps = t(prop.trans), Proportions = t(props)))
}


## Performs t-tests of transformed cell type proportions
.calc_prop_ttest <- function(prop.list=prop.list, design=design,
                             contrasts=contrasts, robust=robust, trend=trend,
                             sort=sort)
{
  prop.trans <- prop.list$TransformedProps
  prop <- prop.list$Proportions

  # Add check for fewer than 3 cell types
  # Robust eBayes doesn't work with fewer than 3 cell types
  if(nrow(prop.trans)<=2){
    message("Setting robust to FALSE for eBayes for less than 3 cell types")
    robust <- FALSE
  }

  fit <- limma::lmFit(prop.trans, design)
  fit.cont <- limma::contrasts.fit(fit, contrasts=contrasts)
  fit.cont <- limma::eBayes(fit.cont, robust=robust, trend=trend)

  # Get mean cell type proportions and relative risk for output
  # If no confounding variable included in design matrix
  if(length(contrasts)==2){
    fit.prop <- limma::lmFit(prop, design)
    z <- apply(fit.prop$coefficients, 1, function(x) x^contrasts)
    RR <- apply(z, 2, prod)
  }
  # If confounding variables included in design matrix exclude them
  else{
    new.des <- design[,contrasts!=0]
    fit.prop <- limma::lmFit(prop,new.des)
    new.cont <- contrasts[contrasts!=0]
    z <- apply(fit.prop$coefficients, 1, function(x) x^new.cont)
    RR <- apply(z, 2, prod)
  }

  fdr <- p.adjust(fit.cont$p.value[,1], method="BH")

  out <- data.frame(PropMean=fit.prop$coefficients, PropRatio=RR,
                    Tstatistic=fit.cont$t[,1], P.Value=fit.cont$p.value[,1],
                    FDR=fdr)
  if(sort){
    o <- order(out$P.Value)
    return(out[o,])
  }else{
    return(out)
  }
}


## Performs F-tests for transformed cell type proportions
.calc_prop_anova <- function(prop.list=prop.list, design=design, coef = coef,
                             robust=robust, trend=trend, sort=sort)
{
  prop.trans <- prop.list$TransformedProps
  prop <- prop.list$Proportions

  # Add check for fewer than 3 cell types
  # Robust eBayes doesn't work with fewer than 3 cell types
  if(nrow(prop.trans)<=2){
    message("Robust eBayes doesn't work with fewer than 3 cell types.
                    Setting robust to FALSE")
    robust <- FALSE
  }

  # get cell type mean proportions ignoring other variables
  # this assumes that the design matrix is not in Intercept format
  fit.prop <- limma::lmFit(prop, design[,coef])

  # Change design matrix to intercept format
  design[,1] <- 1
  colnames(design)[1] <- "Int"

  # Fit linear model taking into account all confounding variables
  fit <- limma::lmFit(prop.trans,design)

  # Get F statistics corresponding to group information only
  # You have to remove the intercept term for this to work
  fit <- limma::eBayes(fit[,coef[-1]], robust=robust, trend=trend)

  # Extract F p-value
  p.value <- fit$F.p.value
  # and perform FDR adjustment
  fdr <- p.adjust(fit$F.p.value, method="BH")

  out <- data.frame(PropMean=fit.prop$coefficients, Fstatistic= fit$F,
                    P.Value=p.value, FDR=fdr)
  if(sort){
    o <- order(out$P.Value)
    return(out[o,])
  }else{
    return(out)
  }
}


test_state_prop_diff <- function(states, sample, group, trend = FALSE, robust = TRUE, transform = "logit") {

  if(is.null(transform)) transform <- "logit"
  # Get transformed proportions
  prop.list <- .trans_props(states, sample, transform)

  # Calculate baseline proportions for each cluster
  baseline.props <- table(states)/sum(table(states))

  # Collapse group information
  group.coll <- table(sample, group)

  design <- matrix(as.integer(group.coll != 0), ncol=ncol(group.coll))
  colnames(design) <- colnames(group.coll)

  if(ncol(design) == 2) {
    message("group variable has 2 levels, t-tests will be performed")
    contrasts <- c(1,-1)
    out <- .calc_prop_ttest(prop.list, design, contrasts=contrasts,
                            robust=robust, trend=trend, sort=FALSE)
    out <- data.frame(BaselineProp=baseline.props,out)
    return(out[order(out$P.Value),])
  }
  else if(ncol(design)>=2){
    message("group variable has > 2 levels, ANOVA will be performed")
    coef <- seq_len(ncol(design))
    out <- .calc_prop_anova(prop.list, design, coef=coef, robust=robust,
                            trend=trend, sort=FALSE)
    out <- data.frame(BaselineProp=as.vector(baseline.props),out)
    return(out[order(out$P.Value),])
  }
}



# General plotting functions ----------------------------------------------


# Emulate ggplot2's default color palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}



elbow_plot <- function(pca_object, ndims = 30) {
  if(!is.list(pca_object)) {stop("pca_object must be list class having the pca matrix + sdev entries, at least")}
  pca_mat <- pca$x
  sdev <- pca$sdev
  if (ndims > length(sdev)) {
    warning("There is information only for ", length(sdev), " dimensions")
    ndims <- length(sdev)
  }
  data <- data.frame(dims = colnames(pca_mat)[1:ndims], sdev = sdev[1:ndims])
  data <- data.frame(dims = factor(colnames(pca_mat)[1:ndims], levels = seq(1:ndims)), sdev = sdev[1:ndims])
  data$dims <- sort(as.integer(gsub("PC", "", levels(data$dims))))
  plot <- ggplot(data, aes(x = dims, y = sdev)) +
    geom_point() + labs(x = "Dimensions", y = "Standard Deviation") + cowplot::theme_cowplot()
  return(plot)
}


co_occur_bars <- function(co_occur_tab, state) {
  stopifnot(state %in% co_occur_tab$spot_type)

  nighbs <- co_occur_tab %>% filter(spot_type == state)
  occur_freqs <- reshape2::melt(sort(table(unlist(nighbs[1:6])), decreasing = TRUE))

  ggplot(data = occur_freqs, aes(x = Var1, y = value)) +
    geom_bar(stat = "identity", fill = "lightsteelblue3") + labs(title = paste0(state, " Neighboring Spots"), x = "", y ="Frequency of co-occurrence") +
    theme_minimal() + geom_text(aes(label = value), position = position_dodge(width = 0.9), vjust = -0.25, color = "firebrick")
}


ggshapes = function(shapeVar,
                    colVar,
                    obs=NULL,
                    dir=c('h','v'),
                    cols = stats::setNames(c('red','blue'), c('G34R','WT')),
                    shapes = stats::setNames(1:3, c('LTC42','LTC77','LTC99')),
                    point.size=3) {

  L = list(i = obs,col = colVar, shape = shapeVar)
  lens = sapply(L, length, simplify = F)
  lens = lens[lens!=0]
  lens = unlist(lens)

  stopifnot(length(lens)>0 & length(unique(lens))==1)

  len = unique(unlist(lens))
  if (is.null(obs)) obs = 1:len
  if (is.null(shapeVar)) shapeVar = rep('',len)
  if (is.null(colVar)) colVar = rep('', len)

  dir = match.arg(dir)
  d = data.frame(id = obs,
                 shapeVar=shapeVar,
                 colVar=colVar,
                 stringsAsFactors=F)

  d = d %>% dplyr::mutate(id = factor(id, levels = unique(id)))

  if (dir=='h') {G = ggplot(d, aes(y=1,x=id,shape=shapeVar,colour=colVar))}
  else {G = ggplot(d, aes(y=id,x=1,shape=shapeVar,colour=colVar))}


  G + geom_point(size=point.size) +
    scale_shape_manual(values=shapes) +
    scale_colour_manual(values=cols) +
    theme_void() +
    theme(legend.position='top',
          plot.margin=margin(0,0,0,0,'cm')) +
    guides(shape=guide_legend(nrow=1),
           color=guide_legend(nrow=3))
}



ggbar = function(colVar,
                 legend_title = NULL,
                 obs=NULL,
                 dir=c('h','v'),
                 cols = c('blue','red','orange','green','magenta'),
                 set_lims = NULL) {

  L = list(i = obs,col = colVar)
  lens = sapply(L, length, simplify = F)
  lens = lens[lens!=0]
  lens = unlist(lens)

  stopifnot(length(lens)>0 & length(unique(lens))==1)

  len = unique(unlist(lens))
  if (is.null(obs)) obs = 1:len
  if (is.null(colVar)) colVar = rep('', len)

  dir = match.arg(dir)
  d = data.frame(id = obs,
                 colVar=colVar,
                 stringsAsFactors=F)

  d = d %>% dplyr::mutate(id = factor(id, levels = unique(id)))

  if (dir=='h') {G = ggplot(d, aes(y=1,x=id,fill=colVar))}
  else {G = ggplot(d, aes(y=id,x=1,fill=colVar))}

  if (is.numeric(d$colVar)) {
    if (min(d$colVar) >= 0) {
      custom_pal <- c(grDevices::colorRampPalette(c(colorspace::lighten("antiquewhite", .75),
                                                    "antiquewhite",
                                                    rev(viridis::magma(30, begin = .25, end = .9))[1]))(5), rev(viridis::magma(30, begin = .25, end = .9)))
    }
    if (min(d$colVar) < 0) {
      custom_pal = c("#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f")
    }
    if (is.null(set_lims)) {lims = c(0, max(d$colVar))} else {lims = set_lims}
    grad_steps <- (max(lims) - min(lims)) / (length(custom_pal) - 1)
    col_vals <- seq(from = min(lims), to = max(lims), by = grad_steps)
    scale_fill <- scale_fill_gradientn(limits = lims, colours = custom_pal, values = scales::rescale(col_vals),
                                       oob = scales::squish, name = legend_title, guide = guide_legend(override.aes = list(size = 10)))
  } else {
    scale_fill <- scale_fill_manual(values=cols, name = legend_title, guide = guide_legend(override.aes = list(size = 10)))
  }

  G + geom_tile() +
    scale_fill +
    theme_void() +
    theme(legend.position='top',
          plot.margin=margin(0,0,0,0,'cm'))
}


### Plot enrichment results as dotplot
enrichment_plot <- function(enricher_res, x = "GeneRatio", color = "p.adjust",
                                 showCategory = 10, size = NULL, split = NULL,
                                 font.size = 12, title = "", orderBy = "x",
                                 label_format = 30, decreasing = TRUE) {

  colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
    if (is.null(size))
      size <- "Count"
  } else if (x == "count" || x == "Count") {
    x <- "Count"
    if (is.null(size))
      size <- "GeneRatio"
  } else if (is(x, "formula")) {
    x <- as.character(x)[2]
    if (is.null(size))
      size <- "Count"
  } else {
    ## message("invalid x, setting to 'GeneRatio' by default")
    ## x <- "GeneRatio"
    ## size <- "Count"
    if (is.null(size))
      size  <- "Count"
  }

  # df <- fortify(enricher_res, showCategory = showCategory, split = split)
  df <- enricher_res[1:showCategory, ]
  df <- df %>% dplyr::mutate(p.adjust = signif(p.adjust, digits = 2))
  ## already parsed in fortify
  ## df$GeneRatio <- parse_ratio(df$GeneRatio)

  if (orderBy !=  'x' && !orderBy %in% colnames(df)) {
    message('wrong orderBy parameter; set to default `orderBy = "x"`')
    orderBy <- "x"
  }

  if (orderBy == "x") {
    df <- dplyr::mutate(df, x = eval(parse(text=x)))
  }

  if (isTRUE(title)) {
    title <- unique(df$name)
  } else {title <- title}

  # label_func <- default_labeller(label_format)
  # if(is.function(label_format)) {
  #   label_func <- label_format
  # }

  idx <- order(df[[orderBy]], decreasing = decreasing)
  df$Description <- factor(df$Description,
                           levels=rev(unique(df$Description[idx])))
  ggplot(df, aes_string(x=x, y="Description", size=size, color=colorBy)) +
    geom_point() +
    scale_color_continuous(low="#d6604d", high="#4393c3", name = color,
                           guide=guide_colorbar(reverse=TRUE)) +
    # scale_y_discrete(labels = label_func) +
    scale_y_discrete() +
    ylab(NULL) + ggtitle(title) + theme_bw(font.size) +
    scale_size(range=c(3, 8)) +
    guides(size  = guide_legend(order = 1),
           color = guide_colorbar(order = 2))
}


### Plot UMAP
plot_umap <- function(metadata, fill_var = NULL, colors = NULL, ...) {
  # check input:
  stopifnot("Error: metadata must contain the UMAP coordinates under the names 'UMAP1' & 'UMAP2'" = all(c("UMAP1", "UMAP2") %in% colnames(metadata)))
  if(!is.null(fill_var)) {
    stopifnot("Error: cannot find the variable to be used for filling UMAP dots, in the metadata columns" = fill_var %in% colnames(metadata))
  }

  if(is.null(colors) && !is.null(fill_var)) {
    colors <- gg_color_hue(length(unique(metadata[[fill_var]])))
  }
  if(is.null(fill_var)) {
    G <- ggplot(data = metadata, aes(x = UMAP1, y = UMAP2)) + geom_point(...) + xlab("UMAP1") + ylab("UMAP2") +
      theme_classic() +
      theme(plot.title = element_text(size = 20, face = "bold"),legend.title=element_blank(), legend.text=element_text(size=16))
  } else {
    G <- ggplot(data = metadata, aes(x = UMAP1, y = UMAP2, color = .data[[fill_var]])) + geom_point(...) + xlab("UMAP1") + ylab("UMAP2") +
      theme_classic() + scale_color_manual(values = colors) +
      theme(plot.title = element_text(size = 20, face = "bold"),legend.title=element_blank(), legend.text=element_text(size=16), legend.position = "top")
  }
  return(G)
}



### Heatmap Functions (Seurat wrapper to enable input from expression matrix, instead of Seurat object):

plot_heatmap <- function(
  mat,
  metadata,
  group.by,
  features = NULL,
  var_genes = 3000,
  cells = NULL,
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  midpoint = 0,
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  sep_line_col = "#000000",
  palette = "dark",
  group.bar.height = 0.02,
  leg.key.size = 1,
  combine = TRUE,
  show_features = FALSE,
  anno_leg_name = NULL,
  expr_leg_name = "Relative\nExpression"
) {

  if(is.null(cells)) {
    cells <- colnames(mat)
  } else {cells <- cells}
  if (is.numeric(cells)) {
    cells <- colnames(mat)[cells]
  }

  if(is.null(features)) {
    s1 <- apply(mat, 1, sd)
    features <- names(tail(sort(s1), n = var_genes))
  } else {features <- features}
  # features <- rev(x = unique(x = features))

  # make sure features are present
  possible.features <- rownames(mat)
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if(length(features) == 0) {
      stop("No requested features found in the expression matrix provided")
    }
    warning("The following features were omitted as they were not found in the expression matrix provided:", paste(bad.features, collapse = ", "))
  }

  data <- as.data.frame(as.matrix(t(mat[features, cells, drop = FALSE])))
  groups.use <- data.frame(metadata[, group.by], row.names = colnames(mat))[cells, , drop = FALSE]
  if(is.null(anno_leg_name)) {anno_leg_name <- "Identity"}

  plots <- vector(mode = 'list', length = ncol(x = groups.use))
  for (i in 1:ncol(x = groups.use)) {
    data.group <- data
    group.use <- groups.use[, i, drop = TRUE]
    if (!is.factor(x = group.use)) {
      group.use <- factor(x = group.use)
    }
    names(x = group.use) <- cells
    if (draw.lines) {
      # create fake cells to serve as the white lines, fill with NAs
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 0.0025)
      placeholder.cells <- sapply(
        X = 1:(length(x = levels(x = group.use)) * lines.width),
        FUN = function(x) {
          return(SeuratObject::RandomName(length = 20))
        }
      )
      placeholder.groups <- rep(x = levels(x = group.use), times = lines.width)
      group.levels <- levels(x = group.use)
      names(x = placeholder.groups) <- placeholder.cells
      group.use <- as.vector(x = group.use)
      names(x = group.use) <- cells
      group.use <- factor(x = c(group.use, placeholder.groups), levels = group.levels)
      na.data.group <- matrix(
        data = NA,
        nrow = length(x = placeholder.cells),
        ncol = ncol(x = data.group),
        dimnames = list(placeholder.cells, colnames(x = data.group))
      )
      data.group <- rbind(data.group, na.data.group)
    }

    lgroup <- length(levels(group.use))
    plot <- singleRasterMap(
      data = data.group,
      raster = raster,
      disp.min = disp.min,
      disp.max = disp.max,
      midpoint = midpoint,
      feature.order = features,
      sep_line_col = sep_line_col,
      palette = palette,
      cell.order = names(x = sort(x = group.use)),
      group.by = group.use,
      expr_leg_name = expr_leg_name,
      leg.key.size = leg.key.size
    )
    if (group.bar) {
      # TODO: Change group.bar to annotation.bar
      default.colors <- c(hue_pal()(length(x = levels(x = group.use))))
      if (!is.null(x = names(x = group.colors))) {
        cols <- unname(obj = group.colors[levels(x = group.use)])
      } else {
        cols <- group.colors[1:length(x = levels(x = group.use))] %||% default.colors
      }
      if (any(is.na(x = cols))) {
        cols[is.na(x = cols)] <- default.colors[is.na(x = cols)]
        cols <- .Col2Hex(cols)
        col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(
          x = cols,
          start = 1,
          stop = 7
        )))))
        through <- length(x = default.colors)
        while (length(x = col.dups) > 0) {
          pal.max <- length(x = col.dups) + through
          cols.extra <- hue_pal()(pal.max)[(through + 1):pal.max]
          cols[col.dups] <- cols.extra
          col.dups <- sort(x = unique(x = which(x = duplicated(x = substr(
            x = cols,
            start = 1,
            stop = 7
          )))))
        }
      }
      group.use2 <- sort(x = group.use)
      if (draw.lines) {
        na.group <- SeuratObject::RandomName(length = 20)
        levels(x = group.use2) <- c(levels(x = group.use2), na.group)
        group.use2[placeholder.cells] <- na.group
        cols <- c(cols, sep_line_col)
      }
      pbuild <- ggplot_build(plot = plot)
      names(x = cols) <- levels(x = group.use2)
      # scale the height of the bar
      y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
      y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
      y.max <- y.pos + group.bar.height * y.range
      x.min <- min(pbuild$layout$panel_params[[1]]$x.range) + 0.1
      x.max <- max(pbuild$layout$panel_params[[1]]$x.range) - 0.1
      plot <- plot +
        annotation_raster(
          raster = t(x = cols[group.use2]),
          xmin = x.min,
          xmax = x.max,
          ymin = y.pos,
          ymax = y.max
        ) +
        coord_cartesian(ylim = c(0, y.max), clip = 'off') +
        scale_color_manual(name = anno_leg_name, na.translate = FALSE, values = cols[-length(cols)])
      # scale_color_discrete(name = anno_leg_name, na.translate = FALSE)
      if (label) {
        x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
        # Attempt to pull xdivs from x.major in ggplot2 < 3.3.0; if NULL, pull from the >= 3.3.0 slot
        x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% attr(x = pbuild$layout$panel_params[[1]]$x$get_breaks(), which = "pos")
        x <- data.frame(group = sort(x = group.use), x = x.divs)
        label.x.pos <- tapply(X = x$x, INDEX = x$group, FUN = function(y) {
          if (isTRUE(x = draw.lines)) {
            mean(x = y[-length(x = y)])
          } else {
            mean(x = y)
          }
        })
        label.x.pos <- data.frame(group = names(x = label.x.pos), label.x.pos)
        plot <- plot + geom_text(
          stat = "identity",
          data = label.x.pos,
          aes_string(label = 'group', x = 'label.x.pos'),
          y = y.max + y.max * 0.03 * 0.5,
          angle = angle,
          hjust = hjust,
          size = size
        )
        plot <- suppressMessages(plot + coord_cartesian(
          ylim = c(0, y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use))) * size),
          clip = 'off')
        )
      }
    }
    plot <- plot + theme(line = element_blank()) + guides(fill = guide_colorbar(frame.colour = "black", ticks.colour = "black"))
    if(show_features) {
      plot <- plot
    } else {plot <- plot + theme(axis.text.y = element_blank())}
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- wrap_plots(plots)
  }
  return(plots)
}



singleRasterMap <- function(
  data,
  raster = TRUE,
  cell.order = NULL,
  feature.order = NULL,
  sep_line_col = "#000000",
  palette = c("dark", "bright"),
  disp.min = -2.5,
  disp.max = 2.5,
  limits = NULL,
  midpoint = 0,
  group.by = NULL,
  expr_leg_name = "Relative\nExpression",
  leg.key.size = 1
) {

  data <- Seurat::MinMax(data = data, min = disp.min, max = disp.max)
  data <- Melt(x = t(x = data))
  colnames(x = data) <- c('Feature', 'Cell', 'Expression')
  if (!is.null(x = feature.order)) {
    data$Feature <- factor(x = data$Feature, levels = unique(x = feature.order))
  }
  if (!is.null(x = cell.order)) {
    data$Cell <- factor(x = data$Cell, levels = unique(x = cell.order))
  }
  if (!is.null(x = group.by)) {
    data$Identity <- group.by[data$Cell]
  }
  limits <- limits %||% c(min(data$Expression), max(data$Expression))
  if (length(x = limits) != 2 || !is.numeric(x = limits)) {
    stop("limits' must be a two-length numeric vector")
  }

  if(palette == "bright") {
    cols = c("#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f")
    # cols = c("#053061", "#2166ac", "#4393c3", "#92c5de", "#f7f7f7", "#f4a582", "#d6604d", "#b2182b", "#67001f")
    if(length(midpoint) == 1 && midpoint == 0) {
      grad_steps <- (disp.max - disp.min) / (length(cols) - 1)
      col_vals <- seq(from = disp.min, to = disp.max, by = grad_steps)
    } else if (length(midpoint) == 1 && midpoint != 0) {
      zero_col_val <- ceiling(length(cols) / 2)
      low_grad_steps <- (midpoint - disp.min) / (zero_col_val - 1)
      high_grad_steps <- (disp.max - midpoint) / (zero_col_val - 1)
      col_vals <- c(seq(from = disp.min, to = midpoint, by = low_grad_steps), seq(from = midpoint, to = disp.max, by = high_grad_steps)[-1])
    } else if (length(midpoint) > 1) {
      zero_col_val <- ceiling(length(cols) / 2)
      cols <- c(cols[1:(zero_col_val - 1)], rep(cols[zero_col_val], 2), cols[(zero_col_val + 1):length(cols)])
      low_grad_steps <- (min(midpoint) - disp.min) / (zero_col_val - 1)
      high_grad_steps <- (disp.max - max(midpoint)) / (zero_col_val - 1)
      col_vals <- c(seq(from = disp.min, to = min(midpoint), by = low_grad_steps),
                    seq(from = max(midpoint), to = disp.max, by = high_grad_steps))
    }
    scale_fill <- ggplot2::scale_fill_gradientn(colors = cols, values = scales::rescale(col_vals), limits = limits, oob = scales::squish, name = expr_leg_name, na.value = sep_line_col)
  } else {
    scale_fill <- scale_fill_gradient2(limits = limits, low = "dodgerblue4", mid = "white", high = "red4", midpoint = 0, na.value = sep_line_col, oob = squish, name = expr_leg_name)
  }

  my_geom <- ifelse(test = raster, yes = geom_raster, no = geom_tile)
  plot <- ggplot(data = data) +
    my_geom(mapping = aes_string(x = 'Cell', y = 'Feature', fill = 'Expression')) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    eval(scale_fill) +
    labs(x = NULL, y = NULL, fill = group.by %iff% 'Expression') +
    Seurat::WhiteBackground() + Seurat::NoAxes(keep.text = TRUE)
  if (!is.null(x = group.by)) {
    plot <- plot + geom_point(
      mapping = aes_string(x = 'Cell', y = 'Feature', color = 'Identity'),
      alpha = 0
    ) +
      guides(color = guide_legend(override.aes = list(alpha = 1, size = leg.key.size)))
  }
  return(plot)
}

.Col2Hex <- function(...) {
  colors <- as.character(x = c(...))
  alpha <- rep.int(x = 255, times = length(x = colors))
  if (sum(sapply(X = colors, FUN = grepl, pattern = '^#')) != 0) {
    hex <- colors[which(x = grepl(pattern = '^#', x = colors))]
    hex.length <- sapply(X = hex, FUN = nchar)
    if (9 %in% hex.length) {
      hex.alpha <- hex[which(x = hex.length == 9)]
      hex.vals <- sapply(X = hex.alpha, FUN = substr, start = 8, stop = 9)
      dec.vals <- sapply(X = hex.vals, FUN = strtoi, base = 16)
      alpha[match(x = hex[which(x = hex.length == 9)], table = colors)] <- dec.vals
    }
  }
  colors <- t(x = col2rgb(col = colors))
  colors <- mapply(
    FUN = function(i, alpha) {
      return(rgb(colors[i, , drop = FALSE], alpha = alpha, maxColorValue = 255))
    },
    i = 1:nrow(x = colors),
    alpha = alpha
  )
  return(colors)
}

`%iff%` <- function(x, y) {
  if (!is_null(x = x)) {
    return(y)
  }
  return(x)
}




## Plot boxplot or violin plot with statistics
plot_t_stats <- function(data, x, y, geom = c("boxplot", "violin"), scale_bar_height = 1) {
  if(length(geom) != 1) {geom <- "boxplot"}
  if(geom == "boxplot") {
    p <- ggplot(data = data, aes_string(x = x, y = y)) +
      geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.25, size = 0.2, alpha = 0.4) +
      scalop::theme_scalop()
  } else if (geom == "violin") {
    p <- ggplot(data = data, aes_string(x = x, y = y)) +
      geom_violin(fill = "lightblue", draw_quantiles = c(0.25, 0.75), linetype = "dashed") + geom_violin(fill = "transparent", draw_quantiles = 0.5, linewidth = 1) +
      scalop::theme_scalop()
  } else {stop("geom must be either boxplot or violin")}

  all_pairs <- combn(as.character(unique(na.omit(data[[x]]))), 2)
  effect_size <- lapply(seq_len(ncol(all_pairs)), function(idx) {
    paste0("Cohen's d: ",
           round(abs(effectsize::cohens_d(data[[y]][data[[x]] == all_pairs[1, idx]], data[[y]][data[[x]] == all_pairs[2, idx]])$Cohens_d), digits = 3))
  })
  names(effect_size) <- lapply(seq_len(ncol(all_pairs)), function(idx) paste0(pmin(all_pairs[1, idx], all_pairs[2, idx]), ".", pmax(all_pairs[1, idx], all_pairs[2, idx])))
  anno_pos <- setNames(attributes(ggplot_build(p)$layout$panel_params[[1]]$x$get_breaks())$pos, ggplot_build(p)$layout$panel_params[[1]]$x$get_labels())

  fmla <- as.formula(paste0(y, " ~", x))
  x <- syms(x)
  y <- syms(y)

  if(ncol(all_pairs) == 1){
    anno_df <- compare_means(fmla, data = data, method = "t.test") %>%
      dplyr::mutate(y.position = data %>% dplyr::summarise(maxi = max(!!!y) * scale_bar_height) %>% pull(maxi),
                    x.position1 = anno_pos[match(.$group1, names(anno_pos))], x.position2 = anno_pos[match(.$group2, names(anno_pos))],
                    p.format = paste0("p < ", gsub(pattern = "<", replacement = "", .$p.format)),
                    effect_size = as.character(effect_size[match(paste0(pmin(.$group1, .$group2), ".", pmax(.$group1, .$group2)), names(effect_size))]),
                    y = y.position - 1.5, x = (x.position1 + x.position2) / 2)
    p_full <- p + ggpubr::stat_pvalue_manual(data = anno_df, label = "p.format") + geom_text(data = anno_df, aes(x = x, y = y, label = effect_size))
  } else {
    anno_df <- compare_means(fmla, data = data, method = "t.test") %>%
      dplyr::mutate( y.position = seq(from = max(data %>% dplyr::pull(!!!y)), by = 7 * scale_bar_height, length.out = nrow(.)),
                     x.position1 = anno_pos[match(.$group1, names(anno_pos))], x.position2 = anno_pos[match(.$group2, names(anno_pos))],
                     p.format = paste0("p < ", gsub(pattern = "<", replacement = "", .$p.format)),
                     effect_size = as.character(effect_size[match(paste0(pmin(.$group1, .$group2), ".", pmax(.$group1, .$group2)), names(effect_size))]),
                     y = y.position - 1.5, x = (x.position1 + x.position2) / 2)
    p_full <- p + ggpubr::stat_pvalue_manual(data = anno_df, label = "p.format") + geom_text(data = anno_df, aes(x = x, y = y, label = effect_size))
  }

  return(p_full)
}

# Nested Pie chart - https://stackoverflow.com/questions/26748069/ggplot2-pie-and-donut-chart-on-same-plot?noredirect=1&lq=1
donuts <- function(x, group = 1, labels = NA, col = NULL, radius = c(.7, 1)) {
  group <- rep_len(group, length(x))
  ug  <- unique(group)
  tbl <- table(group)[order(ug)]

  col <- if (is.null(col))
    seq_along(ug) else rep_len(col, length(ug))
  col.main <- Map(rep, col[seq_along(tbl)], tbl)
  col.sub  <- lapply(col.main, function(x) {
    al <- head(seq(0, 1, length.out = length(x) + 2L)[-1L], -1L)
    Vectorize(adjustcolor)(x, alpha.f = al)
  })

  plot.new()

  par(new = TRUE)
  pie(x, border = NA, radius = radius[2L],
      col = unlist(col.sub), labels = labels)

  par(new = TRUE)
  pie(x, border = NA, radius = radius[1L],
      col = unlist(col.main), labels = NA)
}



# Statistics functions and summaries --------------------------------------

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )

  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}

get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

