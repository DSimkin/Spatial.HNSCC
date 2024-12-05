library(Matrix)
library(tidyverse)
library(Seurat)
library(rjson)
library(cowplot)
library(RColorBrewer)
library(grid)



# Spatial Object Construction ---------------------------------------------


CreateSeuratSpatial <- function(
  dgCMatrix,
  data.dir,
  assay = 'Spatial',
  slice = 'slice1',
  filter.matrix = TRUE,
  to.upper = FALSE,
  ...
) {
  if(class(dgCMatrix)[[1]] == "dgCMatrix"){
    data <- dgCMatrix
  } else if(class(dgCMatrix)[[1]] == "matrix"){
    data <- Matrix::Matrix(dgCMatrix)
  } else {warning("Input has to be of either matrix or dgCMatrix class!")}
  object <- CreateSeuratObject(counts = data, assay = assay)
  image <- Read10X_Image(
    image.dir = file.path(data.dir, 'spatial'),
    filter.matrix = filter.matrix
  )
  image <- image[Cells(x = object)]
  DefaultAssay(object = image) <- assay
  object[[slice]] <- image
  return(object)
}


Read10X_Image_to_Seurat <- function(image.dir, image.name = "tissue_lowres_image.png", filter.matrix = TRUE, ...) {
  dir_files <- list.files(image.dir)
  image <- png::readPNG(file.path(image.dir, dir_files[grep(image.name, dir_files)]))
  scale.factors <- jsonlite::fromJSON(txt = file.path(image.dir, dir_files[grep("scalefactors_json.json", dir_files)]))
  tissue.positions <- read.csv(
    file = file.path(image.dir, dir_files[grep("tissue_positions_list.csv", dir_files)]),
    col.names = c('barcodes', 'tissue', 'row', 'col', 'imagerow', 'imagecol'),
    header = FALSE,
    as.is = TRUE,
    row.names = 1
  )
  if (filter.matrix) {
    tissue.positions <- tissue.positions[which(x = tissue.positions$tissue == 1), , drop = FALSE]
  }
  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
  spot.radius <-  unnormalized.radius / max(dim(x = image))
  return(new(
    Class = 'VisiumV1',
    image = image,
    scale.factors = scalefactors(
      spot = scale.factors$tissue_hires_scalef,
      fiducial = scale.factors$fiducial_diameter_fullres,
      hires = scale.factors$tissue_hires_scalef,
      scale.factors$tissue_lowres_scalef
    ),
    coordinates = tissue.positions,
    spot.radius = spot.radius
  ))
}


load_spatial_image <- function(path_to_image_dir, image_name = "tissue_lowres_image.png", filter = TRUE) {
  dir_files <- list.files(path_to_image_dir)
  image <- png::readPNG(file.path(path_to_image_dir, dir_files[grep(image_name, dir_files)]))
  scale_factors <- jsonlite::fromJSON(file.path(path_to_image_dir, dir_files[grep("scalefactors_json.json", dir_files)]))
  scale_factors <- scale_factors[c("spot_diameter_fullres", "tissue_hires_scalef", "fiducial_diameter_fullres", "tissue_lowres_scalef")]
  names(scale_factors) <- c("spot", "hires", "fiducial", "lowres")
  spot_positions <- read.csv(file.path(path_to_image_dir, dir_files[grep("tissue_positions_list.csv", dir_files)]), header = FALSE)
  colnames(spot_positions) <- c("SpotID", "in_tissue", "array_row", "array_col", "pxl_in_rows", "pxl_in_cols")

  if(filter) {
    spot_positions <- spot_positions[which(spot_positions$in_tissue == 1), , drop = FALSE]
  }

  radius <- scale_factors$fiducial * scale_factors$lowres
  spot_radius <-  radius / max(dim(image))

  return(list(Image = image, Scale = scale_factors, Coordinates = spot_positions, Spot_Radius = spot_radius))
}



# Metaprogram Clustering Algorithm ----------------------------------------

cluster_metaprograms <- function(progs_signatures, progs_intersect, nmf_wBasis, min_init_intersect = 12, min_clust_intersect = 9, min_group_size = 5) {

  # Rank intersections top to bottom
  sorted_intersection <- sort(apply(progs_intersect, 2, function(x) length(which(x >= min_init_intersect)) - 1), decreasing = TRUE)

  # Instantiate objects to be updated in each iteration
  cluster_list <- list()
  k <- 1
  curr_cluster <- c()
  MP_list <- list()

  # Iterate through the ranked intersection vector and break it into clusters
  while(sorted_intersection[1] > min_group_size) {
    curr_cluster <- c(curr_cluster, names(sorted_intersection[1]))

    # Intersection between all remaining gene programs and genes in current metaprogram
    genes_MP <- progs_signatures[, names(sorted_intersection[1])]     # initial genes are those in the first (most intersecting) program. 'genes_MP' always has only 50 genes consisting of the current MP
    progs_signatures <- progs_signatures[, -match(names(sorted_intersection[1]), colnames(progs_signatures))]     # remove selected program
    intersection_with_genes_MP <- sort(apply(progs_signatures, 2, function(x) length(intersect(genes_MP, x))), decreasing = TRUE)     # intersection between all other programs and 'genes_MP'
    program_history <- genes_MP       # has all genes in all programs of the current cluster, for newly defining 'genes_MP' after adding a new program

    # Create gene list - composed of intersecting genes in descending order + genes with highest program rank to add up to 50 genes. Update 'curr_cluster' each time
    while(intersection_with_genes_MP[1] >= min_clust_intersect) {
      curr_cluster <- c(curr_cluster, names(intersection_with_genes_MP)[1])
      genes_MP_temp <- sort(table(c(program_history, progs_signatures[, names(intersection_with_genes_MP)[1]])), decreasing = TRUE)     # 'genes_MP' is newly defined each time according to all programs in the current cluster
      genes_at_border <- genes_MP_temp[which(genes_MP_temp == genes_MP_temp[50])]     # genes with overlap equal to the 50th gene

      if(length(genes_at_border) > 1) {
        # To resolve overlaps with border genes - we use the NMF scores from the NMF programs:
        # Sort last genes in 'genes_at_border' according to maximal NMF gene scores
        # Run over all NMF programs in 'curr_cluster' and extract NMF scores for each gene
        curr_genes_NMF_score <- c()
        for (i in curr_cluster) {
          curr_study <- strsplit(i, "[.]")[[1]][[1]]
          boarder_genes_scores <- nmf_wBasis[[curr_study]][names(genes_at_border)[!is.na(match(names(genes_at_border), rownames(nmf_wBasis[[curr_study]])))], i]
          curr_genes_NMF_score <- c(curr_genes_NMF_score, boarder_genes_scores)
        }
        curr_genes_NMF_score_sort <- sort(curr_genes_NMF_score, decreasing = TRUE)
        curr_genes_NMF_score_sort <- curr_genes_NMF_score_sort[unique(names(curr_genes_NMF_score_sort))]      # take only the maximal score of each gene - which is the first entry after sorting
        genes_MP_temp <- c(names(genes_MP_temp[which(genes_MP_temp > genes_MP_temp[50])]), names(curr_genes_NMF_score_sort))
      } else {
        genes_MP_temp <- names(genes_MP_temp)[1:50]
      }

      program_history <- c(program_history, progs_signatures[, names(intersection_with_genes_MP)[1]])
      genes_MP <- genes_MP_temp[1:50]
      progs_signatures <- progs_signatures[, -match(names(intersection_with_genes_MP)[1], colnames(progs_signatures))]      # remove selected program
      intersection_with_genes_MP <- sort(apply(progs_signatures, 2, function(x) length(intersect(genes_MP, x))), decreasing = TRUE)     # intersection between all other NMFs and Genes_MP
    }

    # Update metaprograms and clusters lists
    cluster_list[[paste0("Cluster_", k)]] <- curr_cluster
    MP_list[[paste0("MP_", k)]] <- genes_MP
    k <- k + 1

    # Remove used programs, sort programs intersection again and reset 'curr_cluster' to empty vector
    progs_intersect <- progs_intersect[-match(curr_cluster, rownames(progs_intersect)), -match(curr_cluster, colnames(progs_intersect))]      # remove current chosen cluster
    sorted_intersection <- sort(apply(progs_intersect, 2, function(x) length(which(x >= min_init_intersect)) - 1), decreasing = TRUE)     # sort intersection of remaining gene programs not included in any of the previous clusters
    curr_cluster <- c()

    print(dim(progs_intersect)[2])
  }
  return(list(Cluster_list = cluster_list, MP_list = MP_list))
}


## Function for selecting robust nonnegative matrix factorization (NMF) programs
robust_nmf_programs <- function(nmf_programs, intra_min = 35, intra_max = 10, inter_filter = TRUE, inter_min = 10) {

  # Select NMF programs based on the minimum overlap with other NMF programs from the same cell line
  intra_intersect <- lapply(nmf_programs, function(z) apply(z, 2, function(x) apply(z, 2, function(y) length(intersect(x, y)))))
  intra_intersect_max <- lapply(intra_intersect, function(x) apply(x, 2, function(y) sort(y, decreasing = TRUE)[2]))
  nmf_sel <- lapply(names(nmf_programs), function(x) nmf_programs[[x]][, intra_intersect_max[[x]] >= intra_min])
  names(nmf_sel) <- names(nmf_programs)

  # Select NMF programs based on i) the maximum overlap with other NMF programs from the same cell line and
  # ii) the minimum overlap with programs from another cell line
  nmf_sel_unlist <- do.call(cbind, nmf_sel)
  inter_intersect <- apply(nmf_sel_unlist, 2, function(x) apply(nmf_sel_unlist, 2, function(y) length(intersect(x, y))))      # calculating intersection between all programs

  final_filter <- NULL
  for(i in names(nmf_sel)) {
    a <- inter_intersect[grep(i, colnames(inter_intersect), invert = TRUE), grep(i, colnames(inter_intersect))]
    b <- sort(apply(a, 2, max), decreasing = TRUE)      # for each cell line, ranks programs based on their maximum overlap with programs of other cell lines
    if(inter_filter == TRUE) {b <- b[b >= inter_min]}     # selects programs with a maximum intersection of at least a minimal number specified by 'inter_min'
    if(length(b) > 1) {
      c <- names(b[1])
      for(y in 2:length(b)) {
        if(max(inter_intersect[c, names(b[y])]) <= intra_max) {c <- c(c, names(b[y]))}     # selects programs iteratively from top-down. Only selects programs that have a intersection smaller than 10 with a previously selected programs
      }
      final_filter <- c(final_filter, c)
    } else {
      final_filter <- c(final_filter, names(b))
    }
  }
  return(final_filter)
}



# Spatial Correlation Analysis --------------------------------------------

spat_corr <- function(expr_mat,
                      sig_scores,
                      prog1,
                      prog2,
                      control,
                      n_rolls,
                      k_spots,
                      lim = c(0.05, 0.95),
                      sm_method = NULL,
                      sm_param = NULL) {

  # Score matrix for epithelial status. use this as a continuous score to sort spots.
  epi_stroma_progs <- readRDS("epi_vs_non_epi_metaprograms.rds")
  epi_stroma_sigs <- scalop::sigScores(expr_mat, epi_stroma_progs, expr.center = TRUE, conserved.genes = 0.5)
  sorted_epi_spots <- epi_stroma_sigs[, control][order(epi_stroma_sigs[, control], decreasing = FALSE)]
  sorted_epi_spots <- setNames(sorted_epi_spots, rownames(epi_stroma_sigs[order(epi_stroma_sigs[, control], decreasing = FALSE), ]))
  ordered_sigs <- sig_scores[order(match(rownames(sig_scores), names(sorted_epi_spots))), ]

  # Create correlation data frame with positions and values and adjust sliding window
  cor_res_tab <- as.data.frame(matrix(nrow = n_rolls + 1, ncol = 2, dimnames = list(c(1:(n_rolls + 1)), c("position", "correlation_coef"))))
  limits <- quantile(epi_stroma_sigs$Epithelial, probs = lim)
  window_size <- (abs(limits[[1]]) + abs(limits[[2]])) / n_rolls
  window_position <- limits[[1]]

  # Compute sliding correlation
  i <- 1
  while(window_position <= limits[[2]]) {
    cor_res_tab$position[[i]] <- window_position
    k_nearest_spots <- sorted_epi_spots[abs(sorted_epi_spots - window_position) %in% sort(abs(sorted_epi_spots - window_position), partial=1:k_spots)[1:k_spots]]
    prog1_vec <- ordered_sigs[rownames(ordered_sigs) %in% names(k_nearest_spots), prog1]
    prog2_vec <- ordered_sigs[rownames(ordered_sigs) %in% names(k_nearest_spots), prog2]
    cor_res_tab$correlation_coef[[i]] <- cor(prog1_vec, prog2_vec, method = "pearson")
    # cor_res_tab$correlation_coef[[i]] <- mean(na.omit((TTR::runCor(x = prog1_vec, y = prog2_vec, n = k_spots/2))))
    window_position <- window_position + window_size
    i <- i + 1
  }

  # Smooth results
  if(is.null(sm_method)) {return(cor_res_tab)}
  if(is.null(sm_param)) {sm_param <- 0.08 * n_rolls}
  if(sm_method == "SMA") {
    smooth_cor_res <- as.data.frame(zoo::rollmean(cor_res_tab, sm_param))
  }
  if(sm_method == "TMA") {
    smooth_cor_res <- as.data.frame(zoo::rollmean(zoo::rollmean(cor_res_tab, sm_param), sm_param))
  }
  return(smooth_cor_res)
}


plot_epi_prog_cor <- function(expr_mat,
                              sig_scores,
                              spec_prog = NULL,
                              n_rolls,
                              k_spots,
                              lim = c(0.05, 0.95),
                              sm_method = NULL,
                              sm_param = NULL) {

  epi_progs <- c("pEMT", "Senescense", "Cell_cycle", "Epithelial", "Hypoxia", "Apoptosis_LowQ")
  existing_epi_progs <- epi_progs[epi_progs %in% colnames(sig_scores)]
  paired_progs <- combn(existing_epi_progs, 2)
  intra_spat_cor <- apply(paired_progs, 2, function(x) spat_corr(expr_mat, sig_scores, x[[1]], x[[2]], "Epithelial", n_rolls = n_rolls, k_spots = k_spots, lim = lim, sm_method = sm_method, sm_param = sm_param))
  position_vec <- intra_spat_cor[[1]][, 1]
  cor_res_tab <- do.call(cbind.data.frame, intra_spat_cor)
  cor_res_tab <- cor_res_tab[, grep("correlation_coef", colnames(cor_res_tab))]
  colnames(cor_res_tab) <- str_c(paired_progs[1, ], paired_progs[2, ], sep = "-")

  plot_res <- reshape2::melt(cor_res_tab)
  plot_res$index <- position_vec
  if(is.null(spec_prog)) {
    plot <- ggplot(plot_res, aes(x = index, y = value, color = variable)) + geom_line() + geom_vline(xintercept = 0) +
      theme(axis.title = element_text()) +
      xlab("Epithelial Score") + ylab("Correlation Coefficient") + scale_color_discrete(name = "Paired\nPrograms")
  }
  if(!is.null(spec_prog)) {
    prog_pairs <- as.character(unique(plot_res$variable[grep(spec_prog, plot_res$variable)]))
    plot_res <- plot_res %>% filter(variable %in% prog_pairs)
    plot <- ggplot(plot_res, aes(x = index, y = value, color = variable)) + geom_line() + geom_vline(xintercept = 0) +
      theme(axis.title = element_text()) +
      xlab("Epithelial Score") + ylab("Correlation Coefficient") + scale_color_discrete(name = "Paired\nPrograms")
  }
  return(list(corr_table = cor_res_tab, plot = plot))
}




# Distance and Neighbor Computations --------------------------------------

### Calculate euclidean distances between pairwise observations in two groups
calc_min_dist <- function(metadata, var, dist_from, dist_to, strict = FALSE, no_info = c("zeor", "inf"), use_decon = FALSE, prop_cut = 0.2) {
  # Check input
  if(isFALSE(use_decon)) {
    stopifnot("Error: metadata must contain at least - spot barcodes under the name SpotID, tissue_position information, and the variable to claculate minimal distance from" =
                sum(c("SpotID", "array_row", "array_col") %in% colnames(metadata)) == 3 & var %in% colnames(metadata))
  }
  if(isTRUE(use_decon)) {
    stopifnot("Error: metadata must contain at least - spot barcodes under the name SpotID, and tissue_position information" =
                sum(c("SpotID", "array_row", "array_col") %in% colnames(metadata)) == 3)
  }

  if(missing(no_info)) {
    no_info <- 10^309
  } else if(no_info == "inf") {
    no_info <- 10^309
  } else {no_info <- 0}

  # Calculate pairwise distances
  pairwise_dist <- as.matrix(dist(as.matrix(metadata[, c("array_row", "array_col")]), method = "euclidean"))
  dimnames(pairwise_dist) <- list(metadata$SpotID, metadata$SpotID)

  if(isFALSE(strict) & isFALSE(use_decon)) {
    # Separate to 2 groups
    resolved_metadata <- metadata[metadata[[var]] == dist_from | metadata[[var]] == dist_to, ]
    interest_groups <- split(resolved_metadata$SpotID, resolved_metadata[[var]] == dist_from)
  }

  if(isTRUE(strict) & isFALSE(use_decon)) {
    # Calculate coherence to set as reference only spots which all neighbors share the same label
    spots_coherence <- check_coherence(metadata = metadata, var = var)
    metadata$Coherence <- spots_coherence$Coherence[match(rownames(spots_coherence), metadata$SpotID)]
    metadata$RefCent <- ifelse(metadata[[var]] == dist_from & metadata$Coherence == "Yes", paste0(dist_from, "_reference"),
                               ifelse(metadata[[var]] == dist_from & metadata$Coherence == "No", paste0(dist_from, "_disregard"), metadata[[var]]))

    # Separate to 2 refined groups
    resolved_metadata <- metadata[metadata$RefCent != paste0(dist_from, "_disregard") & (metadata[[var]] == dist_from | metadata[[var]] == dist_to), ]
    interest_groups <- split(resolved_metadata$SpotID, resolved_metadata$RefCent == paste0(dist_from, "_reference"))
  }

  if(isTRUE(use_decon)) {
    # Separate to groups based on minimum proportion of state in a spot
    resolved_metadata <- metadata[metadata[[dist_from]] > prop_cut | metadata[[dist_to]] > prop_cut, ]
    interest_groups <- split(resolved_metadata$SpotID, resolved_metadata[[dist_from]] > prop_cut)
  }

  # Measure the minimal euclidean distance between pair of spots belonging to different groups
  min_dist <- sapply(interest_groups$"FALSE", function(x) min(pairwise_dist[x, interest_groups$"TRUE"]))

  # Assign the resulting distances back to the metadata
  if(isFALSE(use_decon)) {
    metadata$Dist <- ifelse(metadata$SpotID %in% names(min_dist), yes = unname(min_dist[match(metadata$SpotID, names(min_dist))]), no = no_info)
  } else if(isTRUE(use_decon)) {
    metadata$Dist <- ifelse(metadata$SpotID %in% names(min_dist), yes = unname(min_dist[match(metadata$SpotID, names(min_dist))]),
                            ifelse(metadata$SpotID %in% interest_groups$"TRUE" & metadata[[dist_to]] > prop_cut, yes = 0, no = 10^309))   # consider changing from 0 to 1-proportion to take the state abundancy inside each spot into consideration
  }

  return(metadata)
}



### helper function to describe the 6 neighboring spots surrounding each spot in a sample
.neighbors_table <- function(metadata, var, spot_id = "SpotID") {
  if(is.factor(metadata[[var]])) {metadata[[var]] <- S4Vectors::unfactor(metadata[[var]])}
  neighbors_table <- sapply(metadata[[spot_id]], function(spot){
    spots_row = metadata$array_row[metadata[[spot_id]] == spot]
    spots_col = metadata$array_col[metadata[[spot_id]] == spot]
    c0 = as.character(metadata[metadata[[spot_id]] == spot, var])

    if(spots_col == 0 | spots_row == 0) {
      c1 = NaN
    } else {
      n1 = metadata[[spot_id]][metadata$array_row == spots_row - 1 & metadata$array_col == spots_col - 1]
      if(length(n1) == 0) {
        c1 = NaN
      } else {
        c1 = as.character(metadata[metadata[[spot_id]] == n1, var])
      }
    }

    if(spots_col == 127 | spots_row == 0) {
      c2 = NaN
    } else {
      n2 = metadata[[spot_id]][metadata$array_row == spots_row - 1 & metadata$array_col == spots_col + 1]
      if(length(n2) == 0) {
        c2 = NaN
      } else {
        c2 = as.character(metadata[metadata[[spot_id]] == n2, var])
      }
    }

    if(spots_col == 0 | spots_col == 1) {
      c3 = NaN
    } else {
      n3 = metadata[[spot_id]][metadata$array_row == spots_row & metadata$array_col == spots_col - 2]
      if(length(n3) == 0) {
        c3 = NaN
      } else {
        c3 = as.character(metadata[metadata[[spot_id]] == n3, var])
      }
    }

    if(spots_col == 126 | spots_col == 127) {
      c4 = NaN
    } else {
      n4 = metadata[[spot_id]][metadata$array_row == spots_row & metadata$array_col == spots_col + 2]
      if(length(n4) == 0) {
        c4 = NaN
      } else {
        c4 = as.character(metadata[metadata[[spot_id]] == n4, var])
      }
    }

    if(spots_col == 0 | spots_row == 77) {
      c5 = NaN
    } else {
      n5 = metadata[[spot_id]][metadata$array_row == spots_row + 1 & metadata$array_col == spots_col - 1]
      if(length(n5) == 0) {
        c5 = NaN
      } else {
        c5 = as.character(metadata[metadata[[spot_id]] == n5, var])
      }
    }

    if(spots_col == 127 | spots_row == 77) {
      c6 = NaN
    } else {
      n6 = metadata[[spot_id]][metadata$array_row == spots_row + 1 & metadata$array_col == spots_col + 1]
      if(length(n6) == 0) {
        c6 = NaN
      } else {
        c6 = as.character(metadata[metadata[[spot_id]] == n6, var])
      }
    }

    return(c(c0, c1, c2, c3, c4, c5, c6))
  })
  out <- as.data.frame(t(neighbors_table))
}


### function to check coherence - to define whom all their neighbors are of the same label as them
check_coherence <- function(metadata, var, spot_id = "SpotID", n_neighbs = 6) {
  neighbors_table <- .neighbors_table(metadata = metadata, var = var, spot_id = spot_id)
  dup_obs <- apply(neighbors_table, 1, function(x) sum(duplicated(x)))
  # neighbors_table$Coherence <- ifelse(dup_obs == (ncol(neighbors_table) - 1), "Yes", "No")
  neighbors_table$Coherence <- ifelse(dup_obs >= n_neighbs, "Yes", "No")
  return(neighbors_table)
}


### Function to filter spots that are distant from the bulk of tumor
filter_distant_spots <- function(metadata, neighbor_cut = 4, max_iter = 20) {
  iter <- 1
  repeat {
    check_neighbors <- .check_neighbor_spots_num(metadata = metadata, neighbor_cut = neighbor_cut)
    spots2filter <- rownames(check_neighbors)[check_neighbors$To_Filter == "Yes"]
    if(length(spots2filter) <= 1) {break}
    metadata <- metadata[!metadata$SpotID %in% spots2filter, ]
    if(iter > max_iter) {break}
    iter <- iter + 1
  }
  return(metadata)
}

.check_neighbor_spots_num <- function(metadata, neighbor_cut = 6) {
  neighbors_table <- sapply(metadata$SpotID, function(spot){
    spots_row = metadata$array_row[metadata$SpotID == spot]
    spots_col = metadata$array_col[metadata$SpotID == spot]
    c0 = as.numeric(metadata[metadata$SpotID == spot, "in_tissue"])

    if(spots_col == 0 | spots_row == 0) {
      c1 = 0
    } else {
      n1 = metadata$SpotID[metadata$array_row == spots_row - 1 & metadata$array_col == spots_col - 1]
      if(length(n1) == 0) {
        c1 = 0
      } else {
        c1 = as.numeric(metadata[metadata$SpotID == n1, "in_tissue"])
      }
    }

    if(spots_col == 127 | spots_row == 0) {
      c2 = 0
    } else {
      n2 = metadata$SpotID[metadata$array_row == spots_row - 1 & metadata$array_col == spots_col + 1]
      if(length(n2) == 0) {
        c2 = 0
      } else {
        c2 = as.numeric(metadata[metadata$SpotID == n2, "in_tissue"])
      }
    }

    if(spots_col == 0 | spots_col == 1) {
      c3 = 0
    } else {
      n3 = metadata$SpotID[metadata$array_row == spots_row & metadata$array_col == spots_col - 2]
      if(length(n3) == 0) {
        c3 = 0
      } else {
        c3 = as.numeric(metadata[metadata$SpotID == n3, "in_tissue"])
      }
    }

    if(spots_col == 126 | spots_col == 127) {
      c4 = 0
    } else {
      n4 = metadata$SpotID[metadata$array_row == spots_row & metadata$array_col == spots_col + 2]
      if(length(n4) == 0) {
        c4 = 0
      } else {
        c4 = as.numeric(metadata[metadata$SpotID == n4, "in_tissue"])
      }
    }

    if(spots_col == 0 | spots_row == 77) {
      c5 = 0
    } else {
      n5 = metadata$SpotID[metadata$array_row == spots_row + 1 & metadata$array_col == spots_col - 1]
      if(length(n5) == 0) {
        c5 = 0
      } else {
        c5 = as.numeric(metadata[metadata$SpotID == n5, "in_tissue"])
      }
    }

    if(spots_col == 127 | spots_row == 77) {
      c6 = 0
    } else {
      n6 = metadata$SpotID[metadata$array_row == spots_row + 1 & metadata$array_col == spots_col + 1]
      if(length(n6) == 0) {
        c6 = 0
      } else {
        c6 = as.numeric(metadata[metadata$SpotID == n6, "in_tissue"])
      }
    }

    if(spots_row == 0 | spots_row == 1) {
      c7 = 0
    } else {
      n7 = metadata$SpotID[metadata$array_row == spots_row + 2 & metadata$array_col == spots_col]
      if(length(n7) == 0) {
        c7 = 0
      } else {
        c7 = as.numeric(metadata[metadata$SpotID == n7, "in_tissue"])
      }
    }

    if(spots_row == 76 | spots_row == 77) {
      c8 = 0
    } else {
      n8 = metadata$SpotID[metadata$array_row == spots_row - 2 & metadata$array_col == spots_col]
      if(length(n8) == 0) {
        c8 = 0
      } else {
        c8 = as.numeric(metadata[metadata$SpotID == n8, "in_tissue"])
      }
    }

    return(c(c1, c2, c3, c4, c5, c6, c7, c8))
  })

  out <- as.data.frame(t(neighbors_table))
  # neighbors_sum <- rowSums(out)
  out$Sum <- rowSums(out)
  out$Border <- sapply(rownames(out), function(idx) ifelse(metadata$array_row[metadata$SpotID == idx] %in% c(0, 77) | metadata$array_col[metadata$SpotID == idx] %in% c(0, 127), "TRUE", "FALSE"))
  out$To_Filter <- ifelse(out$Border == "TRUE" & out$Sum < 3, "Yes",
                          ifelse(out$Border == "FALSE" & out$Sum < neighbor_cut, "Yes", "No"))
  return(out)
}


# Zone Classification & Enrichment ----------------------------------------


classify_zones <- function(metadata, decon_mat, only_epi_zones = TRUE, stromal_prop_cutoff = 0.7, strict = FALSE, rm_distant_spots = FALSE, zones_limits = c(3, 6), with_plots = TRUE,
                           pt_size = 3.2, zone_cols = "default") {
  # Check input
  stopifnot("Error: metadata must contain at least - spot barcodes under the name SpotID, and tissue_position information" =
              sum(c("SpotID", "array_row", "array_col") %in% colnames(metadata)) == 3)

  # Define Epithelial vs Non_Epithelial rich spots based on the cell-type proportions of the deconvolution matrix
  if("Skeletal_Muscle" %in% colnames(decon_mat) & !"Secretory_Norm" %in% colnames(decon_mat)) {
    spot_stromal_prop <- apply(decon_mat, 1, function(x) sum(x[c("T_cell", "B_cell", "Macrophage", "Fibroblast", "Endothelial", "Complement", "Skeletal_Muscle")]) / sum(x))
  } else if ("Secretory_Norm" %in% colnames(decon_mat) & !"Skeletal_Muscle" %in% colnames(decon_mat)) {
    spot_stromal_prop <- apply(decon_mat, 1, function(x) sum(x[c("T_cell", "B_cell", "Macrophage", "Fibroblast", "Endothelial", "Complement", "Secretory_Norm")]) / sum(x))
  } else if ("Skeletal_Muscle" %in% colnames(decon_mat) & "Secretory_Norm" %in% colnames(decon_mat)) {
    spot_stromal_prop <- apply(decon_mat, 1, function(x) sum(x[c("T_cell", "B_cell", "Macrophage", "Fibroblast", "Endothelial", "Complement", "Skeletal_Muscle", "Secretory_Norm")]) / sum(x))
  } else if (!"Skeletal_Muscle" %in% colnames(decon_mat) & !"Secretory_Norm" %in% colnames(decon_mat)) {
    spot_stromal_prop <- apply(decon_mat, 1, function(x) sum(x[c("T_cell", "B_cell", "Macrophage", "Fibroblast", "Endothelial", "Complement")]) / sum(x))
  }
  decon_df <- as.data.frame(decon_mat) %>%
    dplyr::mutate(EpiStroma = unname(ifelse(rowSums(decon_mat) == 0, "Filtered_out", ifelse(spot_stromal_prop >= stromal_prop_cutoff, "Non_Epithelial", "Epithelial"))))
  merged_meta <- decon_df %>% rownames_to_column(var = "SpotID") %>% left_join(metadata, ., by = "SpotID")

  if(rm_distant_spots) {
    merged_meta <- filter_distant_spots(merged_meta)
  }

  # Diagnostic plot - spatial distribution of Epithelial vs Stromal/Immune spots
  if("Filtered_out" %in% merged_meta$EpiStroma) {
    set_cols <- setNames(c("#00BFC4", "#F8766D", "#8C8C8C"), c("Epithelial", "Non_Epithelial", "Filtered_out"))
  } else {set_cols <- setNames(c("#00BFC4", "#F8766D"), c("Epithelial", "Non_Epithelial"))}
  EpiStroma_plot <- ggplot(merged_meta, aes(x = pxl_col_in_fullres, y = -pxl_row_in_fullres, color = EpiStroma)) +
    geom_point(size = pt_size) + theme_void() + scale_color_manual(values = set_cols)

  # Calculate pairwise minimal distance from Stromal spots upon which zones will be classified
  merged_meta <- calc_min_dist(merged_meta, "EpiStroma", dist_from = "Non_Epithelial", dist_to = "Epithelial", strict = strict, no_info = "inf")
  if(length(zones_limits) == 1 && zones_limits == "auto") {
    border_2_and_3 <- median(merged_meta$Dist[is.finite(merged_meta$Dist)] %>% .[. >= 2.5])
    zones_limits <- c(2.5, border_2_and_3)
  }
  btw_zones <- sapply(seq_along(zones_limits[-length(zones_limits)]), function(y) ifelse(merged_meta$EpiStroma == "Epithelial" & merged_meta$Dist >= zones_limits[y] & merged_meta$Dist <= zones_limits[y+1], paste0("Zone_", y+1), "Filtered_out"))
  btw_zones <- apply(as.data.frame(btw_zones), 1, function(y) ifelse(all(y == "Filtered_out"), "Filtered_out", y[grep("Zone", y)]))
  merged_meta$Zone <- as.factor(unname(ifelse(merged_meta$EpiStroma == "Non_Epithelial", "Non_Epithelial",
                                              ifelse(merged_meta$EpiStroma == "Epithelial" & merged_meta$Dist < zones_limits[[1]], "Zone_1",
                                                     ifelse(merged_meta$EpiStroma == "Epithelial" & merged_meta$Dist > zones_limits[[length(zones_limits)]], paste0("Zone_", length(zones_limits) + 1),
                                                            btw_zones)))))

  if(!only_epi_zones) {
    merged_meta <- calc_min_dist(merged_meta, "EpiStroma", dist_from = "Epithelial", dist_to = "Non_Epithelial", strict = strict, no_info = "inf")
    btw_zones <- sapply(seq_along(zones_limits[-length(zones_limits)]), function(y) ifelse(merged_meta$EpiStroma == "Non_Epithelial" & merged_meta$Dist >= zones_limits[y] & merged_meta$Dist <= zones_limits[y+1], paste0("NE_Zone_", y+1), "Filtered_out"))
    btw_zones <- apply(as.data.frame(btw_zones), 1, function(y) ifelse(all(y == "Filtered_out"), "Filtered_out", y[grep("NE_Zone", y)]))
    merged_meta$Zone <- as.factor(unname(ifelse(merged_meta$EpiStroma == "Epithelial", levels(merged_meta$Zone)[merged_meta$Zone],
                                                ifelse(merged_meta$EpiStroma == "Non_Epithelial" & merged_meta$Dist < zones_limits[[1]], "NE_Zone_1",
                                                       ifelse(merged_meta$EpiStroma == "Non_Epithelial" & merged_meta$Dist > zones_limits[[length(zones_limits)]], paste0("NE_Zone_", length(zones_limits) + 1),
                                                              btw_zones)))))
    merged_meta$Zone <- factor(merged_meta$Zone, levels = c("NE_Zone_3", "NE_Zone_2", "NE_Zone_1", "Zone_1", "Zone_2", "Zone_3"))
  }

  # Diagnostic plot - spatial distribution of the different zones
  if(zone_cols == "default") {
    set_cols <- setNames(c("#8C8C8C", ggsci::pal_jama("default")(4)), c("Filtered_out", "Non_Epithelial", "Zone_1", "Zone_2", "Zone_3"))
  } else {
    if("Filtered_out" %in% merged_meta$Zone) {
      set_cols <- setNames(c(gg_color_hue(sum(!levels(merged_meta$Zone) %in% "Filtered_out")), "#8C8C8C"), c(levels(merged_meta$Zone)[!levels(merged_meta$Zone) %in% "Filtered_out"], "Filtered_out"))
    } else {set_cols <- setNames(gg_color_hue(length(levels(merged_meta$Zone))), levels(merged_meta$Zone))}
  }
  Zones_plot <- ggplot(merged_meta, aes(x = pxl_col_in_fullres, y = -pxl_row_in_fullres, color = Zone)) +
    geom_point(size = pt_size) + theme_void() + scale_color_manual(values = set_cols)

  if(with_plots) {
    return(list(merged_meta = merged_meta, EpiStroma_plot = EpiStroma_plot, Zones_plot = Zones_plot))
  } else {
    return(merged_meta)
  }
}


zone_enrichment <- function(metadata, decon_mat, max_zone = NULL, min_spots_per_zone = 50, samp_spot_num = 200, samp_itter = 100, set_seed = 5410, with_plot = TRUE, set_colors = NULL) {
  # Check input
  stopifnot("Error: metadata must contain at least - spot barcodes under the name SpotID, and tissue_position information" =
              sum(c("SpotID", "array_row", "array_col", "Zone") %in% colnames(metadata)) == 4)

  # Create a data-frame with abundance of each state per-zone. this is done by summing the per-spot proportion of each state in a zone
  presence_df <- as.data.frame(decon_mat) %>% rownames_to_column(var = "SpotID") %>%
    left_join(metadata[metadata$Zone != "Filtered_out", c("SpotID", "Zone")], ., by = "SpotID") %>% na.omit %>%
    column_to_rownames(var = "SpotID") %>% mutate(Zone = droplevels(Zone), .keep = "unused")
  colnames(presence_df)[grep("Muscle", colnames(presence_df))] <- "Skeletal_Muscle"
  zonation_abund <- as.data.frame(sapply(split(presence_df[, colnames(presence_df) != "Zone"], presence_df$Zone), function(x) round(colSums(x), digits = 0)))

  if(!is.null(max_zone)) {
    zonation_abund <- zonation_abund[, 1:which(colnames(zonation_abund) == max_zone)]
  }

  # Sample X number of spots Y times for each zone
  base_df <- as.data.frame(matrix(data = 0, nrow = nrow(zonation_abund), ncol = 0, dimnames = list(rownames(zonation_abund))))
  spots_per_zone <- table(presence_df$Zone)
  while(min(spots_per_zone) <= min_spots_per_zone) {
    zonation_abund <- zonation_abund[, -which.min(spots_per_zone)]
    spots_per_zone <- spots_per_zone[-which.min(spots_per_zone)]
  }
  if(min(colSums(zonation_abund)) <= samp_spot_num) {
    samp_spot_num <- min(colSums(zonation_abund))
    warning(paste0("Ignoring 'samp_spot_num' argument and setting sampled spots to - ", min(colSums(zonation_abund)), ", which is the maximal number of spots in the smallest zone"))
  } else {samp_spot_num <- samp_spot_num}
  set.seed(set_seed)
  down_samp <- lapply(seq(1, samp_itter), function(loop) {
    convert2vec <- lapply(colnames(zonation_abund), function(y) {
      vec <- unlist(sapply(seq_len(nrow(zonation_abund)), function(x) c(rep(rownames(zonation_abund)[x], zonation_abund[x, y]))))
    })
    samp <- lapply(convert2vec, function(x) sample(x, size = samp_spot_num, replace = FALSE) %>% table %>% as.data.frame %>%
                     left_join(rownames_to_column(base_df), ., by = c("rowname" = ".")) %>% column_to_rownames("rowname"))
    na2zero_samp <- rapply(samp, function(x) ifelse(is.na(x), 0, x), how = "replace")
  })

  merged_samp_trails <- list()
  for(i in seq_len(ncol(zonation_abund))) {
    merged_samp <- list()
    for(j in seq_along(down_samp)){
      merged_samp[[j]] <- down_samp[[j]][[i]]
    }
    merged_samp_trails[[i]] <- do.call(cbind.data.frame, merged_samp)
    names(merged_samp_trails)[i] <- colnames(zonation_abund)[i]
  }

  # Take the mean of all sampled vectors for each state and transform value to its per-zone proportion
  mean_zone_abund <- lapply(merged_samp_trails, function(x) rowMeans(x))
  mean_zone_abund <- do.call(cbind.data.frame, mean_zone_abund)
  zone_prop <- as.data.frame(apply(mean_zone_abund, 2, function(x) x/sum(x)))
  if(!all(levels(presence_df$Zone) %in% colnames(zone_prop))) {
    missing_zones <- setdiff(levels(presence_df$Zone), colnames(zone_prop))
    zone_prop[, missing_zones] <- NA
    mean_zone_abund[, missing_zones] <- NA
  }

  if(!is.null(max_zone)) {
    mean_zone_abund <- mean_zone_abund[, 1:which(colnames(mean_zone_abund) == max_zone)]
  }

  if(isFALSE(with_plot)) {return(zone_prop)}
  if(isTRUE(with_plot)) {
    for_plot <- as.data.frame(zone_prop) %>% rownames_to_column(var = "State") %>%
      dplyr::mutate(States = as.factor(State), .keep = "unused", .before = 1)
    for_plot <- suppressMessages(reshape2::melt(for_plot))
    # Set colors
    default.colors <- c(scales::hue_pal()(length(levels(for_plot$States))))
    if(!is.null(set_colors) & is.null(names(set_colors))) {
      cols <- set_colors[levels(for_plot$States)]; names(cols) <- levels(for_plot$States)
    }
    if(!is.null(set_colors) & !is.null(names(set_colors))) {
      cols <- set_colors[levels(for_plot$States)]
    } else {
      cols <- default.colors; names(cols) <- levels(for_plot$States)
    }
    Zone_enrichment_plot <- ggplot(for_plot, aes(x = variable, y = value, group = States)) +
      geom_line(size = 2, aes(color = States)) + scale_color_manual(values = cols, na.translate = FALSE)
    return(list(zone_prop = zone_prop, mean_zone_abund = mean_zone_abund, enrichment_plot = Zone_enrichment_plot))
  }
}


state_zone_enrichment <- function(state, abundance_tabs, colors = NULL, subset = NULL, z_scale = TRUE,
                                  show_trend = TRUE, with_plot = TRUE, exlude_samp = "none", choose_samp = "none",
                                  trend_style = "linear", line_size = 1, trend_line_size = 3, sample_filter = FALSE, ratio = 1, panel_border_size = 1) {
  if(exlude_samp != "none") {
    abundance_tabs <- abundance_tabs[!names(abundance_tabs) %in% exlude_samp]
  }
  if(all(choose_samp != "none")) {
    abundance_tabs <- abundance_tabs[names(abundance_tabs) %in% choose_samp]
  }
  state_zonation <- lapply(abundance_tabs, function(x) x[state, ]) %>% do.call(plyr::rbind.fill, .) %>%
    dplyr::rename_all(~colnames(abundance_tabs[[1]])) %>% magrittr::set_rownames(names(abundance_tabs)) %>%
    rownames_to_column(var = "Sample") %>% dplyr::mutate(Samples = factor(Sample, levels = unique(Sample)), .keep = "unused", .before = 1)
  if(!is.null(subset)) {
    state_zonation <- state_zonation %>% filter(Samples %in% subset) %>% droplevels
  }
  if(z_scale) {
    state_zonation <- state_zonation %>% column_to_rownames(var = "Samples") %>% t() %>% scale() %>% t() %>%
      as.data.frame() %>% rownames_to_column(var = "Sample") %>% dplyr::mutate(Samples = factor(Sample, levels = unique(Sample)), .keep = "unused", .before = 1)
  }
  if(sample_filter) {
    state_zonation <- na.omit(state_zonation)
  }

  if(isFALSE(with_plot)) {return(state_zonation[, -grep("Samples", colnames(state_zonation))])}
  if(isTRUE(with_plot)) {
    for_plot <- reshape2::melt(state_zonation)
    default.colors <- c(scales::hue_pal()(length(levels(state_zonation$Samples))))
    if(!is.null(colors)) {
      if(length(colors) == length(levels(state_zonation$Samples))) {
        cols <- colors
      } else if(length(colors) != length(levels(state_zonation$Samples)) & !is.null(names(colors))) {
        cols <- colors[names(colors) %in% state_zonation$Samples]
      } else {cols <- c(colors, default.colors[1:(length(levels(state_zonation$Samples)) - length(colors))])}
    } else {cols <- default.colors}
    if(show_trend) {
      if(trend_style == "linear") {
        state_zonation_means <- sapply(state_zonation[, -grep("Samples", colnames(state_zonation))], function(x) mean(na.omit(x)))
        state_zonation_means <- as.data.frame(state_zonation_means) %>% tibble::rownames_to_column("variable") %>% add_column(Samples = NA, .before = 1) %>% dplyr::rename(value = state_zonation_means)
        suppressWarnings(print(ggplot(for_plot, aes(x = variable, y = value, group = Samples)) +
                                 geom_line(size = line_size, alpha = 0.5, aes(color = Samples)) + scale_color_manual(values = cols, na.translate = FALSE) +
                                 geom_line(data = state_zonation_means, aes(x = variable, y = value), color = "#0c2a50ff", size = trend_line_size) +
                                 theme_bw() + theme(axis.title.x = element_blank(), text = element_text(size = 14), legend.title = element_text(size = 16),
                                                    panel.border = element_rect(fill = NA, colour = "black", linewidth = panel_border_size), aspect.ratio = ratio, plot.margin = margin(0.2, , 0.2, 0.4, "cm")) +
                                 ylab("Z-Score") + guides(color = guide_legend(override.aes = list(linewidth = 5))) +
                                 scale_x_discrete(expand = c(0, 0), labels = function(x) gsub("_", " ", x, fixed = TRUE))))
      } else {
        suppressWarnings(print(ggplot(for_plot, aes(x = variable, y = value, group = Samples)) +
                                 geom_line(size = line_size, alpha = 0.5, aes(color = Samples)) +
                                 stat_smooth(aes(group = 1), size = trend_line_size, color = "#0c2a50ff") + scale_color_manual(values = cols, na.translate = FALSE) +
                                 theme_bw() + theme(axis.title.x = element_blank(), text = element_text(size = 14), legend.title = element_text(size = 16),
                                                    panel.border = element_rect(fill = NA, colour = "black", linewidth = panel_border_size), aspect.ratio = ratio, plot.margin = margin(0.2, , 0.2, 0.4, "cm")) +
                                 ylab("Z-Score") + guides(color = guide_legend(override.aes = list(linewidth = 5))) +
                                 scale_x_discrete(expand = c(0, 0), labels = function(x) gsub("_", " ", x, fixed = TRUE))))
      }
    } else if(!show_trend) {
      ggplot(for_plot, aes(x = variable, y = value, group = Samples)) +
        geom_line(size = line_size, aes(color = Samples)) + scale_color_manual(values = cols, na.translate = FALSE) +
        theme_bw() + theme(axis.title.x = element_blank(), text = element_text(size = 14), legend.title = element_text(size = 16),
                           panel.border = element_rect(fill = NA, colour = "black", linewidth = panel_border_size), aspect.ratio = ratio, plot.margin = margin(0.2, , 0.2, 0.4, "cm")) +
        ylab("Z-Score") + guides(color = guide_legend(override.aes = list(linewidth = 5))) +
        scale_x_discrete(expand = c(0, 0), labels = function(x) gsub("_", " ", x, fixed = TRUE))
    }
  }
}





## defining zone via distance from specific TME cell-type
state_relative_enrichment <- function(m, meta_prog, metadata, state, zones_limits = c(3, 6), decon_mat, min_spots_per_zone = 50, samp_spot_num = 200, samp_itter = 100, set_seed = 5410, with_plot = TRUE, set_colors = NULL) {
  spot_MP_score <- metaprog_score(m, meta_prog, score_diff = 0, bin_size = 100, bin_num = 30, conserved.genes = 0.5, center_rows = TRUE, scale = FALSE)
  metadata$MPid <- unname(droplevels(factor(spot_MP_score$meta_program, levels = c(names(meta_prog)))))
  not_state <- paste0("Non_", state)
  metadata$BinClass <- ifelse(metadata$MPid == state, state, not_state)
  metadata$State_Dist <- calc_min_dist(metadata, "BinClass", dist_from = state, dist_to = not_state, strict = FALSE, no_info = "inf")$Dist

  presence_df <- as.data.frame(decon_mat) %>% rownames_to_column(var = "SpotID") %>%
    left_join(metadata[metadata$Zone != "Filtered_out", c("SpotID", "Zone", "State_Dist")], ., by = "SpotID") %>%
    column_to_rownames(var = "SpotID") %>% mutate(Zone = droplevels(Zone), .keep = "unused")
  colnames(presence_df)[grep("Muscle", colnames(presence_df))] <- "Skeletal_Muscle"
  n_zones <- levels(merged_meta$Zone)[grep("Zone", levels(merged_meta$Zone))]

  zone_sum_ls <- list()
  for(i in seq_len(length(zones_limits) + 1)) {
    if(i == 1) {
      zone_sum_ls[[i]] <- presence_df %>% filter(Zone == n_zones[i], State_Dist < zones_limits[[i]]) %>% dplyr::select(!c(Zone, State_Dist)) %>% map_dbl(., sum)
    }
    if(i == length(zones_limits) + 1) {
      zone_sum_ls[[i]] <- presence_df %>% filter(Zone == n_zones[i], State_Dist > zones_limits[[length(zones_limits)]]) %>% dplyr::select(!c(Zone, State_Dist)) %>% map_dbl(., sum)
    }
    else if(i != 1 & i != length(zones_limits) + 1) {
      zone_sum_ls[[i]] <- presence_df %>% filter(Zone == n_zones[i], State_Dist >= zones_limits[[i-1]] & State_Dist <= zones_limits[[i]]) %>% dplyr::select(!c(Zone, State_Dist)) %>% map_dbl(., sum)
    }
  }
  names(zone_sum_ls) <- n_zones
  relative_zone_abund <- do.call(cbind.data.frame, zone_sum_ls)
  Non_Epithelial <- presence_df %>% filter(Zone == "Non_Epithelial", State_Dist < 3) %>% dplyr::select(!c(Zone, State_Dist)) %>% map_dbl(., sum)
  relative_zone_abund <- do.call(cbind.data.frame, list(as.data.frame(Non_Epithelial), relative_zone_abund)) %>% round(digits = 0)

  # Sample X number of spots Y times for each zone
  base_df <- as.data.frame(matrix(data = 0, nrow = nrow(relative_zone_abund), ncol = 0, dimnames = list(rownames(relative_zone_abund))))
  while(min(colSums(relative_zone_abund)) <= min_spots_per_zone) {
    relative_zone_abund <- relative_zone_abund[, -which.min(colSums(relative_zone_abund))]
  }
  if(min(colSums(relative_zone_abund)) <= samp_spot_num) {
    samp_spot_num <- min(colSums(relative_zone_abund))
    warning(paste0("Ignoring 'samp_spot_num' argument and setting sampled spots to - ", min(colSums(relative_zone_abund)), ", which is the maximal number of spots in the smallest zone"))
  } else {samp_spot_num <- samp_spot_num}
  set.seed(set_seed)
  down_samp <- lapply(seq(1, samp_itter), function(loop) {
    convert2vec <- lapply(colnames(relative_zone_abund), function(y) {
      vec <- unlist(sapply(seq_len(nrow(relative_zone_abund)), function(x) c(rep(rownames(relative_zone_abund)[x], relative_zone_abund[x, y]))))
    })
    samp <- lapply(convert2vec, function(x) sample(x, size = samp_spot_num, replace = FALSE) %>% table %>% as.data.frame %>%
                     left_join(rownames_to_column(base_df), ., by = c("rowname" = ".")) %>% column_to_rownames("rowname"))
    na2zero_samp <- rapply(samp, function(x) ifelse(is.na(x), 0, x), how = "replace")
  })

  merged_samp_trails <- list()
  for(i in seq_len(ncol(relative_zone_abund))) {
    merged_samp <- list()
    for(j in seq_along(down_samp)){
      merged_samp[[j]] <- down_samp[[j]][[i]]
    }
    merged_samp_trails[[i]] <- do.call(cbind.data.frame, merged_samp)
    names(merged_samp_trails)[i] <- colnames(relative_zone_abund)[i]
  }

  # Take the mean of all sampled vectors for each state and transform value to its per-zone proportion
  mean_zone_abund <- lapply(merged_samp_trails, function(x) rowMeans(x))
  mean_zone_abund <- do.call(cbind.data.frame, mean_zone_abund)
  zone_prop <- as.data.frame(apply(mean_zone_abund, 2, function(x) x/sum(x)))
  if(!all(n_zones %in% colnames(zone_prop))) {
    missing_zones <- setdiff(n_zones, colnames(zone_prop))
    zone_prop[, missing_zones] <- NA
  }

  if(isFALSE(with_plot)) {return(zone_prop)}
  if(isTRUE(with_plot)) {
    for_plot <- as.data.frame(zone_prop) %>% rownames_to_column(var = "State") %>%
      dplyr::mutate(States = as.factor(State), .keep = "unused", .before = 1)
    for_plot <- suppressMessages(reshape2::melt(for_plot))
    # Set colors
    default.colors <- c(scales::hue_pal()(length(levels(for_plot$States))))
    if(!is.null(set_colors) & is.null(names(set_colors))) {
      cols <- set_colors[levels(for_plot$States)]; names(cols) <- levels(for_plot$States)
    }
    if(!is.null(set_colors) & !is.null(names(set_colors))) {
      cols <- set_colors[levels(for_plot$States)]
    } else {
      cols <- default.colors; names(cols) <- levels(for_plot$States)
    }
    Zone_enrichment_plot <- ggplot(for_plot, aes(x = variable, y = value, group = States)) +
      geom_line(size = 2, aes(color = States)) + scale_color_manual(values = cols, na.translate = FALSE)
    return(list(zone_prop = zone_prop, enrichment_plot = Zone_enrichment_plot))
  }
}


# Get the X most positively or negatively correlating patterns of sifting enrichment of states over zones for each sample
per_sample_top_cor <- function(abundance_tabs, sample, rm_only_TME = TRUE, rm_non_abund_states = FALSE, rm_non_epi_zone = FALSE, top = 10, direction = "positive") {
  if(rm_non_abund_states) {
    samp_abund_tab <- abundance_tabs[[sample]][!rownames(abundance_tabs[[sample]]) %in% c("Skeletal_Muscle", "Secretory"), ]
  } else {samp_abund_tab <- abundance_tabs[[sample]]}
  if(rm_non_epi_zone) {
    samp_abund_tab <- samp_abund_tab[, -grep("Non_Epithelial", colnames(samp_abund_tab))]
  } else {samp_abund_tab <- samp_abund_tab}

  cor_df <- suppressMessages(reshape2::melt(suppressWarnings(cor(t(samp_abund_tab), use = "pairwise.complete.obs"))))
  if(rm_only_TME) {
    tme_states <- c("Macrophage", "B_cell", "T_cell", "Endothelial", "Fibroblast", "Complement_Fibro", "Skeletal_Muscle")
    cor_df <- cor_df[!(cor_df$Var1 %in% tme_states & cor_df$Var2 %in% tme_states), ]
  }
  direction <- ifelse(direction == "positive", TRUE,
                      ifelse(direction == "negative", FALSE, stop("direction must be set to either positive or negative")))
  cor_df <- cor_df[!cor_df$Var1 == cor_df$Var2, ] %>%
    filter(!duplicated(paste0(pmax(as.character(Var1), as.character(Var2)), pmin(as.character(Var1), as.character(Var2))))) %>%
    mutate(State_Pairs = str_c(Var1, Var2, sep = "-"), .keep = "unused", .before = 1) %>% .[order(.$value, decreasing = direction), ]
  if(top == "all") {return(cor_df)}
  else {
    out <- cor_df[1:top, ] %>% magrittr::set_rownames(1:top)
    return(out)
  }
}



# Barplot proportions of states the are spot-level neighbors of a chosen state, either average over all spots or site specific
plot_neighbor_state_abund <- function(metadata, state, site = "all") {
  samples_metadata <- readRDS(file = here("metadata/samples_metadata.rds"))
  if(site != "all") {
    stopifnot("Error: Site argument must be one of the following: 'Laryngeal' / 'Oral' / 'Oropharynx'" = site %in% unique(samples_metadata$Site))
    metadata <- metadata %>% dplyr::filter(Sample %in% samples_metadata$Sample[samples_metadata$Site == site])
  }

  inflam_spots <- nrow(metadata[metadata$MPid == state & !is.na(metadata$MPid), ])
  inflam_spots_enrich <- metadata %>% filter(MPid == state) %>% dplyr::select(Fibroblast:TNF_Signaling, -state) %>% colSums()
  scale_fact <- 1 / sum(inflam_spots_enrich / inflam_spots)
  p_df <- reshape2::melt(inflam_spots_enrich * scale_fact / inflam_spots) %>% rownames_to_column(var = "State")
  ggpubr::ggbarplot(p_df, x = "State", y = "value",
                    order = rev(levels(reorder(p_df$State, p_df$value, FUN = max))), fill = "lightskyblue") +
    ylab("") + xlab("") + ggtitle(paste0("Spot-Level ", state, " Neighboring states proportions"), subtitle = paste0(site, " spots"))
}


# Neighboring state abundance are calculated as spot-level (= spot composition) enrichment, by using the deconvolution matrix
plot_neighbor_state_sum <- function(metadata, state, frac_abund = TRUE) {
  samples_metadata <- readRDS(file = here("metadata/samples_metadata.rds"))
  all_states <- names(readRDS(file = here("results/Generated_Data/Final_Metaprograms_Extended.rds")))
  prop_tab <- lapply(c("All", unique(samples_metadata$Site)), function(site) {
    if(site == "All") {
      metadata <- metadata
    } else {
      metadata <- metadata %>% dplyr::filter(Sample %in% samples_metadata$Sample[samples_metadata$Site == site])
    }
    if(frac_abund) {
      calc_jacc <- sapply(all_states[!all_states %in% state], function(x) {
        neighbs <- sum(metadata[metadata[[state]] != 0 & metadata[[x]] != 0, x])
        uni <- sum(metadata[metadata[[state]] != 0, state]) + sum(metadata[metadata[[x]] != 0, x]) - neighbs
        jac_idx <- neighbs / uni
      })
    } else {
      calc_jacc <- sapply(all_states[!all_states %in% state], function(x) {
        neighbs <- nrow(metadata[metadata[[state]] != 0 & metadata[[x]] != 0, ])
        uni <- nrow(metadata[metadata[[state]] != 0, ]) + nrow(metadata[metadata[[x]] != 0, ]) - neighbs
        jac_idx <- neighbs / uni
      })
    }
  })
  prop_df <- do.call(cbind.data.frame, prop_tab) %>% purrr::set_names("All", "Laryngeal", "Oral", "Oropharynx") %>% tibble::rownames_to_column(var = "State") %>%
    reshape2::melt() %>% mutate(State = as.factor(State))
  prop_df$State <- forcats::fct_reorder(prop_df$State, prop_df$value, max)

  # Plot
  sub_title <- ifelse(frac_abund, "Using per-spot state fraction summation", "Using the number of spots where both states represented")
  ggplot(prop_df, aes(x = State, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.7, color = "black") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1, size = 11), axis.title.x = element_blank(), axis.title.y = element_blank(), ) +
    ggtitle(paste0("Spot-Level ", state, " Neighboring states proportions"), subtitle = sub_title) + scale_y_continuous(expand = c(0, 0), limits = c(0, max(prop_df$value) + 0.01)) +
    guides(fill = guide_legend(title = "Site"))
}



# Presence (vs absence) of neighboring state is calculated between spots (neighboring spots).
# The metadata supplied to the function must contain spatial coordinate information!
neighbor_props <- function(metadata, sample, spot_class = "MPid", spot_name = "Key", plot_states = FALSE) {

  # Subset metadata to include only sample of interest
  metadata <- metadata[metadata$Sample == sample, ]
  all_states <- names(readRDS(file = here("results/Generated_Data/Final_Metaprograms_Extended.rds")))
  existing_states <- all_states[all_states %in% unique(metadata[[spot_class]])]

  # Compute neighboring spots table
  neighbs_tab <- check_coherence(metadata, var = spot_class, spot_id = spot_name)
  neighbs_tab <- neighbs_tab[, -c(1, ncol(neighbs_tab))]
  neighbs_df <- as.data.frame(matrix(data = NA, nrow = length(existing_states), ncol = length(existing_states), dimnames = list(existing_states, existing_states)))
  for(i in colnames(neighbs_df)) {
    state_spots_ids <- dplyr::pull(metadata[metadata[[spot_class]] == i & !is.na(metadata[[spot_class]]), "Key"])
    neighbs_df[[i]] <- sapply(all_states[all_states %in% unique(metadata[[spot_class]])], function(x) sum(apply(neighbs_tab[state_spots_ids, ], 1, function(y) ifelse(x %in% y, 1, 0))))
  }

  # Calculate proportion of neighboring spots
  state_pairs <- t(combn(colnames(neighbs_df), 2))
  states_neighb_prop <- c()
  for(i in 1:nrow(state_pairs)) {
    states_neighb_prop[[i]] <- sum(neighbs_df[state_pairs[i, 1], state_pairs[i, 2]], neighbs_df[state_pairs[i, 2], state_pairs[i, 1]]) / sum(table(metadata[[spot_class]])[c(state_pairs[i, 1], state_pairs[i, 2])])
    names(states_neighb_prop)[i] <- stringr::str_c(state_pairs[i, ], collapse = ".")
  }
  states_neighb_prop <- as.data.frame(t(rbind.data.frame(states_neighb_prop))) %>% tibble::rownames_to_column(var = "State_Pairs") %>% dplyr::rename(Proportion = V1) %>% arrange(desc(Proportion))

  if(plot_states) {
    state_cols <- read.delim(file = here("aux_data/state_col_tab_extend.tsv"), sep = "\t", header = TRUE)
    state_cols <- setNames(state_cols$V2, state_cols$V1)
    p <- ggplot(metadata, aes(x = pxl_in_cols, y = -pxl_in_rows, color = .data[[spot_class]])) +
      geom_point(size = 4) + scale_color_manual(name = "Metaprograms", values = state_cols) +
      theme_classic() + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.line = element_blank())
    return(list(neighbs_df = neighbs_df, neighbs_prop = states_neighb_prop, plot = p))
  } else {
    return(list(neighbs_df = neighbs_df, neighbs_prop = states_neighb_prop))
  }
}


state_neighbor_props <- function(metadata, state, site = "All", plot_res = FALSE, set_colors = NULL) {
  samples_metadata <- readRDS(file = here("metadata/samples_metadata.rds"))
  all_states <- names(readRDS(file = here("results/Generated_Data/Final_Metaprograms_Extended.rds")))

  if(length(state) == 1) {
    if(site == "All") {
      state_props <- as.data.frame(matrix(data = NA, nrow = length(all_states) - 1, ncol = length(samples_metadata$Sample),
                                          dimnames = list(all_states[!all_states %in% state], samples_metadata$Sample)))
    } else {
      stopifnot("Error: Site argument must be one of the following: 'Laryngeal' / 'Oral' / 'Oropharynx'" = site %in% unique(samples_metadata$Site))
      state_props <- as.data.frame(matrix(data = NA, nrow = length(all_states) - 1, ncol = length(samples_metadata$Sample[samples_metadata$Site == site]),
                                          dimnames = list(all_states[!all_states %in% state], samples_metadata$Sample[samples_metadata$Site == site])))
    }
    existing_samps <- colnames(state_props)
  }
  if(length(state) > 1) {
    states_props_ls <- lapply(state, function(selec_states) {
      if(site == "All") {
        state_props <- as.data.frame(matrix(data = NA, nrow = length(all_states) - 1, ncol = length(samples_metadata$Sample),
                                            dimnames = list(all_states[!all_states %in% selec_states], samples_metadata$Sample)))
      } else {
        stopifnot("Error: Site argument must be one of the following: 'Laryngeal' / 'Oral' / 'Oropharynx'" = site %in% unique(samples_metadata$Site))
        state_props <- as.data.frame(matrix(data = NA, nrow = length(all_states) - 1, ncol = length(samples_metadata$Sample[samples_metadata$Site == site]),
                                            dimnames = list(all_states[!all_states %in% selec_states], samples_metadata$Sample[samples_metadata$Site == site])))
      }
    })
    names(states_props_ls) <- state
    existing_samps <- colnames(states_props_ls[[1]])
  }

  for(samp in existing_samps) {
    neighbs <- neighbor_props(metadata, sample = samp, plot_states = FALSE)
    # neighbs <- neighbor_props(extend_metadata, sample = samp, plot_states = FALSE)
    if(length(state) == 1) {
      state_neighbs <- neighbs$neighbs_prop[grep(state, neighbs$neighbs_prop$State_Pairs), ]
      for(selec_state in all_states[!all_states %in% state]) {
        if(S4Vectors::isEmpty(state_neighbs[grep(selec_state, state_neighbs$State_Pairs), "Proportion"])) {
          state_props[selec_state, samp] <- 0
        } else {
          state_props[selec_state, samp] <- state_neighbs[grep(selec_state, state_neighbs$State_Pairs), "Proportion"]
        }
      }
    } else if(length(state) > 1) {
      for(selec_states in state) {
        state_neighbs <- neighbs$neighbs_prop[grep(selec_states, neighbs$neighbs_prop$State_Pairs), ]
        for(selec_state in all_states[!all_states %in% selec_states]) {
          if(S4Vectors::isEmpty(state_neighbs[grep(selec_state, state_neighbs$State_Pairs), "Proportion"])) {
            states_props_ls[[selec_states]][selec_state, samp] <- 0
          } else {
            states_props_ls[[selec_states]][selec_state, samp] <- state_neighbs[grep(selec_state, state_neighbs$State_Pairs), "Proportion"]
          }
        }
      }
    }
  }

  if(length(state) == 1) {
    state_props <- state_props %>% mutate(Mean = rowMeans(.), SE = matrixStats::rowSds(as.matrix(.)) / sqrt(ncol(.)))
    sorted_mean_vals <- state_props[order(state_props$Mean, decreasing = TRUE), "Mean"]
  } else if(length(state) > 1) {
    state_props <- lapply(states_props_ls, function(x) x %>% mutate(Mean = rowMeans(.), SE = matrixStats::rowSds(as.matrix(.)) / sqrt(ncol(.))))
    sorted_mean_vals <- lapply(state_props, function(x) x[order(x$Mean, decreasing = TRUE), "Mean"])
  }

  if(plot_res) {
    # Set colors
    default.colors <- c(scales::hue_pal()(length(all_states)))
    if(!is.null(set_colors) & is.null(names(set_colors))) {
      cols <- set_colors[all_states]; names(cols) <- all_states
    }
    if(!is.null(set_colors) & !is.null(names(set_colors))) {
      cols <- set_colors[all_states]
    } else {
      cols <- default.colors; names(cols) <- all_states
    }

    if(length(state) == 1) {
      plot_df <- state_props %>% dplyr::select(Mean, SE) %>% tibble::rownames_to_column(var = "State") %>% tibble::add_column(Ref_state = state)
    } else if(length(state) > 1) {
      plot_df <- lapply(names(state_props), function(selec_state) state_props[[selec_state]] %>% dplyr::select(Mean, SE) %>% tibble::rownames_to_column(var = "State") %>% tibble::add_column(Ref_state = selec_state))
      plot_df <- do.call(rbind.data.frame, plot_df)
    }

    # Plot
    p <- ggplot(data = plot_df, aes(x = Ref_state, y = Mean, fill = State)) +
      geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.7, color = "black") +
      geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2, position = position_dodge(0.7)) +
      scale_fill_manual(values = cols) + theme_bw() + theme(axis.text.x = element_text(size = 11), axis.title.x = element_blank()) + ylab("Mean Neighbor\nState Proportion") +
      scale_y_continuous(expand = c(0, 0), limits = c(0, max(plot_df$Mean + plot_df$SE) + 0.01))
    return(list(Sample_Neighbors_Props = state_props, Mean_Neighbor_State_Prop = sorted_mean_vals, plot = p))
  } else {
    return(list(Sample_Neighbors_Props = state_props, Mean_Neighbor_State_Prop = sorted_mean_vals))
  }
}








# Deconvolution functions -------------------------------------------------

calc_interactions <- function(x,
                              which = c("norm", "hypergeo"),
                              min_prop = 0) {
  # check validity of input arguments
  which <- match.arg(which)
  stopifnot(
    is.matrix(x), is.numeric(x),
    all(dim(x) > 0), ncol(x) > 1,
    is.numeric(min_prop), length(min_prop) == 1)

  # get interactions table
  if (is.null(colnames(x))) {
    colnames(x) <- seq_len(ncol(x))
  }
  df <- .count_interactions(x, min_prop)

  switch(which,
         norm = .compute_stats(x, df),
         hypergeo = .compute_hypergeometric(x, df))
}

.count_interactions <- function(x, min_prop) {
  # for each pair of groups count how many
  # samples have value above 'min_prop'
  x <- x > min_prop
  ij <- combn(colnames(x), 2)
  y <- apply(ij, 2, function(.) sum(matrixStats::rowAlls(x[, ., drop = FALSE])))

  # construct 'data.frame'
  df <- data.frame(t(ij), y)
  names(df) <- c("from", "to", "n")
  return(df)
}

.compute_stats <- function(x, df) {
  # ensure proper ordering
  y <- colnames(x)
  df$i <- factor(df$from, y)
  df$j <- factor(df$to, rev(y))

  # compute proportion of samples that have all groups
  t <- colSums(x > 0)
  i <- match(df$from, y)
  j <- match(df$to, y)
  df$ti <- t[i]
  df$tj <- t[j]
  df$pi <- df$n / df$ti
  df$pj <- df$n / df$tj
  df$ij_mean <- (df$pi + df$pj) / 2
  df$fi <- df$ti/nrow(x)
  df$fj<-df$tj/nrow(x)
  df$pairs<-paste(df$i,df$j)
  return(df)
}

.compute_hypergeometric <- function(x, df) {
  # ensure proper ordering
  y <- colnames(x)
  df$i <- factor(df$from, y)
  df$j <- factor(df$to, rev(y))

  # compute proportion of samples that have all groups
  t <- colSums(x > 0)
  i <- match(df$from, y)
  j <- match(df$to, y)
  df$ti <- t[i]
  df$tj <- t[j]
  df$pi <- df$n / df$ti
  df$pj <- df$n / df$tj
  df$ij_mean <- (df$pi+df$pj)/2
  df$q.hyper <- df$n
  df$m.hyper <- df$ti
  df$n.hyper <- nrow(x)
  df$k.hyper <- df$tj
  df$pairs <- paste(df$i,df$j)
  return(df)
}

coloc_merge <- function(x, y){
  df <- merge(x, y, by = "row.names", all = TRUE)
  rownames(df) <- df$Row.names
  df$Row.names <- NULL
  return(df)
}




# Spatial Data Plotting Functions -----------------------------------------


geom_spatial <-  function(mapping = NULL,
                          data = NULL,
                          image_obj = image_obj,
                          stat = "identity",
                          position = "identity",
                          na.rm = FALSE,
                          show.legend = NA,
                          inherit.aes = TRUE,
                          image.alpha = image.alpha,
                          crop = crop,
                          ...) {

  GeomCustom <- ggproto("GeomCustom",
                        Geom,
                        required_aes = c("x", "y"),
                        extra_params = c("na.rm", "image_obj", "image.alpha", "crop"),
                        default_aes = aes(
                          shape = 21,
                          colour = "black",
                          point.size.factor = 1.0,
                          fill = NA,
                          alpha = NA,
                          stroke = 0.25
                        ),

                        setup_data = function(self, data, params) {
                          data <- ggproto_parent(Geom, self)$setup_data(data, params)
                          # We need to flip the image as the Y coordinates are reversed
                          data$y = max(data$y) - data$y + min(data$y)
                          data
                        },

                        draw_key = draw_key_point,
                        draw_group = function(data, panel_scales, coord, image_obj, image.alpha, crop) {
                          if (!crop) {
                            y.transform <- c(0, nrow(image_obj$Image)) - panel_scales$y.range
                            data$y <- data$y + sum(y.transform)
                            panel_scales$x$continuous_range <- c(0, ncol(image_obj$Image))
                            panel_scales$y$continuous_range <- c(0, nrow(image_obj$Image))
                            panel_scales$y.range <- c(0, nrow(image_obj$Image))
                            panel_scales$x.range <- c(0, ncol(image_obj$Image))
                          }

                          z <- coord$transform(
                            data.frame(x = c(0, ncol(image_obj$Image)), y = c(0, nrow(image_obj$Image))),
                            panel_scales)
                          z$y <- -rev(z$y) + 1
                          wdth <- z$x[2] - z$x[1]
                          hgth <- z$y[2] - z$y[1]
                          vp <- grid::viewport(x = unit(z$x[1], units = "npc"),
                                               y = unit(z$y[1], units = "npc"),
                                               width = unit(wdth, units = "npc"),
                                               height = unit(hgth, units = "npc"),
                                               just = c("left", "bottom"))

                          i_grob <- grid::rasterGrob(image_obj$Image, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"))
                          img_grob <- grid::editGrob(i_grob, vp = vp)
                          spot.size <- image_obj$Spot_Radius
                          coords <- coord$transform(data, panel_scales)
                          pts <- grid::pointsGrob(x = coords$x,
                                                  y = coords$y,
                                                  pch = data$shape,
                                                  size = unit(spot.size, "npc") * data$point.size.factor,
                                                  gp = gpar(col = alpha(colour = coords$colour, alpha = coords$alpha),
                                                            fill = alpha(colour = coords$fill, alpha = coords$alpha),
                                                            lwd = coords$stroke)
                          )

                          vp <- grid::viewport()
                          gt <- grid::gTree(vp = vp)
                          if (image.alpha > 0) {
                            if (image.alpha != 1) {
                              img_grob$raster = as.raster(matrix(
                                data = alpha(colour = img_grob$raster, alpha = image.alpha),
                                nrow = nrow(img_grob$raster),
                                ncol = ncol(img_grob$raster),
                                byrow = TRUE)
                              )
                            }
                            gt <- grid::addGrob(gTree = gt, child = img_grob)
                          }
                          gt <- grid::addGrob(gTree = gt, child = pts)
                          # Replacement for ggname
                          gt$name <- grobName(grob = gt, prefix = 'geom_spatial')
                          return(gt)
                          # ggplot2:::ggname("geom_spatial", gt)
                        }
  )

  layer(geom = GeomCustom,
        mapping = mapping,
        data = data,
        stat = stat,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(na.rm = na.rm, image_obj = image_obj, image.alpha = image.alpha, crop = crop, ...)
  )
}


.plot_image <- function(x) {
  x <- grid::rasterGrob(x,
                        interpolate = FALSE,
                        width = grid::unit(1, "npc"),
                        height = grid::unit(1, "npc"))

  ggplot() +
    annotation_custom(
      grob = x,
      xmin = 0,
      xmax = ncol(x$raster),
      ymin = 0,
      ymax = nrow(x$raster)) +
    coord_fixed(
      xlim = c(0, ncol(x$raster)),
      ylim = c(0, nrow(x$raster))) +
    theme_void()
}


.get_coords <- function(image_obj, scale = "lowres") {

  if (!is.null(scale)) {
    coords <- image_obj$Coordinates[, c("pxl_in_rows", "pxl_in_cols")]
    rownames(coords) <- image_obj$Coordinates$SpotID
    scale <- match.arg(arg = scale, choices = c("spot", "hires", "fiducial", "lowres"))
    scale_use <- image_obj$Scale[[scale]]
    coords <- coords * scale_use
  } else {
    coords <- image_obj$Coordinates[, c("pxl_in_rows", "pxl_in_cols")]
  }
  return(coords)
}


.plot_hm <- function(dat, program_score, class, size = 2, direction = 1) {
  if(class == "numeric") {
    G <- ggplot(dat, aes(x = pxl_in_cols, y = -pxl_in_rows)) +
      geom_point(aes(col = .data[[program_score]]), size = size) + scale_color_viridis_c(option = "A", na.value = adjustcolor("#0c2a50ff", alpha.f = 0.9), direction = direction) +
      labs(title = element_blank()) +
      theme(legend.position = "right",
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_blank())
  } else {
    G <- ggplot(dat, aes(x = pxl_in_cols, y = -pxl_in_rows)) +
      geom_point(aes(col = .data[[program_score]]), size = size)  +
      labs(title = element_blank()) +
      theme(legend.position = "right",
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_blank())
  }
  return(G)
}

.plot_hm2 <- function(dat, program_score, class, size = 2, stroke = 0) {
  if(class == "numeric") {
    G <- ggplot(dat, aes(x = pxl_in_cols, y = -pxl_in_rows)) +
      geom_point(aes(fill = .data[[program_score]]), shape = 21, size = size, stroke = stroke, color = "black") + scale_color_viridis_c(option = "A", na.value = adjustcolor("#0c2a50ff", alpha.f = 0.9)) +
      labs(title = element_blank()) +
      theme(legend.position = "right",
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_blank())
  } else {
    G <- ggplot(dat, aes(x = pxl_in_cols, y = -pxl_in_rows)) +
      geom_point(aes(fill = .data[[program_score]]), shape = 21, size = size, stroke = stroke, color = "black")  +
      labs(title = element_blank()) +
      theme(legend.position = "right",
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_blank())
  }
  return(G)
}

# Plot Spatial Heatmap of Choosen Program:
spatial_heatmap <- function(sig_scores, image_obj, include_prog = "all", normalize = FALSE, size = 2, direction = 1) {
  # Check input
  if(!is.list(sig_scores)) {stop("Input should be list containing named vectors of signature score. Each program as an entry. Spot IDs as names, for each numeric score")}
  if(all(lengths(sig_scores) != length(sig_scores[[1]]))) {stop("Not all list's entries are of same length!")}
  if(include_prog != "all") {
    if(!all(include_prog %in% names(sig_scores))) {stop("The parameter 'include_prog' must contain name/s of program/s that match the name/s of the 'sig_score' list entries")}
    sig_scores <- sig_scores[include_prog]
  }

  # Construct a dataframe to plot
  if(normalize) {
    spot_scores <- do.call(cbind.data.frame, sig_scores) %>% lapply(., function(sig) minMax_normalize(sig)) %>%
      as.data.frame(.) %>% add_column(SpotID = names(sig_scores[[1]]))
  } else {
    # spot_scores <- do.call(cbind.data.frame, sig_scores) %>% rownames_to_column(var = "SpotID")
    spot_scores <- do.call(cbind.data.frame, sig_scores)
    if(any("SpotID" %in% colnames(spot_scores))) {
      spot_scores <- spot_scores
    } else if(all(grep("[1-9]", rownames(spot_scores)))) {
      rownames(spot_scores) <- names(sig_scores[sapply(sig_scores, function(x) !is.null(names(x)))][[1]])
      spot_scores <- spot_scores %>% rownames_to_column(var = "SpotID")
    } else {spot_scores <- spot_scores %>% rownames_to_column(var = "SpotID")}
  }
  spots_filter <- image_obj$Coordinates[image_obj$Coordinates$SpotID %in% spot_scores$SpotID, ]
  spots_merge <- merge(spots_filter, spot_scores, by = "SpotID")

  # Plot program scores on the tissue spot map
  plots <- map(names(sig_scores)[!names(sig_scores) %in% "SpotID"], ~.plot_hm(dat = spots_merge, program_score = ., class = class(spots_merge[[.]]), size = size, direction = direction))
  names(plots) <- names(sig_scores)[!names(sig_scores) %in% "SpotID"]
  return(plots)
}


# Plot Spatial Scatter-Pie
spatial_scatterpie <- function(image_obj, decon_mtrx, cell_types = colnames(decon_mtrx), img = FALSE, scatterpie_alpha = 1, pie_scale = 0.4, degrees = NULL, axis = NULL, ...) {

  # Check input
  stopifnot(is.matrix(decon_mtrx) | is.data.frame(decon_mtrx), is.character(cell_types) & length(cell_types) <= ncol(decon_mtrx),
            is.numeric(scatterpie_alpha), is.numeric(pie_scale)
  )

  # If image is TRUE extract image from the image_object
  if (isFALSE(img)) {
    p <- ggplot() + coord_fixed()
    ymax <- max(dim(image_obj$Image))
  } else {
    img <- image_obj$Image
    p <- .plot_image(img)
    ymax <- max(p$coordinates$limits$y)
  }

  # Extract coordinate matrix from image_object
  x <- .get_coords(image_obj)
  # Subset rows to match spots present in the deconvolution dataframe
  x <- as.matrix(x[rownames(x) %in% rownames(decon_mtrx), ])
  colnames(x) <- c("coord_y", "coord_x")

  # Create dataframe for plotting
  # df <- merge(x, decon_mtrx, by = "row.names", all = TRUE)
  df <- merge(x, decon_mtrx, by = "row.names", all = FALSE)
  df$coord_y_i <- abs(df$coord_y - ymax)

  # Plot
  p + scatterpie::geom_scatterpie(data = df, aes_string(x = "coord_x", y = "coord_y_i"), cols = cell_types, color = NA, alpha = scatterpie_alpha, pie_scale = pie_scale, ...) +
    theme_void() +
    theme(legend.key.size = unit(0.5, "lines"))
}



# Plot expression state colocalization
plot_colocalization <- function(score_mat, cor_method = "pearson", cor_use = "complete.obs", hc_method = "ward.D", add_hc_plot = FALSE, plot_title = "",
                                palette = "bright", leg_name = "Pearson\nCorrelation", triang_plot = FALSE) {
  state_cor <- cor(score_mat, method = cor_method, use = cor_use)
  state_hc <- hclust(as.dist(1 - state_cor), method = hc_method)
  state_hc <- reorder(as.dendrogram(state_hc), colMeans(state_cor))
  state_cor <- state_cor[order.dendrogram(state_hc), order.dendrogram(state_hc)]
  state_cor_melt <- reshape2::melt(state_cor)

  if(palette == "bright") {
    cols = c("#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f")
    grad_steps <- 2 / (length(cols) - 1)
    col_vals <- seq(from = -1, to = 1, by = grad_steps)
    scale_fill <- ggplot2::scale_fill_gradientn(colors = cols, values = scales::rescale(col_vals), limits = c(-1, 1), oob = scales::squish, name = leg_name, na.value = "#000000")
  } else {
    scale_fill <- scale_fill_gradient2(limits = c(-1, 1), low = "dodgerblue4", mid = "white", high = "red4", midpoint = 0, na.value = "#000000", oob = squish, name = leg_name)
  }

  if(triang_plot) {
    hm_cols = c("#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f")
    coloc_plot <- corrplot::corrplot(state_cor, method = "color", type = "upper", col = hm_cols, tl.col = "black", tl.srt = 45, tl.cex = 1.2, order = "hclust")
  } else {
    coloc_plot <- ggplot(data = state_cor_melt, aes(x = Var1, y = Var2, fill = value)) + geom_raster() + labs(title = plot_title) +
      # scale_fill_gradient2(limits = c(-1, 1), low = "dodgerblue4", mid = "antiquewhite", high = "red4", midpoint = 0 , oob = squish, name = "") +
      eval(scale_fill) +
      theme(axis.ticks = element_blank(), panel.background = element_rect(fill = "white"),  axis.line = element_blank(), legend.text.align = 0.5,
            axis.text.x = element_text(angle = 45, hjust=1), axis.title = element_blank())
  }
  if(add_hc_plot) {
    hc_plot <- plot(state_hc)
    return(list(Colocalization = coloc_plot, Dendrogram = hc_plot))
  } else {
    return(coloc_plot)
  }
}



# Plot spatial features (metaprogram-IDs, tissue zones, etc')
# Numeric (such as, expression values) can also be ploted as a contineous variable bu must be first added to the metadata
plot_spatial_features <- function(metadata,
                                  image_obj,
                                  color_by,
                                  barcodes = "SpotID",
                                  cols = NULL,
                                  magma = FALSE,
                                  magma_rev = FALSE,
                                  limits = NULL,
                                  midpoint = NULL,
                                  image.alpha = 0.6,
                                  crop = TRUE,
                                  scale = "lowres",
                                  pt.size = 1.6,
                                  pt.alpha = 1,
                                  stroke = 0.25,
                                  na.value = 'transparent',
                                  rm_spots = metadata[[barcodes]][is.na(metadata[[color_by]])]) {

  # Spatial image object most be filtered to contain only existing spots
  if(any(image_obj$Coordinates$in_tissue == 0)) {
    stop("The image object provided most be filtered! use - load_spatial_image(path_to_image_dir, filter = TRUE)")
  }

  data <- metadata[, c(barcodes, color_by)]
  if(is.character(data[[color_by]])) {
    data[, color_by] <- factor(x = data[[color_by]])
  }

  coordinates <- .get_coords(image_obj, scale = scale) %>% tibble::rownames_to_column(var = barcodes)
  data <- merge(coordinates, data, by = barcodes, all = TRUE) %>% tibble::column_to_rownames(var = barcodes)

  if(is.numeric(data[[color_by]])) {
    plot <- ggplot(data = data, aes_string(x = colnames(data)[2], y = colnames(data)[1], fill = color_by, alpha = color_by))
  } else {
    plot <- ggplot(data = data, aes_string(x = colnames(data)[2], y = colnames(data)[1], fill = color_by))
  }
  plot <- plot + geom_spatial(data = data,
                              point.size.factor = ifelse(rownames(data) %in% rm_spots, yes = 0, no = pt.size),
                              image_obj = image_obj,
                              image.alpha = image.alpha,
                              crop = crop,
                              stroke = ifelse(rownames(data) %in% rm_spots, yes = 0, no = stroke))
  plot <- plot + geom_spatial(data = data,
                              point.size.factor = ifelse(rownames(data) %in% rm_spots, yes = 0, no = pt.size),
                              image_obj = image_obj,
                              image.alpha = 0,
                              crop = crop,
                              stroke = ifelse(rownames(data) %in% rm_spots, yes = 0, no = stroke),
                              alpha = pt.alpha)

  if(!is.numeric(data[[color_by]])) {
    if(is.null(cols)) {
      cols <- gg_color_hue(length(levels(data[, color_by])))
    } else if (!is.null(cols)) {
      if(length(cols) < length(levels(data[, color_by]))) {
        warning("Not enough colors provided! seting default colors...")
        cols <- gg_color_hue(length(levels(data[, color_by])))
      } else {
        cols <- cols[1:length(levels(data[, color_by]))]
      }
    }
    scale <- scale_fill_manual(name = ,values = cols, na.value = na.value, na.translate = FALSE)
    plot <- plot + scale
  }

  if(is.numeric(data[[color_by]])) {
    if(purrr::is_null(limits)) {
      limits <- range(data[[color_by]], na.rm = TRUE)
    }
    if(purrr::is_null(midpoint)) {
      midpoint <- round(mean(limits), digits = 1)
    }
    if(magma) {
      if(magma_rev) {
        # c_cols <- rev(c(colorRampPalette(c("white", rev(viridis::magma(323, begin = 0.15))[1]))(10), rev(viridis::magma(323, begin = 0.18))))
        c_cols <- rev(viridis::magma(n = 6))
      } else {
        # c_cols <- c(colorRampPalette(c("white", rev(viridis::magma(323, begin = 0.15))[1]))(10), rev(viridis::magma(323, begin = 0.18)))
        c_cols <- viridis::magma(n = 6)
      }
    } else if (isTRUE(min(na.omit(data[[color_by]])) < 0)) {
      c_cols <- grDevices::colorRampPalette(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100)
    } else if(isTRUE(min(na.omit(data[[color_by]])) >= 0)) {
      c_cols <- grDevices::colorRampPalette(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")[1:6]))(50)
    }

    if(length(midpoint) == 1 && midpoint == 0) {
      grad_steps <- (max(limits) - min(limits)) / (length(c_cols) - 1)
      col_vals <- seq(from = min(limits), to = max(limits), by = grad_steps)
    } else if (length(midpoint) == 1 && midpoint != 0) {
      zero_col_val <- round(length(c_cols) / 2)
      low_grad_steps <- (midpoint - min(limits)) / (zero_col_val - 1)
      high_grad_steps <- (max(limits) - midpoint) / (zero_col_val - 1)
      col_vals <- c(seq(from = min(limits), to = midpoint, by = low_grad_steps), seq(from = midpoint, to = max(limits), by = high_grad_steps)[-1])
    } else if (length(midpoint) > 1) {
      zero_col_val <- round(length(c_cols) / 2)
      c_cols <- c(c_cols[1:(zero_col_val - 1)], rep(c_cols[zero_col_val], 2), c_cols[(zero_col_val + 1):length(c_cols)])
      low_grad_steps <- (min(midpoint) - min(limits)) / (zero_col_val - 1)
      high_grad_steps <- (max(limits) - max(midpoint)) / (zero_col_val - 1)
      col_vals <- c(seq(from = min(limits), to = min(midpoint), by = low_grad_steps),
                    seq(from = max(midpoint), to = max(limits), by = high_grad_steps))
    }

    scale <- scale_fill_gradientn(colors = c_cols, values = scales::rescale(col_vals), limits = limits, oob = scales::squish, na.value = na.value)
    plot <- plot + eval(scale) + scale_alpha_continuous(range = c(0, 0.1)) + guides(alpha = "none")
  }

  plot <- plot + coord_fixed() + theme(aspect.ratio = 1) + NoAxes() + theme(panel.background = element_blank())
  return(plot)
}


# Plot croped (or uncroped) spatial image alone
plot_scaled_image <- function(image_obj, crop = TRUE, scale = "lowres") {
  img_coords <- .get_coords(image_obj, scale = scale)
  plot <- ggplot(img_coords, aes_string(x = colnames(img_coords)[2], y = colnames(img_coords)[1]))
  plot <- plot + geom_spatial(image_obj = image_obj,
                              image.alpha = 1,
                              point.size.factor = 0,
                              stroke = 0,
                              crop = crop)
  plot <- plot + coord_fixed() + theme(aspect.ratio = 1) + NoAxes() + theme(panel.background = element_blank())
  return(plot)
}





plot_spatial_feats <- function(
  metadata,
  image_obj,
  color_by,
  barcodes = "SpotID",
  cols = NULL,
  image.alpha = 0.6,
  crop = TRUE,
  pt.size = 1.6,
  pt.alpha = 1,
  stroke = 0.25,
  legend_name = "Vars",
  # na.value = 'grey50'
  # na.value = scales::alpha("grey50", alpha = 0.1)
  # na.value = adjustcolor("transparent", alpha.f = 0.1)
  na.value = "#FFFFFF00",
  rm_spots = metadata[[barcodes]][is.na(metadata[[color_by]])]
) {
  # Spatial image object most be filtered to contain only existing spots
  if(any(image_obj$Coordinates$in_tissue == 0)) {
    image_obj$Coordinates <- image_obj$Coordinates[which(image_obj$Coordinates$in_tissue == 1), , drop = FALSE]
  }

  data <- metadata[, c(barcodes, color_by)]
  data[, color_by] <- factor(data[[color_by]])

  coordinates <- .get_coords(image_obj) %>% tibble::rownames_to_column(var = barcodes)
  data <- merge(coordinates, data, by = barcodes, all = TRUE) %>% tibble::column_to_rownames(var = barcodes)

  if(!is.null(cols)) {
    if(length(cols) < length(levels(data[, color_by]))) {
      warning("Not enough colors provided! seting default colors...")
      cols <- gg_color_hue(length(levels(data[, color_by])))
    } else {
      cols <- cols[1:length(levels(data[, color_by]))]
    }
    names(cols) <- levels(data[, color_by])
  }
  data$fill <- cols[match(data[[color_by]], names(cols))]

  plot <- ggplot(data = data, aes_string(x = colnames(data)[2], y = colnames(data)[1]))   # , fill = colnames(data)[3]
  plot <- plot + geom_point(data = data, aes_string(fill = colnames(data)[4])) +
    scale_fill_manual(name = legend_name, values = cols)
  plot <- plot + geom_spatial(data = data,
                              fill = data$fill,
                              point.size.factor = ifelse(rownames(data) %in% rm_spots, yes = 0, no = pt.size),
                              image_obj = image_obj,
                              image.alpha = image.alpha,
                              crop = crop,
                              stroke = ifelse(rownames(data) %in% rm_spots, yes = 0, no = stroke),
                              alpha = pt.alpha)

  plot <- plot + coord_fixed() + theme(aspect.ratio = 1) + NoAxes() + theme(panel.background = element_blank())
  return(plot)
}
