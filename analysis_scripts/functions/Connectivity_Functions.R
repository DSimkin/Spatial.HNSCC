
# Connectivity Helper Functions -------------------------------------------


.neighbors_table_func <- function(metadata, spot_class = "MPid", spot_name = "Key") {
  if(is.factor(metadata[[spot_class]])) {metadata[[spot_class]] <- as.character(metadata[[spot_class]])}
  neighbors_table <- sapply(metadata[[spot_name]], function(spot) {
    spots_row = metadata$array_row[metadata[[spot_name]] == spot]
    spots_col = metadata$array_col[metadata[[spot_name]] == spot]

    n1_temp = metadata[[spot_name]][metadata$array_row == as.numeric(spots_row) - 1 & metadata$array_col == as.numeric(spots_col) - 1]
    if(length(n1_temp) == 0) {
      n1 = NA
    } else {
      n1 = as.character(metadata[[spot_class]][metadata[[spot_name]] == n1_temp])
    }

    n2_temp = metadata[[spot_name]][metadata$array_row == as.numeric(spots_row) - 1 & metadata$array_col == as.numeric(spots_col) + 1]
    if (length(n2_temp) == 0) {
      n2 = NA
    } else {
      n2 = as.character(metadata[[spot_class]][metadata[[spot_name]] == n2_temp])
    }

    n3_temp = metadata[[spot_name]][metadata$array_row == as.numeric(spots_row) & metadata$array_col == as.numeric(spots_col) - 2]
    if (length(n3_temp) == 0) {
      n3 = NA
    } else {
      n3 = as.character(metadata[[spot_class]][metadata[[spot_name]] == n3_temp])
    }

    n4_temp = metadata[[spot_name]][metadata$array_row == as.numeric(spots_row) & metadata$array_col == as.numeric(spots_col) + 2]
    if (length(n4_temp) == 0) {
      n4 = NA
    } else {
      n4 = as.character(metadata[[spot_class]][metadata[[spot_name]] == n4_temp])
    }

    n5_temp = metadata[[spot_name]][metadata$array_row == as.numeric(spots_row) + 1 & metadata$array_col == as.numeric(spots_col) - 1]
    if (length(n5_temp) == 0) {
      n5 = NA
    } else {
      n5 = as.character(metadata[[spot_class]][metadata[[spot_name]] == n5_temp])
    }

    n6_temp = metadata[[spot_name]][metadata$array_row == as.numeric(spots_row) + 1 & metadata$array_col == as.numeric(spots_col) + 1]
    if (length(n6_temp) == 0) {
      n6 = NA
    } else {
      n6 = as.character(metadata[[spot_class]][metadata[[spot_name]] == n6_temp])
    }

    return(c(n1, n2, n3, n4, n5, n6))
  })

  neighbors_table = t(neighbors_table)
  rownames(neighbors_table) = metadata[[spot_name]]
  return(neighbors_table)
}


.prog_connectivity_score <- function(program_neighbors, state) {
  state_neighbors_bin <- ifelse(program_neighbors == state, 0, 1)
  if(is.null(dim(state_neighbors_bin))) {
    prog_connect <- 1
  } else {
    prog_connect <- length(which(apply(state_neighbors_bin, 1, function(x) {sum(na.omit(x))}) > 0))
  }
  return(prog_connect)
}


.calc_adj_mat <- function(neighbored_state, spots_states, state_neighbors_table, spot_class = "MPid") {
  if(!(neighbored_state %in% spots_states[[spot_class]])) {
    return(0)
  } else {
    state_neighbors_bin <- ifelse(state_neighbors_table == neighbored_state, 1, 0)
    if(is.null(dim(state_neighbors_bin))) {
      state_neighbors_sum <- sum(na.omit(state_neighbors_bin))
    } else {
      state_neighbors_sum <- sum(apply(state_neighbors_bin, 1, function(x) {sum(na.omit(x))}))
    }
    return(state_neighbors_sum)
  }
}


.obs_coherence <- function(program_neighbors, state) {
  state_neighbors_bin <- ifelse(program_neighbors == state, 1, 0)
  if(is.null(dim(state_neighbors_bin))) {
    state_neighbors_sum <- sum(state_neighbors_bin)
  } else {
    state_neighbors_sum <- apply(state_neighbors_bin, 1, function(x) {sum(na.omit(x))})
  }
  obs <- mean(state_neighbors_sum)
  return(obs)
}


.max_expected_val <- function(spots_num){
  a <- sqrt((4 * spots_num) / (6 * sqrt(3)))
  maxval <- (6 * spots_num - 24 * a + 12) / spots_num
  return(maxval)
}


.min_expected_val <- function(rand_table, spots_states, state, spot_class = "MPid", spot_name = "Key"){
  all_minvals <- sapply(rand_table, function(neighbors_rand_table) {
    program_rand_neighbors_table <- neighbors_rand_table[rownames(neighbors_rand_table) %in% spots_states[[spot_name]][spots_states[[spot_class]] == state & !is.na(spots_states[[spot_class]])], ]
    rand_obs <- .obs_coherence(program_rand_neighbors_table, state)
    return(rand_obs)
  })
  minval <- mean(all_minvals)
  return(minval)
}



# Connectivity Bootstrapping Method ---------------------------------------


# The metadata supplied to the function must contain spatial coordinate information!
calc_spatial_neighborhood <- function(metadata,
                                      spot_class = "MPid",
                                      spot_name = "Key",
                                      samples = "all",
                                      iter = 20) {

  # Check inputs
  if(is.factor(metadata[[spot_class]])) {
    all_labels <- na.omit(levels(metadata[[spot_class]]))
  } else if(is.character(metadata[[spot_class]])) {
    all_labels <- na.omit(unique(metadata[[spot_class]]))
    metadata[[spot_class]] <- as.factor(metadata[[spot_class]])
  } else {stop("spot_class should either be a factor or character vector!")}
  stopifnot("Error: metadata must contain a column specifying sample name." =
              any(grepl("sample*", colnames(metadata), ignore.case = TRUE)))
  metadata$enumerate <- metadata[[grep("sample*", colnames(metadata), ignore.case = TRUE)]]
  if(length(samples) == 1 && samples == "all") {
    samples <- unique(metadata$enumerate)
  } else if((length(samples) > 1 || samples != "all") & all(samples %in% unique(metadata$enumerate))) {
    samples <- samples
  } else {stop("samples supplied do not match the samples present in the metadata!")}

  message("This proccess may take a while...")
  connectivity <- lapply(seq_along(samples), function(i) {
    message(paste0("Start proccessing sample: ", samples[[i]]))
    samp_meta <- metadata[metadata$enumerate == samples[[i]], ]

    ## Construct Neighbors tables
    # Randomized spot ids & position table
    rand_neighbors_table <- lapply(seq_len(iter), function(j) {
      shuffled_spots <- sample(samp_meta[[spot_name]], length(samp_meta[[spot_name]]), replace = FALSE)
      shuffled_meta <- samp_meta %>% dplyr::mutate({{spot_name}} := shuffled_spots)
      neighbors_table <- .neighbors_table_func(shuffled_meta, spot_class = spot_class, spot_name = spot_name)
      return(neighbors_table)
    })

    # Bootstrapped spots class identity
    boots_neighbors_table <- lapply(seq_len(iter), function(j) {
      new_spots <- unique(sample(samp_meta[[spot_name]], length(samp_meta[[spot_name]]), replace = TRUE))
      bootsteped_meta <- samp_meta[samp_meta[[spot_name]] %in% new_spots, ]
      neighbors_table <- .neighbors_table_func(bootsteped_meta, spot_class = spot_class, spot_name = spot_name)
      return(neighbors_table)
    })

    # Actual (observed) neighborhood table
    neighbors_table <- .neighbors_table_func(samp_meta, spot_class = spot_class, spot_name = spot_name)

    # Calculate connectivity scores (spatial coherence)
    programs_connectivity_score <- sapply(sort(all_labels), function(cluster) {
      if (!(cluster %in% samp_meta[[spot_class]])) {
        prog_score <- NaN
      } else {
        program_neighbors_table = neighbors_table[rownames(neighbors_table) %in% samp_meta[[spot_name]][samp_meta[[spot_class]] == cluster & !is.na(samp_meta[[spot_class]])], ]
        prog_score <- .prog_connectivity_score(program_neighbors_table, cluster)
      }
      return(prog_score)
    })


    ## Connectivity
    boots_adj_mat <- lapply(boots_neighbors_table, function(b_table) {
      obs_adj_mat <- sapply(sort(all_labels), function(cluster) {
        if (!(cluster %in% samp_meta[[spot_class]])) {
          zero_neigh <- rep(0, length(all_labels))
          names(zero_neigh) <- all_labels
          return(zero_neigh)
        } else {
          cluster_neighbors_table = b_table[rownames(b_table) %in% samp_meta[[spot_name]][samp_meta[[spot_class]] == cluster & !is.na(samp_meta[[spot_class]])], ]
          num_of_neighbores = sapply(sort(as.character(all_labels)), function(neighbored_cluster) {
            num <- .calc_adj_mat(neighbored_cluster, samp_meta, cluster_neighbors_table)
            return(num)
          })
          return(num_of_neighbores)
        }})

      diag(obs_adj_mat) <- 0
      weighted_adj_mat <- apply(obs_adj_mat, 2, function(x) {x / sum(x)})

      comp4mat <- programs_connectivity_score[colnames(weighted_adj_mat)]
      comp4mat[is.na(comp4mat)] <- 0
      weighted_denominator_v2 <- sapply(c(names(comp4mat)), function(prog){
        new_comp <- comp4mat
        new_comp[prog] <- 0
        new_comp <- new_comp / sum(new_comp)
        return(new_comp)
      })

      norm_adj_mat <- weighted_adj_mat / weighted_denominator_v2

      upper_mat <- (norm_adj_mat[upper.tri(norm_adj_mat)] + t(norm_adj_mat)[upper.tri(t(norm_adj_mat))]) / 2
      lower_mat <- rep(NaN, length(upper_mat))
      avg_mat <- norm_adj_mat
      avg_mat[upper.tri(avg_mat)] <- upper_mat
      avg_mat[lower.tri(avg_mat)] <- lower_mat

      rownames(avg_mat) <- colnames(avg_mat)
      avg_mat <- t(avg_mat)
      avg_mat <- as.data.frame(avg_mat)
      avg_mat$pair2 <- rownames(avg_mat)
      long <- reshape2::melt(data.table::setDT(avg_mat), id.vars = c("pair2"), variable.name = "pair1")
      return(long)
    })


    # Random connectivity
    rand_adj_mat <- lapply(rand_neighbors_table, function(b_table) {
      obs_adj_mat <- sapply(sort(all_labels), function(cluster) {
        if (!(cluster %in% samp_meta[[spot_class]])) {
          return(rep(0, length(all_labels)))
        } else {
          cluster_neighbors_table = b_table[row.names(b_table) %in% samp_meta[[spot_name]][samp_meta[[spot_class]] == cluster & !is.na(samp_meta[[spot_class]])], ]
          num_of_neighbores = sapply(sort(as.character(all_labels)), function(neighbored_cluster) {
            num <- .calc_adj_mat(neighbored_cluster, samp_meta, cluster_neighbors_table)
            return(num)
          })
          return(num_of_neighbores)
        }})

      diag(obs_adj_mat) <- 0
      weighted_adj_mat <- apply(obs_adj_mat, 2, function(x) {x / sum(x)})

      comp4mat <- programs_connectivity_score[colnames(weighted_adj_mat)]
      comp4mat[is.na(comp4mat)] <- 0
      weighted_denominator_v2 <- sapply(c(names(comp4mat)), function(prog) {
        new_comp <- comp4mat
        new_comp[prog] <- 0
        new_comp <- new_comp / sum(new_comp)
        return(new_comp)
      })

      norm_adj_mat <- weighted_adj_mat / weighted_denominator_v2

      upper_mat <- (norm_adj_mat[upper.tri(norm_adj_mat)] + t(norm_adj_mat)[upper.tri(t(norm_adj_mat))]) / 2
      lower_mat <- rep(NaN, length(upper_mat))
      avg_mat <- norm_adj_mat
      avg_mat[upper.tri(avg_mat)] <- upper_mat
      avg_mat[lower.tri(avg_mat)] <- lower_mat

      rownames(avg_mat) <- colnames(avg_mat)
      avg_mat <- t(avg_mat)
      avg_mat <- as.data.frame(avg_mat)
      avg_mat$pair2 <- rownames(avg_mat)
      long <- reshape2::melt(data.table::setDT(avg_mat), id.vars = c("pair2"), variable.name = "pair1")
      return(long)
    })

    final_adj_mat <- data.frame(pair1 = boots_adj_mat[[1]]$pair1,
                                pair2 = boots_adj_mat[[1]]$pair2,
                                connectivity = apply(sapply(boots_adj_mat, function(x) {return(x$value)}), 1, function(k) {mean(na.omit(k))}),
                                effect_size = apply(sapply(boots_adj_mat, function(x) {return(x$value)}), 1, function(k) {mean(na.omit(k))}) / apply(sapply(rand_adj_mat, function(x) {return(x$value)}), 1, function(k) {mean(na.omit(k))}),
                                sd = apply(sapply(boots_adj_mat, function(x) {return(x$value)}), 1, function(k) {sd(na.omit(k))}))


    ## Add p-value
    pval <- sapply(seq_len(nrow(final_adj_mat)), function(j1) {
      obs <- na.omit(as.numeric(sapply(seq_len(iter), function(j2) {
        return(boots_adj_mat[[j2]][j1, "value"])
      })))
      exp <- na.omit(as.numeric(sapply(seq_len(iter), function(j2) {
        return(rand_adj_mat[[j2]][j1, "value"])
      })))
      if (length(obs) == 0) {
        return(NA)
      } else {
        t.res <- t.test(obs, exp, alternative = "two.sided", var.equal = FALSE)
        return(t.res$p.value)
      }
    })

    final_adj_mat$pval <- pval
    return(final_adj_mat)
  })
  names(connectivity) <- samples
  return(connectivity)
}





# Connectivity Permutation Method ----------------------------------------------------

# The metadata supplied to the function must contain spatial coordinate information!
neighbor_spot_props <- function(metadata,
                                zone = "All",
                                site = "All",
                                samples = "All",
                                spot_class = "MPid",
                                spot_name = "Key",
                                zone_by = "EpiStroma",
                                n_cores = 10,
                                n_perm = 1000,
                                signif_val = 0.001,
                                plot_perm_distr = FALSE,
                                filter_signif = TRUE,
                                zscore_thresh = 1) {
  # Load variables
  all_states <- names(readRDS(file = here("results/Generated_Data/Final_Metaprograms_Extended.rds")))
  samples_metadata <- readRDS(file = here("metadata/samples_metadata.rds"))
  neighbs_stats_ls <- list()
  signif_neighbs_ls <- list()
  distr_plots_ls <- list()

  # Filter metadata input by selected site and zone
  stopifnot("Error: metadata must contain a column specifying sample name." =
              any(grepl("sample*", colnames(metadata), ignore.case = TRUE)))
  if(site != "All") {
    stopifnot("Error: Site argument must be one of the following: 'Laryngeal' / 'Oral' / 'Oropharynx'." =
                site %in% unique(samples_metadata$Site))
    metadata <- metadata[metadata$Sample %in% samples_metadata$Sample[samples_metadata$Site == site], ]
  }
  if(zone != "All") {
    stopifnot("Error: Argument `zone_by` must be either 'EpiStroma' (for separating Epithelial, Stromal or Mixed spots) or 'Zone' (for separating Epithelial zonation)." =
                zone_by %in% c("EpiStroma", "Zone") & length(zone_by) == 1)
    if(zone_by == "EpiStroma") {
      stopifnot("Error: Argument `zone` should specify on which tumor region neighboring states will be computed - 'Epithelial' / 'Stroma' / 'Mixed'." =
                  zone %in% unique(metadata$EpiStroma) & length(zone) == 1)
      metadata <- metadata[metadata$EpiStroma == zone, ]
    }
    if(zone_by == "Zone") {
      stopifnot("Error: The variable `Zone` is not found in the metadata. Run classify_zones function first." =
                  any(colnames(metadata) %in% "Zone"))
      stopifnot("Error: Argument `zone` should specify on which Epithelial zone neighboring states will be computed - 'Zone_1' / 'Zone_2' / 'Zone_3'." =
                  zone %in% unique(metadata$Zone) & length(zone) == 1)
      metadata <- metadata[metadata$Zone == zone, ]
    }
  }
  metadata$enumerate <- metadata[[grep("sample*", colnames(metadata), ignore.case = TRUE)]]
  if(length(samples) == 1 && samples == "All") {
    samples <- unique(metadata$enumerate)
  } else if((length(samples) > 1 || samples != "All") & all(samples %in% unique(metadata$enumerate))) {
    samples <- samples
  } else {stop("samples supplied do not match the samples present in the metadata!")}

  for(samp in samples) {
    message(paste0("Processing sample: ", samp))

    ### ============ Actual state pairs connectivity values ============
    # Construct neighboring table for the sample
    samp_meta <- metadata[metadata$Sample == samp, ]
    neighbors_table <- .neighbors_table_func(samp_meta, spot_class = spot_class, spot_name = spot_name)

    # Calculate connectivity scores (spatial coherence)
    programs_connectivity_score <- sapply(sort(all_states), function(state) {
      if (!(state %in% samp_meta[[spot_class]])) {
        prog_score <- NaN
      } else {
        program_neighbors_table = neighbors_table[rownames(neighbors_table) %in% samp_meta[[spot_name]][samp_meta[[spot_class]] == state & !is.na(samp_meta[[spot_class]])], ]
        prog_score <- .prog_connectivity_score(program_neighbors_table, state)
      }
      return(prog_score)
    })

    # Count neighboring states for each spot (number of free (non-coherent) X state classified spots that neighbor free Y reference-state)
    obs_adj_mat <- sapply(sort(all_states), function(state) {
      if (!(state %in% samp_meta[[spot_class]])) {
        zero_neigh <- rep(0, length(all_states))
        names(zero_neigh) <- all_states
        return(zero_neigh)
      } else {
        state_neighbors_table = neighbors_table[rownames(neighbors_table) %in% samp_meta[[spot_name]][samp_meta[[spot_class]] == state & !is.na(samp_meta[[spot_class]])], ]
        num_of_neighbors = sapply(sort(as.character(all_states)), function(neighbored_state) {
          num <- .calc_adj_mat(neighbored_state, samp_meta, state_neighbors_table)
          return(num)
        })
        return(num_of_neighbors)
      }})
    diag(obs_adj_mat) <- 0
    weighted_adj_mat <- apply(obs_adj_mat, 2, function(x) {x / sum(x)})

    # Calculate corrected proportion of neighboring spots
    comp4mat <- programs_connectivity_score[colnames(weighted_adj_mat)]
    comp4mat[is.na(comp4mat)] <- 0
    weighted_denominator_v2 <- sapply(c(names(comp4mat)), function(prog){
      new_comp <- comp4mat
      new_comp[prog] <- 0
      new_comp <- new_comp / sum(new_comp)
      return(new_comp)
    })

    norm_adj_mat <- weighted_adj_mat / weighted_denominator_v2
    baseline_stat <- Melt(norm_adj_mat)

    ### ============ Permuted state pairs connectivity values ============
    # Permute neighbor table `n_perm` times, to create a permuted sampling distribution which will serve as the null distribution
    permute_neighbs <- parallel::mclapply(1:n_perm, function(x) {
      shuffled_spots <- sample(samp_meta[[spot_name]], length(samp_meta[[spot_name]]), replace = FALSE)
      shuffled_meta <- samp_meta %>% dplyr::mutate({{spot_name}} := shuffled_spots)
      neighbors_table <- .neighbors_table_func(shuffled_meta, spot_class = spot_class, spot_name = spot_name)
      neighbors_table <- neighbors_table[match(samp_meta[[spot_name]], rownames(neighbors_table)), ]

      obs_adj_mat <- sapply(sort(all_states), function(state) {
        if (!(state %in% samp_meta[[spot_class]])) {
          zero_neigh <- rep(0, length(all_states))
          names(zero_neigh) <- all_states
          return(zero_neigh)
        } else {
          state_neighbors_table = neighbors_table[rownames(neighbors_table) %in% samp_meta[[spot_name]][samp_meta[[spot_class]] == state & !is.na(samp_meta[[spot_class]])], ]
          num_of_neighbors = sapply(sort(as.character(all_states)), function(neighbored_state) {
            num <- .calc_adj_mat(neighbored_state, samp_meta, state_neighbors_table)
            return(num)
          })
          return(num_of_neighbors)
        }})
      diag(obs_adj_mat) <- 0
      weighted_adj_mat <- apply(obs_adj_mat, 2, function(x) {x / sum(x)})

      comp4mat <- programs_connectivity_score[colnames(weighted_adj_mat)]
      comp4mat[is.na(comp4mat)] <- 0
      weighted_denominator_v2 <- sapply(c(names(comp4mat)), function(prog){
        new_comp <- comp4mat
        new_comp[prog] <- 0
        new_comp <- new_comp / sum(new_comp)
        return(new_comp)
      })

      norm_adj_mat <- weighted_adj_mat / weighted_denominator_v2
      perm_stat <- Melt(norm_adj_mat)
    }, mc.cores = n_cores)

    perm_stats <- cbind.data.frame(permute_neighbs[[1]]$rows, permute_neighbs[[1]]$cols, do.call(cbind.data.frame, lapply(permute_neighbs, function(y) y$vals))) %>%
      magrittr::set_colnames(c("rows", "cols", paste0("perm_", seq(n_perm))))



    # Compare observed neighboring state proportion to the permuted null distribution - extract Neighbor-proportion, P-value & Permuted distribution statistics
    # merged_df <- merge(baseline_stat, perm_stats, by = c("rows", "cols"))
    merged_df <- na.omit(merge(baseline_stat, perm_stats, by = c("rows", "cols")))
    scaled_vals <- merged_df %>% dplyr::select(-c("rows", "cols")) %>% t() %>% scale() %>% t() %>%
      as.data.frame() %>% dplyr::mutate(rows = merged_df$rows, cols = merged_df$cols, .before = 1) %>% dplyr::pull(vals)
    neighbs_stats <- lapply(seq_len(nrow(merged_df)), function(i) {
      summary_stats <- summary(as.numeric(merged_df[i, grep("perm_", colnames(merged_df))]))
      if(merged_df[i, "vals"] == 0) {
        pvals <- 1
        tail_direction <- NA
      } else {
        Rtail_pvals <- (sum(as.numeric(merged_df[i, grep("perm_", colnames(merged_df))]) >= merged_df[i, "vals"]) + 1) / (n_perm + 1)
        Ltail_pvals <- (n_perm - sum(as.numeric(merged_df[i, grep("perm_", colnames(merged_df))]) > merged_df[i, "vals"]) + 1) / (n_perm + 1)
        tail_direction <- Rtail_pvals < Ltail_pvals
        pvals <- Rtail_pvals * tail_direction + Ltail_pvals * !tail_direction
      }
      tail_dir <- ifelse(is.na(tail_direction), NA,
                         ifelse(tail_direction, "Drawn", "Repelled"))
      sds <- sd(as.numeric(merged_df[i, grep("perm_", colnames(merged_df))]))
      out <- data.frame(Neighb_State = merged_df$rows[i], Ref_State = merged_df$cols[i], Prop_Neighb = merged_df$vals[i], Z_Score = scaled_vals[i],
                        Perm_Min = summary_stats[["Min."]], Perm_Max = summary_stats[["Max."]], Perm_Mean = summary_stats[["Mean"]], Perm_Median = summary_stats[["Median"]],
                        P_val = pvals, Interaction_Type = tail_dir, SD = sds)
    })
    neighbs_stats <- do.call(rbind.data.frame, neighbs_stats)
    neighbs_stats$Significant <- neighbs_stats$P_val < signif_val
    neighbs_stats_ls[[samp]] <- neighbs_stats

    if(filter_signif) {
      top_neighbs_df <- neighbs_stats %>%
        dplyr::filter(Neighb_State != Ref_State) %>%
        dplyr::mutate(State_Pair = paste0(Ref_State, ".", Neighb_State)) %>%
        dplyr::filter(Significant == TRUE) %>%
        mutate(Signif_Pair = unname(sapply(.$State_Pair, function(x) {
          ifelse(paste0(scalop::substri(x, pos = 2, sep = "\\."), ".", scalop::substri(x, pos = 1, sep = "\\.")) %in% .$State_Pair,
                 yes = "Yes", no = "No")
        }))) %>% dplyr::filter(abs(Z_Score) >= zscore_thresh, Signif_Pair == "Yes") %>% dplyr::arrange(desc(Z_Score))
      signif_neighbs_ls[[samp]] <- top_neighbs_df
    }

    # Generate plots
    if(plot_perm_distr) {
      plot_df <- reshape2::melt(perm_stats)
      distr_plot <- ggplot(plot_df, aes(x = value)) +
        facet_grid(rows ~ cols) +
        geom_histogram() +
        geom_vline(data = baseline_stat, aes(xintercept = vals), color = "red")
      distr_plots_ls[[samp]] <- distr_plot
    }
  }

  if(isTRUE(plot_perm_distr) & isFALSE(filter_signif)) {return(list(Neighbor_Stats = neighbs_stats_ls, Distr_plots = distr_plots_ls))}
  else if(isTRUE(filter_signif) & isFALSE(plot_perm_distr)) {return(list(Neighbs_Stats = neighbs_stats_ls, Top_Neighbs = signif_neighbs_ls))}
  else if(isTRUE(plot_perm_distr) & isTRUE(filter_signif)) {return(list(Neighbs_Stats = neighbs_stats_ls, Top_Neighbs = signif_neighbs_ls, Distr_plots = distr_plots_ls))}
  else {return(neighbs_stats_ls)}
}





# Spatial Coherence Score Calculation Function ----------------------------

# The metadata supplied to the function must contain spatial coordinate information!
calc_spatial_coherence <- function(metadata,
                                zone = "All",
                                site = "All",
                                samples = "All",
                                spot_class = "MPid",
                                spot_name = "Key",
                                zone_by = "EpiStroma",
                                n_cores = 10,
                                n_perm = 100) {
  # Load variables
  samples_metadata <- readRDS(file = here("metadata/samples_metadata.rds"))
  coherence_scores_ls <- list()
  if(spot_class == "MPid") {
    all_states <- names(readRDS(file = here("results/Generated_Data/Final_Metaprograms_Extended.rds")))
  } else if(spot_class == "Subclone") {
    all_states <- seq(1, max(metadata$Subclone[!is.na(metadata$Subclone)]))
  } else if(spot_class == "binCNAstatus") {
    all_states <- c("Malignant", "Mixed", "Non_Malignant")
  } else if(spot_class == "TLS") {
    all_states <- c("Not TLS", "Posible TLS")
  }

  # Filter metadata input by selected site and zone
  stopifnot("Error: metadata must contain a column specifying sample name." =
              any(grepl("sample*", colnames(metadata), ignore.case = TRUE)))
  if(site != "All") {
    stopifnot("Error: Site argument must be one of the following: 'Laryngeal' / 'Oral' / 'Oropharynx'." =
                site %in% unique(samples_metadata$Site))
    metadata <- metadata[metadata$Sample %in% samples_metadata$Sample[samples_metadata$Site == site], ]
  }
  if(zone != "All" || length(zone) > 1) {
    stopifnot("Error: Argument `zone_by` must be either 'EpiStroma' (for separating Epithelial, Stromal or Mixed spots) or 'Zone' (for separating Epithelial zonation)." =
                zone_by %in% c("EpiStroma", "Zone") & length(zone_by) == 1)
    if(zone_by == "EpiStroma") {
      stopifnot("Error: Argument `zone` should specify on which tumor region neighboring states will be computed - 'Epithelial' / 'Stroma' / 'Mixed'." =
                  zone %in% unique(metadata$EpiStroma) & length(zone) == 1)
      metadata <- metadata[metadata$EpiStroma == zone, ]
    }
    if(zone_by == "Zone") {
      stopifnot("Error: The variable `Zone` is not found in the metadata. Run classify_zones function first." =
                  any(colnames(metadata) %in% "Zone"))
      if(length(zone) == 1) {
        stopifnot("Error: Argument `zone` should specify on which Epithelial zone neighboring states will be computed - 'Zone_1' / 'Zone_2' / 'Zone_3'." =
                    zone %in% unique(metadata$Zone) & length(zone) == 1)
        metadata <- metadata[metadata$Zone == zone, ]
      }
      if(length(zone) > 1) {
        stopifnot("Error: Argument `zone` should specify on which Epithelial zone, or combination of zones, neighboring states will be computed - 'Zone_1', 'Zone_2', 'Zone_3'." =
                    zone %in% unique(metadata$Zone) & length(zone) > 1)
        metadata <- metadata[metadata$Zone %in% zone, ]
      }
    }
  }
  metadata$enumerate <- metadata[[grep("sample*", colnames(metadata), ignore.case = TRUE)]]
  if(length(samples) == 1 && samples == "All") {
    samples <- unique(metadata$enumerate)
  } else if((length(samples) > 1 || samples != "All") & all(samples %in% unique(metadata$enumerate))) {
    samples <- samples
  } else {stop("samples supplied do not match the samples present in the metadata!")}


  for(samp in samples) {
    message(paste0("Processing sample: ", samp))

    # Construct neighboring table for the sample
    samp_meta <- metadata[metadata$Sample == samp, ] %>% dplyr::filter(!is.na(MPid))
    neighbors_table <- .neighbors_table_func(samp_meta, spot_class = spot_class, spot_name = spot_name)
    neighbors_table[is.na(neighbors_table)] <- NaN

    # Construct several random (shuffled) neighboring tables for the sample
    permute_neighbs <- parallel::mclapply(1:n_perm, function(x) {
      shuffled_spots <- sample(samp_meta[[spot_name]], length(samp_meta[[spot_name]]), replace = FALSE)
      shuffled_meta <- samp_meta %>% dplyr::mutate({{spot_name}} := shuffled_spots)
      neighbors_table <- .neighbors_table_func(shuffled_meta, spot_class = spot_class, spot_name = spot_name)
      neighbors_table[is.na(neighbors_table)] <- NaN
      neighbors_table <- neighbors_table[match(samp_meta[[spot_name]], rownames(neighbors_table)), ]
    })

    # Calculate spatial coherence scores
    programs_coherence_score <- sapply(sort(all_states), function(state) {
      if(!(state %in% samp_meta[[spot_class]]) || sum(samp_meta[[spot_class]] == state & !is.na(samp_meta[[spot_class]])) < 2) {
        prog_score <- NaN
      } else {
        program_neighbors_table = neighbors_table[rownames(neighbors_table) %in% samp_meta[[spot_name]][samp_meta[[spot_class]] == state & !is.na(samp_meta[[spot_class]])], ]
        obs_coherence <- .obs_coherence(program_neighbors_table, state)
        max_exp_score <- .max_expected_val(nrow(program_neighbors_table))
        min_exp_score <- .min_expected_val(permute_neighbs, samp_meta, state, spot_class = spot_class, spot_name = spot_name)
        if(obs_coherence > max_exp_score) {
          obs_coherence <- max_exp_score
          }
        if(obs_coherence < min_exp_score) {
          obs_coherence <- min_exp_score
        }

        prog_score <- (obs_coherence - min_exp_score)/(max_exp_score - min_exp_score)
      }
      return(prog_score)
    })
    coherence_scores_ls[[samp]] <- programs_coherence_score
  }
  return(coherence_scores_ls)
}



# Calculating total-tumor, MPs and subclones coherence for malignant spots
calc_malig_spatial_coherence <- function(metadata,
                                   samples = "All",
                                   spot_name = "Key",
                                   n_cores = 10,
                                   n_perm = 100) {
  # Load variables
  samples_metadata <- readRDS(file = here("metadata/samples_metadata.rds"))
  coherence_scores_ls <- list()
  all_states <- names(readRDS(file = here("results/Generated_Data/Final_Metaprograms.rds")))
  all_clones <- seq(1, max(metadata$Subclone[!is.na(metadata$Subclone)]))

  # Filter metadata input by selected site and zone
  stopifnot("Error: metadata must contain a column specifying sample name." =
              any(grepl("sample*", colnames(metadata), ignore.case = TRUE)))

  metadata$enumerate <- metadata[[grep("sample*", colnames(metadata), ignore.case = TRUE)]]
  if(length(samples) == 1 && samples == "All") {
    samples <- unique(metadata$enumerate)
  } else if((length(samples) > 1 || samples != "All") & all(samples %in% unique(metadata$enumerate))) {
    samples <- samples
  } else {stop("samples supplied do not match the samples present in the metadata!")}

  # metadata <- metadata[metadata$binCNAstatus == "Malignant", ]

  for(samp in samples) {
    message(paste0("Processing sample: ", samp))
    samp_df <- data.frame(Tumor = NA, MPs_Mean = NA, MPs_Sd = NA, Subclones_Mean = NA, Subclones_Sd = NA, row.names = samp)

    ### ============ Total-Tumor Coherence score ============
    # Construct neighboring table for the sample
    samp_meta <- metadata[metadata$Sample == samp, ] %>% dplyr::filter(!is.na(MPid))
    neighbors_table <- .neighbors_table_func(samp_meta, spot_class = "binCNAstatus", spot_name = spot_name)
    neighbors_table[is.na(neighbors_table)] <- NaN

    # Construct several random (shuffled) neighboring tables for the sample
    permute_neighbs <- parallel::mclapply(1:n_perm, function(x) {
      shuffled_spots <- sample(samp_meta[[spot_name]], length(samp_meta[[spot_name]]), replace = FALSE)
      shuffled_meta <- samp_meta %>% dplyr::mutate({{spot_name}} := shuffled_spots)
      neighbors_table <- .neighbors_table_func(shuffled_meta, spot_class = "binCNAstatus", spot_name = spot_name)
      neighbors_table[is.na(neighbors_table)] <- NaN
      neighbors_table <- neighbors_table[match(samp_meta[[spot_name]], rownames(neighbors_table)), ]
    })

    # Calculate spatial coherence scores
    if(!("Malignant" %in% samp_meta$binCNAstatus) || sum(samp_meta$binCNAstatus == "Malignant" & !is.na(samp_meta$binCNAstatus)) < 2) {
      prog_score <- NaN
    } else {
      program_neighbors_table = neighbors_table[rownames(neighbors_table) %in% samp_meta[[spot_name]][samp_meta$binCNAstatus == "Malignant" & !is.na(samp_meta$binCNAstatus)], ]
      obs_coherence <- .obs_coherence(program_neighbors_table, "Malignant")
      max_exp_score <- .max_expected_val(nrow(program_neighbors_table))
      min_exp_score <- .min_expected_val(permute_neighbs, samp_meta, "Malignant", spot_class = "binCNAstatus", spot_name = spot_name)
      if(obs_coherence > max_exp_score) {
        obs_coherence <- max_exp_score
      }
      if(obs_coherence < min_exp_score) {
        obs_coherence <- min_exp_score
      }

      prog_score <- (obs_coherence - min_exp_score)/(max_exp_score - min_exp_score)
    }
    samp_df$Tumor <- prog_score

    ### ============ MPs Coherence score ============
    # Filter sample metadata to include only malignant spots
    filt_samp_meta <- samp_meta[samp_meta$binCNAstatus == "Malignant", ]
    # Construct neighboring table for the sample
    neighbors_table <- .neighbors_table_func(filt_samp_meta, spot_class = "MPid", spot_name = spot_name)
    neighbors_table[is.na(neighbors_table)] <- NaN

    # Construct several random (shuffled) neighboring tables for the sample
    permute_neighbs <- parallel::mclapply(1:n_perm, function(x) {
      shuffled_spots <- sample(filt_samp_meta[[spot_name]], length(filt_samp_meta[[spot_name]]), replace = FALSE)
      shuffled_meta <- filt_samp_meta %>% dplyr::mutate({{spot_name}} := shuffled_spots)
      neighbors_table <- .neighbors_table_func(shuffled_meta, spot_class = "MPid", spot_name = spot_name)
      neighbors_table[is.na(neighbors_table)] <- NaN
      neighbors_table <- neighbors_table[match(filt_samp_meta[[spot_name]], rownames(neighbors_table)), ]
    })

    # Calculate spatial coherence scores
    programs_coherence_score <- sapply(sort(all_states), function(state) {
      if(!(state %in% filt_samp_meta$MPid) || sum(filt_samp_meta$MPid == state & !is.na(filt_samp_meta$MPid)) < 2) {
        prog_score <- NaN
      } else {
        program_neighbors_table = neighbors_table[rownames(neighbors_table) %in% filt_samp_meta[[spot_name]][filt_samp_meta$MPid == state & !is.na(filt_samp_meta$MPid)], ]
        obs_coherence <- .obs_coherence(program_neighbors_table, state)
        max_exp_score <- .max_expected_val(nrow(program_neighbors_table))
        min_exp_score <- .min_expected_val(permute_neighbs, filt_samp_meta, state, spot_class = "MPid", spot_name = spot_name)
        if(obs_coherence > max_exp_score) {
          obs_coherence <- max_exp_score
        }
        if(obs_coherence < min_exp_score) {
          obs_coherence <- min_exp_score
        }

        prog_score <- (obs_coherence - min_exp_score)/(max_exp_score - min_exp_score)
      }
      return(prog_score)
    })
    # samp_df$MPs_Mean <- mean(programs_coherence_score[c("Cell_Cycle", "Epithelial", "Hypoxia", "LowQ", "pEMT", "Secretory", "Senescence", "TNF_Signaling")], na.rm = TRUE)
    # samp_df$MPs_Sd <- sd(programs_coherence_score[c("Cell_Cycle", "Epithelial", "Hypoxia", "LowQ", "pEMT", "Secretory", "Senescence", "TNF_Signaling")], na.rm = TRUE)
    samp_df$MPs_Mean <- mean(programs_coherence_score, na.rm = TRUE)
    samp_df$MPs_Sd <- sd(programs_coherence_score, na.rm = TRUE)

    ### ============ Subclones Coherence score ============
    # Construct neighboring table for the sample
    neighbors_table <- .neighbors_table_func(filt_samp_meta, spot_class = "Subclone", spot_name = spot_name)
    neighbors_table[is.na(neighbors_table)] <- NaN

    # Construct several random (shuffled) neighboring tables for the sample
    permute_neighbs <- parallel::mclapply(1:n_perm, function(x) {
      shuffled_spots <- sample(filt_samp_meta[[spot_name]], length(filt_samp_meta[[spot_name]]), replace = FALSE)
      shuffled_meta <- filt_samp_meta %>% dplyr::mutate({{spot_name}} := shuffled_spots)
      neighbors_table <- .neighbors_table_func(shuffled_meta, spot_class = "Subclone", spot_name = spot_name)
      neighbors_table[is.na(neighbors_table)] <- NaN
      neighbors_table <- neighbors_table[match(filt_samp_meta[[spot_name]], rownames(neighbors_table)), ]
    })

    # Calculate spatial coherence scores
    clones_coherence_score <- sapply(sort(all_clones), function(clone) {
      if(!(clone %in% filt_samp_meta$Subclone) || sum(filt_samp_meta$Subclone == clone & !is.na(filt_samp_meta$Subclone)) < 2) {
        prog_score <- NaN
      } else {
        program_neighbors_table = neighbors_table[rownames(neighbors_table) %in% filt_samp_meta[[spot_name]][filt_samp_meta$Subclone == clone & !is.na(filt_samp_meta$Subclone)], ]
        obs_coherence <- .obs_coherence(program_neighbors_table, clone)
        max_exp_score <- .max_expected_val(nrow(program_neighbors_table))
        min_exp_score <- .min_expected_val(permute_neighbs, filt_samp_meta, clone, spot_class = "Subclone", spot_name = spot_name)
        if(obs_coherence > max_exp_score) {
          obs_coherence <- max_exp_score
        }
        if(obs_coherence < min_exp_score) {
          obs_coherence <- min_exp_score
        }

        prog_score <- (obs_coherence - min_exp_score)/(max_exp_score - min_exp_score)
      }
      return(prog_score)
    })
    samp_df$Subclones_Mean <- mean(clones_coherence_score, na.rm = TRUE)
    samp_df$Subclones_Sd <- sd(clones_coherence_score, na.rm = TRUE)

    # Output to a list
    coherence_scores_ls[[samp]] <- samp_df
  }
  return(coherence_scores_ls)
}
