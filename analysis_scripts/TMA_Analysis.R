library(tidyverse)
library(scalop)
library(ggpubr)
library(here)


# pEMT Pattern Intra-tumor Consistancy ------------------------------------

labels_tab <- read.csv(here("aux_data/TMA_LAMC2_patterns.csv"))
labels_tab <- labels_tab %>% dplyr::filter(grepl("HNOC", ID)) %>% 
  dplyr::mutate(LAMC2_pattern = dplyr::case_when(grepl("diffuse", LAMC2_pattern, ignore.case = TRUE) ~ "Diffuse",
                                                 grepl("edge", LAMC2_pattern, ignore.case = TRUE) ~ "Edge", .default = LAMC2_pattern)) %>% 
  dplyr::filter(!is.na(LAMC2_pattern), !LAMC2_pattern %in% c("Normal", "Stromal"))
grouping_tab <- labels_tab %>% dplyr::count(ID, LAMC2_pattern) %>% dplyr::group_by(ID) %>% dplyr::mutate(Replicates = sum(n)) %>% dplyr::ungroup() %>% 
  dplyr::group_by(ID) %>% dplyr::slice(which.max(n)) %>% dplyr::mutate(Consistency = n/Replicates) %>% ungroup()
order_patients <- grouping_tab %>% dplyr::group_by(LAMC2_pattern) %>% dplyr::arrange(desc(Consistency), .by_group = TRUE) %>% dplyr::pull(ID)
plot_df <- as.data.frame(table(labels_tab$ID, labels_tab$LAMC2_pattern)) %>% 
  dplyr::mutate(Var1 = factor(Var1, levels = order_patients),
                Group = ifelse(grouping_tab$Consistency[match(Var1, grouping_tab$ID)] >= 0.75, grouping_tab$LAMC2_pattern[match(Var1, grouping_tab$ID)], "Non-consistant"))
p <- ggplot(plot_df, aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "LAMC2 Pattern", values = setNames(c("#B24745FF", "#DF8F44FF", "#374E55FF", "#79AF97FF"), c("Diffuse", "Edge", "Low", "Mixed"))) +
  theme_classic() + scale_y_continuous(expand = c(0, 0)) + labs(x = "", y = "Core Number") +
  theme(text = element_text(size = 20), title = element_text(size = 24), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.text = element_text(size = 24),
        strip.background = element_blank(), strip.text.x = element_text(size = 24, face = "bold")) + 
  facet_grid(.~Group, scales = 'free', space = 'free', drop = F)
# ggplot2::ggsave(filename = here("results/Paper_Figures/Main_figs/Fig_5c.pdf"), plot = p, device = "pdf", dpi = 300, width = 25, height = 6)


# Percentage of samples which show intra-sample core agreement greater than 3/4
# Calculate observed proportions
consist_thresh <- 0.75
consist_props <- labels_tab %>% dplyr::count(ID, LAMC2_pattern) %>% dplyr::group_by(ID) %>% dplyr::mutate(Replicates = sum(n)) %>% dplyr::ungroup() %>% 
  dplyr::group_by(ID) %>% dplyr::slice(which.max(n)) %>% dplyr::mutate(Consistency = n/Replicates) %>% ungroup() %>% 
  dplyr::mutate(Group = ifelse(Consistency >= consist_thresh, "Consistant", "Non-consistant"))
obs_prop <-  consist_props %>% dplyr::group_by(Group) %>% dplyr::summarize(Freq = n()) %>%
  dplyr::summarize(sum(Freq[Group == "Consistant"]) / sum(Freq)) %>% as.numeric()

# Bootstrapping for confidence intervals
set.seed(5410) 
n_boot <- 10000
bootstrap_props <- numeric(n_boot)

for (i in 1:n_boot) {
  bootstrap_prop <- consist_props %>% dplyr::mutate(Group = sample(Group, replace = TRUE)) %>% dplyr::group_by(Group) %>% dplyr::summarize(Freq = n()) %>%
    dplyr::summarize(sum(Freq[Group == "Consistant"]) / sum(Freq)) %>% as.numeric()
  bootstrap_props[i] <- bootstrap_prop
}

# Confidence intervals
ci_lower <- quantile(bootstrap_props, 0.025)
ci_upper <- quantile(bootstrap_props, 0.975)
cat("95% Confidence Interval for Percentage Agreement:", ci_lower, "-", ci_upper, "\n")

# Permutation testing for p-value
n_perm <- 10000
permuted_props <- numeric(n_perm)

for (i in 1:n_perm) {
  permuted_prop <- labels_tab %>% dplyr::select(ID, LAMC2_pattern) %>% dplyr::mutate(LAMC2_pattern = sample(LAMC2_pattern, replace = FALSE)) %>% 
    dplyr::count(ID, LAMC2_pattern) %>% dplyr::group_by(ID) %>% dplyr::mutate(Replicates = sum(n)) %>% dplyr::ungroup() %>% 
    dplyr::group_by(ID) %>% dplyr::slice(which.max(n)) %>% dplyr::mutate(Consistency = n/Replicates) %>% ungroup() %>% 
    dplyr::mutate(Group = ifelse(Consistency >= consist_thresh, "Consistant", "Non-consistant")) %>% dplyr::group_by(Group) %>% dplyr::summarize(Freq = n()) %>%
    dplyr::summarize(sum(Freq[Group == "Consistant"]) / sum(Freq)) %>% as.numeric()
  permuted_props[i] <- permuted_prop
}

# Calculate p-value (two-tailed)
Rtail_pvals <- (sum(as.numeric(permuted_props) >= obs_prop) + 1) / (n_perm + 1)
Ltail_pvals <- (n_perm - sum(as.numeric(permuted_props) > obs_prop) + 1) / (n_perm + 1)
tail_direction <- Rtail_pvals < Ltail_pvals
pvals <- Rtail_pvals * tail_direction + Ltail_pvals * !tail_direction
cat("P-value for Percentage Agreement:", pvals, "\n")

# Plot the bootstrapped and permuted distributions
plot_df <- data.frame(Permuted = permuted_props,
                      Bootstrapped = bootstrap_props) %>% reshape2::melt(.)
# saveRDS(plot_df, file = here("results/Generated_data/TMA_label_agreement_of_0.75_bootstrap_and_permute_vals.rds"))

plot_df <- readRDS(file = here("results/Generated_data/TMA_label_agreement_of_0.75_bootstrap_and_permute_vals.rds"))
p <- ggplot(plot_df, aes(x = value)) + geom_histogram(bins = 100, fill = "skyblue", color = "black") +
  annotate("segment", x = obs_prop, xend = obs_prop, y = 0, yend = 850, color = "red", linewidth = 1, linetype = "dashed") +
  annotate("text", label = "Permuted Distribution", x = median(plot_df$value[plot_df$variable == "Permuted"]), y = 1100, size = 6) +
  annotate("text", label = "Bootstrap Distribution", x = obs_prop, y = 970, size = 6) +
  annotate("text", label = paste("Observed Proportion =", round(obs_prop, 3)), x = obs_prop, y = 900, size = 4, color = "red") +
  annotate("text", label = paste0("P-value = ", round(pvals, 5)), x = 0.1, y = 1150, size = 5) +
  theme_scalop() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "Proporation of samples with at least 3/4 core agreement", y = "Frequency") + 
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1180))
ggplot2::ggsave(filename = here("results/Paper_Figures/Supp_figs/SFig_4e.pdf"), device = "pdf", plot = p, dpi = 300)



# Krippendorff’s Alpha coefficients for labels agreement (variation of Fleiss Kappa that better handle NAs generated by the unequal number of cores per tumor)
library(irr)
per_tumor_labs <- labels_tab %>% dplyr::select(ID, LAMC2_pattern) %>% dplyr::arrange(ID) %>%
  dplyr::mutate(idx = ave(ID, ID, FUN = seq_along)) %>% reshape2::dcast(ID ~ idx, value.var = "LAMC2_pattern") %>% dplyr::rename_at(.vars = c(2:5), ~paste0("Core_", seq(1,4)))
per_tumor_labs[, -1] <- lapply(per_tumor_labs[, -1], function(x) as.numeric(factor(x, levels = c("Diffuse", "Edge", "Low", "Mixed"))))
agr_mat <- as.matrix(per_tumor_labs[, -1])
rownames(agr_mat) <- per_tumor_labs$ID
observed_alpha <- irr::kripp.alpha(t(agr_mat), method = "nominal")

# Assess the precision of the alpha estimate
# Step 1: Bootstrapping for confidence intervals
set.seed(5410) 
n_boot <- 10000
bootstrap_alphas <- numeric(n_boot)

for (i in 1:n_boot) {
  sample_indices <- sample(1:nrow(agr_mat), replace = TRUE)
  bootstrap_sample <- agr_mat[sample_indices, ]
  bootstrap_alpha <- kripp.alpha(t(bootstrap_sample), method = "nominal")$value
  bootstrap_alphas[i] <- bootstrap_alpha
}

# Confidence intervals
ci_lower <- quantile(bootstrap_alphas, 0.025)
ci_upper <- quantile(bootstrap_alphas, 0.975)
cat("95% Confidence Interval for Krippendorff's Alpha:", ci_lower, "-", ci_upper, "\n")

# Step 2: Permutation testing for p-value
n_perm <- 10000
permuted_alphas <- numeric(n_perm)

for (i in 1:n_perm) {
  permuted_sample <- apply(agr_mat, 2, sample, replace = FALSE)  # Shuffle the labels
  permuted_alpha <- kripp.alpha(t(permuted_sample), method = "nominal")$value
  permuted_alphas[i] <- permuted_alpha
}

# Calculate p-value (two-tailed)
Rtail_pvals <- (sum(as.numeric(permuted_alphas) >= observed_alpha$value) + 1) / (n_perm + 1)
Ltail_pvals <- (n_perm - sum(as.numeric(permuted_alphas) > observed_alpha$value) + 1) / (n_perm + 1)
tail_direction <- Rtail_pvals < Ltail_pvals
pvals <- Rtail_pvals * tail_direction + Ltail_pvals * !tail_direction
cat("P-value for Krippendorff's Alpha:", pvals, "\n")

# Plot the bootstrapped and permuted distributions
plot_df <- data.frame(Permuted = permuted_alphas,
                      Bootstrapped = bootstrap_alphas) %>% reshape2::melt(.)
p <- ggplot(plot_df, aes(x = value)) + geom_histogram(bins = 100, fill = "skyblue", color = "black") +
  annotate("segment", x = observed_alpha$value, xend = observed_alpha$value, y = 0, yend = 740, color = "red", linewidth = 1, linetype = "dashed") +
  annotate("text", label = "Permuted Distribution", x = 0, y = 1100, size = 6) +
  annotate("text", label = "Bootstrap Distribution", x = observed_alpha$value, y = 810, size = 6) +
  annotate("text", label = paste("Observed Alpha =", round(observed_alpha$value, 3)), x = observed_alpha$value + 0.09, y = 745, size = 4, color = "red") +
  annotate("text", label = paste0("P-value = ", round(pvals, 5)), x = 0.58, y = 1100, size = 5) +
  theme_scalop() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "Krippendorff's Alpha", y = "Frequency") + 
  scale_x_continuous(expand = c(0, 0), limits = c(min(plot_df$value) - 0.01, max(plot_df$value) + 0.01)) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1150))
ggplot2::ggsave(filename = here("results/Diagnostic_Plots/TMA_Per_Tumor_LAMC2_Pattern_Label_Agreement.png"), plot = p)
