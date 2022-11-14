#### note ####
# This script implements Poisson based model with adaptations, including 
# "Allowing for Uncertainty in the Pause Site" and "Allowing for Collisions"
# Phi is estimated after EM

#### log file ####
log <- file(snakemake@log[[1]], open="wt")
sink(file = log, type = "output")
sink(file = log, type = "message")

#### load packages ####
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(ggpubr)
library(scales)
library(ggpointdensity)
library(viridis)
library(cowplot)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(ggseqlogo)
library(genomation)

#### snakemake files ####
grng_in <- snakemake@input[["grng"]]
bwp1_p3_in <- snakemake@input[["bwp1_p3"]]
bwm1_p3_in <- snakemake@input[["bwm1_p3"]]

helper_in <- snakemake@params[["helper"]]
em_in <- snakemake@params[["em"]]

result_dir <- snakemake@params[["result_dir"]]

rate_tbl_out <- snakemake@output[["rate_tbl"]]

# #### testing files ####
# root_dir <- "~/Desktop/github/unimod_human"
# 
# grng_in <- file.path(root_dir, "outputs/read_dt/granges_for_read_counting.RData")
# 
# bwp1_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-dukler-control-SE_plus.bw")
# bwm1_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-dukler-control-SE_minus.bw")
# 
# # bwp1_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-vihervaara-control-SE_plus.bw")
# # bwm1_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-vihervaara-control-SE_minus.bw")
# 
# # bwp1_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-chivu-treated-PE_plus.bw")
# # bwm1_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-chivu-treated-PE_minus.bw")
# 
# helper_in <- file.path(root_dir, "scripts/unimod/helper_function.R")
# em_in <- file.path(root_dir, "scripts/unimod/helper_function_em_pause_release.R")
# 
# result_dir <-
#   file.path(root_dir, "outputs/within_sample", "PROseq-K562-dukler-control-SE", "adaptation")
# # result_dir <-
# #   file.path(root_dir, "outputs/within_sample", "PROseq-K562-vihervaara-control-SE", "adaptation")
# 
# rate_tbl_out <- file.path(result_dir, "rate.csv")
# 
# # load saved files to save some time
# # load(file.path(result_dir, "model_adaptation.RData"))

#### end of parsing arguments ####
# set up parameters
k <- 50
kmin <- 1
kmax <- 200 # also used as k on the poisson case

rnap_size <- 50
zeta <- 2000
med_init <- 1 # median initiation rate

phi_cap <- 1 # cap phi values since some of them may be larger than 1 

rc_cutoff <- 20 # read count cut-off for both gene body and pause peak

threads <- 12

# set theme for plotting
theme_set(theme_cowplot())

# create dirs
walk(
  file.path(result_dir,
            "post_Yk", c("random", "broad", "narrow",
                         "sd5", "sd10", "sd15")),
  dir.create, recursive = TRUE, showWarnings = FALSE
)

source(helper_in)
source(em_in)
# load granges for read counting
load(grng_in)

# import and process bigwigs for 3' end
bwp1_p3 <- import.bw(bwp1_p3_in)
bwm1_p3 <- import.bw(bwm1_p3_in)

bwp1_p3 <- process_bw(bw = bwp1_p3, strand = "+")
bwm1_p3 <- process_bw(bw = bwm1_p3, strand = "-")
bw1_p3 <- c(bwp1_p3, bwm1_p3)
rm(bwp1_p3, bwm1_p3)

# resize pause region same as kmax used in EM
bw_pause_filtered <-
  promoters(bw_pause_filtered, upstream = 0, downstream = kmax)

# summarize read counts
rc1_pause <- summarise_bw(bw1_p3, bw_pause_filtered, "sp1")
rc1_pause$pause_length <- kmax  

rc1_gb <- summarise_bw(bw1_p3, bw_gb_filtered, "sb1")
rc1_gb$gb_length <-
  width(bw_gb_filtered)[match(rc1_gb$gene_id, bw_gb_filtered$gene_id)]

# prepare read count table
rc1 <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE),
              list(rc1_pause, rc1_gb))

# clean up some genes with missing values in tss length or gene body length
rc1 <- rc1[!(is.na(rc1$pause_length) | is.na(rc1$gb_length)), ]
rc1 <- rc1[(rc1$sp1 > rc_cutoff) & (rc1$sb1 > rc_cutoff), ]

#### Original model: Poisson-based Maximum Likelihood Estimation ####
# Use all analyzed genes to calculate lambda
lambda_1 <- sum(rc1$sb1, na.rm = TRUE) / sum(rc1$gb_length, na.rm = TRUE)

# analytical solution based on Poisson assumption
omega_1 <- rc1$sb1 / (rc1$gb_length * lambda_1)
beta_1 <- (rc1$sb1 / rc1$gb_length) / (rc1$sp1 / rc1$pause_length)

analytical_rate_tbl <-
  tibble(gene_id = rc1$gene_id,
         omega = omega_1,
         beta = beta_1)

#### Adapted model: allow uncertainty in the pause site and steric hindrance ####
# prepare data for running EM
em_rate <-
  DataFrame(gene_id = rc1$gene_id,
            s = rc1$sb1,
            N = rc1$gb_length)

# use read count within gene body to pre-estimate chi hat
em_rate$chi = em_rate$s / em_rate$N

# get read counts on each position within pause peak (from kmin to kmax)
Xk <-
  BRGenomics::getCountsByPositions(bw1_p3, bw_pause_filtered, melt = TRUE)
Xk <- splitAsList(Xk$signal, Xk$region)
names(Xk) <- bw_pause_filtered$gene_id

em_rate$Xk <- Xk[em_rate$gene_id]

# initialize beta using sum of read counts within pause peak
em_rate$Xk_sum <- sapply(em_rate$Xk, sum)
em_rate$beta_int <- em_rate$chi / em_rate$Xk_sum 

# initialize fk with some reasonable values based on heuristic
fk_int <- dnorm(kmin:kmax, mean = 50, sd = 100)
fk_int <- fk_int / sum(fk_int)

# estimate rates using EM
em_ls <- list()
for (i in 1:NROW(em_rate)) {
  # message("Dealing with the ", i, " row...")
  rc <- em_rate[i, ]
  # message("This is gene ", rc$gene_id)
  try(em_ls[[i]] <- main_EM(fk_int = fk_int, Xk = rc$Xk[[1]], kmin = kmin, kmax = kmax,
          beta_int = rc$beta_int[[1]], chi_hat = rc$chi, max_itr = 100, tor = 1e-10),
      silent = TRUE)

  }

names(em_ls) <- em_rate$gene_id

# get rate estimates and posterior distribution of pause sites
em_rate$beta <- map_dbl(em_ls, "beta", .default = NA)
em_rate$Yk <- map(em_ls, "Yk", .default = NA)
em_rate$fk <- map(em_ls, "fk", .default = NA)
em_rate$fk_mean <- map_dbl(em_ls, "fk_mean", .default = NA)
em_rate$fk_var <- map_dbl(em_ls, "fk_var", .default = NA)
# calculate Yk / Xk
em_rate$t <- sapply(em_rate$Yk, sum)
em_rate$proportion_Yk <- em_rate$t / sapply(em_rate$Xk, sum)

# calculate phi
em_rate <- em_rate %>% as_tibble() 

## Method 1: match alpha with median initiation rate
omega_scale <- med_init / median(em_rate$chi)

em_rate_matched <- em_rate %>%
  mutate(omega_zeta = chi * omega_scale,
         beta_zeta = beta * zeta,
         phi = omega_zeta / beta_zeta + omega_zeta / zeta * (k - 1),
         phi_cap = if_else(phi > .env$phi_cap, .env$phi_cap, phi),
         alpha = omega_zeta / (1 - phi_cap)) 

## Method 2: Use t to normalize alpha
# quantile of t to explore its effect 
t_quantile <- quantile(em_rate$t, probs = c(0.8, 0.85, 0.9, 0.95, 0.99, 1))

calculate_phi <- function(em_rate, t_val) {
  em_rate <- em_rate %>% 
    mutate(omega = chi / t_val,
           phi = omega / beta + omega * (k - 1),
           # cap some edge cases with phi slightly higher than 1
           # https://rlang.r-lib.org/reference/topic-data-mask-ambiguity.html#:~:text=Data%20masking%20is%20an%20R,defined%20in%20the%20current%20environment.
           phi_cap = if_else(phi > .env$phi_cap, .env$phi_cap, phi),
           alpha = omega / (1 - phi_cap))
  return(em_rate)
}

em_rate_t_max <- calculate_phi(em_rate, t_val = t_quantile[["100%"]])
em_rate_t_99 <- calculate_phi(em_rate, t_val = t_quantile[["99%"]])
em_rate_t_95 <- calculate_phi(em_rate, t_val = t_quantile[["95%"]])
em_rate_t_90 <- calculate_phi(em_rate, t_val = t_quantile[["90%"]])
em_rate_t_85 <- calculate_phi(em_rate, t_val = t_quantile[["85%"]])
em_rate_t_80 <- calculate_phi(em_rate, t_val = t_quantile[["80%"]])

#### visualize results ####
#### rate estimate comparisons ####
# scatter plot with margin
# https://stackoverflow.com/questions/8545035/scatterplot-with-marginal-histograms-in-ggplot2
rate_tbl <- analytical_rate_tbl %>%
  left_join(as_tibble(em_rate), by = "gene_id", suffix = c("_org", "_adp"))

rate_tbl <- rate_tbl %>%
  mutate(beta_zeta_org = beta_org * zeta,
         beta_zeta_adp = beta_adp * zeta)

# beta estimates, original model vs. adapted model
p <- rate_tbl %>%
  ggscatterhist(x = "beta_zeta_org", y = "beta_zeta_adp",
                cor.coef = FALSE,
                alpha = 0.1,
                xlab =  expression("Estimated "*beta*zeta*" (averaging)"),
                ylab = expression("Estimated "*beta*zeta*" (model)")
  )

sp <- p$sp +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(10^-0.5, 10^5)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(10^-2.5, 10^3)) +
  geom_abline(slope = 1, linetype = "dashed", color = "gray", size = 1) +
  theme(strip.text.x = element_blank(),
        # axis.text.x = element_blank(),
        # legend.position = c(0.83, 0.8),
        text = element_text(size = 25),
        axis.text = element_text(size = 25))

p <- rate_tbl %>%
  mutate(across(where(is.numeric), ~log10(.))) %>%
  ggscatterhist(x = "beta_zeta_org", y = "beta_zeta_adp",
                cor.coef = FALSE,
                cor.coef.coord = c(-12, 0),
                cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
                alpha = 0.1,
                # add = "reg.line",  # Add regression line
                # add.params = list(color = "blue", linetype = "dashed"), # Customize reg. line
                # conf.int = FALSE, # Add confidence interval
                xlim = c(-0.5, 5), ylim = c(-2.5, 3),
                xlab =  expression("Estimated "*beta*zeta*" (averaging)"),
                ylab = expression("Estimated "*beta*zeta*" (model)")
  )

p$sp <- sp
# convert the plot list into a single plot
p <- print(p)
ggsave(filename = file.path(result_dir, "analytical_vs_em_beta_scatter_with_margin.png"), plot = p,
       width = 8, height = 8)

sp <- sp + cowplot::theme_cowplot() +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15))
ggsave(filename = file.path(result_dir, "analytical_vs_em_beta_scatter.png"), plot = sp,
       width = 4, height = 4)

# alpha, omega and beta
em_rate_ls <-
  list(em_rate_t_max, em_rate_t_99, em_rate_t_95,
       em_rate_t_90, em_rate_t_85, em_rate_t_80)

scatter_plot <- function(tbl) {
  tbl %>% ggplot(aes(x = beta, y = omega)) +
  geom_pointdensity() +
  scale_color_viridis() +
  geom_abline(slope = 1, intercept = log10(c(0.1, 0.2, 0.5, 0.8, 0.9, 1)),
              color = RColorBrewer::brewer.pal(6, "Reds"), linetype = "dashed") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  cowplot::theme_cowplot() +
  labs(x = expression("Estimated "*beta), y = expression("Estimated "*omega))
  }

beta_omega_plt_ls <- map(em_rate_ls, scatter_plot)

p <-
  cowplot::plot_grid(plotlist = beta_omega_plt_ls, align = "hv", nrow = 3,
                     labels = paste0(rev(names(t_quantile)), " percentile of t"))

ggsave(filename = file.path(result_dir, "omega_vs_beta.png"), plot = p,
       width = 16, height = 18)

# visualize distribution of phi
histogram <- function(tbl, sel_col) {
  tbl %>%
    ggplot(aes(x = {{sel_col}})) +
    geom_histogram(aes(y = (..count..)/sum(..count..)), bins = 50) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = expression("Estimated "*phi), y = "density") +
    cowplot::theme_cowplot()
}

phi_hist_ls <- map(em_rate_ls, histogram, sel_col = phi_cap)
p <- cowplot::plot_grid(plotlist = phi_hist_ls, nrow = 3,
                        labels = paste0(rev(names(t_quantile)), " percentile of t"))

ggsave(filename = file.path(result_dir, "phi_cap.png"), plot = p,
       width = 10, height = 9)

phi_hist_ls <- map(em_rate_ls, histogram, sel_col = phi)
p <- cowplot::plot_grid(plotlist = phi_hist_ls, nrow = 3,
                        labels = paste0(rev(names(t_quantile)), " percentile of t"))

ggsave(filename = file.path(result_dir, "phi.png"), plot = p,
       width = 10, height = 9)

p <- histogram(em_rate_matched, phi)
ggsave(filename = file.path(result_dir, "phi_matched_initiation.png"), plot = p,
       width = 6, height = 3)

p <- histogram(em_rate_matched, phi_cap)
ggsave(filename = file.path(result_dir, "phi_cap_matched_initiation.png"), plot = p,
       width = 6, height = 3)

# # try a subset of genes
# long_gn <- bw_gb_filtered[width(bw_gb_filtered) > 1e4, ]$gene_id
# 
# phi_hist_ls <- map(em_rate_ls, function(x) {
#   x %>%
#     filter(gene_id %in% long_gn) %>%
#     slice_min(order_by = beta, prop = 0.3) %>%
#     histogram(sel_col = phi_cap)
# })
# 
# p <- cowplot::plot_grid(plotlist = phi_hist_ls, nrow = 3,
#                         labels = paste0(rev(names(t_quantile)), " percentile of t"))
# 
# ggsave(filename = file.path(result_dir, "phi_cap_filtered.png"), plot = p,
#        width = 10, height = 9)
# 
# map(em_rate_ls, function(x) {
#   x %>%
#     filter(gene_id %in% long_gn) %>%
#     slice_min(order_by = beta, prop = 0.3) %>%
#     pull(phi_cap) %>%
#     summary()
# }) 

# p <- em_rate %>%
#   ggplot(aes(x = t)) +
#   geom_histogram(bins = 50) +
#   geom_vline(xintercept = max(em_rate$t)) + 
#   labs(x = "Estimated t") +
#   cowplot::theme_cowplot()
# 
# ggsave(filename = file.path(result_dir, "t.png"), plot = p,
#        width = 10, height = 3)

# p <- rate_tbl %>% ggplot(aes(x = s, y = t)) +
#   geom_pointdensity() +
#   scale_color_viridis() +
#   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                 labels = trans_format("log10", math_format(10^.x))) +
#   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                 labels = trans_format("log10", math_format(10^.x))) +
#   cowplot::theme_cowplot() +
#   labs(x = "Sum of reads count in gene body", y = "Sum of Yk")
# 
# ggsave(filename = file.path(result_dir, "s_vs_t.png"), plot = p,
#        width = 8, height = 6)
# 
# p <- rate_tbl %>% ggplot(aes(x = chi, y = t)) +
#   geom_pointdensity() +
#   scale_color_viridis() +
#   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                 labels = trans_format("log10", math_format(10^.x))) +
#   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                 labels = trans_format("log10", math_format(10^.x))) +
#   cowplot::theme_cowplot() +
#   labs(x = expression(chi), y = "Sum of Yk")
# 
# ggsave(filename = file.path(result_dir, "chi_vs_t.png"), plot = p,
#        width = 8, height = 6)
# 
# p <- rate_tbl %>% ggplot(aes(x = alpha_adp, y = 1 / beta_adp)) +
#   geom_point(alpha = 0.1) +
#   cowplot::theme_cowplot() +
#   coord_cartesian(xlim = c(0, 1e-4), ylim = c(0, 1e5)) +
#   labs(x = expression(omega), y = expression(1/beta))
# 
# ggsave(filename = file.path(result_dir, "alpha_vs_1_over_beta.png"), plot = p,
#        width = 8, height = 6)

# visualize pause duration (half lives)
em_rate <- em_rate %>%
  mutate(halflife = log(2) * (k - 1 + 1 / beta) / zeta)

hl_mean <- mean(em_rate$halflife)
hl_med <- median(em_rate$halflife)

p <- em_rate %>%
  ggplot(aes(x = halflife)) +
  geom_histogram(aes(y = (..count..)/sum(..count..)), binwidth = 0.5) +
  geom_vline(xintercept = hl_mean, linetype = "dashed", color = "red", lwd = 1) +
  geom_vline(xintercept = hl_med, linetype = "dashed", color = "blue", lwd = 1) +
  # https://stackoverflow.com/questions/47916307/specify-position-of-geom-text-by-keywords-like-top-bottom-left-right
  annotate(geom = 'text', label = paste("mean =", round(hl_mean, 2)), color = "red",
           x = 7, y = Inf, vjust = 3) +
  annotate(geom = 'text', label = paste("median =", round(hl_med, 2)), color = "blue",
           x = 7, y = Inf, vjust = 5) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(xlim = c(0, 10)) +
  labs(x = "Paused Pol II half-life (min)", y = "Proportion") 

ggsave(filename = file.path(result_dir, "pause_halflife.png"), plot = p,
       width = 6, height = 3)
ggsave(filename = file.path(result_dir, "pause_halflife.pdf"), plot = p,
       width = 6, height = 3)

# visualize Yk / Xk
p <- rate_tbl %>%
  ggplot(aes(x = proportion_Yk)) +
  geom_histogram(bins = 50) +
  labs(x = "Yk / Xk") +
  cowplot::theme_cowplot()

ggsave(filename = file.path(result_dir, "yk_proportion.png"), plot = p,
       width = 6, height = 3)

# # Quintile of T
# high_t_gene <- rate_tbl %>% slice_max(order_by = t, n = 20) %>% pull(gene_id)
# 
# p <- rate_tbl %>%
#   mutate(group = ifelse(gene_id %in% high_t_gene, T, F)) %>%
#   ggplot(aes(x = group, y = log2(alpha_adp))) +
#   geom_boxplot() +
#   labs(x = "high t", y = expression(log[2]*alpha)) +
#   cowplot::theme_cowplot()
# 
# ggsave(filename = file.path(result_dir, "alpha_for_high_t_gene.png"), plot = p,
#        width = 5, height = 5)
# 
# p <- rate_tbl %>%
#   mutate(group = ifelse(gene_id %in% high_t_gene, T, F)) %>%
#   ggplot(aes(x = group, y = log2(beta_adp))) +
#   geom_boxplot() +
#   labs(x = "high t", y = expression(log[2]*beta)) +
#   cowplot::theme_cowplot()
# 
# ggsave(filename = file.path(result_dir, "beta_for_high_t_gene.png"), plot = p,
#        width = 5, height = 5)

# rate_tbl %>%
#   mutate(t_quantile = factor(ntile(t, 5))) %>%
#   ggplot(aes(x = t_quantile, y = log2(alpha_adp))) +
#   geom_boxplot() +
#   labs(x = "Quintile of t", y = expression(log[2]*alpha)) +
#   cowplot::theme_cowplot()
# 
# rate_tbl %>%
#   mutate(t_quantile = factor(ntile(t, 5))) %>%
#   ggplot(aes(x = t_quantile, y = log2(beta_adp))) +
#   geom_boxplot() +
#   labs(x = "Quintile of t", y = expression(log[2]*beta)) +
#   cowplot::theme_cowplot()
# 
# rate_tbl %>%
#   mutate(t_quantile = factor(ntile(t, 5))) %>%
#   ggplot(aes(x = t_quantile, y = t)) +
#   geom_boxplot() +
#   labs(x = "Quintile of t", y = "t") +
#   cowplot::theme_cowplot()

# visualize Xk, Yk and fk for some genes
p <- em_rate %>%
  mutate(fk_std = fk_var ^ 0.5) %>% 
  ggplot(aes(x = fk_mean, y = fk_std)) +
  geom_pointdensity() +
  scale_color_viridis() +
  geom_vline(xintercept  = 50, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = 25, linetype = "dashed", color = "gray") +
  labs(x = "Mean of fk", y = "SD of fk")

ggsave(filename = file.path(result_dir, "fk_mean_vs_std.png"), plot = p,
       width = 8, height = 5)

p <- em_rate %>% 
  ggplot(aes(x = fk_mean)) +
  geom_histogram() +
  labs(x = "Mean of fk", y = "Number of genes")

ggsave(filename = file.path(result_dir, "fk_mean_histogram.png"), plot = p,
       width = 8, height = 5)

fk_sd_med <- median(em_rate$fk_var ^ 0.5)

p <- em_rate %>% 
  ggplot(aes(x = fk_var ^ 0.5)) +
  geom_histogram() +
  geom_vline(xintercept = fk_sd_med, linetype = "dashed", color = "blue", lwd = 1) +
  labs(x = "SD of fk", y = "Number of genes")

ggsave(filename = file.path(result_dir, "fk_sd_histogram.png"), plot = p,
       width = 8, height = 5)

# plot pause peak
em_rate_plt_df <- em_rate[!is.na(em_rate$beta), ]
# set seed to ensure results are reproducible
set.seed(1234)
em_rate_plt_df <- em_rate_plt_df[sample(nrow(em_rate_plt_df),), ]

plot_pause_peak <- function(em_rate_plt_df, i, path, format = ".png") {
  df <- em_rate_plt_df[i, c("gene_id", "Xk", "Yk", "fk")]
  file_name <- file.path(result_dir, "post_Yk", path, paste0(df$gene_id, format))
  
  df <- data.frame(Xk = df$Xk[[1]], Yk = df$Yk[[1]], fk = df$fk[[1]])
  scale_yaxis <- sum(df$Xk)
  # plot first 100bp for clarity
  kmax_p <- min(100, kmax)
  
  p <- df[kmin:kmax_p,] %>% 
    mutate(Xk = Xk / scale_yaxis) %>%  ggplot() +
    geom_col(aes(x = seq(kmin:kmax_p), y = Xk)) +
    geom_line(aes(x = seq(kmin:kmax_p), y = fk), color = "red") +
    scale_y_continuous(
      "Xk",
      breaks = seq(0, 500, 5) / scale_yaxis,
      labels = ~ . * scale_yaxis,
      sec.axis = sec_axis(~ ., name = "fk")
    ) +
    labs(x = "", y = "Xk") +
    cowplot::theme_cowplot()
  
  # p1 <- df %>% ggplot() +
  #   geom_col(aes(x = seq(kmin:kmax), y = Xk)) +
  #   labs(x = "", y = "Xk") +
  #   cowplot::theme_cowplot()
  # p2 <- df %>% ggplot() +
  #   geom_col(aes(x = seq(kmin:kmax), y = Yk)) +
  #   labs(x = "", y = "Yk") +
  #   cowplot::theme_cowplot()
  # p3 <- df %>% ggplot() +
  #   geom_col(aes(x = seq(kmin:kmax), y = fk)) +
  #   labs(x = "", y = "fk") +
  #   cowplot::theme_cowplot()
  # p <- cowplot::plot_grid(p1, p2, p3, nrow = 3, align = "v")
  
  ggsave(file_name, plot = p, width = 6, height = 2)
}

 if (NROW(em_rate_plt_df) > 200){
   # pick 200 genes randomly and see model's behavior
   for (i in 1:200) {
     plot_pause_peak(em_rate_plt_df, i, "random")
   }
   
   # sd5_tbl <- em_rate_plt_df %>% filter(fk_var < 7)
   # sd10_tbl <- em_rate_plt_df %>% filter(fk_var > 8, fk_var < 12)
   # sd15_tbl <- em_rate_plt_df %>% filter(fk_var > 13, fk_var < 17)
   # 
   # for (i in 1:NROW(sd5_tbl)) {
   #   plot_pause_peak(sd5_tbl, i, "sd5")
   # }
   # 
   # for (i in 1:NROW(sd10_tbl)) {
   #   plot_pause_peak(sd10_tbl, i, "sd10")
   # }
   # 
   # for (i in 1:NROW(sd15_tbl)) {
   #   plot_pause_peak(sd15_tbl, i, "sd15")
   # }

   # plot_pause_peak(em_rate, which(em_rate$gene_id == "ENSG00000086061"), "random")
   
   # examine some genes with broad or narrow peaks
   em_rate_high_Xk <- em_rate_plt_df %>%
     filter(Xk_sum >= 50)
   
   em_rate_broad_peak <-
     em_rate_high_Xk %>% slice_max(order_by = fk_var, n = 100)
   em_rate_narrow_peak <-
     em_rate_high_Xk %>% slice_min(order_by = fk_var, n = 100)
   
   for (i in 1:100) {
     plot_pause_peak(em_rate_broad_peak, i, "broad")
   }
   
   for (i in 1:100) {
     plot_pause_peak(em_rate_narrow_peak, i, "narrow")
   }
 }

# sequence logo
pause_max <- map_dbl(em_rate$Xk, which.max)
pause_model <- round(em_rate$fk_mean)

get_seqlogo <- function(em_rate, sites, len = 21) {
  pause_grng <- bw_pause_filtered[match(em_rate$gene_id, bw_pause_filtered$gene_id), ]
  pause_grng <- promoters(pause_grng, upstream = 0, downstream = sites)
  
  pause_grng <- as.data.frame(pause_grng) %>%
    mutate(start = ifelse(strand == "+", end, start),
           end = ifelse(strand == "-", start, end)) %>%
    select(-width) %>%
    GRanges() %>% 
    resize(width = len, fix = "center")
  
  seq <- getSeq(BSgenome.Hsapiens.NCBI.GRCh38, pause_grng)
  seq_vec <- as.character(seq)
  
  p <- ggseqlogo(seq_vec) +
    scale_x_continuous(breaks=1:len, labels=c(1:len - 11))
  return(list("seq" = seq, "p" = p))
} 

p <- get_seqlogo(em_rate, pause_max)
ggsave(file.path(result_dir, "seqlogo_max_pause_all_genes.png"), plot = p$p, width = 5, height = 2)

p <- get_seqlogo(em_rate, pause_model)
ggsave(file.path(result_dir, "seqlogo_mean_k_all_genes.png"), plot = p$p, width = 5, height = 2)

fk_var_head <- em_rate %>% arrange(fk_var) %>% slice_head(prop = 0.1)
fk_var_tail <- em_rate %>% arrange(fk_var) %>% slice_tail(prop = 0.1)

p <- get_seqlogo(fk_var_head, map_dbl(fk_var_head$Xk, which.max))
ggsave(file.path(result_dir, "seqlogo_small_fk_var_max_pause.png"), plot = p$p, width = 5, height = 2)

p <- get_seqlogo(fk_var_tail, map_dbl(fk_var_tail$Xk, which.max))
ggsave(file.path(result_dir, "seqlogo_large_fk_var_max_pause.png"), plot = p$p, width = 5, height = 2)

p <- get_seqlogo(fk_var_head, map_dbl(fk_var_head$Xk, ~ order(.x, decreasing = T)[2]))
ggsave(file.path(result_dir, "seqlogo_small_fk_var_second_max.png"), plot = p$p, width = 5, height = 2)
p <- get_seqlogo(fk_var_tail, map_dbl(fk_var_tail$Xk, ~ order(.x, decreasing = T)[2]))
ggsave(file.path(result_dir, "seqlogo_large_fk_var_second_max.png"), plot = p$p, width = 5, height = 2)

# p <- get_seqlogo(fk_var_head, map_dbl(fk_var_head$Xk, ~ order(.x, decreasing = T)[3]))
# ggsave(file.path(result_dir, "seqlogo_small_fk_var_third_max.png"), plot = p$p, width = 5, height = 2)
# p <- get_seqlogo(fk_var_tail, map_dbl(fk_var_tail$Xk, ~ order(.x, decreasing = T)[3]))
# ggsave(file.path(result_dir, "seqlogo_large_fk_var_third_max.png"), plot = p$p, width = 5, height = 2)

p <- get_seqlogo(fk_var_head, map_dbl(fk_var_head$Xk, which.min))
ggsave(file.path(result_dir, "seqlogo_small_fk_var_min_pause.png"), plot = p$p, width = 5, height = 2)
p <- get_seqlogo(fk_var_tail, map_dbl(fk_var_tail$Xk, which.min))
ggsave(file.path(result_dir, "seqlogo_large_fk_var_min_pause.png"), plot = p$p, width = 5, height = 2)

# extract local GC content
p1 <- get_seqlogo(fk_var_head, map_dbl(fk_var_head$Xk, which.max), len = 101)
p2 <- get_seqlogo(fk_var_tail, map_dbl(fk_var_tail$Xk, which.max), len = 101)

p <- tibble(
  position = -50:50,
  narrow = colSums(consensusMatrix(p1$seq, as.prob = TRUE)[c("C", "G"), ]),
  broad = colSums(consensusMatrix(p2$seq, as.prob = TRUE)[c("C", "G"), ])
) %>% pivot_longer(cols = c(narrow, broad)) %>%
  ggplot(aes(x = position, y = value, color = name)) +
  geom_line() +
  labs(y = "GC content") +
  theme_cowplot()

ggsave(file.path(result_dir, "local_GC_content_around_pause_site.png"), plot = p, width = 8, height = 3)

# TAF1 Chip-seq
read_chipseq <- function(bw_in) {
  grng <- import.bw(bw_in)
  keepStandardChromosomes(grng, pruning.mode="coarse")
  seqlevelsStyle(grng) <- "Ensembl"
  return(grng)
}

# metaplot for TAF1 Chip-seq downstream of TSS
# use TAF1 Chip-seq from ENCODE
# https://www.encodeproject.org/experiments/ENCSR000BKS/
taf1_grng <- read_chipseq("ext_data/encode/chipseq/ENCFF101GBL.bigWig")

narrow_sm <- ScoreMatrix(taf1_grng, bw_pause_filtered[fk_var_head$gene_id, ], weight.col="score")
broad_sm <- ScoreMatrix(taf1_grng, bw_pause_filtered[fk_var_tail$gene_id, ], weight.col="score")

taf_sm_ls <- as(list(narrow_sm, broad_sm), "ScoreMatrixList")
taf_sm_ls@names <- c("Narrow peak", "Broad peak")

set_color <- RColorBrewer::brewer.pal(name = "Set1", n = max(length(taf_sm_ls), 3))

# add alpha to the colors
# https://gotellilab.github.io/GotelliLabMeetingHacks/NickGotelli/ColorPalettes.html
addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

meta_plot <- function(sm, smlcolors, meta.rescale = FALSE, xcoords = c(1, 200),
                      dispersion = "se", ylab = "Average score", legend = "topright") {
  plotMeta(sm, xcoords = xcoords, meta.rescale = meta.rescale,
           line.col = smlcolors, dispersion = dispersion,
           xlab = "Bases downstream of TSS", ylab = ylab,
           dispersion.col = addalpha(smlcolors, alpha = 0.5))
  legend(legend, names(sm), lty=c(1,1), lwd=c(2.5,2.5), col=smlcolors, cex = 1.5)
}

save_png <- function(file_name, plot_fun, width = 800, height = 500) {
  png(filename = file_name, width = width, height = height, units = "px")
  plot_fun
  invisible(dev.off())
}

save_png(file_name = paste0(file.path(result_dir, "taf1_chipseq_for_narrow_and_broad_peaks"), ".png"),
         plot_fun = meta_plot(taf_sm_ls, set_color))

# write results
em_rate_matched %>% 
  select(-c(Xk, Yk, fk)) %>%
  write_csv(file = rate_tbl_out)

# x <- as.data.frame(em_rate_matched)
# saveRDS(x, file.path(result_dir, "rate.RDS"))
# readRDS(file.path(result_dir, "rate.RDS"))

# save image for easier access
save.image(file = file.path(result_dir, "model_adaptation.RData"))
