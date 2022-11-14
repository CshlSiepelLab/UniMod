#### note ####
# This script implements Poisson based model with adaptations, including 
# "Allowing for Uncertainty in the Pause Site" and "Allowing for Collisions"
# Phi estimated is included in EM

#### log file ####
log <- file(snakemake@log[[1]], open="wt")
sink(file = log, type = "output")
sink(file = log, type = "message")

#### load packages ####
library(tidyverse)
library(rtracklayer)
library(ggpubr)
library(scales)
library(ggpointdensity)
library(viridis)
library(cowplot)
# library(BSgenome.Hsapiens.NCBI.GRCh38)
# library(ggseqlogo)
# library(genomation)

#### snakemake files ####
grng_in <- snakemake@input[["grng"]]
bwp1_p3_in <- snakemake@input[["bwp1_p3"]]
bwm1_p3_in <- snakemake@input[["bwm1_p3"]]
scale_in <- snakemake@input[["scale"]]

helper_in <- snakemake@params[["helper"]]
em_in <- snakemake@params[["em"]]

sample_id <- snakemake@params[["sample_id"]]
result_dir <- snakemake@params[["result_dir"]]

rate_tbl_out <- snakemake@output[["rate_tbl"]]

# #### testing files ####
# root_dir <- "~/Desktop/github/unimod_human"
# 
# grng_in <- file.path(root_dir, "outputs/read_dt/granges_for_read_counting.RData")
# 
# # bwp1_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-dukler-control-SE_plus.bw")
# # bwm1_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-dukler-control-SE_minus.bw")
# 
# bwp1_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-vihervaara-control-SE_plus.bw")
# bwm1_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-vihervaara-control-SE_minus.bw")
# 
# # bwp1_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-chivu-treated-PE_plus.bw")
# # bwm1_p3_in <- file.path(root_dir, "outputs/bigwig/p3/PROseq-K562-chivu-treated-PE_minus.bw")
# 
# sample_id <- "PROseq-K562-vihervaara-control-SE"
# 
# em_in <- file.path(root_dir, "scripts/unimod/helper_function_em_steric_hindrance.R")
# helper_in <- file.path(root_dir, "scripts/unimod/helper_function.R")
# 
# scale_in <- file.path(root_dir, "outputs/between_samples/table/scale_factor.csv")
# 
# # result_dir <-
# #   file.path(root_dir, "outputs/within_sample", "PROseq-K562-dukler-control-SE", "pause_release")
# result_dir <-
#   file.path(root_dir, "outputs/within_sample", "PROseq-K562-vihervaara-control-SE", "steric_hindrance")
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

phi_cap <- 1 # cap phi values since some of them may be larger than 1 

rc_cutoff <- 20 # read count cut-off for both gene body and pause peak

threads <- 12

# set theme for plotting
theme_set(theme_cowplot())

# create dirs
walk(
  file.path(result_dir),
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
Xk <- BRGenomics::getCountsByPositions(bw1_p3, bw_pause_filtered, melt = TRUE)
Xk <- splitAsList(Xk$signal, Xk$region)
names(Xk) <- bw_pause_filtered$gene_id

em_rate$Xk <- Xk[em_rate$gene_id]

# initialize beta using sum of read counts within pause peak
em_rate$Xk_sum <- sapply(em_rate$Xk, sum)
em_rate$beta_int <- em_rate$chi / em_rate$Xk_sum 

# initialize fk with some reasonable values based on heuristic
fk_int <- dnorm(kmin:kmax, mean = 50, sd = 100)
fk_int <- fk_int / sum(fk_int)

run_EM <- function(em_rate, omega_scale) {
  em_rate$omega_zeta <- em_rate$chi * omega_scale
  em_rate$omega <- em_rate$omega_zeta / zeta
  
  # compute a scaling factor lambda for the purpose of using the same parameters as simulations
  lambda <- zeta ^ 2 / omega_scale 
  
  # estimate rates using EM
  em_ls <- list()
  
  main_EM <- possibly(main_EM, otherwise = NA)
  
  for (i in 1:NROW(em_rate)) {
    # message("Dealing with the ", i, " row...")
    rc <- em_rate[i, ]
    
    em_ls[[i]] <- main_EM(Xk = rc$Xk[[1]], kmin = kmin, kmax = kmax, f1 = 0.517, f2 = 0.024,
                          fk_int = fk_int, beta_int = rc$beta_int[[1]], phi_int = 0.5,
                          chi_hat = rc$chi, lambda = lambda, zeta = zeta,
                          max_itr = 500, tor = 1e-3)
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
  em_rate$phi <- map_dbl(em_ls, "phi", .default = NA)
  
  # calculate alpha
  em_rate <- em_rate %>% as_tibble() 
  
  em_rate <- em_rate %>%
    mutate(beta_zeta = beta * zeta,
           alpha_zeta = omega_zeta / (1 - phi)) 
  
  return(em_rate)
}

scale_tbl <- read_csv(scale_in, show_col_types = FALSE)
scale_tbl <- scale_tbl %>% filter(sample_id == {{sample_id}})

# scale chi to get omega by assuming median omega * zeta is 0.2 event / min.
em_rate_os1 <- run_EM(em_rate, scale_tbl$omega_scale_l)
# show top 75% genes by chi
em_rate_os1_plot <- em_rate_os1 %>% slice_max(order_by = chi, prop = 0.75)

# scale chi to get omega by assuming house keeping genes have same median omega
em_rate_os2 <- run_EM(em_rate, scale_tbl$omega_scale_h)
# use the same set of gene for scaling 1
em_rate_os2_plot <- em_rate_os2 %>% filter(gene_id %in% em_rate_os1_plot$gene_id)

#### visualize results ####
# visualize distribution of phi
histogram <- function(tbl, sel_col, xlab, bin_n = 20) {
  
  x <- summary(tbl$phi)
  
  # https://stackoverflow.com/questions/32123288/position-ggplot-text-in-each-corner
  annotations <- data.frame(
    xpos = rep(0.4, 3), ypos =  c(Inf, Inf, Inf),
    annotateText = c(paste0("Median: ", round(x[["Median"]], digits = 2)),
                     paste0("Mean: ", round(x[["Mean"]], digits = 2)),
                     paste0("N: ", NROW(tbl))),
    hjustvar = rep(0, 3),
    vjustvar = seq(4, 8, 2))
  
  tbl %>%
    ggplot(aes(x = {{sel_col}})) +
    # geom_histogram(aes(y = (..count..)/sum(..count..)), bins = bin_n) +
    geom_histogram(bins = bin_n, color="gray", fill="lightblue") + 
    geom_text(data = annotations,
              aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar,
                  label = annotateText), size = 5) +
    # scale_y_continuous(labels = scales::percent) +
    labs(x = xlab, y = "Frequency")
}

# phi estimates
p1 <- histogram(em_rate_os1_plot, phi, expression("Estimated "*phi))
ggsave(file.path(result_dir, "phi_estimates.png"), width = 4, height = 3)

# histogram(em_rate_os1_plot %>% filter(phi < 0.999), phi, expression("Estimated "*phi))

x <- summary(em_rate_os1_plot$omega_zeta)
y <- summary(em_rate_os1_plot$alpha_zeta)

annotation_x <- data.frame(
  xpos = c(0.4, 0.4), ypos =  c(0.15, 0.135),
  annotateText = c(paste0(""),
                   paste0("Median: ", round(x[["Median"]], digits = 2))
                   # paste0("Mean: ", round(x[["Mean"]], digits = 2)),
                   ),
  hjustvar = c(0, 0),
  vjustvar = c(0, 0))

annotation_y <- data.frame(
  xpos = c(0.8, 0.8), ypos =  c(0.045, 0.03),
  annotateText = c(paste0(""),
                   paste0("Median: ", round(y[["Median"]], digits = 2))
                   # paste0("Mean: ", round(y[["Mean"]], digits = 2)),
                   ),
  hjustvar = c(0, 0),
  vjustvar = c(0, 0))

x_lab <- seq(0, 2, 0.5)
names(x_lab) <- x_lab  
x_lab[length(x_lab)] <- "≥2.0"

p2 <- em_rate_os1_plot %>%
  select(omega_zeta, alpha_zeta) %>% 
  pivot_longer(cols = c(omega_zeta, alpha_zeta)) %>%
  mutate(value = ifelse(value > 2, 2, value)) %>% 
  ggplot(aes(x = value)) +
  # geom_histogram(aes(y =), bins = bin_n) +
  geom_histogram(aes(y =  (..count..)/sum(..count..), fill = name), position="identity",
                 color = "gray", binwidth = 0.1, alpha=0.3) +
  geom_text(data = annotation_x,
            aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar,
                label = annotateText), size = 5) +
  geom_text(data = annotation_y,
            aes(x = xpos, y = ypos, hjust = hjustvar, vjust = vjustvar,
                label = annotateText), size = 5) +
  geom_text(aes(x = 1.2, y = 0.095, label = "Steric Hindrance"), fontface = 'italic', size = 4) +
  geom_segment(aes(x = 1.5, y = 0.08, xend = 0.5, yend = 0.08),
               lineend = "round", linejoin = "round", color = "gray30",
               arrow = arrow(length = unit(0.2, "cm")), size = 1.5) +
  geom_vline(xintercept = 1.95, linetype = "dashed") +
  coord_cartesian(xlim = c(0, 2)) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks=seq(0, 2, 0.5), labels = x_lab) +
  scale_fill_discrete(labels=c(expression(alpha*zeta), expression(omega*zeta))) +
  labs(x = "Initiation Rate (events/min.)", y = "Percentage", fill = "")

ggsave(file.path(result_dir, "alpha_vs_omega_estimates.png"), width = 4, height = 3)

# full occupancy proportion 
proportion_df <-
  tibble(group = factor(c("L", "H"), levels = c("L", "H")),
         value = c(mean(em_rate_os1_plot$phi == 0.999),
                   mean(em_rate_os2_plot$phi == 0.999)))

p3 <- proportion_df %>%
  ggplot(aes(x = group, y = value)) +
  geom_col(fill = "darksalmon", color = "gray") +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Percentage") +
  theme(axis.title.x = element_blank())

ggsave(file.path(result_dir, "phi_estimates_by_different_scaling.png"), width = 2, height = 3)

p4 <- em_rate_os1 %>%
  # mutate(across(.cols = c(omega_zeta, beta_zeta), log2),
  #        group = ifelse(phi < 0.999, "others", "0.999")) %>%
  # filter(phi != 0.999) %>% 
  ggplot(aes(omega_zeta, beta_zeta)) +
  geom_pointdensity() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(10^-3.5, 10^1)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(10^-3.5, 10^3)) +
  scale_color_viridis() +
  geom_abline(slope = 1, intercept = 0, color = "gray", lwd = 1, linetype = "dashed") +
  stat_cor(method = "pearson") +
  labs(x = expression(omega*zeta), y = expression(beta*zeta), color = "# of neighbors")

ggsave(file.path(result_dir, "omega_vs_beta_estimates.png"), width = 5, height = 4)

p5 <- em_rate_os1 %>%
  # mutate(across(.cols = c(alpha_zeta, beta_zeta), log2),
  #        group = ifelse(phi < 0.999, "others", "0.999")) %>%
  filter(phi != 0.999) %>% 
  ggplot(aes(alpha_zeta, beta_zeta)) +
  geom_pointdensity() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(10^-3.5, 10^1)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(10^-3.5, 10^3)) +
  scale_color_viridis() +
  geom_abline(slope = 1, intercept = 0, color = "gray", lwd = 1, linetype = "dashed") +
  stat_cor(method = "pearson", label.x.npc = "left") +
  labs(x = expression(alpha*zeta), y = expression(beta*zeta), color = "# of neighbors")

ggsave(file.path(result_dir, "alpha_vs_beta_estimates.png"), width = 5, height = 4)

g1 <- plot_grid(p1, p2, nrow = 1, rel_widths = c(2, 2.5), align = "hv", axis = "l")
ggsave(file.path(result_dir, "steric_hindrance_2.png"), width = 9, height = 4, plot = g1,
       device = png)

# g2 <- plot_grid(p4, p5, nrow = 1, align = "hv", axis = "l")
# p <- plot_grid(g1, g2, align = "hv", axis = "tblr", nrow = 2)

# write results
em_rate_os1 %>% 
  select(-c(Xk, Yk, fk)) %>%
  write_csv(file = rate_tbl_out)

em_rate_os2 %>% 
  select(-c(Xk, Yk, fk)) %>%
  write_csv(file = file.path(result_dir, "rate_matched_gressel.csv"))

# proportion of genes with fully occupied landing pad
proportion_df %>%
  write_csv(file = file.path(result_dir, "proportion_of_landing_pad_with_full_occupancy.csv"))

# save image for easier access
save.image(file = file.path(result_dir, "model_adaptation.RData"))
