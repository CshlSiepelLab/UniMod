#!/usr/bin/env Rscript

#### load optparse ####
library(optparse)

#### parse command line arguments ####
parser <-
  OptionParser(prog = "./estimate_rates_simulation.R",
               description = "Estimate transcription rates based on simulated data")

parser <- add_option(parser, c("-r", "--rds"), action="store", type="character",
                     default=NULL, help="Input file produced by SimPol [default %default]",
                     metavar="character")

parser <- add_option(parser, c("-s", "--steric"), action="store", type="logical",
                     default=FALSE, help="Infer landing-pad occupancy or not [default %default]",
                     metavar="logical")

parser <- add_option(parser, c("-d", "--outputDir"), action="store", type="character",
                     default=".", help="Directory for saving results [default %default]",
                     metavar="character")

cmd_args <- parse_args(parser)

rds_in <- cmd_args$rds
steric_hindrance <- cmd_args$steric
output_dir <- cmd_args$outputDir

if (!steric_hindrance) {
  helper_in <- "helper_function_em_pause_escape.R"
} else {
  helper_in <- "helper_function_em_steric_hindrance.R"
}

#### load packages ####
suppressPackageStartupMessages({
  library(tidyverse)
  library(plyranges)
})

# #### testing files ####
# root_dir <- "~/Desktop/github/unimod"
#
# bw_dir <- file.path(root_dir, "data/simulation")
# output_dir <- file.path(root_dir, "outputs/simulation")
#
# rds_in <-
#   file.path(bw_dir,
#             paste0("k50ksd25kmin17kmax200l1950a1b1z2000zsd1000zmin1500zmax2500t20n20000s33h17",
#                    "_pos.RDS"))
#
# steric_hindrance <- FALSE
#
# if (!steric_hindrance) {
#   helper_in <- file.path(root_dir, "scripts/helper_function_em_pause_escape.R")
# } else {
#   helper_in <- file.path(root_dir, "scripts/helper_function_em_steric_hindrance.R")
# }
# # rnap_out <- file.path(output_dir, paste0(sel_sample, "_rnap.csv"))

#### set parameters ####
sel_sample <- str_remove(basename(rds_in), "_pos.RDS")
rate_out <- file.path(output_dir, paste0(sel_sample, "_rate.csv"))
walk(c(output_dir), dir.create, showWarnings = FALSE, recursive = TRUE)

sample_cell <- 5000 # number of cells sampled each time
sample_n <- 50 # number of times to sample
matched_len <- 2e4 # gene length to be matched

kmin <- 1
kmax <- 200
matched_gb_len <- matched_len - kmax
count_rnap <- FALSE

# calculate total spacing from sample name
spacing <-
  sum(as.integer(
      str_split(str_split(sel_sample, "[:digit:]s",  simplify = TRUE)[, 2],
                "h")[[1]]))

k <- as.integer(str_remove(str_split(sel_sample, "ksd",  simplify = TRUE)[, 1], "k"))

start_point <- 0.99 * 1e6 # set a start coordinate for the simulated gene

# set lambda according to the read coverage of PRO-seq in the control samples
# from Dukler et al. 2017
lambda <- 102.1
# lambda <- NULL

#### end of parsing arguments ####
source(helper_in)

# retrieve RNAP positions (the last step) for re-sampling
rds <- readRDS(rds_in)
rnap_pos <- rds[[length(rds)]]$pos
# total number of cells
total_cell <- NCOL(rnap_pos)
# gene length of the simulated gene
gene_len <- NROW(rnap_pos) - 1

# count number of RNAP before the pause site
if (count_rnap) {
  # get probability vector
  prob <- readRDS(str_replace(rds_in, ".RDS", "_prob_init.RDS"))
  alpha <- as.double(gsub(".*a([0-9].*)b.*", "\\1",sel_sample))
  beta <- as.double(gsub(".*b([0-9].*)g.*", "\\1",sel_sample))
  zeta <- as.double(gsub(".*z([0-9].*)zsd.*", "\\1",sel_sample))

  # calculate time slice first then get the corresponding probability for beta
  beta_prob <- prob[1, 1] / alpha * beta
  # get pause position for every cell
  idx <- which(prob == beta_prob, arr.ind = TRUE)
  idx <- idx[idx[, 1] != 1, ]
}

#### initiation and pause release rate estimates ####
# generate regions for read counting
gn_rng <-
  GRanges(seqnames = rep("chr1", 3),
          IRanges(start = c(1, kmax + 1, 1),
                  end = c(kmax, gene_len, spacing)))

gn_rng <- shift(gn_rng, shift = start_point)

region_names <- c("tss", "gb", "landing")
names(gn_rng) <- region_names

len <- as.list(width(gn_rng))
names(len) <- region_names

# set seeds for random sampling
seeds <- seq(from = 2013, by = 1, length.out = sample_n)
# a list to Granges for rnap positions
rnap_grng <- list()
# a list recording number of RNAPs at or before the pause site
if (count_rnap) rnap_n_ls <- list()

for (i in 1:sample_n) {
  set.seed(seeds[[i]])
  sel_cells <- sample(1:total_cell, size = sample_cell, replace = TRUE)
  res_pos <- rnap_pos[, sel_cells]
  # get rid of position 1, which is always 1
  res_pos <- res_pos[-1, ]
  if (count_rnap) {
    # get pause sites
    pause_site <- idx[sel_cells, 1] - 1
    # generate data mask
    # inspired by https://stackoverflow.com/questions/47732085/sum-of-some-positions-in-a-row-r
    res_shape <- dim(res_pos)
    after_pause_len <- res_shape[1] - pause_site
    mask_mx <- map2(pause_site, after_pause_len,
         function(x, y) c(rep(TRUE, x), rep(FALSE, y)))
    mask_mx <- Matrix::Matrix(unlist(mask_mx), nrow = res_shape[1], ncol = res_shape[2])
    # calculate rnap number before pause site for every cell
    rnap_n_ls[[i]] <- colSums(res_pos * mask_mx)
  }
  # combine rnap positions across all cells
  res_all <- rowSums(res_pos)
  # generate bigwig for positive strand
  rnap_grng[[i]] <-
    GRanges(seqnames = "chr1",
            IRanges(start = (1 + start_point) : (gene_len + start_point),
                    width = 1),
            score = res_all,
            strand = "+",
            seqlengths = c("chr1" = gene_len * 10)  + start_point)

  rm(res_pos, res_all)
}

# read counting
summarise_bw <-
  function(bw, grng) {
    rc <- grng %>%
      plyranges::group_by_overlaps(bw) %>%
      plyranges::group_by(query) %>%
      plyranges::summarise(score = sum(score))
    if (!1 %in% rc$query) {
      rc <- rbind(DataFrame(list(query = 1, score = 0)), rc)
    }
    rc <- as.list(rc$score)
    names(rc) <- region_names
    return(rc)
  }

bw_dfs <- tibble(trial = 1:sample_n)
bw_dfs$rc_region <- map(rnap_grng, ~ summarise_bw(.x, gn_rng))

bw_dfs$rc_tss <- map_dbl(bw_dfs$rc_region, "tss")
bw_dfs$rc_gb <-map_dbl(bw_dfs$rc_region, "gb")
bw_dfs$rc_landing <-map_dbl(bw_dfs$rc_region, "landing")

#### empirical way to calculate steric hindrance at pause site ####
bw_dfs <- bw_dfs %>%
  mutate(# number of RNAPs per cell per gene
         R = (rc_tss + rc_gb) / sample_cell,
         # number of RNAPs in the pause peak per cell per gene
         R_pause = rc_tss / sample_cell,
         # proportion of landing pad being occupied by RNAP, i.e., empirical phi
         rnap_prop = rc_landing / sample_cell
         )

# whether to match the simulated number of RNAPs to read coverage in experimental data or not
# here match RNAP number within kmin to kmax, RNAP in gene body will be taken care
# of afterwards
if (!is.null(lambda)) {
  rnap_grng <- map(rnap_grng, function(grng) {
    grng$score[kmin:kmax] <-
      rpois(length(kmin:kmax), grng$score[kmin:kmax] / sample_cell * lambda)
    # first 20bp get removed because they are usually not seen in sequencing
    grng$score[1:20] <- 0
    return(grng)
  })
  bw_dfs$rc_region <- map(rnap_grng, ~ summarise_bw(.x, gn_rng))
  bw_dfs$rc_tss <- map_dbl(bw_dfs$rc_region, "tss")
}

# match RNAP number within gene bodies to desired read coverage
if (!is.null(lambda)) {
  pois_mean <- (lambda * bw_dfs$rc_gb / sample_cell) * (matched_gb_len / len$gb)
  bw_dfs$rc_gb <- rpois(length(pois_mean), pois_mean)
  # assign matched gene body length as gene body length
  len$gb <- matched_gb_len
}

# get read counts on each position within pause peak (from kmin to kmax)
bw_dfs$Xk <- map(rnap_grng,
                 ~ .x[(start(.x) >= 990000 + kmin) & (start(.x) <= 990000 + kmax), ]$score)

## Initial model: Poisson based MLEs ##
# use read count within gene body to pre-estimate chi hat
bw_dfs$chi <- bw_dfs$rc_gb / len$gb

# take care of single pause site or variable pause sites (pause peak)
if (str_detect(sel_sample, "ksd0")) {
  bw_dfs$beta_org <- bw_dfs$chi / map_dbl(bw_dfs$Xk, k)
} else {
  bw_dfs$beta_org <- bw_dfs$chi / (bw_dfs$rc_tss / len$tss)
}

bw_dfs$beta_max_rc <- bw_dfs$chi / map_dbl(bw_dfs$Xk, max)

## Adapted model: allows uncertainty in the pause site and steric hindrance ##
# initialize beta using sum of read counts within pause peak
bw_dfs$Xk_sum <- sapply(bw_dfs$Xk, sum)
bw_dfs$beta_int <- bw_dfs$chi / bw_dfs$Xk_sum

# initialize fk with some reasonable values based on heuristic
fk_int <- dnorm(kmin:kmax, mean = 50, sd = 100)
fk_int <- fk_int / sum(fk_int)

# estimate rates using EM
em_ls <- list()
# wrap the EM function in case there is an error
main_EM <- possibly(main_EM, otherwise = NA)

if (steric_hindrance) {
  f <- calculate_f(s = spacing, k = k)
  phi_int <- 0.5
  zeta <- 2000 # elongation rate
  # lambda used for scaling in EM, different from the one used to match coverage
  lambda1 <- 0.0505 * zeta ^ 2
}

for (i in 1:NROW(bw_dfs)) {
  rc <- bw_dfs[i, ]

  if (!steric_hindrance) {

    em_ls[[i]] <- main_EM(Xk = rc$Xk[[1]], kmin = kmin, kmax = kmax,
                          fk_int = fk_int, beta_int = rc$beta_int[[1]], chi_hat = rc$chi,
                          max_itr = 500, tor = 1e-3)

  } else {

    em_ls[[i]] <- main_EM(Xk = rc$Xk[[1]], kmin = kmin, kmax = kmax, f1 = f[["f1"]], f2 = f[["f2"]],
                          fk_int = fk_int, beta_int = rc$beta_int[[1]], chi_hat = rc$chi,
                          phi_int = phi_int, lambda = lambda1, zeta = zeta,
                          max_itr = 500, tor = 1e-4)
  }
}

# get rate estimates and posterior distribution of pause sites
bw_dfs$beta_adp <- map_dbl(em_ls, "beta", .default = NA)
bw_dfs$Yk <- map(em_ls, "Yk", .default = NA)
bw_dfs$fk <- map(em_ls, "fk", .default = NA)
bw_dfs$fk_mean <- map_dbl(em_ls, "fk_mean", .default = NA)
bw_dfs$fk_var <- map_dbl(em_ls, "fk_var", .default = NA)
# calculate Yk / Xk
# bw_dfs$proportion_Yk <- sapply(bw_dfs$Yk, sum) / sapply(bw_dfs$Xk, sum)
bw_dfs$flag <- map_chr(em_ls, "flag", .default = NA)

if (steric_hindrance) bw_dfs$phi <- map_dbl(em_ls, "phi", .default = NA)

# #### empirical "effective" initiation rate ####
# # get idx for the first positions on each gene (cell)
# init_site <- rds[[1]]$pos
# init_site[2:NROW(init_site), ] <- 0
# init_site <- which(init_site == 1)
#
# get_moved_rnap_num <- function(x, init_site) {
#   init_site_idx <- x$pos_pending %in% init_site
#   # RNAP potentially can move
#   c1 <- x$c1[init_site_idx]
#   c2 <- x$c2[init_site_idx]
#   # RNAP actually move
#   y <- c1 & c2
#   # sum RNAP number
#   c1_sum <- sum(c1)
#   y_sum <- sum(y)
#
#   return(c("c1" = c1_sum, "y" = y_sum))
# }
#
# # We are recording 100 time slices, concat the results to get a more robust estimate
# init_rnap <- bind_rows(map(rds, get_moved_rnap_num, init_site))
# init_rnap %>%
#   summarise(
#     # potential initiation rate
#     potential = sum(c1) / total_cell,
#     # effective initiation rate
#     actual = sum(y) / total_cell) %>%
#   write_csv(file = rnap_out)

# add number of RNAP before pause site to output if it exists
if (count_rnap) bw_dfs$rnap_n <- rnap_n_ls

if (steric_hindrance) {
  bw_dfs %>%
    select(trial, chi, beta_org, beta_adp, phi) %>%
    write_csv(rate_out)
} else {
  bw_dfs %>%
    select(trial, chi, beta_org, beta_adp) %>%
    write_csv(rate_out)
}
