#!/usr/bin/env Rscript

#### load optparse ####
library(optparse)

#### parse command line arguments ####
parser <-
  OptionParser(prog = "./estimate_rates_experiment.R",
               description = "Estimate transcription rates based on experimental data")

parser <- add_option(parser, c("-v", "--verbose"), action="store_true",
                     default=TRUE, help="Print messages [default]")

parser <- add_option(parser, c("-q", "--quietly"), action="store_false",
                     dest="verbose", help="Print no messages")

parser <- add_option(parser, c("--bwp"), action="store", type="character",
                     default=NULL, help="Input bigwig file from the plus strand [default %default]",
                     metavar="character")

parser <- add_option(parser, c("--bwm"), action="store", type="character",
                     default=NULL, help="Input bigwig file from the minus strand [default %default]",
                     metavar="character")

parser <- add_option(parser, c("--grng"), action="store", type="character",
                     default=NULL, help="Gene regions for read counting [default %default]",
                     metavar="character")

parser <- add_option(parser, c("-s", "--steric"), action="store", type="logical",
                     default=FALSE, help="Infer landing-pad occupancy or not [default %default]",
                     metavar="logical")

parser <- add_option(parser, c("--scale"), action="store", type="character",
                     default=NULL, help="A file provides scaling factors for omega [default %default]",
                     metavar="character")

parser <- add_option(parser, c("--type"), action="store", type="character",
                     default="L", help="Scale omega based on [L]ow or [H]igh initiation rate [default %default]",
                     metavar="character")

parser <- add_option(parser, c("-d", "--outputDir"), action="store", type="character",
                     default=".", help="Directory for saving results [default %default]",
                     metavar="character")

cmd_args <- parse_args(parser)

bwp1_p3_in <- cmd_args$bwp
bwm1_p3_in <- cmd_args$bwm
grng_in <- cmd_args$grng
steric_hindrance <- cmd_args$steric
scale_in <- cmd_args$scale
scale_type <- cmd_args$type
result_dir <- cmd_args$outputDir

if (!steric_hindrance) {
  em_in <- "helper_function_em_pause_escape.R"
} else {
  em_in <- "helper_function_em_steric_hindrance.R"
}

helper_in <- "helper_function.R"

#### load packages ####
suppressPackageStartupMessages(
  {
    library(tidyverse)
    library(GenomicRanges)
    library(rtracklayer)
    library(BRGenomics)
    library(plyranges)
  }
)

# #### testing files ####
# root_dir <- "~/Desktop/github/unimod"
# 
# grng_in <- file.path(root_dir, "data/granges_for_read_counting.RData")
# 
# bwp1_p3_in <- file.path(root_dir, "data/PROseq-K562-vihervaara-control-SE_plus.bw")
# bwm1_p3_in <- file.path(root_dir, "data/PROseq-K562-vihervaara-control-SE_minus.bw")
# 
# helper_in <- file.path(root_dir, "scripts/helper_function.R")
# 
# # steric_hindrance <- TRUE
# steric_hindrance <- FALSE
# scale_type <- "L"
# 
# if (!steric_hindrance) {
#   em_in <- file.path(root_dir, "scripts/helper_function_em_pause_escape.R")
# } else {
#   em_in <- file.path(root_dir, "scripts/helper_function_em_steric_hindrance.R")
#   scale_in <- file.path(root_dir, "data/scale_factor.csv")
# }
# 
# if (!steric_hindrance) {
#   result_dir <-
#     file.path(root_dir, "outputs/experiment", sample_id, "pause_escape")
# } else {
#   result_dir <-
#     file.path(root_dir, "outputs/experiment", sample_id, "steric_hindrance")
# }

#### end of parsing arguments ####
# create dir to save output
walk(result_dir, dir.create, recursive = TRUE, showWarnings = FALSE)
rate_tbl_out <- file.path(result_dir, "rate.csv")
sample_id <- str_remove(basename(bwp1_p3_in), "_plus.bw")

source(helper_in)
source(em_in)

# set up parameters
k <- 50
kmin <- 1
kmax <- 200 # also used as k on the poisson case

rnap_size <- 50
zeta <- 2000

rc_cutoff <- 20 # read count cut-off for both gene body and pause peak

if(cmd_args$verbose) message("loading data...")

# load granges for read counting
load(grng_in)

# import and process bigwigs for 3' end
bwp1_p3 <- import.bw(bwp1_p3_in)
bwm1_p3 <- import.bw(bwm1_p3_in)

bwp1_p3 <- process_bw(bw = bwp1_p3, strand = "+")
bwm1_p3 <- process_bw(bw = bwm1_p3, strand = "-")
bw1_p3 <- c(bwp1_p3, bwm1_p3)
rm(bwp1_p3, bwm1_p3)

# make sure pause region is the same as kmax used in EM
bw_pause_filtered <- promoters(bw_pause_filtered, upstream = 0, downstream = kmax)

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

if(cmd_args$verbose) message("estimating rates...")

#### Initial model: Poisson-based Maximum Likelihood Estimation ####
analytical_rate_tbl <-
  tibble(gene_id = rc1$gene_id,
         beta_org =  (rc1$sb1 / rc1$gb_length) / (rc1$sp1 / rc1$pause_length))

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
main_EM <- possibly(main_EM, otherwise = NA)

if (steric_hindrance) {
  scale_tbl <- read_csv(scale_in, show_col_types = FALSE)
  scale_tbl <- scale_tbl %>% filter(sample_id == {{sample_id}})
  
  if (scale_type == "L") {
    omega_scale <- scale_tbl$omega_scale_l
  } else {
    omega_scale <- scale_tbl$omega_scale_h  
    }
  
  em_rate$omega_zeta <- em_rate$chi * omega_scale
  em_rate$omega <- em_rate$omega_zeta / zeta
  
  # compute a scaling factor lambda for the purpose of using the same parameters as simulations
  lambda <- zeta ^ 2 / omega_scale 
  
}

for (i in 1:NROW(em_rate)) {
  rc <- em_rate[i, ]

  if(!steric_hindrance) {
    em_ls[[i]] <- main_EM(Xk = rc$Xk[[1]], kmin = kmin, kmax = kmax,
                          fk_int = fk_int, beta_int = rc$beta_int[[1]],
                          chi_hat = rc$chi, max_itr = 500, tor = 1e-4)
  } else {
    em_ls[[i]] <- main_EM(Xk = rc$Xk[[1]], kmin = kmin, kmax = kmax, f1 = 0.517, f2 = 0.024,
                          fk_int = fk_int, beta_int = rc$beta_int[[1]], phi_int = 0.5,
                          chi_hat = rc$chi, lambda = lambda, zeta = zeta,
                          max_itr = 500, tor = 1e-4)
  }
}

names(em_ls) <- em_rate$gene_id

# get rate estimates and posterior distribution of pause sites
em_rate$beta_adp <- map_dbl(em_ls, "beta", .default = NA)
em_rate$Yk <- map(em_ls, "Yk", .default = NA)
em_rate$fk <- map(em_ls, "fk", .default = NA)
em_rate$fk_mean <- map_dbl(em_ls, "fk_mean", .default = NA)
em_rate$fk_var <- map_dbl(em_ls, "fk_var", .default = NA)
# calculate Yk / Xk
em_rate$t <- sapply(em_rate$Yk, sum)
em_rate$proportion_Yk <- em_rate$t / sapply(em_rate$Xk, sum)

em_rate <- em_rate %>% as_tibble()

if (steric_hindrance) {
  em_rate$phi <- map_dbl(em_ls, "phi", .default = NA)
  em_rate <- em_rate %>%
    mutate(alpha_zeta = omega_zeta / (1 - phi)) 
}

em_rate <- em_rate %>% left_join(analytical_rate_tbl, by = "gene_id")

print(em_rate)

if (!steric_hindrance) {
  em_rate %>%
    select(gene_id, chi, beta_org, beta_adp, fk_mean, fk_var) %>%
    write_csv(rate_tbl_out)
} else {
  em_rate %>%
    select(gene_id, chi, beta_org, beta_adp, fk_mean, fk_var, phi, omega_zeta) %>%
    mutate(beta_zeta = beta_adp * zeta,
           alpha_zeta = omega_zeta / (1 - phi)) %>% 
    write_csv(rate_tbl_out)
}

