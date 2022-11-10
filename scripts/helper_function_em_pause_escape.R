#### function for adapted model ####
# model allows pause sites to vary across cells
# EM doesn't include phi estimates
# functions for EM based on Gaussian distributed k
get_expectation <- function(fk, Xk, beta) {

  Yk <- Xk / (1 - beta + beta / fk)

  return(Yk)
}

get_maximization <- function(chi_hat, Xk, Yk, fk, kmin, kmax) {

  t <- sum(Yk)
  u <- sum(Yk * seq(kmin, kmax))
  v <- sum(Yk * seq(kmin, kmax) ^ 2)

  w <- sum(fk / (1- fk) * (Xk - Yk))
  z <- sum(fk / (1- fk) * (Xk - Yk) *  seq(kmin, kmax))
  r <- sum(fk / (1- fk) * (Xk - Yk) *  seq(kmin, kmax) ^ 2)

  fk_mean <- (u - z) / (t - w)
  fk_var <- (v - r) / (t - w) - fk_mean ^ 2

  # avoid small and negative values
  if (fk_var < 1e-10) {
    fk[1:length(fk)] <- 0
    # sometimes it looks like an integer but actually it's not
    fk[round(fk_mean)] <- 1
  } else {
    fk <- dnorm(kmin:kmax, mean = fk_mean, sd = fk_var ^ 0.5)
    fk <- fk / sum(fk)
  }

  beta <- chi_hat / t

  return(list("beta" = beta, "fk" = fk, "fk_mean" = fk_mean, "fk_var" = fk_var))
}

get_likelihood <- function(beta, chi, Xk, Yk, fk) {
  # part of the original likelihood function associated with beta, Xk and Yk
  # used as criteria to terminate EM
  t <- sum(Yk)
  # take care of the 0s
  idx_1 <- fk != 0
  idx_2 <- 1 - fk != 0
  likelihood <- -t * log(beta) - chi / beta +
    sum(Yk[idx_1] * log(fk[idx_1])) + sum((Xk - Yk)[idx_2] * log(1 - fk[idx_2]))
  return(likelihood)
}

main_EM <- function(fk_int, Xk, kmin, kmax, beta_int, chi_hat, max_itr = 100,
                    tor = 1e-3) {
  # lists to record changes of likelihood and betas in iterations
  betas <- list()
  likelihoods <- list()
  # default flag is normal
  flag <- "normal"

  for (i in 1:max_itr) {
    if (i == 1) {
      Yk <- get_expectation(fk_int, Xk, beta_int)
      hats <- get_maximization(chi_hat, Xk, Yk, fk_int, kmin, kmax)
      beta <- beta_int
    }
    if (i != 1) {
      Yk <- get_expectation(hats$fk, Xk, hats$beta)
      hats <- get_maximization(chi_hat, Xk, Yk, hats$fk, kmin, kmax)
    }

    likelihoods[[i]] <-
      get_likelihood(beta = hats$beta, chi = chi_hat, Xk = Xk, Yk = Yk, fk = hats$fk)

    betas[[i]] <- hats$beta

    if (any(hats$fk == 1)) {
      hats$beta <- chi_hat / Xk[which(hats$fk == 1)]
      flag <- "single_site"
      break}

    if (i > 1) {
      diff <- likelihoods[[i]] - likelihoods[[i-1]]
      if (diff <= tor) break
      }
    }

  if (i == max_itr) flag <- "max_iteration"

  # message("Done!")

  return(list("beta" = hats$beta, "Yk" = Yk, "fk" = hats$fk,
              "fk_mean" = hats$fk_mean, "fk_var" = hats$fk_var,
              "betas" = betas, "likelihoods" = likelihoods, "flag" = flag))
}
