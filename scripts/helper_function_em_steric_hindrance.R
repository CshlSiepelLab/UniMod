#### function for adapted model ####
# model allows both varied pause sites and steric hindrance, EM contains phi estimations
# functions for EM based on Gaussian distributed k
get_expectation <- function(fk, Xk, beta) {

  Yk <- Xk / (1 - beta + beta / fk)

  return(Yk)
}

mult.RNAP.phi <-
  function(alpha, beta, f1, f2) {
    return(
      (1 - f1 -f2) * alpha / (alpha + beta) +
        f1 * alpha^2 / (alpha^2 + beta^2 + alpha*beta) +
        f2 * alpha^3 / (beta^2 * alpha + beta^3 + alpha^2 * beta + alpha^3)
    )
   }

phi.polynom <- function(phi, beta, omega, f1, f2) {
  
  alpha <- omega / (1-phi)  
  
  return(mult.RNAP.phi(alpha, beta, f1, f2) - phi)	
}

# find phi corresponding to omega and beta by solving polynomial that results from
# substituting alpha = omega/(1-phi)into the equation that defines phi in terms of alpha and beta 
mult.RNAP.phi.omega <- function(omega, beta, f1, f2) {
  # set bounds for solution
  lb = 1e-6
  ub = 1-1e-6
  epsilon = 1e-3
  
  # make sure opposite signs at bounds; if not treat as an edge case
  phi1 = phi.polynom(lb, beta, omega, f1, f2)
  phi2 = phi.polynom(ub, beta, omega, f1, f2)
  
  if ((phi1 > 0 & phi2 > 0) | (phi1 < 0 & phi2 < 0)) {
    # in this case phi is almost certainly close to 0 or 1
    # pick the closer case
    if (abs(phi1) < epsilon)
      return(epsilon)
    else if (abs(phi2) < epsilon)
      return(1-epsilon)
    else
      simpleError("Edge case fail")
  }
  
  try(phi.root <- uniroot(phi.polynom, c(lb,ub), beta, omega, f1, f2), silent = F)	
  
  return(phi.root$root)
}

# beta.ecll.omega <- function(args, omega, chi, t, f1, f2) {
#   beta = args[1]
#   
#   # fix Beta prior with a=b=2
#   a=2
#   b=2
#   
#   phi = mult.RNAP.phi.omega(omega, beta, f1, f2)
#   
#   if (phi <= 0 | phi >= 1) { 
#     retval = -Inf 
#   }
#   else {
#     retval = -t*log(beta) - chi/beta + (a-1)*log(phi) + (b-1)*log(1-phi)
#   }
#   
#   return(-retval)  # minimization!
# }
# 
# beta.M.step.omega = function(chi, t, f1, f2, oldphi, oldbeta, lambda, zeta) {
#   omega = chi*zeta / lambda  # this will be const; could just be passed in
#   ret = list("par" = NA_integer_)	
#   try(ret<-optim(c(oldbeta), beta.ecll.omega, gr=NULL, omega=omega, chi=chi, t=t, f1=f1, f2=f2, method="L-BFGS-B", lower=1e-6),silent=T)
#   beta <- ret$par
#   # set a ceiling for phi
#   phi <- mult.RNAP.phi.omega(omega, beta, f1, f2)
#   
#   return(list("beta" = beta, "phi" = phi))
# }

# version of above that uses log parameterization of beta
beta.ecll.omega.log <- function(args, omega, chi, t, f1, f2) {
  beta = exp(args[1])
  
  # fix Beta prior with a=b=2
  a = 2
  b = 2
  
  phi = mult.RNAP.phi.omega(omega, beta, f1, f2)
  
  if (phi <= 0 | phi >= 1) { 
    retval = -Inf 
  }
  else {
    retval = -t*log(beta) - chi/beta + (a-1)*log(phi) + (b-1)*log(1-phi)
  }
  return(-retval) # minimization!
}

beta.M.step.omega <- function(chi, t, f1, f2, oldphi, oldbeta, lambda, zeta) {
  
  omega = chi*zeta / lambda # this will be const; could just be passed in
  ret = list("par" = NA_integer_)	

  try(ret <-
        optim(c(log(oldbeta)), beta.ecll.omega.log, gr=NULL, omega=omega, chi=chi,
              t=t, f1=f1, f2=f2, method="L-BFGS-B", lower=log(omega)), silent = F)

  beta <- exp(ret$par)
  # set a ceiling for phi
  phi <- mult.RNAP.phi.omega(omega, beta, f1, f2)
  
  return(list("beta" = beta, "phi" = phi))
}

get_maximization <- function(chi_hat, Xk, Yk, fk, kmin, kmax, f1, f2,
                             phi, beta, lambda, zeta) {

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

  param <- beta.M.step.omega(chi = chi_hat, t = t, f1 = f1, f2 = f2, oldphi = phi,
                      oldbeta = beta, lambda = lambda, zeta = zeta)
  
  return(list("beta" = param$beta, "phi" = param$phi, "fk" = fk,
              "fk_mean" = fk_mean, "fk_var" = fk_var))
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

main_EM <- function(Xk, kmin, kmax, f1, f2, fk_int, beta_int, phi_int, chi_hat,
                    lambda, zeta, max_itr = 100, tor = 1e-3) {
  # lists to record changes of likelihood and betas in iterations
  betas <- list()
  likelihoods <- list()
  # default flag is normal
  flag <- "normal"

  for (i in 1:max_itr) {
    if (i == 1) {
      Yk <- get_expectation(fk_int, Xk, beta_int)
      hats <- get_maximization(chi_hat, Xk, Yk, fk_int, kmin, kmax, f1, f2,
                               phi_int, beta_int, lambda, zeta)
      beta <- beta_int
    }
    if (i != 1) {
      Yk <- get_expectation(hats$fk, Xk, hats$beta)
      hats <- get_maximization(chi_hat, Xk, Yk, hats$fk, kmin, kmax, f1, f2,
                               hats$phi, hats$beta, lambda, zeta)
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
              "betas" = betas, "likelihoods" = likelihoods,
              "phi"= hats$phi,  "flag" = flag))
}

calculate_f <- function(s, k) {
  # sd is set as 25 here
  x <- round(rnorm(1e7, mean = k, sd = 25))
  x <- x[x >= 17 & x <= 200]
  f <- mean(x > s) 
  f1 <- mean((x > s) & (x <= 2 * s))
  f2 <- mean(x > 2 * s)
  return(c("f" = f, "f1" = f1, "f2" = f2))
}
