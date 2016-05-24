### This file contains functions related to ASL or ASL*.
### Reference: Kotz, Kozubowski, and Podgorski (2001).

rasl <- function(n, theta = 0, mu = 0, sigma = 1){
  if(sigma <= 0){
    stop("sigma > 0 is required.")
  }

  if(length(n) > 1){
    n <- length(n)
  }

  ### Ref: Kotz, et. al. (2001), p145, eqn (3.2.1)
  ### Y ~ theta + mu * W + sigma  sqrt(W) Z
  ### where Z ~ N(0, 1) and W ~ Exp(1) are independent.
  Z <- rnorm(n)
  W <- rexp(n)
  ret <- theta + mu * W + sigma * sqrt(W) * Z
  ret
} # End of rasl().

rasla <- function(n, theta = 0 , kappa = 1, sigma = 1){
  if(sigma <= 0){
    stop("sigma > 0 is required.")
  }
  if(kappa <= 0){
    stop("kappa > 0 is required.")
  }

  ### Ref: Kotz, et. al. (2001), p136, eqn (3.1.4) & (3.1.5).
  ### kappa = (sqrt(2 * sigma^2 + mu^2) - mu) / sqrt(2 * sigma)
  ### mu = sigma / sqrt(2) * (1 / kappa - kappa)
  mu <- sigma / sqrt(2) * (1 / kappa - kappa)
  ret <- rasl(n, theta, mu, sigma)
  ret
} # End of rasla().

