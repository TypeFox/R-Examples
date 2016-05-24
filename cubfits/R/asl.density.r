### This file contains functions related to ASL or ASL*.
### Reference: Kotz, Kozubowski, and Podgorski (2001).

dasl <- function(x, theta = 0, mu = 0, sigma = 1, log = FALSE){
  if(sigma <= 0){
    stop("sigma > 0 is required.")
  }

  ### Ref: Kotz, et. al. (2001), p136, eqn (3.1.4) & (3.1.5).
  ### kappa = (sqrt(2 * sigma^2 + mu^2) - mu) / sqrt(2 * sigma)
  ### mu = sigma / sqrt(2) * (1 / kappa - kappa)
  kappa <- (sqrt(2 * sigma^2 + mu^2) - mu) / sqrt(2 * sigma)
  dasla(x, theta, kappa, sigma, log = log)
} # End of dasl().

dasla <- function(x, theta = 0 , kappa = 1, sigma = 1, log = FALSE){
  if(sigma <= 0){
    stop("sigma > 0 is required.")
  }
  if(kappa <= 0){
    stop("kappa > 0 is required.")
  }

  ret <- do.call("c", lapply(x, dasla.one.log, theta, kappa, sigma))
  if(! log){
    ret <- exp(ret)
  }
  ret
} # End of dasla().

dasla.one.log <- function(x, theta, kappa, sigma){
  ### Ref: Kotz, et. al. (2001), p137, eqn (3.1.10).
  ### if log = FALSE
  # ret <- (sqrt(2) / sigma) * (kappa / (1 + kappa^2))
  # if(x >= theta){
  #   ret <- ret * exp(-sqrt(2) / sigma * kappa * abs(x - theta))
  # } else{
  #   ret <- ret * exp(-sqrt(2) / sigma / kappa * abs(x - theta))
  # }

  if(x == -Inf || x == Inf){
    ret <- 0
  } else{
    ### Ref: Kotz, et. al. (2001), p137, eqn (3.1.10).
    ret <- log((sqrt(2) / sigma) * (kappa / (1 + kappa^2)))
    if(x >= theta){
      ret <- ret + (-sqrt(2) / sigma * kappa * (x - theta))
    } else{
      ret <- ret + (-sqrt(2) / sigma / kappa * (theta - x))
    }
  }

  ret
} # End of dasla.one.log().

