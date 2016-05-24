### This file contains functions related to ASL or ASL*.
### Reference: Kotz, Kozubowski, and Podgorski (2001).

qasl <- function(p, theta = 0, mu = 0, sigma = 1, lower.tail = TRUE,
    log.p = FALSE){
  if(sigma <= 0){
    stop("sigma > 0 is required.")
  }

  ### Ref: Kotz, et. al. (2001), p136, eqn (3.1.4) & (3.1.5).
  ### kappa = (sqrt(2 * sigma^2 + mu^2) - mu) / sqrt(2 * sigma)
  ### mu = sigma / sqrt(2) * (1 / kappa - kappa)
  kappa <- (sqrt(2 * sigma^2 + mu^2) - mu) / sqrt(2 * sigma)
  qasla(p, theta, kappa, sigma, lower.tail = lower.tail, log.p = log.p)
} # End of qasl().

qasla <- function(p, theta = 0, kappa = 1, sigma = 1, lower.tail = TRUE,
    log.p = FALSE){
  if(sigma <= 0){
    stop("sigma > 0 is required.")
  }
  if(kappa <= 0){
    stop("kappa > 0 is required.")
  }

  ret <- do.call("c", lapply(p, qasla.one.log, theta, kappa, sigma,
                             lower.tail = lower.tail, log.p = log.p))
  ret
} # End of qasla().

qasla.one.log <- function(p, theta, kappa, sigma, lower.tail = TRUE,
    log.p = FALSE){
  if(log.p){
    if(p > 0){
      return(NaN)
    }
  } else{
    if(p < 0 | p > 1){
      return(NaN)
    }
  }

  ### convert to p any way.
  if(log.p){
    p <- exp(p)
  }
  if(! lower.tail){
    p <- 1 - p
  }

  if(p == 0){
    ret <- -Inf
  } else if(p == 1){
    ret <- Inf
  } else{
    ### Ref: Kotz, et. al. (2001), p143, eqn (3.1.32).
    ret <- theta
    if(p <= kappa^2 / (1 + kappa^2)){
      ret <- ret + sigma * kappa / sqrt(2) * log((1 + kappa^2) / kappa^2 * p)
    } else{
      ret <- ret - sigma / (sqrt(2) * kappa) * log((1 + kappa^2)* (1- p))
    }
  }
  ret
} # End of qasla.one.log().
