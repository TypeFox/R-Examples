### This file contains functions related to ASL or ASL*.
### Reference: Kotz, Kozubowski, and Podgorski (2001).

pasl <- function(q, theta = 0, mu = 0, sigma = 1, lower.tail = TRUE,
    log.p = FALSE){
  if(sigma <= 0){
    stop("sigma > 0 is required.")
  }

  ### Ref: Kotz, et. al. (2001), p136, eqn (3.1.4) & (3.1.5).
  ### kappa = (sqrt(2 * sigma^2 + mu^2) - mu) / sqrt(2 * sigma)
  ### mu = sigma / sqrt(2) * (1 / kappa - kappa)
  kappa <- (sqrt(2 * sigma^2 + mu^2) - mu) / sqrt(2 * sigma)
  pasla(q, theta, kappa, sigma, lower.tail = lower.tail, log.p = log.p)
} # End of pasl().

pasla <- function(q, theta = 0, kappa = 1, sigma = 1, lower.tail = TRUE,
    log.p = FALSE){
  if(sigma <= 0){
    stop("sigma > 0 is required.")
  }
  if(kappa <= 0){
    stop("kappa > 0 is required.")
  }

  ret <- do.call("c", lapply(q, pasla.one.log, theta, kappa, sigma))
  if(! log.p){
    ret <- exp(ret)
  }
  if(! lower.tail){
    if(! log.p){
      ret <- 1 - ret
    } else{
      ret <- log(1 - exp(ret))
    }
  }
  ret
} # End of pasla().

pasla.one.log <- function(q, theta, kappa, sigma){
  ### Ref: Kotz, et. al. (2001), p138, eqn (3.1.11).
  ### if log = FALSE
  # if(q >= theta){
  #   ret <- 1 - 1 / (1 + kappa^2) * 
  #          exp(-sqrt(2) / sigma * kappa * (q - theta))
  # } else{
  #   ret <- kappa^2 / (1 + kappa^2) *
  #          exp(-sqrt(2) / sigma / kappa * (theta - q))
  # }

  if(q == -Inf){
    ret <- -Inf
  } else if(q == Inf){
    ret <- 0 
  } else{
    ### Ref: Kotz, et. al. (2001), p138, eqn (3.1.11).
    if(q >= theta){
      ret <- log(1 / (1 + kappa^2)) + 
             -sqrt(2) / sigma * kappa * (q - theta)
      ret <- log(1 - exp(ret))
    } else{
      ret <- log(kappa^2 / (1 + kappa^2)) +
             -sqrt(2) / sigma / kappa * (theta - q)
    }
  }

  ret
} # End of pasla.one.log().
