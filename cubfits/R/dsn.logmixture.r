### For log mixture normal distribution.
### param is stored in the format of
###   c(eta_1, ..., eta_K, mu_1, ..., mu_K, sd_1, ..., sd_K)
### with sum_k eta_k = 1 and mu_1 < ... < mu_K.

### Density of mixture normal distribution in K components.
dmixnorm <- function(x, param, log = FALSE){
  K <- length(param) / 3

  ret <- list()
  for(i.k in 1:K){
    log.eta.k <- log(param[i.k])
    mu.k <- param[i.k + K]
    sd.k <- param[i.k + K + K]
    ret[[i.k]] <- log.eta.k + dnorm(x, mean = mu.k, sd = sd.k, log = TRUE)
  }

  ### Return.
  ret <- rowSums(exp(do.call("cbind", ret)))
  if(log){
    ret <- log(ret)
  }
  ret
} # End of dmixnorm().

### Density of log mixture normal distribution in K components.
dlmixnorm <- function(x, paramlog, log = FALSE, x.in.log = TRUE){
  if(!x.in.log){
    log.x <- log(x)
  } else{
    log.x <- x
  }
  ret <- dmixnorm(log.x, paramlog, log = TRUE) - log.x

  ### Return.
  if(!log){
    ret <- exp(ret)
  }
  ret
} # End of dlmixnorm().

### Compute eta_k * f(x_i; mu_k, sigma_k)
### Return a N by K matrix.
proplmixnorm <- function(x, paramlog, x.in.log = TRUE){
  K <- length(paramlog) / 3

  if(!x.in.log){
    log.x <- log(x)
  } else{
    log.x <- x
  }

  ret <- list()
  for(i.k in 1:K){
    log.eta.k <- log(paramlog[i.k])
    mu.k <- paramlog[i.k + K]
    sigma.k <- paramlog[i.k + K + K]
    ret[[i.k]] <- log.eta.k + dnorm(x, mean = mu.k, sd = sigma.k, log = TRUE)
  }

  ### Return.
  ret <- exp(do.call("cbind", ret))
  ret
} # End of proplmixnorm().

