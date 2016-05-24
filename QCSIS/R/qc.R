qc <- function(x, y, tau) {
  n   <- length(y)
  k   <- length(tau)
  rho <- NULL
  for(i in 1:k) {
    ytau   <- as.numeric(stats::quantile(y, probs=tau[i]))
    psiy   <- rep(tau[i], n) - I(y - ytau<0)
    rho[i] <- (x-mean(x)) %*% psiy / (n * sqrt(tau[i] - tau[i]^2) * stats::sd(x)) 
  }
  return(list(tau = tau, rho = rho))
}

