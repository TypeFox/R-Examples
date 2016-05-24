cqc <- function(x, y) {
  n   <- length(y)
  tau <- 1:(n-1) / n
  k   <- length(tau)
  rho <- NULL
  for(i in 1:k){
    ytau   <- as.numeric(stats::quantile(y, probs = tau[i]))
    psiy   <- rep(tau[i], n) - I(y - ytau < 0)
    rho[i] <- (x - mean(x)) %*% psiy / (n * sqrt(tau[i] - tau[i] ^ 2) * stats::sd(x)) 
  }
  cqc <- mean(rho)
  return(cqc)
}
