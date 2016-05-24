CQCSIS <- function(x, y, d) {
  p    <-  dim(x)[2] 
  n    <-  dim(x)[1]
  tau  <- 1:(n-1) / n
  k    <- length(tau)
  psi  <- function(x, tau) (tau - I(x < 0))
  MCQC <- matrix(0, k, p)
  for(i in 1:k) {
    ytau     <- as.numeric(stats::quantile(y, probs = tau[i]))
    psiy     <- psi(x = (y - ytau), tau = tau[i])
    w        <- t(scale(x, center = apply(x, 2, mean), scale = F)) %*% psiy / (n * sqrt(tau[i] - tau[i] ^ 2) * apply(x, 2, stats::sd)) 
    MCQC[i,] <- w
  }
  w.cqc <- abs(apply(MCQC, 2, mean))
  M.cqc <- order(w.cqc, decreasing = T)[1:d]
  return(list(w = w.cqc, M = M.cqc))
}