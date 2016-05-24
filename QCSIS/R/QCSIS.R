QCSIS <- function(x, y, tau = 1:(n-1)/n, d) {
  p   <- dim(x)[2] 
  n   <- dim(x)[1]
  k   <- length(tau)
  psi <- function(x, taus) (taus-I(x<0))
  MQC <- matrix(0, k, p)
  for(i in 1:k) {
    ytaus   <- as.numeric(stats::quantile(y, probs = tau[i]))
    psiy    <- psi(x = (y-ytaus), taus = tau[i])
    w       <- t(scale(x, center=apply(x, 2, mean), scale = F)) %*% psiy / (n * sqrt(tau[i] - tau[i] ^ 2) * apply(x, 2, stats::sd)) 
    MQC[i,] <- w^2
  }
  w.qc <- apply(MQC, 2, mean)
  M.qc <- order(w.qc, decreasing=T)[1:d]
  return(list(w = w.qc, M = M.qc))
}