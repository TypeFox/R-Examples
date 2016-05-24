expectedSignatureLifetimeExp <- function(s, rate=1) {
  n <- length(s)
  xin <- rep(0, n)
  for(i in 1:n) {
    xin[i] <- rate * sum(1/(n-1:i+1))
  }
  sum(s * xin)
}

expectedSystemLifetimeExp <- function(g, rate=1) {
  s <- computeSystemSignature(g)
  expectedSignatureLifetimeExp(s, rate)
}

expectedNetworkLifetimeExp <- function(g, rate=1) {
  s <- computeNetworkSignature(g)
  expectedSignatureLifetimeExp(s, rate)
}
