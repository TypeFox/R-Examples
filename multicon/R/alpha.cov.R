alpha.cov <-
function(sigma) {
  p <- dim(sigma)[1]
  alpha <- (p / (p - 1)) * (1 - sum(diag(sigma))/sum(sigma))
  return(alpha)
}
