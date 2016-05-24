##====================================================================
## Prefer corr and sd to covariance?
##====================================================================

"cov2corr" <- function(cov) {
  s <- sqrt(diag(cov))
  res <- cov / outer(X = s, Y = s, FUN = "*")
  res
}
