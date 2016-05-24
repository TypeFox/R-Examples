mpinv <- function(M,eps=1e-13) {
  s <- svd(M)
  e <- s$d
  e[e>eps] <- 1/e[e>eps]
  return(s$v%*%diag(e)%*%t(s$u))
}
