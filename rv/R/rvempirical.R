
rvempirical <- function (n, data) {
  n.sims <- getnsims()
  n.all <- (n.sims * n)
  s <- sample(data, size=n.all, replace=TRUE)
  m <- matrix(s, nrow=n.sims, ncol=n)
  x <- rvsims(m, permute=FALSE)
  return(x)
}
