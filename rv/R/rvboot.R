


rvboot <- function (data) {
#  empirical (bootstrap) distribution
#
  n.sims <- getnsims()
  n <- n.sims*length(data)
  s <- matrix(sample(data, size=n, replace=TRUE), nrow=n.sims)
  r <- rvsims(s)
  dim(r) <- dim(data)
  return(r)
}

