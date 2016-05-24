

rvdiscrete <- function (n=1, x, prob=NULL) {
  n.sims <- getnsims()
  rvsims(matrix(sample(x=x, size=n * n.sims, prob=prob, replace=TRUE), nrow=n.sims))
}

