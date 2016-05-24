

rvdirichlet <- function (n = 1, alpha)  {
  x <- NULL
  for (i in 1:n) {
    g <- rvgamma(n = 1, shape = alpha, scale = 1)
    x <- cbind.rv(x, g/sum(g))
  }
  return(x)
}

