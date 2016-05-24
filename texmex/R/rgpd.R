rgpd <- function(n, sigma, xi, u = 0) {
  ## a standard GPD is (exp(xi * Exp(1)) - 1) / xi
  ## and the rest follows
  exponentials <- rexp(n)

  ## expand mu, sigma, and xi to be n long
  ## this is necessary to ensure that we get
  ## exactly n random numbers if mu, sigma, xi
  ## are greater than n long
  n     <- length(exponentials)
  sigma <- rep(sigma, length.out=n)
  xi    <- rep(xi, length.out=n)
  u     <- rep(u, length.out=n)

  u + sigma * .exprel(exponentials * xi) * exponentials
}


