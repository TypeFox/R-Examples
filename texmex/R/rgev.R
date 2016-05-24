rgev <- function(n, mu, sigma, xi){
  ## use standard GEV ~ (exp(-xi log E) - 1) / xi

  ## where E is a standard exponential

  neg.log.exp <- -log(rexp(n))

  ## expand mu, sigma, and xi to be n long
  ## this is necessary to ensure that we get
  ## exactly n random numbers if mu, sigma, xi
  ## are greater than n long
  n     <- length(neg.log.exp)
  mu    <- rep(mu, length.out=n)
  sigma <- rep(sigma, length.out=n)
  xi    <- rep(xi, length.out=n)

  ## and here we go
  standard.gev <- .exprel(neg.log.exp * xi) * neg.log.exp

  mu + sigma * standard.gev
}
