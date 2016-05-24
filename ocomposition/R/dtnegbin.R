dtnegbin <- function(x, mu, dispersion, l.bound){
  N <- length(x)
  mu <- rep(mu,length=N)
  l.bound <- rep(l.bound,length=N)
  pmf <- dnbinom(x, mu = mu, size = dispersion)
  pmf <- pmf/pnbinom(l.bound - .1, mu = mu, size = dispersion, lower.tail = FALSE)
  pmf[x <= l.bound - .1] <- 0
  pmf
}

