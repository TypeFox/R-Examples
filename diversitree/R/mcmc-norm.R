## Updates for MCMC from a homogeneous normal distribution, with
## single parameter updates.  The 'w' parameter now represents the
## standard deviation of the normal updates.
sampler.norm <- function(lik, x.init, y.init, w, lower, upper, control) {
  for ( i in seq_along(x.init) ) {
    xy <- mcmc.norm.1d(make.unipar(lik, x.init, i),
                       x.init[i], y.init, w[i], lower[i], upper[i])
    x.init[i] <- xy[1]
    y.init    <- xy[2]
  }

  list(x.init, y.init)
}

mcmc.norm.1d <- function(f, x.init, y.init, w, lower, upper, control) {
  x.new <- rnorm(1, x.init, w)

  if ( x.new < lower || x.new > upper )
    y.new <- -Inf
  else
    y.new <- f(x.new)

  alpha <- exp(y.new - y.init)
  if ( alpha >= 1 || runif(1) < alpha )
    c(x.new, y.new)
  else
    c(x.init, y.init)
}
