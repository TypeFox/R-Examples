fmodel2pl <-
function(zeta, y, apar, bpar, prior = dnorm, ...) {
  m <- length(bpar)
  return(fmodel4pl(zeta, y, apar = apar, bpar = bpar,
    cpar = rep(0,m), dpar = rep(1,m), prior = prior, ...))
}
