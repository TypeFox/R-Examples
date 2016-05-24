fmodel3pl <-
function(zeta, y, apar, bpar, cpar, prior = dnorm, ...) {
  return(fmodel4pl(zeta, y, apar, bpar, cpar, 
    dpar = rep(1, length(bpar)), prior = prior, ...))
}