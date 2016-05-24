fmodel3pp <-
function(zeta, y, apar, bpar, cpar, prior = dnorm, ...) {
  return(fmodel4pp(zeta, y, apar, bpar, cpar, 
    dpar = rep(1, length(bpar)), prior = prior, ...))
}