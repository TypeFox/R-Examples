fmodel1pp <-
function(zeta, y, bpar, prior = dnorm, ...) {
  m <- length(bpar)
  return(fmodel4pp(zeta, y, apar = rep(1,m), bpar = bpar,
    cpar = rep(0,m), dpar = rep(1,m), prior = prior, ...))
}
