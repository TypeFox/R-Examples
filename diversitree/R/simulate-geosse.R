## There is no make.tree.geosse().  Mooch make.tree.classe() instead.

tree.geosse <- function(pars, max.taxa=Inf, max.t=Inf,
                        include.extinct=FALSE, x0=NA) {
  check.pars.nonnegative(pars, 7)
  pars.cl <- pars.ge.to.cl(pars)
  k <- 3
  x0 <- x0 + 1

  if ( is.na(x0) )
    x0 <- sample(k, 1, FALSE, stationary.freq.classe(pars.cl, k))
  if ( length(x0) != 1 || is.na(x0) || x0 < 1 || x0 > k )
    stop(paste("Invalid root state", x0-1))

  info <- make.tree.classe(pars.cl, 3, max.taxa, max.t, x0)
  info$state <- info$state - 1

  phy <- me.to.ape.bisse(info[-1,], info$state[1])
  if ( include.extinct || is.null(phy) )
    phy
  else
    prune(phy)
}

## Frame geosse parameters in classe terms
pars.ge.to.cl <- function(pars.ge)
{
  if (is.null(names(pars.ge)))
    names(pars.ge) <- default.argnames.geosse()
  pars.cl <- rep(0, 27)
  names(pars.cl) <- default.argnames.classe(3)
  pars.cl['lambda222'] <- pars.cl['lambda112'] <- pars.ge['sA']
  pars.cl['lambda333'] <- pars.cl['lambda113'] <- pars.ge['sB']
  pars.cl['lambda123'] <-  pars.ge['sAB']
  pars.cl['mu2'] <- pars.cl['q13'] <- pars.ge['xA']
  pars.cl['mu3'] <- pars.cl['q12'] <- pars.ge['xB']
  pars.cl['q21'] <- pars.ge['dA']
  pars.cl['q31'] <- pars.ge['dB']
  pars.cl
}
