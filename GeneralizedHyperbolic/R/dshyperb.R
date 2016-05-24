dshyperb <- function(x, rho = 0, zeta = 1, param = c(rho,zeta)){
  ## Purpose: Calculate the density of the standardised hyperbolic
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date:  5 Jan 2011, 05:16
  rho <- param[1]
  zeta <- param[2]

  paramStar <- hyperbStandPars(rho, zeta)
  dens <- dhyperb(x, param = paramStar)
  return(dens)
}
