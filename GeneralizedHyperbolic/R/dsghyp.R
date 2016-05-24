dsghyp <- function(x, rho = 0, zeta = 1, lambda = 1,
                   param = c(rho,zeta,lambda)){
  ## Purpose: Calculate the density of the standardised generalized hyperbolic
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date:  5 Jan 2011, 05:16
  rho <- param[1]
  zeta <- param[2]
  lambda <- param[3]

  paramStar <- ghypStandPars(rho, zeta)
  dens <- dghyp(x, param = paramStar)
  return(density)
}
