hyperbStandPars <- function(rho, zeta)
{
  ## get standardized parameters in (alpha, beta) parameterization
  out <- ghypStandPars(rho, zeta, 1)[1:4]
  return(out)
}
