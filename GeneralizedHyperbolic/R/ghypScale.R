ghypScale <- function(newMean, newSD,
                      mu = 0, delta = 1, alpha = 1, beta = 0, lambda = 1,
                      param = c(mu,delta,alpha,beta,lambda))
{
  ## Purpose: Rescale generalized hyperbolic to specified mean and sd
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date:  7 Jan 2011, 22:49

  ## Lambda defaults to one if omitted from param vector
  if (length(param) == 4){
    param <- c(param,1)
  }

  ## check parameters
  parResult <- ghypCheckPars(param)
  case <- parResult$case
  errMessage <- parResult$errMessage

  if (case == "error")
    stop(errMessage)

  param <- as.numeric(param)
  mu <- param[1]
  delta <- param[2]
  alpha <- param[3]
  beta <- param[4]
  lambda <- param[5]
  gamma <- sqrt(alpha^2 - beta^2)
  rho <- beta/alpha
  zeta <- delta*gamma

  varFactor <- besselRatio(zeta, lambda, 1)/zeta + (rho^2/(1 - rho^2))*
                 (besselRatio(zeta, lambda, 2) - besselRatio(zeta, lambda, 1)^2)
  newDelta <- newSD/sqrt(varFactor)
  newMu <- newMean - (newDelta*rho/sqrt(1 - rho^2))*besselRatio(zeta, lambda, 1)
  newParam <- ghypChangePars(2, 1, c(newMu,newDelta,rho,zeta,lambda))
  return(newParam)
}
