kappaLambda <- function(lambda, zeta){
  kL <- besselRatio(zeta, lambda, 1)/zeta
  return(kL)
}

alphaStand <- function(rho, zeta, lambda = 1){
  DK <- kappaLambda(lambda + 1, zeta) - kappaLambda(lambda, zeta)
  kL <- kappaLambda(lambda, zeta)
  alphaStar <- (zeta^2/(1 - rho^2)*kL*(1 + rho^2*zeta^2/(1 - rho^2)*DK))^(1/2)
  return(alphaStar)
}

betaStand <- function(rho, zeta, lambda = 1){
  DK <- kappaLambda(lambda + 1, zeta) - kappaLambda(lambda, zeta)
  kL <- kappaLambda(lambda, zeta)
  alphaStar <- (zeta^2/(1 - rho^2)*kL*(1 + rho^2*zeta^2/(1 - rho^2)*DK))^(1/2)
  betaStar <- alphaStar*rho
  return(betaStar)
}

ghypStandPars <- function(rho, zeta, lambda = 1)
{
  DK <- kappaLambda(lambda + 1, zeta) - kappaLambda(lambda, zeta)
  kL <- kappaLambda(lambda, zeta)
  alphaStar <- (zeta^2/(1 - rho^2)*kL*(1 + rho^2*zeta^2/(1 - rho^2)*DK))^(1/2)
  betaStar <- alphaStar*rho
  deltaStar <- (1/alphaStar)*zeta/sqrt(1 - rho^2)
  muStar <- -(1/alphaStar)*rho*zeta^2*kL/(1 - rho^2)

  out <- c(muStar,deltaStar,alphaStar, betaStar,lambda)
  names(out) <- c("muStar","deltaStar","alphaStar"," betaStar","lambda")
  return(out)
}
