ghypMomBNS <- function(order,
                          mu = 0, delta = 1, alpha = 1, beta = 0, lambda = 1,
                          param = c(mu, delta, alpha, beta, lambda),
                          tol = 10^(-12))
{
  ## Purpose: Calculate moments using Barndorf-Nielsen and Steltzer (2005)
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date:  9 Nov 2010, 15:16

  ## unpack parameters
  param <- as.numeric(param)

  if (length(param) == 4)
    param <- c(param, 1)

  ## check parameters
  parResult <- ghypCheckPars(param)
  case <- parResult$case
  errMessage <- parResult$errMessage

  if (case == "error")
    stop(errMessage)

  mu <- param[1]
  delta <- param[2]
  alpha <- param[3]
  beta <- param[4]
  lambda <- param[5]

  alphabar <- delta*alpha
  betabar <- delta*beta
  gamma <- sqrt(alpha^2 - beta^2)
  gammabar <- delta * gamma

  m <- order%%2
  ord2 <- ceiling(order/2)

  ## initialize sum
  const <- 2^ord2*gammabar^lambda*delta^(2*ord2)*beta^m/
           (sqrt(pi)*besselK(gammabar, lambda)*alphabar^(lambda + ord2))
  ak <- const*gamma(ord2 + 0.5)*besselK(alphabar, lambda + ord2)/gamma(m + 1)
  mom <- ak
  k <- 0
  cont <- TRUE
  while (cont){
    ak <- ak*2*betabar^2*(k  + ord2 + 0.5)/
          (alphabar*(2*k + 2 + m)*(2*k + 1 + m))*
          besselRatio(alphabar, lambda + k  + ord2, 1)
    if (ak/mom < tol){
      cont <- FALSE
    }
    mom <- mom + ak
    ## cat(k, 2*betabar^2*(k  + ord2 + 0.5),
    ##     besselRatio(alphabar, lambda + k  + ord2, 1), ak, mom,"\n")
    k <- k + 1
  }

  return(mom)
}
