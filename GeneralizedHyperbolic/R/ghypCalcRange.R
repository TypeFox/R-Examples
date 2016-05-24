### Function to calculate an effective range for the distribution function
### or for the density function
### DJS 19/12/06
ghypCalcRange <- function(mu = 0, delta = 1, alpha = 1, beta = 0, lambda = 1,
                          param = c(mu, delta, alpha, beta, lambda),
                          tol = 10^(-5), density = TRUE, ...) {

  param <- as.numeric(param)

  if (length(param) == 4)
    param <- c(param, 1)

  mu <- param[1]
  delta <- param[2]
  alpha <- param[3]
  beta <- param[4]
  lambda <- param[5]

  if (!density) {
    ## bounds are for distribution function
    stop("Distribution function bounds not yet implemented")
  } else {
    ## bounds are for the density function
    mode <- ghypMode(param = param)
    xHigh <- mode + sqrt(ghypVar(param = param))

    while (dghyp(xHigh, param = param) > tol) {
      xHigh <- xHigh + sqrt(ghypVar(param = param))
    }

    zeroFun <- function(x) {
      dghyp(x, param = param) - tol
    }

    xUpper <- uniroot(zeroFun, interval = c(mode, xHigh), ...)$root
    xLow <- mode - sqrt(ghypVar(param = param))

    while (dghyp(xLow, param = param) > tol) {
      xLow <- xLow - sqrt(ghypVar(param = param))
    }

    xLower <- uniroot(zeroFun, interval = c(xLow, mode), ...)$root
    range <- c(xLower, xUpper)
  }
  return(range)
} ## End of ghypCalcRange()
