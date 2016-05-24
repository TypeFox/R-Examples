### Function to calculate an effective range for the distribution function
### or for the density function
### DJS 10/01/07
gigCalcRange <- function(chi = 1, psi = 1, lambda = 1,
                         param = c(chi, psi, lambda),
                         tol = 10^(-5), density = TRUE, ...) {

  param <- as.numeric(param)

  ## check parameters
  parResult <- gigCheckPars(param)
  case <- parResult$case
  errMessage <- parResult$errMessage

  if (case == "error")
    stop(errMessage)

  chi <- param[1]
  psi <- param[2]
  lambda <- param[3]
  KLambda <- besselK(x = sqrt(chi * psi), nu = lambda)

  if (!density) {
    ## bounds are for distribution function
    stop("Distribution function bounds not yet implemented")
  } else {
    ## bounds are for the density function
    mode <- gigMode(param = param)
    xHigh <- mode + sqrt(gigVar(param = param))

    while (dgig(xHigh, param = param) > tol) {
      xHigh <- xHigh + sqrt(gigVar(param = param))
    }

    zeroFun <- function(x) {
       dgig(x, param = param) - tol
    }

    xUpper <- uniroot(zeroFun, interval = c(mode, xHigh), ...)$root
    xLow <- 0
    xLower <- uniroot(zeroFun, interval = c(xLow, mode), ...)$root
    range <- c(xLower, xUpper)
  }

  return(range)
} ## End of gigCalcRange()
