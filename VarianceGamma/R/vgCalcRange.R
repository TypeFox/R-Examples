vgCalcRange <- function (vgC = 0, sigma = 1, theta = 0, nu = 1,
  param = c(vgC,sigma,theta,nu), tol = 10^(-5), density = TRUE, ...) {

  #check parameters
  parResult <- vgCheckPars(param = param)
  case <- parResult$case
  errMessage <- parResult$errMessage
  if (case == "error"){
    stop(errMessage)
  }
  vgC <- param[1]
  sigma <- param[2]
  theta <- param[3]
  nu <- param[4]

  if (density == FALSE) {
    stop("Distribution function bounds not yet implemented")
  }

  if (nu < 2) {
    modeDist <- vgMode(param = param)
    sDev <- sqrt(sigma^2 + theta^2*nu)
    xHigh <- modeDist + sDev
    while (dvg(x = xHigh, param = param, log = FALSE) > tol) {
      xHigh <- xHigh + sDev
    }
    zeroFun <- function (x) {
      dvg(x = x, param = param, log = FALSE) - tol
    }
    xUpper <- uniroot(zeroFun, interval = c(modeDist,xHigh), ...)$root

    xLow <- modeDist - sDev
    while (dvg(x = xLow, param = param, log = FALSE) > tol) {
      xLow <- xLow - sDev
    }
    zeroFun <- function (x) {
      dvg(x = x, param = param, log = FALSE) - tol
    }
    xLower <- uniroot(zeroFun, interval = c(xLow,modeDist), ...)$root
    range <- c(xLower,xUpper)
  }

  if (nu >= 2) {
    modeDist <- vgMode(param = param)
    sDev <- sqrt(sigma^2 + theta^2*nu)
    xHigh <- modeDist + sDev
    while (dvg(x = xHigh, param = param, log = FALSE) > tol) {
      xHigh <- xHigh + sDev
    }
    zeroFun <- function (x) {
      dvg(x = x, param = param, log = FALSE) - tol
    }
    xUpper <- uniroot(zeroFun, interval = c(modeDist + 0.001,xHigh), ...)$root

    xLow <- modeDist - sDev
    while (dvg(x = xLow, param = param, log = FALSE) > tol) {
      xLow <- xLow - sDev
    }
    zeroFun <- function (x) {
      dvg(x = x, param = param, log = FALSE) - tol
    }
    xLower <- uniroot(zeroFun, interval = c(xLow,modeDist - 0.001), ...)$root
    range <- c(xLower,xUpper)
    }
  return(range)
}
