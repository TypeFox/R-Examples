### Function to calculate the theoretical raw moments (moments about 0)
### of a generalized inverse Gaussian distribution given its parameters.
gigRawMom <- function(order, chi = 1, psi = 1, lambda = 1,
                      param = c(chi, psi, lambda)) {

  param <- as.numeric(param)

  parResult <- gigCheckPars(param)
  case <- parResult$case
  errMessage <- parResult$errMessage

  if (case == "error")
    stop(errMessage)

  if (case =="normal") {
    chi <- param[1]
    psi <- param[2]
    lambda <- param[3]
    omega <- sqrt(chi * psi)
    eta <- sqrt(chi / psi)
    mom <- eta^order * besselRatio(omega, lambda, order)
  }

  if (case == "gamma") {
    shape <- param[3]
    scale <- param[2] / 2
    mom <- gammaRawMom(order, shape, scale)
  }

  if (case == "invgamma") {
    shape <- -param[3]
    scale <- 2 / param[1]
    mom <- gammaRawMom(-order, shape, scale)
  }

  return(mom)
} ## End of gigRawMom()


gigMom <- function(order, chi = 1, psi = 1, lambda = 1,
                   param = c(chi, psi, lambda), about = 0) {

  if (about != 0 & !is.wholenumber(order) & length(order) == 1)
    stop("Order must be a whole number except for moments about 0")

  if (order < 0 & about != 0)
    stop("Order must be positive except for moments about 0")

  param <- as.numeric(param)
  chi <- param[1]
  psi <- param[2]
  lambda <- param[3]
  omega <- sqrt(chi * psi)
  eta <- sqrt(chi / psi)

  if (about == 0) {
    mom <- gigRawMom(order, param = param)
  } else {
    if (order == 0) {
      mom <- 1
    } else {
      rawMoments <- sapply(1:order, gigRawMom, param = param)
      mom <- momChangeAbout(order, rawMoments, 0, about)
    }
  }

  return(mom)
} ## End of gigMom()


### Calculate moments of special cases
### chi = 0, gamma; psi = 0, inverse gamma
### Inverse gamma raw moments can be obtained from gamma raw moments
### An alternative would be to use mgamma and minvgamma from actuar
gammaRawMom <- function(order, shape = 1, rate = 1, scale = 1/rate) {

  if (order <= -shape) {
    mom <- Inf
  } else {
    mom <- (scale^order) * gamma(shape + order)/gamma(shape)
  }

  return(mom)
} ## End of gammaRawMom()
