### Function to calculate the theoretical mean of a
### generalized inverse Gaussian distribution given its parameters.
gigMean <- function(chi = 1, psi = 1, lambda = 1,
                    param = c(chi, psi, lambda)) {

  ## check parameters
  parResult <- gigCheckPars(param)
  case <- parResult$case
  errMessage <- parResult$errMessage

  if (case == "error")
    stop(errMessage)

  param <- as.numeric(param)
  chi <- param[1]
  psi <- param[2]
  lambda <- param[3]
  omega <- sqrt(chi * psi)
  eta <- sqrt(chi / psi)
  eta * besselRatio(omega, lambda, 1)
} ## End of gigMean()

### Function to calculate the theoretical variance of a
### generalized inverse Gaussian distribution given its parameters.
gigVar <- function(chi = 1, psi = 1, lambda = 1,
                   param = c(chi, psi, lambda)) {

  ## check parameters
  parResult <- gigCheckPars(param)
  case <- parResult$case
  errMessage <- parResult$errMessage

  if (case == "error")
    stop(errMessage)

  m1 <- gigMean(param = param)
  var <- gigMom(2, param = param, about = m1)
  return(var)
} ## End of gigVar()

### Function to calculate the theoretical skewness of a
### generalized inverse Gaussian distribution given its parameters.
gigSkew <- function(chi = 1, psi = 1, lambda = 1,
                    param = c(chi, psi, lambda)) {

  ## check parameters
  parResult <- gigCheckPars(param)
  case <- parResult$case
  errMessage <- parResult$errMessage

  if (case == "error")
    stop(errMessage)

  m1 <- gigMean(param = param)
  skew <- gigMom(3, param = param, about = m1) / (gigVar(param = param)^(3 / 2))
  return(skew)
} ## End of gigSkew()

### Function to calculate the theoretical kurtosis of a
### generalized inverse Gaussian distribution given its parameters.
gigKurt <- function(chi = 1, psi = 1, lambda = 1,
                    param = c(chi, psi, lambda)) {

  ## check parameters
  parResult <- gigCheckPars(param)
  case <- parResult$case
  errMessage <- parResult$errMessage

  if (case == "error")
    stop(errMessage)

  m1 <- gigMean(param = param)
  kurt <- gigMom(4, param = param, about = m1) / (gigVar(param = param)^2) - 3
  return(kurt)
} ## End of gigKurt()


### Function to calculate the theoretical mode point of a
### generalized inverse Gaussian distribution given its parameters.
gigMode <- function(chi = 1, psi = 1, lambda = 1,
                    param = c(chi, psi, lambda)) {

  ## check parameters
  parResult <- gigCheckPars(param)
  case <- parResult$case
  errMessage <- parResult$errMessage

  if (case == "error")
    stop(errMessage)

  param <- as.numeric(param)
  chi <- param[1]
  psi <- param[2]
  lambda <- param[3]
  (lambda - 1 + sqrt((lambda - 1)^2 + chi * psi)) / psi
} ## End of gigMode()
