### Function to calculate the theoretical mean of a
### hyperbolic distribution given its parameters.
hyperbMean <- function(mu = 0, delta = 1, alpha = 1, beta = 0,
                       param = c(mu, delta, alpha, beta)) {

  param <- as.numeric(param)
  ghypMean(param = c(param, 1))
} ## End of hyperbMean()

### Function to calculate the theoretical variance of a
### hyperbolic distribution given its parameters.
hyperbVar <- function(mu = 0, delta = 1, alpha = 1, beta = 0,
                      param = c(mu, delta, alpha, beta)) {

  param <- as.numeric(param)
  ghypVar(param = c(param, 1))
} ## End of hyperbVar()

### Function to calculate the theoretical skewness of a
### hyperbolic distribution given its parameters.
hyperbSkew <- function(mu = 0, delta = 1, alpha = 1, beta = 0,
                       param = c(mu, delta, alpha, beta)) {

  param <- as.numeric(param)
  ghypSkew(param = c(param, 1))
} ## End of hyperbSkew()

### Function to calculate the theoretical kurtosis of a
### hyperbolic distribution given its parameters.
hyperbKurt <- function(mu = 0, delta = 1, alpha = 1, beta = 0,
                       param = c(mu, delta, alpha, beta)) {

  param <- as.numeric(param)
  ghypKurt(param = c(param, 1))
} ## End of hyperbKurt()


### Function to calculate the theoretical mode point of a
### hyperbolic distribution given its parameters.
hyperbMode <- function(mu = 0, delta = 1, alpha = 1, beta = 0,
                       param = c(mu, delta, alpha, beta)) {

  param <- as.numeric(hyperbChangePars(2, 1, param = param))
  mu <- param[1]
  hyperbPi <- param[3]
  nu <- mu + delta*hyperbPi
  return(nu)
} ## End of hyperbMode()
