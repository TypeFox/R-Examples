### Function to calculate the theoretical mean of a 
### normal inverse Gaussian distribution given its parameters.
nigMean <- function(mu = 0, delta = 1, alpha = 1, beta = 0,
                       param = c(mu, delta, alpha, beta)) {

  param <- as.numeric(param)
  ghypMean(param = c(param, -1/2))
} ## End of nigMean() 

### Function to calculate the theoretical variance of a 
### normal inverse Gaussian distribution given its parameters.
nigVar <- function(mu = 0, delta = 1, alpha = 1, beta = 0,
                      param = c(mu, delta, alpha, beta)) {

  param <- as.numeric(param)
  ghypVar(param = c(param, -1/2))
} ## End of nigVar()

### Function to calculate the theoretical skewness of a 
### normal inverse Gaussian distribution given its parameters.
nigSkew <- function(mu = 0, delta = 1, alpha = 1, beta = 0,
                       param = c(mu, delta, alpha, beta)) {

  param <- as.numeric(param)
  ghypSkew(param = c(param, -1/2))
} ## End of nigSkew()

### Function to calculate the theoretical kurtosis of a 
### normal inverse Gaussian distribution given its parameters.
nigKurt <- function(mu = 0, delta = 1, alpha = 1, beta = 0,  
                       param = c(mu, delta, alpha, beta)) {

  param <- as.numeric(param)
  ghypKurt(param = c(param, -1/2))
} ## End of nigKurt()


### Function to calculate the theoretical mode point of a 
### normal inverse Gaussian distribution given its parameters.
nigMode <- function(mu = 0, delta = 1, alpha = 1, beta = 0,
                       param = c(mu, delta, alpha, beta)) {

  param <- as.numeric(param)
  ghypMode(param = c(param, -1/2))
} ## End of nigMode()
