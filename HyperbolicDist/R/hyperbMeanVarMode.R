### Function to calculate the theoretical mean of a 
### hyperbolic distribution given its parameters.
hyperbMean <- function(Theta){
  Theta <- as.numeric(Theta)
  hyperbPi <- Theta[1]
  zeta <- Theta[2]
  delta <- Theta[3]
  mu <- Theta[4]
  mu + delta*hyperbPi*RLambda(zeta, lambda = 1)
} ## End of hyperbMean() 

### Function to calculate the theoretical variance of a 
### hyperbolic distribution given its parameters.
hyperbVar <- function(Theta){
  Theta <- as.numeric(Theta)
  hyperbPi <- Theta[1]
  zeta <- Theta[2]
  delta <- Theta[3]
  mu <- Theta[4]
  delta^2*(1/zeta*RLambda(zeta) + hyperbPi^2*SLambda(zeta))
} ## End of hyperbVar()

### Function to calculate the theoretical skewness of a 
### hyperbolic distribution given its parameters.
hyperbSkew <- function(Theta){
  Theta <- as.numeric(Theta)
  hyperbPi <- Theta[1]
  zeta  <- Theta[2]
  gammaLambda1(hyperbPi, zeta)
} ## End of hyperbSkew()

### Function to calculate the theoretical kurtosis of a 
### hyperbolic distribution given its parameters.
hyperbKurt <- function(Theta){
  Theta <- as.numeric(Theta)
  hyperbPi <- Theta[1]
  zeta  <- Theta[2]
  gammaLambda2(hyperbPi, zeta)
} ## End of hyperbKurt()


### Function to calculate the theoretical mode point of a 
### hyperbolic distribution given its parameters.
hyperbMode <- function(Theta){
  Theta <- as.numeric(Theta)
  hyperbPi <- Theta[1]
  zeta <- Theta[2]
  delta <- Theta[3]
  mu <- Theta[4]
  nu <- mu + delta*hyperbPi
  nu
} ## End of hyperbMode()
