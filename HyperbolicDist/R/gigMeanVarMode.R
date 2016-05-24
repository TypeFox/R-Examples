### Function to calculate the theoretical mean of a 
### generalized inverse Gaussian distribution given its parameters.
gigMean <- function (Theta) {
  Theta <- as.numeric(Theta)
  lambda <- Theta[1]
  chi <- Theta[2]
  psi <- Theta[3]
  omega <- sqrt(chi*psi)
  eta <- sqrt(chi/psi)
  eta*besselRatio(omega, lambda, 1)
}## End of gigMean() 

### Function to calculate the theoretical variance of a 
### generalized inverse Gaussian distribution given its parameters.
gigVar <- function(Theta){
  m1 <- gigMean(Theta)
  var <- gigMom(2, Theta, about = m1)
  return(var)
} ## End of gigVar()

### Function to calculate the theoretical skewness of a 
### generalized inverse Gaussian distribution given its parameters.
gigSkew <- function(Theta){
  m1 <- gigMean(Theta)
  skew <- gigMom(3, Theta, about = m1)/(gigVar(Theta)^(3/2))
  return(skew)
} ## End of gigSkew()

### Function to calculate the theoretical kurtosis of a 
### generalized inverse Gaussian distribution given its parameters.
gigKurt <- function(Theta){
  m1 <- gigMean(Theta)
  kurt <- gigMom(4, Theta, about = m1)/(gigVar(Theta)^2) - 3
  return(kurt)
} ## End of gigKurt()


### Function to calculate the theoretical mode point of a 
### generalized inverse Gaussian distribution given its parameters.
gigMode <- function(Theta) {
  Theta <- as.numeric(Theta)
  lambda <- Theta[1]
  chi <- Theta[2]
  psi <- Theta[3]
  (lambda - 1 + sqrt((lambda - 1)^2 + chi*psi))/psi
} ## End of gigMode()
