### Function to calculate the theoretical mean of a 
### generalized hyperbolic distribution given its parameters.
ghypMean <- function(Theta){
  Theta <- as.numeric(Theta)
  if(length(Theta)==4) Theta <- c(1,Theta)
  lambda <- Theta[1]
  alpha <- Theta[2]
  beta <- Theta[3]
  delta <- Theta[4]
  mu <- Theta[5]

  gamma <- sqrt(alpha^2 - beta^2)
  
  mu + delta*beta*besselRatio(delta*gamma, lambda, 1)/gamma
} ## End of ghypMean() 

### Function to calculate the theoretical variance of a 
### generalized hyperbolic distribution given its parameters.
ghypVar <- function(Theta){
  var <- ghypMom(2, Theta, momType = "central")
  return(var)
} ## End of ghypVar()

### Function to calculate the theoretical skewness of a 
### generalized hyperbolic distribution given its parameters.
ghypSkew <- function(Theta){
  skew <- ghypMom(3, Theta, momType = "central")/(ghypVar(Theta)^(3/2))
  return(skew)
} ## End of ghypSkew()

### Function to calculate the theoretical kurtosis of a 
### generalized hyperbolic distribution given its parameters.
ghypKurt <- function(Theta){
  kurt <- ghypMom(4, Theta, momType = "central")/(ghypVar(Theta)^2) - 3
  return(kurt)
} ## End of ghypKurt()

### Function to calculate the theoretical mode point of a 
### generalized hyperbolic distribution given its parameters.
ghypMode <- function(Theta){
  Theta <- as.numeric(Theta)
  if(length(Theta)==4) Theta <- c(1,Theta)
  lambda <- Theta[1]
  alpha <- Theta[2]
  beta <- Theta[3]
  delta <- Theta[4]
  mu <- Theta[5]

  modeFun <- function(x){
    log(dghyp(x, Theta))
  }
  start <- ghypMean(Theta)
  optResult <- optim(start, modeFun,
                     control = list(fnscale = -1, maxit = 1000),
                     method = "BFGS")
  if (optResult$convergence == 0){
    mode <- optResult$par
  }else{
    mode <- NA
  }
  mode
} ## End of ghypMode()

