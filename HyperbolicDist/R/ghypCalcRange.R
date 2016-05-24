### Function to calculate an effective range for the distribution function
### or for the density function
### DJS 19/12/06
ghypCalcRange <- function(Theta, tol = 10^(-5), density = TRUE, ...) 
{
  Theta <- as.numeric(Theta)
  if(length(Theta)==4) Theta <- c(1,Theta)
  lambda <- Theta[1]
  alpha <- Theta[2]
  beta <- Theta[3]
  delta <- Theta[4]
  mu <- Theta[5]

  if (density == FALSE){
    ## bounds are for distribution function
    stop("Distribution function bounds not yet implemented")
  }else{
    ## bounds are for the density function
    mode <- ghypMode(Theta)
    xHigh <- mode + sqrt(ghypVar(Theta))
    while (dghyp(xHigh, Theta) > tol){
      xHigh <- xHigh + sqrt(ghypVar(Theta))
    }
    zeroFun<-function(x){ 
       dghyp(x, Theta) - tol 
     } 
    xUpper <- uniroot(zeroFun, interval = c(mode, xHigh), ...)$root
    xLow <- mode - sqrt(ghypVar(Theta))
    while (dghyp(xLow, Theta) > tol){
      xLow <- xLow - sqrt(ghypVar(Theta))
    }
    zeroFun<-function(x){ 
       dghyp(x, Theta) - tol 
     } 
    xLower <- uniroot(zeroFun, interval = c(xLow, mode), ...)$root
    range <- c(xLower, xUpper)
  }
  return(range)
} ## End of hyperbCalcRange()
