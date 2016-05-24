### Function to calculate an effective range for the distribution function
### or for the density function
### DJS 10/01/07
gigCalcRange <- function (Theta, tol = 10^(-5), density = TRUE, ...) 
{
  Theta <- as.numeric(Theta)
  lambda <- Theta[1]
  chi <- Theta[2]
  psi <- Theta[3]
  KLambda <- besselK(x = sqrt(chi*psi), nu = lambda)

  if (density == FALSE){
    ## bounds are for distribution function
    stop("Distribution function bounds not yet implemented")
  }else{
    ## bounds are for the density function
    mode <- gigMode(Theta)
    xHigh <- mode + sqrt(gigVar(Theta))
    while (dgig(xHigh, Theta) > tol){
      xHigh <- xHigh + sqrt(gigVar(Theta))
    }
    zeroFun<-function(x){ 
       dgig(x, Theta) - tol 
     } 
    xUpper <- uniroot(zeroFun, interval = c(mode, xHigh), ...)$root
    xLow <- 0
    zeroFun<-function(x){ 
       dgig(x, Theta) - tol 
     } 
    xLower <- uniroot(zeroFun, interval = c(xLow, mode), ...)$root
    range <- c(xLower, xUpper)
  }
  return(range)
} ## End of gigCalcRange()
