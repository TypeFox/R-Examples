### Function to calculate an effective range for the distribution function
### or for the density function of the Hyperbolic Distribution
### DJS 8/09/06
hyperbCalcRange <- function(Theta, tol = 10^(-5), density = FALSE) 
{
  Theta <- as.numeric(Theta)
  hyperbPi <- Theta[1]
  zeta <- Theta[2]
  delta <- Theta[3]
  mu <- Theta[4]
  KNu <- besselK(zeta, nu = 1)
  phi <- as.numeric(hyperbChangePars(1, 3, Theta)[1])
  gamma <- as.numeric(hyperbChangePars(1, 3, Theta)[2])
  alpha <- (phi + gamma)/2
  const <- 1/(2*delta*(1 + hyperbPi^2)^(1/2)*KNu)
  ## this shouldn't make a difference but is put in
  ## to follow the theory
  tol <- min(tol, const/gamma, const/phi)
  if (density == FALSE){
    ## bounds are for distribution function
    xLower <- min(mu - (1/phi)*log(const/(tol*phi)), mu - delta)
    xUpper <- max(mu + (1/gamma)*log(const/(tol*gamma)), mu + delta) 
    range <- c(xLower, xUpper)
  }else{
    ## bounds are for the density function
    xLower <- min(mu - 1/phi*log(const/tol), mu - delta) 
    xUpper <- max(mu + 1/gamma*log(const/tol), mu + delta)
    range <- c(xLower, xUpper)
  }
  return(range)
} ## End of hyperbCalcRange()
