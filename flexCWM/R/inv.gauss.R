####################################################
## Re-parameterized Inverse Gaussian distribution ##
####################################################

library(statmod)

dig <- function(x,mu,var,log=FALSE){
  
  dinvgauss(x, mean=mu, dispersion=var/mu^3, log=log)
  
}

rig <- function(n,mu,var){
  
  statmod::rinvgauss(n, mean=mu, dispersion=var/mu^3)
  
}