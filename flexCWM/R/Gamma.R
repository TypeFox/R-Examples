#########################################
## Re-parameterized Gamma distribution ##
#########################################

dgam <- function(x,mu,nu,log=FALSE){
  dgamma(x, shape=nu, scale = mu/nu, log = log)
}

rgam <- function(x,mu,nu){
  rgamma(x, shape=nu, scale = mu/nu)
}