###
### pimom.R
###

pimom <- function(q,V1=1,tau=1,nu=1) {

  z <- (.5*pgamma(tau*V1/q^2,nu/2,1))
  return(z*(q<=0)+(1-z)*(q>0))
}

