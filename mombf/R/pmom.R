###
### pmom.R
###

pmom <- function(q,V1=1,tau=1) {

  z <- .5-(pnorm(abs(q)/sqrt(V1*tau)) - abs(q)/sqrt(2*pi*V1*tau) * exp(-.5*q^2/(tau*V1)) - .5)
  return(z*(q<=0)+(1-z)*(q>0))
}

