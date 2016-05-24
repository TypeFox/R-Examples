##### For BPSR to work, we need to specify a JAGS model first

modelpsr <- function(){
  for (i in 1:N) {
    logit(p[i]) <- X[i,] %*% theta
    phat[i] <- max(0.00001, min(0.99999, p[i]))
    t[i] ~ dbern(phat[i])
    
    
    Y[i] ~ dnorm(mu[i], tau)
    mu[i] <-  beta * t[i] +  zeta * phat[i]
  }
  
  for (d in 1:D) { 
    theta[d] ~ dnorm(ET[d], 1E-5)
  }
  
  beta ~ dnorm(ET[D+1], 1E-5)
  zeta ~ dnorm(ET[D+2], 1E-5)
  tau <- pow(sigma, -2)
  sigma ~ dunif(0, 100)
}
