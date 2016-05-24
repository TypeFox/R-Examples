##### For BPSM to work, we need to specify a JAGS model first

### Declare global varaibles so that the JAGS model can pass the package checking 
utils::globalVariables(c("M", "N", "N1", "N0", "Y", "Y0", "Y1", "X", "logit<-", "theta", "bet", "pow", "sigma", "zeta" ))

modelpsm <- function(){
  for ( i in 1 : N ) {
    logit ( p[i] ) <- X[i,] %*% theta
    phat[i] <- max ( .00001, min(.99999, p[i]) )
    t[i] ~ dbern ( phat[i] )
  }
  
  for ( i in 1 : N1 ) {
    d[i, 1:N0] <- abs ( phat[(N1+1):N] - phat[i] )
  }
  
  for ( i in 1 : N1 ) {
    r[i, 1:N0] <- rank ( d[i, 1:N0] ) <= M
    summ[i] <- r[i, 1:N0] %*% Y0
    Y0h[i] <- summ[i] / M  
    Y[i] ~ dnorm( mu[i], tau[i] )
    mu[i] <- Y0h[i] + bet[i]      
  }
  
  for ( i in ( N1 + 1 ) : N ) {
    r[i, 1:N1] <- rank ( d[1:N1, i - N1] ) <= M
    summ[i] <- r[i, 1:N1] %*% Y1
    Y1h[i] <- summ[i] / M  
    Y[i] ~ dnorm ( mu[i], tau[i] )
    mu[i] <- Y1h[i] - bet[i]
  }
  
  for (d in 1:D) { 
    theta[d] ~ dnorm ( ET[d], 1E-5 )
  }
  
  for (i in 1:N) { 
    bet[i] ~ dnorm ( ET[D+1], 1E-5 )
    tau[i] <- pow(sigma[i], -2)
    sigma[i] ~ dunif(0, 10)
  }
  
  ate <- mean(bet)
  att <- mean(bet[1:N1])
  atc <- mean(bet[(N1+1):N])
}