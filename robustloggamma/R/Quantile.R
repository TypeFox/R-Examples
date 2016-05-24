### DERIVATA DELLA FUNZIONE QUANTILE RISPETTO A MU
d.Q1 <- function(p, mu, sigma, lambda) {
  return(1)  
}

### DERIVATA DELLA FUNZIONE QUANTILE RISPETTO A SIGMA
d.Q2 <- function(p, mu, sigma, lambda) {
  zero <- 0.001
  if (lambda < -zero) 
    p <- 1 - p
  if (abs(lambda) < zero) {
    res <- qnorm(p)
  } else {
    res <- (log(qgamma(p=p, shape=lambda^(-2), rate=1)) + log(lambda^2))/lambda
  }
  return(res)
}

### DERIVATA DELLA FUNZIONE QUANTILE RISPETTO A LAMBDA
d.Q3  <- function(p, lambda) {
  Q <- qloggamma(p, mu=0, sigma=1, lambda=lambda)
  if (!is.nan(Q))
    dQ <- -d.FL(Q, lambda=lambda)/dloggamma(Q, mu=0, sigma=1, lambda=lambda)
  else
    dQ <- NaN
  return(dQ)
}

## Questa funzione presenta instabilita' numerica per lambda piccoli!!
## d.Q3 <- function(p, mu, sigma, lambda) {
##   zero <- 1e-04
##   if (lambda < -zero) 
##     p <- 1 - p
##   if (abs(lambda) < zero) {
##     res <- NA
##   } else {
##     qg <- qgamma(p=p, shape=lambda^(-2), rate=1)
##      myenv <- new.env()    
##     temp <- function(x) qgamma(p=p, shape=x^(-2), rate=1)
##     assign("x", lambda, envir = myenv)
##     dqg <- drop(attr(numericDeriv(quote(temp(x)), "x", myenv), "gradient"))
##     res <- sigma*(dqg/qg*lambda - log(qg) + 2 - log(lambda^2))/lambda^2    
##   }
##   return(res)
## }

######### Per vedere i problemi
#####plot(Vectorize(function(x) d.Q3(0.25, 0, 1, x)), from=-0.01, to=0.01, n=10000)
