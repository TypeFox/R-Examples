lsdnorm <- function(x, mean=0, sd=1, tau=0.2) {
  if (tau >= 1 | tau <= 0)
    stop("quantile order of 'tau' must be in the range (0,1). For order equal to 1 use function 'sdnorm'")
  tau <- 1-(1-tau)/2
  tau <- qnorm(tau, 0, sqrt(2)*sd)
  ldnorm.int <- function(x, tau, mean, sd) {
    pdf <- function(x) dnorm(x,mean,sd)
    cdf <- function(x) pnorm(x,mean,sd)
    acca <- function(x) cdf(x+tau)+cdf(x-tau)-2*cdf(x)    
    int1 <- function(t) cdf(t+tau)*pdf(t)
    integr1 <- function(x) integrate(int1,x-tau,x)[1]
    int2 <- function(t) cdf(t-tau)*pdf(t)
    integr2 <- function(x) integrate(int2,x,x+tau)[1]
    a <- integr1(x); b <- integr2(x) 
    res <- a$value-b$value+cdf(x)*acca(x)
    return(res)
  }
  y <- sapply(X=x, FUN=ldnorm.int, simplify=TRUE, tau=tau, mean=mean, sd=sd)
  return(y)
}

hsnorm <- function(x, mean=0, sd=1, tau=0.2) {
  if (tau >= 1 | tau <= 0)
    stop("quantile order of 'tau' must be in the range (0,1)")
  tau <- 1-(1-tau)/2
  tau <- qnorm(tau, 0, sqrt(2)*sd)
  acca <- pnorm(x+tau, mean, sd)+pnorm(x-tau, mean, sd)-2*pnorm(x, mean, sd)  
  return(acca)
}

sdnorm <- function(x, mean=0, sd=1) {
  y <- 2*pnorm(x, mean, sd)*pnorm(x, mean, sd, lower.tail=FALSE)
  return(y)
}
