# Utility functions that are not exported
optimalCScalar2 <- function(lambda, alpha, min.loc, n, p, k, var.known) {
  out <- uniroot(f=eqnForMin, interval=c(0,500), lambda=lambda,
    alpha=alpha, min.loc=min.loc, n=n, p=p, k=k, var.known=var.known)$root

  return(out)
}

# For finding optimal lambda
h2 <- function(lambda, alpha, n, p, k, var.known) {
abs(optimalC(lambda, alpha, 'zero', n, p, k, var.known) - 
    optimalC(lambda, alpha, 'infty', n, p, k, var.known))
}
h2.v <- Vectorize(h2, vectorize.args="lambda")

####### Unknown variance case, k=1 #########################
# For the unknown variance, single p case
component1a <- function(x, c.val, lambda, p)
  pnorm(c.val*x)^p - pnorm(-lambda*c.val*x)^p
component1b <- function(x, c.val, lambda)
  pnorm(c.val*x) - pnorm(-lambda*c.val*x)

component2 <- function(x, c.val, lambda, n, p) {
  dgamma(x, shape=0.5*p*(n-1), scale=2/(p*(n-1)))
}

component3 <- function(x, c.val, lambda, n, p) {
  log.scale <- 0.5*p*(n-1)*log(x) + 0.5*x*p*(n-1) - 0.5*p*(n-1)*x^2
  2*exp(log.scale)
}

#2a corresponds to 'zero'
#2b corresponds to 'infty'
integrand2a <- function(x, c.val, lambda, n, p) {
  component1a(x,c.val,lambda,p)*
  component2(x,c.val,lambda,n,p)*
  component3(x,c.val,lambda,n,p)
}
integrand2b <- function(x, c.val, lambda, n, p) {
  component1b(x,c.val,lambda)*
  component2(x,c.val,lambda,n,p)*
  component3(x,c.val,lambda,n,p)
}

####### Unknown variance case, k>1 #########################
# For the unknown variance, multiple k case

#4a corresponds to 'zero'
#4b corresponds to 'infty'
integrand4a <- function(x, c.val, lambda, n, p, k) {
  component1a(x,c.val,lambda,p-k+1)*
  component2(x,c.val,lambda,n,p)*
  component3(x,c.val,lambda,n,p)*
  component1b(x,c.val,lambda)^(k-1)
}
integrand4b <- function(x, c.val, lambda, n, p, k) {
  (component1b(x,c.val,lambda)^k)*
  component2(x,c.val,lambda,n,p)*
  component3(x,c.val,lambda,n,p)
}

###### integration functions ###########################

## Case k=1, min. at 0
integrate2a <- function(lambda, c.val, n, p) {
  integrate(integrand2a, 0, Inf, c.val=c.val, lambda=lambda, n=n, p=p)$value 
}

## Case k=1, min. at infty 
integrate2b <- function(lambda, c.val, n, p) {
  integrate(integrand2b, 0, Inf, c.val=c.val, lambda=lambda, n=n, p=p)$value 
}

## Case k>1, min. at zero
integrate4a <- function(lambda, c.val, n, p, k) {
  integrate(integrand4a, 0, Inf, c.val=c.val, lambda=lambda, n=n, p=p, k=k)$value 
}

## Case k>1, min. at infty
integrate4b <- function(lambda, c.val, n, p, k) {
  integrate(integrand4b, 0, Inf, c.val=c.val, lambda=lambda, n=n, p=p, k=k)$value 
}

####### statistic() function for bootstrapping ############
computeMax <- function(d, i, k=1) {
  d.2 <- d$values[i]
  sample.means <- tapply(d.2, d$ind, mean)
  sample.means <- sort(sample.means, decreasing=TRUE)
  sample.means[1:k]
}

