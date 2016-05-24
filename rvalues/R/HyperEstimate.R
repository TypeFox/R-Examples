HyperEstimate <- function(estimate, nuisance, family) {
    ##  function to estimate the hyperparameters for the parametric models.

    if(family=="gaussian")  {
        mu <- mean(estimate)
        vv <- nuisance^2
        fit <- nlminb( start=mean(vv) , objective=nnnll,gradient=nnnll.g,hessian=nnnll.hess, 
                     lower=0, x=(estimate-mu), sigma2=vv )
    
        tau2 <- fit$par  ## MLE of prior variance, pre-scaling
        hypers <- c(mu,tau2)
    }
    else if(family=="binomial") {
        ## nuisance = nn
        ## esitmate = xx
        
        ok <- (nuisance>0)
        tmp1 <- mean( (estimate/nuisance)[ok] )
        tmp2 <- var( (estimate/nuisance)[ok] )
        den <- tmp1*(1-tmp1)/tmp2
    
        a0 <- den*tmp1
        b0 <- den - a0
        ## get MLE
        fit <- nlminb( start=c(a0,b0), objective=bbnll, 
                       x=estimate[ok], n=nuisance[ok] , lower=c(1,1) )
    
        hypers <- c(fit$par[1],fit$par[2])
    }
    else if(family=="poisson")  {
        ## nuisance = eta
        ## estimate = xx
        
        a0 <- 2
        b0 <- a0/( sum(estimate)/sum(nuisance) )   ## simple method of moments starter
        ## get MLE
        fit <- nlminb( start=c(a0,b0), objective=pgnll, gradient=pgnll.g,
                       x=estimate, eta=nuisance , lower=c(.1,.1) )  ## a constraint, could be lifted.
        ### maybe think about one-dimensional optimization with profile likelihood.
        hypers <- c(fit$par[1],fit$par[2])    
    }
    else if(family=="Gamma") {
        a0 <- 4
        b0 <- mean(nuisance)/((a0 - 1)*mean(estimate))
        b0 <- 3
        fit <- nlminb( start=c(a0,b0), objective=ggnll, gradient=ggnll.g,
                       x = estimate, shapes = nuisance, lower=c(.2,.2))
        hypers <- c(fit$par[1], fit$par[2])
    }
    return(hypers)
}


bbnll <- function(shapes, x, n)
{
  ## neg log likelihood; beta binomial
  ## log bb density; not counting choose(n,x) term
  a <- shapes[1]
  b <- shapes[2]
  tt <- lgamma(a+b) - lgamma(a+b+n)+lgamma(a+x)+lgamma(b+n-x)-
    lgamma(a) - lgamma(b)
  return( -sum(tt) )
}


bbnll.g <- function(shapes,x,n) {
  a <- shapes[1]
  b <- shapes[2]
  tta <- digamma(a+b) - digamma(a+b+n)+digamma(a+x) -digamma(a)
  ttb <- digamma(a+b) - digamma(a+b+n)+digamma(b+n-x)-digamma(b)
  return(c(tta,ttb))
}

pgnll <- function(shapes, x, eta)
{
  ## neg log likelihood; poisson gamma
  a <- shapes[1]
  b <- shapes[2]
  u <- b/(b+eta)
  tt <-  a*log(u) + x*log(1-u) + lgamma(a+x) -lgamma(x+1) -lgamma(a)
  return( -sum(tt) )
}


pgnll.g <- function(shapes, x, eta)
{
  ## neg log likelihood; poisson gamma
  a <- shapes[1]
  b <- shapes[2]
  u <- b/(b+eta)
  dudb <- eta/((eta + b)^2)
  tta <-  log(u) + digamma(a+x) - digamma(a)
  ttb <- (a/u - x/(1-u))*dudb
  return( c(-sum(tta),-sum(ttb)) )
}


nnnll <- function( tau2, x, sigma2 )
{
  ## neg log likelihood for prior variance tau2;  normal/normal
  ss <- tau2 + sigma2
  tmp <- sum( log(ss) ) + sum( (x^2)/ss )
  nll <- (1/2)*tmp
  nll
}

nnnll.g <- function(tau2,x,sigma2) {
  ### gradient function associated with nnnll
  ss <- tau2 + sigma2
  ans <- sum(1/ss) - sum((x/ss)^2)
  return(ans/2)
}

nnnll.hess <- function(tau2,x,sigma2) {
  ### hessian associated with nnnll
  ss <- tau2 + sigma2
  ans <- -sum(1/(ss^2)) + 2*sum((x^2)/(ss^3))
  return(as.matrix(ans/2))
}

ggnll <- function(pars, x, shapes) {
  ### negative log-likelihood gamma-gamma
  ### assumes that X|theta,shapes ~ Gamma(shapes, theta)
  ### and theta ~ InvGamma(a, b)
 
  a <- pars[1]
  b <- pars[2]
  n <- length(x)
  
  ans <- sum((a + shapes)*log(b + x)) + n*lgamma(a) - n*a*log(b) - sum(lgamma(shapes + a))
  return(ans)
}

ggnll.g <- function(pars, x, shapes) {
  ### gradient associated with ggnll
 
  a <- pars[1]
  b <- pars[2]
  n <- length(x)

  a.partial <- sum(log(b + x)) + n*digamma(a) - n*log(b) - sum(digamma(shapes + a))
  b.partial <- sum((a + shapes)/(b + x)) - (n*a)/b 
   
  return(c(a.partial, b.partial))
}






