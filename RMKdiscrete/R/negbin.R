dnegbin <- function(x,nu,p,mu,log=FALSE){
  if( sum(c(missing(nu),missing(p),missing(mu)))!=1 ){stop("exactly 2 of arguments 'nu', 'p', and 'mu' must be provided")}
  if(missing(p)){return(dnbinom(x=x,size=nu,mu=mu,log=log))}
  if(missing(nu)){nu <- mu*p/(1-p)}  #<--Might not be a smart way to compute nu, if p is very close to 0 or 1.
  return(dnbinom(x=x,size=nu,prob=p,log=log))
}

negbinMVP <- function(nu, p, mu, sigma2){
  if( sum(missing(mu),missing(sigma2),missing(nu),missing(p))!=2 ){
    stop("exactly 2 of 4 arguments must be provided")
  }
  if(missing(nu) & missing(p)){
    p <- mu/sigma2
    nu <- mu*p/(1-p)
    out <- cbind(nu=nu, p=p)
  }
  if(missing(mu) & missing(sigma2)){
    mu <- nu*(1-p)/p
    sigma2 <- nu*(1-p)/p/p
    out <- cbind(mu=mu, sigma2=sigma2)
  }
  if(missing(p) & missing(sigma2)){
    sigma2 <- mu+((mu^2)/nu)
    p <- nu/(nu+mu)
    out <- cbind(sigma2=sigma2, p=p)
  }
  if(missing(mu) & missing(nu)){
    mu <- p*sigma2
    nu <- mu*p/(1-p)
    out <- cbind(mu=mu, nu=nu)
  }
  if(missing(mu) & missing(p)){stop("no solution when 'mu' and 'p' are both missing")}
  if(missing(sigma2) & missing(nu)){
    sigma2 <- mu/p
    nu <- mu*p/(1-p)
    out <- cbind(sigma2=sigma2, nu=nu)
  }
  if(any(nu<0) | any(mu<0) | any(p<=0) | any(p>1) | any(!is.finite(out))){
    warning("some input values appear to be outside parameter space")
  }
  return(out)
}
#
