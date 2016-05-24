# R code of package 'mnormpow'

pmnormpow <- function(x, varcov, ...)
    sadmvng(lower=rep(-Inf, length(x)), upper=x, varcov, ...) 

imnormpow <- function(lower,upper, varcov, ...)
    sadmvng(lower=lower, upper=upper, varcov, ...)

sadmvng <- function(lower, upper, varcov,
                    ipuiss=1, puiss=0,
                    maxpts=2000*d, abseps=1e-6, releps=0)
{
  if(any(lower > upper)) stop("lower>upper integration limits")
  if(any(lower == upper)) return(0)
  d <- as.integer(if(is.matrix(varcov)) ncol(varcov) else 1)
  varcov <- matrix(varcov, d, d)
  sd  <- sqrt(diag(varcov))
  rho <- cov2cor(varcov)
  lower <- lower/sd
  upper <- upper/sd
  coeff<-(sd[ipuiss])^puiss
  if(d == 1) return(pnorm(upper)-pnorm(lower))
  infin <- rep(2,d)
  infin <- replace(infin, (upper == Inf) & (lower > -Inf), 1)
  infin <- replace(infin, (upper < Inf) & (lower == -Inf), 0)
  infin <- replace(infin, (upper == Inf) & (lower == -Inf), -1)
  infin <- as.integer(infin)
  if(any(infin == -1)) {
    if(all(infin == -1)) return(1)
    k <- which(infin != -1)
    d <- length(k)
    lower <- lower[k]
    upper <- upper[k]
    if(d == 1) return(pnorm(upper) - pnorm(lower))
    rho <- rho[k, k]
    infin <- infin[k]
    }
  lower <- replace(lower, lower == -Inf, 0)
  upper <- replace(upper, upper == Inf, 0)
  lower <- as.double(lower)
  upper <- as.double(upper)
  correl <- as.double(rho[upper.tri(rho, diag=FALSE)])
  maxpts <- as.integer(maxpts)
  abseps <- as.double(abseps/coeff)
  releps <- as.double(releps)
  ipuiss <- as.integer(ipuiss)
  puiss <- as.integer(puiss)
  error  <- as.double(0)
  value  <- as.double(0)
  inform <- as.integer(0)
  result <- .Fortran("sadmvng", d, ipuiss, puiss, lower, upper, infin, correl, maxpts,
               abseps, releps, error, value, inform, PACKAGE="mnormpow")
  prob <- coeff*result[[12]]
  attr(prob,"error")  <- coeff*result[[11]]
  attr(prob,"status") <- switch(1+result[[13]], 
                                "normal completion", "accuracy non achieved", "oversize")
  return(prob)
}

#----

.onLoad <- .First.lib <- function(library, pkg)
{ 
   Rv <- R.Version()
   library.dynam("mnormpow", pkg, library)
   invisible()
}



