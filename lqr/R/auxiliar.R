# require(ghyp)
# require(LaplacesDemon)
# require(spatstat)

# GENERATE FROM SKD

genSLK = function(n,mu=0,sigma=1,p=0.5,dist = "normal",nu="",gama="")
{
  if(dist == ""){dist = "normal"}
  if(nu=="" && dist == "t"){nu=4}
  if(nu=="" && dist == "slash"){nu=2}
  if(nu=="" && dist == "cont"){nu=0.1}
  if(gama=="" && dist == "cont"){gama=0.1}
  
  rhn = abs(rnorm(n,mean = 0,sd = 1))
  
  if(dist == "normal"){K = 1}
  if(dist == "t"){K = rinvgamma(n,shape = nu/2,scale = nu/2)}
  if(dist == "laplace"){K = rexp(n,rate = 1/2)}
  if(dist == "slash"){K = 1/rbeta(n,nu,1)}
  if(dist == "cont"){K = 1/sample(x = c(gama,1),size = n,replace = T,prob = c(nu,1-nu))}
  
  I = ifelse(runif(n)<p,-(1/(2*(1-p))),1/(2*p))
  y = I*rhn
  return(mu + sigma*sqrt(K)*y)
}

# GENERATE FROM SKD DISTANCE
gendistSLK = function(n,dist = "normal",nu="",gama="")
{
if(dist == ""){dist = "normal"}
if(nu=="" && dist == "t"){nu=4}
if(nu=="" && dist == "slash"){nu=2}
if(nu=="" && dist == "cont"){nu=0.1}
if(gama=="" && dist == "cont"){gama=0.1}

if(dist == "normal"){K = 1}
if(dist == "t"){K = rinvgamma(n,shape = nu/2,scale = nu/2)}
if(dist == "laplace"){K = rexp(n,rate = 1/2)}
if(dist == "slash"){K = 1/rbeta(n,nu,1)}
if(dist == "cont"){K = 1/sample(x = c(gama,1),size = n,replace = T,prob = c(nu,1-nu))}

return(0.5*sqrt(K)*abs(rnorm(n,mean = 0,sd = 1)))
}

# GENERATE SIGNIFICANCE INDEX

defast = function(x){
  if(x>0.1){ast = " "}else
  {
    if(x>0.05){ast = "."}else
    {
      if(x>0.01){ast = "*"}else
      {
        if(x>0.001){ast = "**"}else
        {
{
  ast = "***"
}
        }
      }
    }
  }
return(ast)
}

#DISTANCE DENSITY

densdist = function(di,dens,...)
{
  mu=0;sigma=1;p=0.5
  a1 = -(1-p)/sigma
  b1 = (1-p)*mu/sigma
  a2 = p/sigma
  b2 = -p*mu/sigma
  y = seq(from = 0,to = 6,length.out = 1000)
  return(abs(1/a1)*dens(x = (di-b1)/a1,...) + abs(1/a2)*dens(x = (di-b2)/a2,...))
}

########################################################################
#LOG LIKELIHOOD FUNCTIONS
########################################################################

loglikN = function(x,mu,sigma,p)
{
  return(sum(log(ifelse(x-mu<=0,p*dnorm(x,mu,sqrt((sigma^2)/(4*(1-p)^2))),(1-p)*dnorm(x,mu,sqrt((sigma^2)/(4*p^2)))))))
}

loglikL = function(x,mu,sigma,p)
{
  return(sum(log(densL(x,mu,sigma,p))))
}

loglikT = function(x,mu,sigma,nu,p)
{
  return(sum(log(
    ifelse(test=x<mu,
           yes=(4*p*(1-p)*gamma((nu+1)/2)/(gamma(nu/2)*sigma*sqrt(nu*pi)))*((4*((x-mu)^2)/(nu*sigma^2))*(1-p)^2 +1)^(-(nu+1)/2),
           no=(4*p*(1-p)*gamma((nu+1)/2)/(gamma(nu/2)*sigma*sqrt(nu*pi)))*((4*((x-mu)^2)/(nu*sigma^2))*(p)^2 +1)^(-(nu+1)/2))
  )))
}

loglikSl = function(x,mu,sigma,nu,p)
{
  return(sum(log(
    densSl(x,mu,sigma,nu,p)
  )))
}

loglikNC = function(x,mu,sigma,nu,gama,p)
{
  return(sum(log(
    densNC(x,mu,sigma,nu,gama,p)
  )))
}

AUXloglikNC = function(x,mu,sigma,par,p)
{
  return(sum(log(
    densNC(x,mu,sigma,par[1],par[2],p)
  )))
}

#APPENDIX

rinvgamma <- function(n, shape=1, scale=1)
{return(1 / rgamma(n=n, shape=shape, rate=scale))}

###########################################################################
# Truncated Distribution                                                  #                          #
###########################################################################

###########################################################################
# Truncated Distribution                                                  #
#                                                                         #
# These functions are similar to those from Nadarajah, S. and Kotz, S.    #
# (2006). R Programs for Computing Truncated Distributions. Journal of    #
# Statistical Software, 16, Code Snippet 2, 1-8. These functions have     #
# corrected to work with log-densities.                                   #
###########################################################################

dtrunc <- function(x, spec, a=-Inf, b=Inf, log=FALSE, ...)
{
  if(a >= b) stop("Lower bound a is not less than upper bound b.")
  if(any(x < a) | any(x > b))
    stop("At least one instance of (x < a) or (x > b) found.")
  dens <- rep(0, length(x))
  g <- get(paste("d", spec, sep=""), mode="function")
  G <- get(paste("p", spec, sep=""), mode="function")
  if(log == TRUE) {
    dens <- g(x, log=TRUE, ...) - log(G(b, ...) - G(a, ...))
  }
  else {
    dens <- g(x, ...) / (G(b, ...) - G(a, ...))}
  return(dens)
}
extrunc <- function(spec, a=-Inf, b=Inf, ...)
{
  f <- function(x) x * dtrunc(x, spec, a=a, b=b, log=FALSE, ...)
  return(integrate(f, lower=a, upper=b)$value)
}
ptrunc <- function(x, spec, a=-Inf, b=Inf, ...)
{
  if(a >= b) stop("Lower bound a is not less than upper bound b.")
  if(any(x < a) | any(x > b))
    stop("At least one instance of (x < a) or (x > b) found.")
  p <- x
  aa <- rep(a, length(x))
  bb <- rep(b, length(x))
  G <- get(paste("p", spec, sep=""), mode="function")
  p <- G(apply(cbind(apply(cbind(x, bb), 1, min), aa), 1, max), ...)
  p <- p - G(aa, ...)
  p <- p / {G(bb, ...) - G(aa, ...)}
  return(p)
}
qtrunc <- function(p, spec, a=-Inf, b=Inf, ...)
{
  if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
  if(a >= b) stop("Lower bound a is not less than upper bound b.")
  q <- p
  G <- get(paste("p", spec, sep=""), mode="function")
  Gin <- get(paste("q", spec, sep=""), mode="function")
  q <- Gin(G(a, ...) + p*{G(b, ...) - G(a, ...)}, ...)
  return(q)
}
rtrunc <- function(n, spec, a=-Inf, b=Inf, ...)
{
  if(a >= b) stop("Lower bound a is not less than upper bound b.")
  x <- u <- runif(n)
  x <- qtrunc(u, spec, a=a, b=b,...)
  return(x)
}
vartrunc <- function(spec, a=-Inf, b=Inf, ...)
{
  ex <- extrunc(spec, a=a, b=b, ...)
  f <- function(x) {
    {x - ex}^2 * dtrunc(x, spec, a=a, b=b, log=FALSE, ...)}
  sigma2 <- integrate(f, lower=a, upper=b)$value
  return(sigma2)
}


###########################################################################
###########################################################################
###########################################################################


genSLK = function(n,mu=0,sigma=1,p=0.5,dist = "normal",nu="",gama="")
{
  if(dist == ""){dist = "normal"}
  if(nu=="" && dist == "t"){nu=4}
  if(nu=="" && dist == "slash"){nu=2}
  if(nu=="" && dist == "cont"){nu=0.1}
  if(gama=="" && dist == "cont"){gama=0.1}
  
  rhn = abs(rnorm(n,mean = 0,sd = 1))
  
  if(dist == "normal"){K = 1}
  if(dist == "t"){K = rinvgamma(n,shape = nu/2,scale = nu/2)}
  if(dist == "laplace"){K = rexp(n,rate = 1/2)}
  if(dist == "slash"){K = 1/rbeta(n,nu,1)}
  if(dist == "cont"){K = 1/sample(x = c(gama,1),size = n,replace = T,prob = c(nu,1-nu))}
  
  I = ifelse(runif(n)<p,-(1/(2*(1-p))),1/(2*p))
  y = I*rhn
  return(mu + sigma*sqrt(K)*y)
}

#DENSITIES

densN = function(x,mu=0,sigma=1,p=0.5)
{
  return(2*ifelse(x-mu<=0,p*dnorm(x,mu,sqrt((sigma^2)/(4*(1-p)^2))),(1-p)*dnorm(x,mu,sqrt((sigma^2)/(4*p^2)))))
}

densT = function(x,mu=0,sigma=1,nu=4,p=0.5)
{
  return(
    ifelse(test=x<mu,
           yes=(4*p*(1-p)*gamma((nu+1)/2)/(gamma(nu/2)*sigma*sqrt(nu*pi)))*((4*((x-mu)^2)/(nu*sigma^2))*(1-p)^2 +1)^(-(nu+1)/2),
           no=(4*p*(1-p)*gamma((nu+1)/2)/(gamma(nu/2)*sigma*sqrt(nu*pi)))*((4*((x-mu)^2)/(nu*sigma^2))*(p)^2 +1)^(-(nu+1)/2))
  )
}

densT(0,mu=0,sigma=1,nu=4,p=0.75)



densL = function(x,mu=0,sigma=1,p=0.5)
{
  return(ifelse(test=x<mu,yes=(2*p*(1-p)/sigma)*exp(2*(1-p)*(x-mu)/sigma),no=(2*p*(1-p)/sigma)*exp(-(2*p)*(x-mu)/sigma)))
}

densSl = function(x,mu,sigma,nu,p)
{
  u   = ifelse(x<=mu,yes = rtrunc(n = 1,spec = "gamma", a=0, b=1,shape = nu + (1/2),rate = 2*(((x-mu)/sigma)^2)*(1-p)^2),no = rtrunc(n = 1,spec = "gamma", a=0, b=1,shape = nu + (1/2),rate = 2*(((x-mu)/sigma)^2)*(p)^2))
  gu   = densN(x,mu,sigma*(u^(-1/2)),p)
  fu   = dbeta(x = u,shape1 = nu,shape2 = 1)
  qu   = ifelse(x<=mu,yes = dtrunc(x = u,spec = "gamma", a=0, b=1,shape = nu + (1/2),rate = 2*(((x-mu)/sigma)^2)*(1-p)^2),no = dtrunc(x = u,spec = "gamma", a=0, b=1,shape = nu + (1/2),rate = 2*(((x-mu)/sigma)^2)*(p)^2))
  return(gu*fu/qu)
}

densNC = function(x,mu=0,sigma=1,nu=0.1,gama=0.1,p=0.5)
{
  return(nu*densN(x,mu,sigma/sqrt(gama),p) + (1-nu)*densN(x,mu,sigma,p))
}

inverse <- function (f, lower=-Inf, upper=Inf,...) {
  # DESCRIPTION: Function to calculate the inverse function of a function
  #
  # INPUTS:
  #       f: function for which we want to obtain its inverse
  #       lower: lower limit of f domain (support of the random variable)
  #       upper: upper limit of f domain (support of the random variable)
  #
  # OUTPUT:
  #       A function, the inverse of f
  if (!is.numeric(lower) || !is.numeric(upper) || lower >=upper)
    stop("lower < upper is not fulfilled")
  
  if (!is.finite(f(lower)) & is.finite(lower)) lower <- lower+.Machine$double.eps
  if (!is.finite(f(upper)) & is.finite(upper)) upper <- upper-.Machine$double.eps
  
  if (is.infinite(lower) & is.infinite(upper)) {
    function(y) optim(0,(function (x) abs(f(x) - y)), method="L-BFGS-B")$par
  } else if (is.infinite(lower) & is.finite(upper)) {
    function(y) optim(upper,(function (x) abs(f(x) - y)), lower = lower, upper = upper, method="L-BFGS-B")$par
  } else if (is.finite(lower) & is.infinite(upper)) {
    function(y) optim(lower,(function (x) abs(f(x) - y)), lower = lower, upper = upper, method="L-BFGS-B")$par
  } else {
    if (f(lower)<f(upper)) {
      function(y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)$root
    } else {
      function(y) optim((lower+upper)/2,(function (x) (f(x) - y)^2), lower = lower, upper = upper, method="L-BFGS-B")$par
    }
  }
}