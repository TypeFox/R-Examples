# amended 1/12/2007
#source('/Volumes/Data/Users/stasinom/Documents/gamlss/projects/DISTRIBUTIONS/testingDistributions/testContDist.R')
#testContDist("GP", y.range=c(0,Inf), mu.range = c(0,Inf), sigma.range=c(0, Inf),   mu.val = c(1,5), sigma.val=c(1,2,10))
GP <- function (mu.link="log", sigma.link="log")
{
    mstats <- checklink(   "mu.link", "Generalized Pareto", substitute(mu.link), 
                           c("1/mu^2", "log", "identity", "own"))
    dstats <- checklink(  "sigma.link", "Generalized Pareto", substitute(sigma.link),   
                           c("1/sigma^2", "log", "identity", "own")) 
    structure(
          list(family = c("GP", "Generalized Pareto"),
           parameters = list(mu=TRUE, sigma=TRUE), 
                nopar = 2, 
                 type = "Continuous",
              mu.link = as.character(substitute(mu.link)),  
           sigma.link = as.character(substitute(sigma.link)), 
           mu.linkfun = mstats$linkfun, 
        sigma.linkfun = dstats$linkfun, 
           mu.linkinv = mstats$linkinv, 
        sigma.linkinv = dstats$linkinv,
                mu.dr = mstats$mu.eta, 
             sigma.dr = dstats$mu.eta,  
    dldm = function(y,mu,sigma) 
                      { 
             z <- (y/mu)
          dldm <- -1/mu + (1/mu)*(1+sigma)*z/(1+z)
          dldm
                      },
   d2ldm2 = function(y,mu,sigma)
                      {
        d2ldm2 <- -(sigma)/((2+sigma)*mu^2)
        d2ldm2
                      },     
   dldd = function(y,mu,sigma)  
                      {  
              z <- (y/mu)
           dldd <- -log(1+z)-digamma(sigma)+digamma(1+sigma)
           dldd                     
                      },
   d2ldd2 = function(y,mu,sigma,nu,tau)
                      {
          d2ldd2 <- trigamma(1+sigma) - trigamma(sigma) 
          d2ldd2
                      },   
   d2ldmdd = function(mu,sigma) 
                      {
       d2ldmdd <- 1/(mu*(1+sigma))
       d2ldmdd
                      },
 G.dev.incr  = function(y,mu,sigma,...) 
                      { 
             -2*dGP(y,mu,sigma,log=TRUE)
                      } ,                     
         rqres = expression(   
               rqres(pfun="pGP", type="Continuous", y=y, mu=mu, sigma=sigma)) ,
    mu.initial = expression(mu <- (y+mean(y))/2),    
 sigma.initial = expression(sigma<- rep(1, length(y))),  
      mu.valid = function(mu) all(mu > 0), 
   sigma.valid = function(sigma)  all(sigma > 0),
       y.valid = function(y)  all(y > 0)
          ),
            class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------
dGP <- function(x, mu = 1, sigma = 1, log = FALSE)
 {
          if (any(mu < 0))  stop(paste("mu must be positive", "\n", ""))  
          if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))  
       z <- (x/mu)
       loglik <- log(z)-log(x)-lgamma(1)-lgamma(sigma)+
       lgamma(1+sigma)-(1+sigma)*log(1+z)
       if(log==FALSE) ft  <- exp(loglik) else ft <- loglik 
       ft
  }    
#-----------------------------------------------------------------  
pGP <- function(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
 {  
      if (any(mu < 0))  stop(paste("mu must be positive", "\n", "")) 
      if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))           
      z <- (q/mu)
      w <- (sigma)*z 
      p <- pf(w,2,2*sigma)
    #if (length(sigma)>1)  p <- ifelse(sigma<0,1-p,p)
    #else p <- if (sigma<0) 1-p else p
      if(lower.tail==TRUE) p  <- p else  p <- 1-p 
      if(log.p==FALSE) p  <- p else  p <- log(p) 
      p
 }
#-----------------------------------------------------------------  

qGP <-  function(p, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE)
 {   
    if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))  
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
    if (lower.tail==TRUE) p <- p else p <- 1-p
    #if (length(sigma)>1)  p <- ifelse(sigma<0,1-p,p)
    #else p <- if (sigma<0) 1-p else p
    w <- qf(p,2,2*sigma)   
    q <- mu*(((1/sigma)*w))
    q
 }
#-----------------------------------------------------------------  
rGP <- function(n, mu=1, sigma=1)
  {
    if (any(mu < 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))  
    n <- ceiling(n)
    p <- runif(n)
    r <- qGB2(p,mu=mu,sigma=sigma)
    r
  }
#-----------------------------------------------------------------  
