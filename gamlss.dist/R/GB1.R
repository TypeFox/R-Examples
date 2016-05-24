#   27_11_2007
GB1 <- function (mu.link="logit", sigma.link="logit", nu.link ="log", tau.link="log")
{
    mstats <- checklink("mu.link", "Generalized beta type 1", substitute(mu.link),
                           c("logit", "probit", "cloglog", "log", "own"))
    dstats <- checklink("sigma.link", "Generalized beta type 1", substitute(sigma.link), 
                           c("logit", "probit", "cloglog", "log", "own"))    
    vstats <- checklink(   "nu.link", "Generalized beta type 1", substitute(nu.link),    
                           c("inverse", "log", "identity", "own"))
    tstats <- checklink(  "tau.link", "Generalized beta type 1", substitute(tau.link),   
                           c("inverse", "log", "identity", "own")) 
    structure(
          list(family = c("GB1", "Generalized beta type 1"),
           parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE, tau=TRUE), 
                nopar = 4, 
                 type = "Continuous",
              mu.link = as.character(substitute(mu.link)),  
           sigma.link = as.character(substitute(sigma.link)), 
              nu.link = as.character(substitute(nu.link)), 
             tau.link = as.character(substitute(tau.link)), 
           mu.linkfun = mstats$linkfun, 
        sigma.linkfun = dstats$linkfun, 
           nu.linkfun = vstats$linkfun,
          tau.linkfun = tstats$linkfun,  
           mu.linkinv = mstats$linkinv, 
        sigma.linkinv = dstats$linkinv,
           nu.linkinv = vstats$linkinv,
          tau.linkinv = tstats$linkinv, 
                mu.dr = mstats$mu.eta, 
             sigma.dr = dstats$mu.eta, 
                nu.dr = vstats$mu.eta,
               tau.dr = tstats$mu.eta, 
    dldm = function(y,mu,sigma,nu,tau) 
                   { 
                        a <- mu*((1/sigma^2)-1)
                        b <- a*(1-mu)/mu
                     dldm <- ((1/sigma^2)-1)*( -digamma(a) +
                                  digamma(b) + tau*log(y) - log(1-y^tau) - log(nu) )
                     dldm
                   },
   d2ldm2 = function(y,mu,sigma,nu,tau)
                   {
                        a <- mu*((1/sigma^2)-1)
                        b <- a*(1-mu)/mu
                   d2ldm2 <- -(((1-sigma^2)^2)/(sigma^4))*(trigamma(a) +trigamma(b))
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)  
                   d2ldm2
                    },     
   dldd = function(y,mu,sigma,nu,tau) 
                    {  
                        a <- mu*((1/sigma^2)-1)
                        b <- a*(1-mu)/mu
                    dldd <- -(2/(sigma^3))*( mu*(-digamma(a)+digamma(a+b)+tau*log(y)) +
                              (1-mu)*(log(nu)-digamma(b)+digamma(a+b)+log(1-y^tau)) )
                    dldd <- dldd + (2/(sigma^3))*log(nu+(1-nu)*(y^tau))
                    dldd
                    } ,
   d2ldd2 = function(y,mu,sigma,nu,tau)
                    {
                      a <- mu*((1/sigma^2)-1)
                      b <- a*(1-mu)/mu
                 d2ldd2 <- -(4/(sigma^6))*((mu^2)*trigamma(a) +((1-mu)^2)*trigamma(b) -
                                    trigamma(a+b))
                 d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
                 d2ldd2
                    },   
     dldv = function(y,mu,sigma,nu,tau) 
                    { 
                      a <- mu*((1/sigma^2)-1)
                      b <- a*(1-mu)/mu
                   dldv <- (b/nu) - (a+b)*(1-y^tau)/(nu+(1-nu)*(y^tau))
                   dldv 
                     } ,
    d2ldv2 = function(y,mu,sigma,nu,tau) 
                     { 
                      a <- mu*((1/sigma^2)-1)
                      b <- a*(1-mu)/mu
                 d2ldv2 <- -b + (a+b)*(1-2*mu+(sigma^2)*mu*(1-mu)+(mu^2))
                 d2ldv2 <- d2ldv2/(nu^2)
                 d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)  
                 d2ldv2
                      },
      dldt = function(y,mu,sigma,nu,tau) 
                      {
                      a <- mu*((1/sigma^2)-1)
                      b <- a*(1-mu)/mu
                   dldt <- (1/tau)+a*log(y) - (b-1)*(y^tau)*log(y)/(1-y^tau)
                   dldt <- dldt -(a+b)*(1-nu)*(y^tau)*log(y)/(nu+(1-nu)*(y^tau))
                   dldt
                      } ,
      d2ldt2 = function(y,mu,sigma,nu,tau) 
                      { 
                      a <- mu*((1/sigma^2)-1)
                      b <- a*(1-mu)/mu
                   dldt <- (1/tau)+a*log(y) - (b-1)*(y^tau)*log(y)/(1-y^tau)
                   dldt <- dldt -(a+b)*(1-nu)*(y^tau)*log(y)/(nu+(1-nu)*(y^tau))
                 d2ldt2 <- -dldt*dldt
                 d2ldt2 <- ifelse(d2ldt2 < -1e-15, d2ldt2,-1e-15)                                    
                 d2ldt2
                      } ,
  d2ldmdd = function(y,mu,sigma,nu,tau) 
                      {
                      a <- mu*((1/sigma^2)-1)
                      b <- a*(1-mu)/mu
                d2ldmdd <- (2*(1-sigma^2)/(sigma^5))*(mu*trigamma(a)-(1-mu)*trigamma(b))
                d2ldmdd 
                       },
  d2ldmdv = function(y,mu,sigma,nu,tau) 
                       {
                d2ldmdv <- -(1-(sigma^2))/(nu*(sigma^2))
                d2ldmdv
                       },
  d2ldddv = function(y){
                d2ldddv <- rep(0,length(y))
                d2ldddv
                       },
  d2ldmdt = function(y,mu,sigma,nu,tau) 
                       {
                      a <- mu*((1/sigma^2)-1)
                      b <- a*(1-mu)/mu
                   dldm <- ((1/sigma^2)-1)*( -digamma(a) +
                                  digamma(b) + tau*log(y) - log(1-y^tau) - log(nu) ) 
                   dldt <- (1/tau)+a*log(y) - (b-1)*(y^tau)*log(y)/(1-y^tau)
                   dldt <- dldt -(a+b)*(1-nu)*(y^tau)*log(y)/(nu+(1-nu)*(y^tau))
                d2ldmdt <- -dldm*dldt
                d2ldmdt
                       },
  d2ldddt = function(y,mu,sigma,nu,tau) 
                       {
                      a <- mu*((1/sigma^2)-1)
                      b <- a*(1-mu)/mu
                    dldd <- -(2/(sigma^3))*( mu*(-digamma(a)+digamma(a+b)+tau*log(y)) +
                              (1-mu)*(log(nu)-digamma(b)+digamma(a+b)+log(1-y^tau)) )
                    dldd <- dldd + (2/(sigma^3))*log(nu+(1-nu)*(y^tau))
                    dldt <- (1/tau)+a*log(y) - (b-1)*(y^tau)*log(y)/(1-y^tau)
                    dldt <- dldt -(a+b)*(1-nu)*(y^tau)*log(y)/(nu+(1-nu)*(y^tau))
                 d2ldddt <- -dldd*dldt
                 d2ldddt
                       },
  d2ldvdt = function(y,mu,sigma,nu,tau) 
                       {
                      a <- mu*((1/sigma^2)-1)
                      b <- a*(1-mu)/mu
                   dldv <- (b/nu) - (a+b)*(1-y^tau)/(nu+(1-nu)*(y^tau))
                   dldt <- (1/tau)+a*log(y) - (b-1)*(y^tau)*log(y)/(1-y^tau)
                   dldt <- dldt -(a+b)*(1-nu)*(y^tau)*log(y)/(nu+(1-nu)*(y^tau))
                d2ldvdt <- -dldv*dldt
                d2ldvdt
                       },
 G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
                       { 
                       -2*dGB1(y,mu,sigma,nu,tau,log=TRUE)
                       } ,                     
         rqres = expression(rqres(pfun="pGB1", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)), 
    mu.initial = expression(mu <- (y+mean(y))/2),  
 sigma.initial = expression(sigma<- rep(0.4, length(y))), 
    nu.initial = expression(nu <- rep(1, length(y))), 
   tau.initial = expression(tau <-rep(1, length(y))),
      mu.valid = function(mu) all(mu > 0 & mu < 1) ,  
   sigma.valid = function(sigma)   all(sigma > 0 & sigma < 1), 
      nu.valid = function(nu) all(nu > 0), 
     tau.valid = function(tau) all(tau > 0), 
       y.valid = function(y)  all(y > 0 & y < 1)
          ),
            class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------
dGB1 <- function(x,  mu = 0.5, sigma = 0.4, nu = 1, tau = 1, log = FALSE)
 {
          if (any(mu <= 0) | any(mu >= 1) )  stop(paste("mu must be between 0 and 1", "\n", "")) 
          if (any(sigma <= 0) | any(sigma >= 1))  stop(paste("sigma must be between 0 and 1", "\n", "")) 
          if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))  
          if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
          if (any(x <= 0) | any(x >= 1))  stop(paste("x must be between 0 and 1", "\n", ""))  
                      a <- mu*((1/sigma^2)-1)
                      b <- a*(1-mu)/mu
       loglik <- log(tau) + b*log(nu) + (tau*a-1)*log(x) + (b-1)*log(1-x^tau)
       loglik <- loglik -lgamma(a)-lgamma(b)+lgamma(a+b) - (a+b)*log(nu+(1-nu)*(x^tau))
       if(log==FALSE) ft  <- exp(loglik) else ft <- loglik 
       ft
  }    
#-----------------------------------------------------------------  
pGB1 <- function(q, mu = 0.5, sigma = 0.4, nu = 1, tau = 1, lower.tail = TRUE, log.p = FALSE)
 {  
          if (any(mu <= 0) | any(mu >= 1) )  stop(paste("mu must be between 0 and 1", "\n", "")) 
          if (any(sigma <= 0) | any(sigma >= 1))  stop(paste("sigma must be between 0 and 1", "\n", "")) 
          if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))  
          if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
          if (any(q <= 0) | any(q >= 1))  stop(paste("y must be between 0 and 1", "\n", ""))  
                      a <- mu*((1/sigma^2)-1)
                      b <- a*(1-mu)/mu
                      z <- (q^tau)/(nu+(1-nu)*(q^tau))
                      p <- pbeta(z,a,b)
      if(lower.tail==TRUE) p  <- p else  p <- 1-p 
      if(log.p==FALSE) p  <- p else  p <- log(p) 
      p
 }
#-----------------------------------------------------------------  

qGB1 <-  function(p,mu = 0.5, sigma = 0.4, nu = 1, tau = 1, lower.tail = TRUE, log.p = FALSE)
 {   
          if (any(mu <= 0) | any(mu >= 1) )  stop(paste("mu must be between 0 and 1", "\n", "")) 
          if (any(sigma <= 0) | any(sigma >= 1))  stop(paste("sigma must be between 0 and 1", "\n", "")) 
          if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))  
          if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
    if (lower.tail==TRUE) p <- p else p <- 1-p
                      a <- mu*((1/sigma^2)-1)
                      b <- a*(1-mu)/mu
                      z <- qbeta(p,a,b)
                      q <- (nu/((1/z)-(1-nu)))^(1/tau)
                      q
 }
#-----------------------------------------------------------------  
rGB1 <- function(n, mu = 0.5, sigma = 0.4, nu = 1, tau = 1)
  {
          if (any(mu <= 0) | any(mu >= 1) )  stop(paste("mu must be between 0 and 1", "\n", "")) 
          if (any(sigma <= 0) | any(sigma >= 1))  stop(paste("sigma must be between 0 and 1", "\n", "")) 
          if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))  
          if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
    n <- ceiling(n)
    p <- runif(n)
    r <- qGB1(p,mu=mu,sigma=sigma,nu=nu,tau=tau)
    r
  }
#-----------------------------------------------------------------  
