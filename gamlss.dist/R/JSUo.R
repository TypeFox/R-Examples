# amended 27_11_2007
JSUo <- function (mu.link="identity", sigma.link="log", nu.link ="identity", tau.link="log")
{
    mstats <- checklink(   "mu.link", "Johnson SU", substitute(mu.link), 
                           c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Johnson SU", substitute(sigma.link), 
                           c("inverse", "log", "identity", "own"))
    vstats <- checklink(   "nu.link", "Johnson SU", substitute(nu.link),    
                           c("inverse", "log", "identity","own"))
    tstats <- checklink(  "tau.link", "Johnson SU", substitute(tau.link),   
                           c("inverse", "log", "identity", "own")) 
    structure(
          list(family = c("JSUo", "Johnson SU original"),
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
    dldm = function(y,mu,sigma,nu,tau) { 
      z <- (y-mu)/sigma
      r <- nu + tau*asinh(z)
   dldm <- (z/(sigma*(z*z+1)))+((r*tau)/(sigma*(z*z+1)^(0.5)))
   dldm
                      },
   d2ldm2 = function(y,mu,sigma,nu,tau){
      z <- (y-mu)/sigma
      r <- nu + tau*asinh(z)
   dldm <- (z/(sigma*(z*z+1)))+((r*tau)/(sigma*(z*z+1)^(0.5)))
   d2ldm2 <- -dldm*dldm
    d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)                                    
    d2ldm2
                      },
                       
   dldd = function(y,mu,sigma,nu,tau) {  
      z <- (y-mu)/sigma
      r <- nu + tau*asinh(z)
      dldd <- (-1/(sigma*(z*z+1)))+((r*tau*z)/(sigma*(z*z+1)^(0.5)))
      dldd 
                      } ,
   d2ldd2 = function(y,mu,sigma,nu,tau){
      z <- (y-mu)/sigma
      r <- nu + tau*asinh(z)
      dldd <- (-1/(sigma*(z*z+1)))+((r*tau*z)/(sigma*(z*z+1)^(0.5))) 
      d2ldd2 <- -dldd*dldd
   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
   d2ldd2
                      },   
     dldv = function(y,mu,sigma,nu,tau) { 
      z <- (y-mu)/sigma
      r <- nu + tau*asinh(z)
      dldv <- -r
      dldv
                        } ,
    d2ldv2 = function(y,mu,sigma,nu,tau) { 
      z <- (y-mu)/sigma
      r <- nu + tau*asinh(z)
      d2ldv2 <-  -r^2                
   d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)  
   d2ldv2
                        },
      dldt = function(y,mu,sigma,nu,tau) {
      z <- (y-mu)/sigma
      r <- nu + tau*asinh(z)                   
      dldt <- (1+r*nu-r*r)/tau
      dldt
                        } ,
      d2ldt2 = function(y,mu,sigma,nu,tau) { 
       z <- (y-mu)/sigma
      r <- nu + tau*asinh(z)                   
      dldt <- (1+r*nu-r*r)/tau
      d2ldt2 <-   -dldt*dldt  
    d2ldt2 <- ifelse(d2ldt2 < -1e-15, d2ldt2,-1e-15)                                    
    d2ldt2
                            } ,
  d2ldmdd = function(y,mu,sigma,nu,tau) {
        z <- (y-mu)/sigma
        r <- nu + tau*asinh(z)
     dldm <- (z/(sigma*(z*z+1)))+((r*tau)/(sigma*(z*z+1)^(0.5)))
     dldd <- (-1/(sigma*(z*z+1)))+((r*tau*z)/(sigma*(z*z+1)^(0.5))) 
  d2ldmdd <- -(dldm*dldd)
  d2ldmdd
                       },
  d2ldmdv = function(y,mu,sigma,nu,tau) {
        z <- (y-mu)/sigma
        r <- nu + tau*asinh(z)
     dldm <- (z/(sigma*(z*z+1)))+((r*tau)/(sigma*(z*z+1)^(0.5)))
     dldv <- -r 
  d2ldmdv <- -(dldm*dldv)
  d2ldmdv
                       },
  d2ldmdt = function(y,mu,sigma,nu,tau) {
        z <- (y-mu)/sigma
        r <- nu + tau*asinh(z)
     dldm <- (z/(sigma*(z*z+1)))+((r*tau)/(sigma*(z*z+1)^(0.5)))
     dldt <- (1+r*nu-r*r)/tau
  d2ldmdt <- -(dldm*dldt)
  d2ldmdt
                       },
  d2ldddv = function(y,mu,sigma,nu,tau) {
        z <- (y-mu)/sigma
        r <- nu + tau*asinh(z)
     dldd <- (-1/(sigma*(z*z+1)))+((r*tau*z)/(sigma*(z*z+1)^(0.5))) 
     dldv <- -r 
  d2ldddv <- -(dldd*dldv)
  d2ldddv
                       },
  d2ldddt = function(y,mu,sigma,nu,tau) {
        z <- (y-mu)/sigma
        r <- nu + tau*asinh(z)
     dldd <- (-1/(sigma*(z*z+1)))+((r*tau*z)/(sigma*(z*z+1)^(0.5))) 
     dldt <- (1+r*nu-r*r)/tau
  d2ldddt <- -(dldd*dldt)
  d2ldddt  
                       },
  d2ldvdt = function(y,mu,sigma,nu,tau) {
        z <- (y-mu)/sigma
        r <- nu + tau*asinh(z)
     dldv <- -r 
     dldt <- (1+r*nu-r*r)/tau
  d2ldvdt <- -(dldv*dldt)
  d2ldvdt
                       },
 G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
                       { 
       -2*dJSUo(y,mu,sigma,nu,tau,log=TRUE)
                        } ,                     
         rqres = expression(   
                rqres(pfun="pJSUo", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)),
    mu.initial = expression(mu <- (y+mean(y))/2), 
 sigma.initial = expression(sigma<- rep(.1, length(y))),
    nu.initial = expression(nu <- rep(0, length(y))), 
   tau.initial = expression(tau <-rep(0.5, length(y))), 
      mu.valid = function(mu) TRUE, 
   sigma.valid = function(sigma)  all(sigma > 0),
      nu.valid = function(nu) TRUE , 
     tau.valid = function(tau) all(tau > 0), 
       y.valid = function(y)  TRUE
          ),
            class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------
dJSUo <- function(x, mu=0, sigma=1, nu=0, tau=1, log=FALSE)
 {
          if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
       z <- (x-mu)/sigma
       r <- nu + tau*asinh(z)
       loglik <- -log(sigma)+ log(tau)- .5*log(z*z+1) -.5*log(2*pi)-.5*r*r
       if(log==FALSE) ft  <- exp(loglik) else ft <- loglik 
       ft
  }    
#-----------------------------------------------------------------  
pJSUo <- function(q, mu=0, sigma=1, nu=0, tau=1, lower.tail = TRUE, log.p = FALSE)
 {  
         if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
         if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))           
      z <- (q-mu)/sigma
      r <- nu + tau * asinh(z)    
      p <- pNO(r,0,1)
      if(lower.tail==TRUE) p  <- p else  p <- 1-p 
      if(log.p==FALSE) p  <- p else  p <- log(p) 
      p
 }
#-----------------------------------------------------------------  

qJSUo <-  function(p, mu=0, sigma=1, nu=0, tau=1, lower.tail = TRUE, log.p = FALSE)
 {   
    if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
    if (lower.tail==TRUE) p <- p else p <- 1-p
    r <- qNO(p,0,1)
    z <- sinh((r-nu)/tau)     
    q <- mu+sigma*z   
    q
 }
#-----------------------------------------------------------------  
rJSUo <- function(n, mu=0, sigma=1, nu=0, tau=1)
  {
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    n <- ceiling(n)
    p <- runif(n)
    r <- qJSUo(p,mu=mu,sigma=sigma,nu=nu,tau=tau)
    r
  }
#-----------------------------------------------------------------  
