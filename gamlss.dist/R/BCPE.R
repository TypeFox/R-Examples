#----------------------------------------------------------------------------------------
# last change MS Wednesday, September 17, 2003 at 08:43 
BCPE <- function (mu.link="identity", sigma.link="log", nu.link ="identity", tau.link="log")
{
    mstats <- checklink(   "mu.link", "Box Cox Power Exponential", substitute(mu.link), 
                           c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Box Cox Power Exponential", substitute(sigma.link), 
                           c("inverse", "log", "identity", "own"))
    vstats <- checklink(   "nu.link", "Box Cox Power Exponential", substitute(nu.link),    
                           c("inverse", "log", "identity", "own"))
    tstats <- checklink(  "tau.link", "Box Cox Power Exponential", substitute(tau.link),   
                           c("logshiftto1", "log", "identity", "own")) 
    structure(
          list(family = c("BCPE", "Box-Cox Power Exponential"),
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
      z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
  log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
      c <- exp(log.c)
   dldm <- (tau/(2*mu*sigma*c^2))*(z+sigma*nu*z^2)*
                        ((abs(z/c))^(tau-2))-(nu/mu)
   dldm
                      },
   d2ldm2 = function(y,mu,sigma,nu,tau){
       z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
   log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
       c <- exp(log.c)
    dldm <- (tau/(2*mu*sigma*c^2))*(z+sigma*nu*z^2)*
                        ((abs(z/c))^(tau-2))-(nu/mu)
   d2ldm2 <- -((tau*tau)*gamma(2-(1/tau))*gamma(3/tau))/(mu^2*sigma^2*(gamma(1/tau))^2)
   d2ldm2 <- d2ldm2-(tau*nu^2)/mu^2 # BR MS Wednesday, February 4, 2004 at 15:25
   d2ldm2 <- if (any(tau<1.05)) -dldm*dldm else d2ldm2
                      },
                       
   dldd = function(y,mu,sigma,nu,tau) {  
   f.T <- function(t,log=FALSE)
        {
      log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
          c <- exp(log.c)
    loglik <- log(tau)-log.c-(0.5*(abs(t/c)^tau))-(1+(1/tau))*log(2)-lgamma(1/tau)
      if(log==FALSE) fT  <- exp(loglik) else fT <- loglik
        fT     
        }
    F.T <- function(t)
        {
     log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
         c <- exp(log.c) 
         s <- 0.5*((abs(t/c))^tau)
       F.s <- pgamma(s,shape = 1/tau, scale = 1)
       cdf <- 0.5*(1+F.s*sign(t))        
       cdf   
        } 
       z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
   log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
       c <- exp(log.c) 
       h <- f.T(1/(sigma*abs(nu)))/F.T(1/(sigma*abs(nu)))
    dldd <- (1/sigma)*((tau/2)*(abs(z/c))^tau-1)+h/(sigma^2*abs(nu))
    dldd
                      } ,
   d2ldd2 = function(sigma,tau) -tau/sigma^2,
     dldv = function(y,mu,sigma,nu,tau) { 
    f.T <- function(t,log=FALSE)
       {
     log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
         c <- exp(log.c)
   loglik <- log(tau)-log.c-(0.5*(abs(t/c)^tau))-(1+(1/tau))*log(2)-lgamma(1/tau)
    if(log==FALSE) fT  <- exp(loglik) else fT <- loglik
        fT     
        }   
    F.T <- function(t)
        {
     log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
         c <- exp(log.c) 
         s <- 0.5*((abs(t/c))^tau)
       F.s <- pgamma(s,shape = 1/tau, scale = 1)
       cdf <- 0.5*(1+F.s*sign(t))        
       cdf   
       } 
      z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
  log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
      c <- exp(log.c)
      h <- f.T(1/(sigma*abs(nu)))/F.T(1/(sigma*abs(nu)))                                      
   dldv <-  -(tau/(2*nu*c^2))*((abs(z/c))^(tau-2))*z*((nu*z+1/sigma)
                       *log(y/mu)-z)
   dldv <- dldv+log(y/mu)+sign(nu)*h/(sigma*nu^2)
                        } ,
    d2ldv2 = function(sigma,tau) { 
                      -sigma^2*(3*tau+1)/4  
                        },
      dldt = function(y,mu,sigma,nu,tau) {
      F.T <- function(t,tau)
          {
    log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
        c <- exp(log.c) 
        s <- 0.5*((abs(t/c))^tau)
      F.s <- pgamma(s,shape = 1/tau, scale = 1)
      cdf <- 0.5*(1+F.s*sign(t))        
      cdf   
          } 
        z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
    log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
        c <- exp(log.c)
 dlogc.dt <- (1/(2*tau^2))*(2*log(2)-digamma(1/tau)+3*digamma(3/tau))   
        j <- (log(F.T(1/(sigma*abs(nu)),tau+0.001))-log(F.T(1/(sigma*abs(nu)),tau)))/0.001
     dldt <- (1/tau)-0.5*(log(abs(z/c)))*(abs(z/c))^tau+
                      (1/tau^2)*(log(2)+digamma(1/tau))+
                      ((tau/2)*((abs(z/c))^tau)-1)*dlogc.dt-j               
     dldt                                       
                        } ,
      d2ldt2 = function(y,mu,sigma,nu,tau) {
        p <- (tau+1)/tau
 dlogc.dt <- (1/(2*tau^2))*(2*log(2)-digamma(1/tau)+3*digamma(3/tau))         
    part1 <- p*trigamma(p)+2*(digamma(p))^2
    part2 <- digamma(p)*(log(2)+3-3*digamma(3/tau)-tau)
    part3 <- -3*(digamma(3/tau))*(1+log(2))    
    part4 <- -(tau+log(2))*log(2)
    part5 <- -tau+(tau^4)*(dlogc.dt)^2
   d2ldt2 <- part1+part2+part3+part4+part5
   d2ldt2 <- -d2ldt2/tau^3    
   d2ldt2 <- ifelse(d2ldt2 < -1e-15, d2ldt2,-1e-15)                                    
   d2ldt2
                            } ,
  d2ldmdd = function(mu,sigma,nu,tau) -(nu*tau)/(mu*sigma),
  d2ldmdv = function(mu,sigma,nu,tau) (2*(tau-1)-(tau+1)*(sigma^2)*(nu^2))/(4*mu),
  d2ldmdt = function(mu,sigma,nu,tau) (nu/(mu*tau))*(1+tau+(3/2)*(digamma(1/tau)-digamma(3/tau))),
  d2ldddv = function(sigma,nu,tau) -(sigma*nu*tau)/2,
  d2ldddt = function(sigma,tau) (1/(sigma*tau))*(1+tau+(3/2)*(digamma(1/tau)-digamma(3/tau))),
  d2ldvdt = function(mu,sigma,nu,tau) {
    d2ldvdt  <- (((sigma^2)*nu)/(2*tau))*(1+(tau/3)+0.5*(digamma(1/tau)-digamma(3/tau)))
    d2ldvdt
                        },
 G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
                       -2*dBCPE(y,mu,sigma,nu,tau,log=TRUE),                     
         rqres = expression(rqres(pfun="pBCPE", type="Continuous", y=y, mu=mu, 
                                              sigma=sigma, nu=nu, tau=tau) ),
    mu.initial = expression(mu <- (y+mean(y))/2), #
 sigma.initial = expression(sigma<- rep(0.1, length(y))),
    nu.initial = expression(nu <- rep(1, length(y))), 
   tau.initial = expression(tau <-rep(2, length(y))), 
      mu.valid = function(mu) all(mu > 0), 
   sigma.valid = function(sigma)  all(sigma > 0),
      nu.valid = function(nu) TRUE , 
     tau.valid = function(tau) all(tau > 0), 
       y.valid = function(y)  all(y > 0)
          ),
            class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------

#-----------------------------------------------------------------
dBCPE <- dBCPEo <- function(x, mu=5, sigma=0.1, nu=1, tau=2, log=FALSE)
 {
     f.T <- function(t,log=FALSE){
         log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
             c <- exp(log.c)
       log.lik <- log(tau)-log.c-(0.5*(abs(t/c)^tau))-(1+(1/tau))*log(2)-lgamma(1/tau)
       if(log==FALSE) fT  <- exp(log.lik) else fT <- log.lik
           fT                    }
     F.T <- function(t){
        log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
            c <- exp(log.c) 
            s <- 0.5*((abs(t/c))^tau)
          F.s <- pgamma(s,shape = 1/tau, scale = 1)
          cdf <- 0.5*(1+F.s*sign(t))        
                cdf   }      
          if (any(mu < 0))  stop(paste("mu must be positive", "\n", ""))  
          if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
          if (any(x < 0))  stop(paste("x must be positive", "\n", ""))  
          if(length(nu)>1)  z <- ifelse(nu != 0,(((x/mu)^nu-1)/(nu*sigma)),log(x/mu)/sigma)
          else   if (nu != 0) z <- (((x/mu)^nu-1)/(nu*sigma)) else z <- log(x/mu)/sigma
        logfZ <- f.T(z,log=TRUE)-log(F.T(1/(sigma*abs(nu))))
       logder <- (nu-1)*log(x)-nu*log(mu)-log(sigma)
       loglik <- logder+logfZ
       if(log==FALSE) ft  <- exp(loglik) else ft <- loglik 
       ft
  }    
#-----------------------------------------------------------------  
pBCPE <- pBCPEo <- function(q, mu=5, sigma=0.1, nu=1, tau=2, lower.tail = TRUE, log.p = FALSE)
 {  
        F.T <- function(t,tau){
             log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
                 c <- exp(log.c) 
                 s <- 0.5*((abs(t/c))^tau)
               F.s <- pgamma(s,shape = 1/tau, scale = 1)
               cdf <- 0.5*(1+F.s*sign(t))        
               cdf       }  
         if (any(mu < 0))  stop(paste("mu must be positive", "\n", "")) 
         if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
         if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
         if (any(q < 0))  stop(paste("y must be positive", "\n", ""))  
         
         if(length(nu)>1)  z <- ifelse(nu != 0,(((q/mu)^nu-1)/(nu*sigma)),log(q/mu)/sigma)
         else   if (nu != 0) z <- (((q/mu)^nu-1)/(nu*sigma)) else z <- log(q/mu)/sigma
          FYy1 <- F.T(z,tau)
         if(length(nu)>1)  FYy2 <- ifelse(nu > 0, F.T( -1/(sigma*abs(nu)),tau),0)
         else   if (nu>0)  FYy2 <-  F.T(-1/(sigma*abs(nu)),tau) else FYy2 <- 0
         FYy3 <- F.T(1/(sigma*abs(nu)),tau)
         FYy  <- (FYy1-FYy2)/FYy3
         if(lower.tail==TRUE) FYy  <- FYy else  FYy <- 1-FYy 
         if(log.p==FALSE) FYy  <- FYy else  FYy<- log(FYy) 
         FYy     
 }
#-----------------------------------------------------------------  

qBCPE <- qBCPEo <-  function(p, mu=5, sigma=0.1, nu=1, tau=2, lower.tail = TRUE, log.p = FALSE)
 {   F.T <- function(t,tau) # cdf of PE(0,1,tau)
        {
     log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
         c <- exp(log.c) 
         s <- 0.5*((abs(t/c))^tau)
       F.s <- pgamma(s,shape = 1/tau, scale = 1)
       cdf <- 0.5*(1+F.s*sign(t))        
      cdf   
       } 
   q.T <- function(p, tau, lower.tail = TRUE, log.p = FALSE)# inverse of PE(0,1,tau)
       {  
   log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
       c <- exp(log.c)
       s <- qgamma((2*p-1)*sign(p-0.5),shape=(1/tau),scale=1)
       z <- sign(p-0.5)*((2*s)^(1/tau))*c
       z  
      }   
    if (any(mu < 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
    if (lower.tail==TRUE) p <- p else p <- 1-p
    if(length(nu)>1){ 
                za <- ifelse(nu<0, 
                             q.T(p*F.T(1/(sigma*abs(nu)),tau),tau),
                             q.T((1-(1-p)*F.T(1/(sigma*abs(nu)),tau)),tau))
                za <- ifelse(nu==0, q.T(p,tau), za)                  
                    } 
               else { if (nu<0)  {za <- q.T(p*F.T(1/(sigma*abs(nu)),tau),tau)}
                      if (nu==0) {za <- q.T(p,tau)} 
                      if (nu>0)  {za <- q.T((1-(1-p)*F.T(1/(sigma*abs(nu)),tau)),tau)}
                    }    
     if(length(nu)>1)  ya <- ifelse(nu != 0,mu*(nu*sigma*za+1)^(1/nu),mu*exp(sigma*za))
       else   if (nu != 0) ya <- mu*(nu*sigma*za+1)^(1/nu) else ya <- mu*exp(sigma*za)
     ya   
 }
#-----------------------------------------------------------------  
rBCPE <- rBCPEo <- function(n, mu=5, sigma=0.1, nu=1, tau=2)
  {
    if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))  
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))    
    n <- ceiling(n)
    p <- runif(n)
    r <- qBCPE(p,mu=mu,sigma=sigma,nu=nu,tau=tau)
    r
  }

#-----------------------------------------------------------------  
BCPEuntr <- function (mu.link="identity", sigma.link="log", nu.link ="identity", tau.link="log")
{
    mstats <- checklink(   "mu.link", "Box.Cox.Power.Exponential", substitute(mu.link),    
                            c("1/mu^2", "log", "identity"))
    dstats <- checklink("sigma.link", "Box.Cox.Power.Exponential", substitute(sigma.link), 
                            c("inverse", "log", "identity"))
    vstats <- checklink(   "nu.link", "Box.Cox.Power.Exponential", substitute(nu.link),    
                            c("1/nu^2", "log", "identity"))
    tstats <- checklink(  "tau.link", "Box.Cox.Power.Exponential", substitute(tau.link),   
                            c("1/tau^2", "log", "identity"))   
    structure(
          list(family = c("BCPEuntr", "Box-Cox Power Exponential untrucated"),
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
                       z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
                       log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
                           c <- exp(log.c)
                        dldm <- (tau/(2*mu*sigma*c^2))*(z+sigma*nu*z^2)*
                                    ((abs(z/c))^(tau-2))-(nu/mu)
                        dldm
                        },
   d2ldm2 = function(mu,sigma,nu,tau){
             d2ldm2 <- -((tau*tau)*gamma(2-(1/tau))*gamma(3/tau))/(mu^2*sigma^2*(gamma(1/tau))^2)
             d2ldm2 <- d2ldm2-(tau*nu^2)/mu^2
                       },
   dldd = function(y,mu,sigma,nu,tau) {    z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
                      log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
                          c <- exp(log.c) 
                       dldd <- (1/sigma)*((tau/2)*(abs(z/c))^tau-1)
                       dldd
                        } ,
   d2ldd2 = function(sigma,tau) -tau/sigma^2,
     dldv = function(y,mu,sigma,nu,tau) {    
                         z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
                     log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
                         c <- exp(log.c)
                 # dlogc.dt <- (1/(2*tau^2))*(2*log(2)-digamma(1/tau)+3*digamma(3/tau))                                          
                      dldv <- -(tau/(2*nu*c^2))*((abs(z/c))^(tau-2))*z*((nu*z+1/sigma)
                               *log(y/mu)-z)
                      dldv <- dldv+log(y/mu)
                        } ,
    d2ldv2 = function(sigma,tau) { 
                         -sigma^2*(3*tau+1)/4  
                        },
      dldt = function(y,mu,sigma,nu,tau) {
                      z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
                  log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
                      c <- exp(log.c)
               dlogc.dt <- (1/(2*tau^2))*(2*log(2)-digamma(1/tau)+3*digamma(3/tau))   
                   dldt <- 1/tau-0.5*(log(abs(z/c)))*(abs(z/c))^tau+(1/tau^2)*log(2)
                   dldt <- dldt+(1/tau^2)*digamma(1/tau)
                   dldt <- dldt+((tau/2)*((abs(z/c))^tau)-1)*dlogc.dt                                
                        } ,
      d2ldt2 = function(y,mu,sigma,nu,tau) { 
        p <- (tau+1)/tau
 dlogc.dt <- (1/(2*tau^2))*(2*log(2)-digamma(1/tau)+3*digamma(3/tau))         
    part1 <- p*trigamma(p)+2*(digamma(p))^2
    part2 <- digamma(p)*(log(2)+3-3*digamma(3/tau)-tau)
    part3 <- -3*(digamma(3/tau))*(1+log(2))    
    part4 <- -(tau+log(2))*log(2)
    part5 <- -tau+(tau^4)*(dlogc.dt)^2
   d2ldt2 <- part1+part2+part3+part4+part5
   d2ldt2 <- -d2ldt2/tau^3    
   d2ldt2 <- ifelse(d2ldt2 < -1e-15, d2ldt2,-1e-15)                                    
   d2ldt2
                            } ,
 d2ldmdd = function(mu,sigma,nu,tau) -(nu*tau)/(mu*sigma),
  d2ldmdv = function(mu,sigma,nu,tau) (2*(tau-1)-(tau+1)*(sigma^2)*(nu^2))/(4*mu),
  d2ldmdt = function(mu,sigma,nu,tau) (nu/(mu*tau))*(1+tau+(3/2)*(digamma(1/tau)-digamma(3/tau))),
  d2ldddv = function(sigma,nu,tau) -(sigma*nu*tau)/2,
  d2ldddt = function(sigma,tau) (1/(sigma*tau))*(1+tau+(3/2)*(digamma(1/tau)-digamma(3/tau))),
  d2ldvdt = function(y,mu,sigma,nu,tau) {
    d2ldvdt  <- (((sigma^2)*nu)/(2*tau))*(1+(tau/3)+0.5*(digamma(1/tau)-digamma(3/tau)))
    d2ldvdt
                        },
          G.dev.incr  = function(y,mu,sigma,nu,tau,...) {
                 z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
             log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
                 c <- exp(log.c)
                 t <- 1/(sigma*abs(nu))
                 s <- 0.5*((abs(t/c))^tau)
                tf <- pgamma(s,shape = 1/tau, scale = 1)  
                 if(any(tf<.99))
                     { 
                      warning(paste("The truncation factor F(1/(sigma*|nu|)) \n", 
                             "of BCPE is less than 0.99 for some observations",
                                "\n"))       
                     }
               lik <- (nu-1)*log(y)-nu*log(mu)-log(sigma)+log(tau)-log.c
               lik <-  lik-.5*((abs(z/c))^tau)-(1+(1/tau))*log(2)-lgamma(1/tau)
        G.dev.incr <- -2*lik
        G.dev.incr
                                                        } , 
                   cdf = function(y,mu,sigma,nu,tau,...){
               z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
             log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
                 c <- exp(log.c) 
                 s <- 0.5*((abs(z/c))^tau)
               F.s <- pgamma(s,shape = 1/tau, scale = 1)
               cdf <- 0.5*(1+F.s*sign(z)) 
               cdf         
                                                    },                                            
                rqres =expression(   {                                                               
                 z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
             log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
                 c <- exp(log.c) 
                 s <- 0.5*((abs(z/c))^tau)
               F.s <- pgamma(s,shape = 1/tau, scale = 1)
               cdf <- 0.5*(1+F.s*sign(z))      
             rqres <- qnorm(cdf)
                                      }) ,
            mu.initial = expression(mu <- y+0.0001), 
         sigma.initial = expression(sigma<- rep(0.1, length(y))),#rep(sd(y)/mean(y),length(y))),
            nu.initial = expression(nu <- rep(1, length(y))), 
           tau.initial = expression(tau <-rep(2, length(y))), 
              mu.valid = function(mu) all(mu > 0), # MS Friday, September 5, 2003 at 19:42
           sigma.valid = function(sigma)  all(sigma > 0),
              nu.valid = function(nu) TRUE , 
             tau.valid = function(tau) all(tau > 0), 
               y.valid = function(y)  all(y > 0)
          ),
            class = c("gamlss.family","family"))
}
#---------------------------------------------------------------------
#------------------------------------------------------------------------
BCPEo <- function (mu.link="log", sigma.link="log", nu.link ="identity", tau.link="log")
{
    mstats <- checklink(   "mu.link", "Box Cox Power Exponential", substitute(mu.link), 
                           c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Box Cox Power Exponential", substitute(sigma.link), 
                           c("inverse", "log", "identity", "own"))
    vstats <- checklink(   "nu.link", "Box Cox Power Exponential", substitute(nu.link),    
                           c("inverse", "log", "identity", "own"))
    tstats <- checklink(  "tau.link", "Box Cox Power Exponential", substitute(tau.link),   
                           c("logshiftto1", "log", "identity", "own")) 
    structure(
          list(family = c("BCPEo", "Box-Cox Power Exponential-orig."),
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
      z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
  log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
      c <- exp(log.c)
   dldm <- (tau/(2*mu*sigma*c^2))*(z+sigma*nu*z^2)*
                        ((abs(z/c))^(tau-2))-(nu/mu)
   dldm
                      },
   d2ldm2 = function(y,mu,sigma,nu,tau){
       z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
   log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
       c <- exp(log.c)
    dldm <- (tau/(2*mu*sigma*c^2))*(z+sigma*nu*z^2)*
                        ((abs(z/c))^(tau-2))-(nu/mu)
   d2ldm2 <- -((tau*tau)*gamma(2-(1/tau))*gamma(3/tau))/(mu^2*sigma^2*(gamma(1/tau))^2)
   d2ldm2 <- d2ldm2-(tau*nu^2)/mu^2 # BR MS Wednesday, February 4, 2004 at 15:25
   d2ldm2 <- if (any(tau<1.05)) -dldm*dldm else d2ldm2
                      },
                       
   dldd = function(y,mu,sigma,nu,tau) {  
   f.T <- function(t,log=FALSE)
        {
      log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
          c <- exp(log.c)
    loglik <- log(tau)-log.c-(0.5*(abs(t/c)^tau))-(1+(1/tau))*log(2)-lgamma(1/tau)
      if(log==FALSE) fT  <- exp(loglik) else fT <- loglik
        fT     
        }
    F.T <- function(t)
        {
     log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
         c <- exp(log.c) 
         s <- 0.5*((abs(t/c))^tau)
       F.s <- pgamma(s,shape = 1/tau, scale = 1)
       cdf <- 0.5*(1+F.s*sign(t))        
       cdf   
        } 
       z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
   log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
       c <- exp(log.c) 
       h <- f.T(1/(sigma*abs(nu)))/F.T(1/(sigma*abs(nu)))
    dldd <- (1/sigma)*((tau/2)*(abs(z/c))^tau-1)+h/(sigma^2*abs(nu))
    dldd
                      } ,
   d2ldd2 = function(sigma,tau) -tau/sigma^2,
     dldv = function(y,mu,sigma,nu,tau) { 
    f.T <- function(t,log=FALSE)
       {
     log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
         c <- exp(log.c)
   loglik <- log(tau)-log.c-(0.5*(abs(t/c)^tau))-(1+(1/tau))*log(2)-lgamma(1/tau)
    if(log==FALSE) fT  <- exp(loglik) else fT <- loglik
        fT     
        }   
    F.T <- function(t)
        {
     log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
         c <- exp(log.c) 
         s <- 0.5*((abs(t/c))^tau)
       F.s <- pgamma(s,shape = 1/tau, scale = 1)
       cdf <- 0.5*(1+F.s*sign(t))        
       cdf   
       } 
      z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
  log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
      c <- exp(log.c)
      h <- f.T(1/(sigma*abs(nu)))/F.T(1/(sigma*abs(nu)))                                      
   dldv <-  -(tau/(2*nu*c^2))*((abs(z/c))^(tau-2))*z*((nu*z+1/sigma)
                       *log(y/mu)-z)
   dldv <- dldv+log(y/mu)+sign(nu)*h/(sigma*nu^2)
                        } ,
    d2ldv2 = function(sigma,tau) { 
                      -sigma^2*(3*tau+1)/4  
                        },
      dldt = function(y,mu,sigma,nu,tau) {
      F.T <- function(t,tau)
          {
    log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
        c <- exp(log.c) 
        s <- 0.5*((abs(t/c))^tau)
      F.s <- pgamma(s,shape = 1/tau, scale = 1)
      cdf <- 0.5*(1+F.s*sign(t))        
      cdf   
          } 
        z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
    log.c <- 0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
        c <- exp(log.c)
 dlogc.dt <- (1/(2*tau^2))*(2*log(2)-digamma(1/tau)+3*digamma(3/tau))   
        j <- (log(F.T(1/(sigma*abs(nu)),tau+0.001))-log(F.T(1/(sigma*abs(nu)),tau)))/0.001
     dldt <- (1/tau)-0.5*(log(abs(z/c)))*(abs(z/c))^tau+
                      (1/tau^2)*(log(2)+digamma(1/tau))+
                      ((tau/2)*((abs(z/c))^tau)-1)*dlogc.dt-j               
     dldt                                       
                        } ,
      d2ldt2 = function(y,mu,sigma,nu,tau) {
        p <- (tau+1)/tau
 dlogc.dt <- (1/(2*tau^2))*(2*log(2)-digamma(1/tau)+3*digamma(3/tau))         
    part1 <- p*trigamma(p)+2*(digamma(p))^2
    part2 <- digamma(p)*(log(2)+3-3*digamma(3/tau)-tau)
    part3 <- -3*(digamma(3/tau))*(1+log(2))    
    part4 <- -(tau+log(2))*log(2)
    part5 <- -tau+(tau^4)*(dlogc.dt)^2
   d2ldt2 <- part1+part2+part3+part4+part5
   d2ldt2 <- -d2ldt2/tau^3    
   d2ldt2 <- ifelse(d2ldt2 < -1e-15, d2ldt2,-1e-15)                                    
   d2ldt2
                            } ,
  d2ldmdd = function(mu,sigma,nu,tau) -(nu*tau)/(mu*sigma),
  d2ldmdv = function(mu,sigma,nu,tau) (2*(tau-1)-(tau+1)*(sigma^2)*(nu^2))/(4*mu),
  d2ldmdt = function(mu,sigma,nu,tau) (nu/(mu*tau))*(1+tau+(3/2)*(digamma(1/tau)-digamma(3/tau))),
  d2ldddv = function(sigma,nu,tau) -(sigma*nu*tau)/2,
  d2ldddt = function(sigma,tau) (1/(sigma*tau))*(1+tau+(3/2)*(digamma(1/tau)-digamma(3/tau))),
  d2ldvdt = function(mu,sigma,nu,tau) {
    d2ldvdt  <- (((sigma^2)*nu)/(2*tau))*(1+(tau/3)+0.5*(digamma(1/tau)-digamma(3/tau)))
    d2ldvdt
                        },
 G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
                       -2*dBCPEo(y,mu,sigma,nu,tau,log=TRUE),                     
         rqres = expression(rqres(pfun="pBCPEo", type="Continuous", y=y, mu=mu, 
                                              sigma=sigma, nu=nu, tau=tau) ),
    mu.initial = expression(mu <- (y+mean(y))/2), #
 sigma.initial = expression(sigma<- rep(0.1, length(y))),
    nu.initial = expression(nu <- rep(1, length(y))), 
   tau.initial = expression(tau <-rep(2, length(y))), 
      mu.valid = function(mu) all(mu > 0), 
   sigma.valid = function(sigma)  all(sigma > 0),
      nu.valid = function(nu) TRUE , 
     tau.valid = function(tau) all(tau > 0), 
       y.valid = function(y)  all(y > 0)
          ),
            class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------