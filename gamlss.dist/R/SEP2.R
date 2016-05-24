# 27_11_2007  
# note this distribution is the same as SEP
# the first derivatives squares have been used here
SEP2 <- function (mu.link="identity", sigma.link="log", nu.link ="identity", tau.link="log")
{
    mstats <- checklink(   "mu.link", "Skew Exponential Power", substitute(mu.link), 
                           c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Skew Exponential Power", substitute(sigma.link), 
                           c("inverse", "log", "identity", "own"))
    vstats <- checklink(   "nu.link", "Skew Exponential Power", substitute(nu.link),    
                           c("inverse", "log", "identity", "own"))
    tstats <- checklink(  "tau.link", "Skew Exponential Power", substitute(tau.link),   
                           c("inverse", "log", "identity", "own")) 
    structure(
          list(family = c("SEP2", "Skew Exponential Power type 2"),
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
      w <- sign(z)*(abs(z)^(tau/2))*nu*(sqrt(2/tau))
   dwdz <- (abs(z)^((tau/2)-1))*nu*(sqrt(tau/2)) 
   dldm <- -(dnorm(w)/pnorm(w))*dwdz/sigma + sign(z)*(abs(z)^(tau-1))/sigma
   dldm
                    },
   d2ldm2 = function(y,mu,sigma,nu,tau){
      z <- (y-mu)/sigma
      w <- sign(z)*(abs(z)^(tau/2))*nu*(sqrt(2/tau))
   dwdz <- (abs(z)^((tau/2)-1))*nu*(sqrt(tau/2)) 
   dldm <- -(dnorm(w)/pnorm(w))*dwdz/sigma + sign(z)*(abs(z)^(tau-1))/sigma
   d2ldm2 <- -dldm*dldm
   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)                                    
   d2ldm2
                      }, 
   dldd = function(y,mu,sigma,nu,tau) {  
      z <- (y-mu)/sigma
      w <- sign(z)*(abs(z)^(tau/2))*nu*(sqrt(2/tau))
   dwdz <- (abs(z)^((tau/2)-1))*nu*(sqrt(tau/2)) 
   dldd <- -(dnorm(w)/pnorm(w))*dwdz*z/sigma + ((abs(z)^(tau))-1)/sigma
   dldd
      } ,
   d2ldd2 = function(y,mu,sigma,nu,tau){
      z <- (y-mu)/sigma
      w <- sign(z)*(abs(z)^(tau/2))*nu*(sqrt(2/tau))
   dwdz <- (abs(z)^((tau/2)-1))*nu*(sqrt(tau/2)) 
   dldd <- -(dnorm(w)/pnorm(w))*dwdz*z/sigma + ((abs(z)^(tau))-1)/sigma
      d2ldd2 <- -dldd*dldd
   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
   d2ldd2
                      },   
     dldv = function(y,mu,sigma,nu,tau) { 
      z <- (y-mu)/sigma
      w <- sign(z)*(abs(z)^(tau/2))*nu*(sqrt(2/tau))
   dwdv <- w/nu 
   dldv <- (dnorm(w)/pnorm(w))*dwdv
   dldv
                        } ,
    d2ldv2 = function(y,mu,sigma,nu,tau) { 
      z <- (y-mu)/sigma
      w <- sign(z)*(abs(z)^(tau/2))*nu*(sqrt(2/tau))
   dwdv <- w/nu 
   dldv <- (dnorm(w)/pnorm(w))*dwdv
      d2ldv2 <-  -dldv*dldv             
   d2ldv2 <- ifelse(d2ldv2 < -1e-4, d2ldv2,-1e-4)  
   d2ldv2
                        },
      dldt = function(y,mu,sigma,nu,tau) {
      z <- (y-mu)/sigma
      w <- sign(z)*(abs(z)^(tau/2))*nu*(sqrt(2/tau))
   dwdt <- (log(abs(z))-1/tau)*w/2 
   dldt <- (dnorm(w)/pnorm(w))*dwdt 
   dldt <- dldt+(log(tau)+tau-1+digamma(1/tau)-tau*((abs(z))^tau)*log(abs(z))+((abs(z))^tau))/(tau*tau)
   dldt
                        } ,
      d2ldt2 = function(y,mu,sigma,nu,tau) 
                       { 
      z <- (y-mu)/sigma
      w <- sign(z)*(abs(z)^(tau/2))*nu*(sqrt(2/tau))
   dwdt <- (log(abs(z))-1/tau)*w/2 
   dldt <- (dnorm(w)/pnorm(w))*dwdt 
   dldt <- dldt+(log(tau)+tau-1+digamma(1/tau)-tau*((abs(z))^tau)*log(abs(z))+((abs(z))^tau))/(tau*tau)
      d2ldt2 <-   -dldt*dldt   
    d2ldt2 <- ifelse(d2ldt2 < -1e-4, d2ldt2,-1e-4)                                    
    d2ldt2
                       },
  d2ldmdd = function(y,mu,sigma,nu,tau) {
         z <- (y-mu)/sigma
         w <- sign(z)*(abs(z)^(tau/2))*nu*(sqrt(2/tau))
      dwdz <- (abs(z)^((tau/2)-1))*nu*(sqrt(tau/2)) 
      dldm <- -(dnorm(w)/pnorm(w))*dwdz/sigma + sign(z)*(abs(z)^(tau-1))/sigma 
      dldd <- -(dnorm(w)/pnorm(w))*dwdz*z/sigma + ((abs(z)^(tau))-1)/sigma
   d2ldmdd <- -(dldm*dldd)
   d2ldmdd
   
                       },
  d2ldmdv = function(y,mu,sigma,nu,tau) {
        z <- (y-mu)/sigma
        w <- sign(z)*(abs(z)^(tau/2))*nu*(sqrt(2/tau))
     dwdz <- (abs(z)^((tau/2)-1))*nu*(sqrt(tau/2)) 
     dldm <- -(dnorm(w)/pnorm(w))*dwdz/sigma + sign(z)*(abs(z)^(tau-1))/sigma 
     dwdv <- w/nu 
     dldv <- (dnorm(w)/pnorm(w))*dwdv
  d2ldmdv <- -(dldm*dldv)
  d2ldmdv
                       },    
  d2ldmdt = function(y,mu,sigma,nu,tau) {
        z <- (y-mu)/sigma
        w <- sign(z)*(abs(z)^(tau/2))*nu*(sqrt(2/tau))
     dwdz <- (abs(z)^((tau/2)-1))*nu*(sqrt(tau/2)) 
     dldm <- -(dnorm(w)/pnorm(w))*dwdz/sigma + sign(z)*(abs(z)^(tau-1))/sigma 
     dwdt <- (log(abs(z))-1/tau)*w/2 
     dldt <- (dnorm(w)/pnorm(w))*dwdt 
     dldt <- dldt+(log(tau)+tau-1+digamma(1/tau)-tau*((abs(z))^tau)*log(abs(z))+((abs(z))^tau))/(tau*tau)   
  d2ldmdt <- -(dldm*dldt)
  d2ldmdt
                       },
  d2ldddv = function(y,mu,sigma,nu,tau) {
       z <- (y-mu)/sigma
       w <- sign(z)*(abs(z)^(tau/2))*nu*(sqrt(2/tau))
    dwdz <- (abs(z)^((tau/2)-1))*nu*(sqrt(tau/2)) 
    dldd <- -(dnorm(w)/pnorm(w))*dwdz*z/sigma + ((abs(z)^(tau))-1)/sigma
    dwdv <- w/nu 
    dldv <- (dnorm(w)/pnorm(w))*dwdv  
 d2ldddv <- -(dldd*dldv)
 d2ldddv
                       },
  d2ldddt = function(y,mu,sigma,nu,tau) {
      z <- (y-mu)/sigma
      w <- sign(z)*(abs(z)^(tau/2))*nu*(sqrt(2/tau))
   dwdz <- (abs(z)^((tau/2)-1))*nu*(sqrt(tau/2)) 
   dldd <- -(dnorm(w)/pnorm(w))*dwdz*z/sigma + ((abs(z)^(tau))-1)/sigma
   dwdt <- (log(abs(z))-1/tau)*w/2 
   dldt <- (dnorm(w)/pnorm(w))*dwdt 
   dldt <- dldt+(log(tau)+tau-1+digamma(1/tau)-tau*((abs(z))^tau)*log(abs(z))+((abs(z))^tau))/(tau*tau)
d2ldddt <- -(dldd*dldt)
d2ldddt  
                       },
  d2ldvdt = function(y,mu,sigma,nu,tau) {
      z <- (y-mu)/sigma
      w <- sign(z)*(abs(z)^(tau/2))*nu*(sqrt(2/tau))
   dwdv <- w/nu 
   dldv <- (dnorm(w)/pnorm(w))*dwdv
   dwdt <- (log(abs(z))-1/tau)*w/2 
   dldt <- (dnorm(w)/pnorm(w))*dwdt 
   dldt <- dldt+(log(tau)+tau-1+digamma(1/tau)-tau*((abs(z))^tau)*log(abs(z))+((abs(z))^tau))/(tau*tau)
d2ldvdt <- -(dldv*dldt)
d2ldvdt 
                       },
 G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
                       { 
   -2*dSEP2(y,mu,sigma,nu,tau,log=TRUE)
                        } ,                     
         rqres = expression(   
               rqres(pfun="pSEP2", type="Continuous", y=y, mu=mu, 
                                sigma=sigma, nu=nu, tau=tau)) ,
    mu.initial = expression(mu <- (y+mean(y))/2),    
 sigma.initial = expression(sigma <- rep(sd(y)/4, length(y))),
    nu.initial = expression(nu <- rep(0.1, length(y))), 
   tau.initial = expression(tau <-rep(1.6, length(y))), 
      mu.valid = function(mu) TRUE, 
   sigma.valid = function(sigma)  all(sigma > 0),
      nu.valid = function(nu) TRUE , 
     tau.valid = function(tau) all(tau > 0), 
       y.valid = function(y)  TRUE
          ),
            class = c("gamlss.family","family"))
}
#------------------------------------------------------------------------------------------
dSEP2 <- function(x, mu = 0, sigma = 1, nu = 0, tau = 2, log = FALSE)
 {
          if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
      z <- (x-mu)/sigma
      w <- sign(z)*(abs(z)^(tau/2))*nu*(sqrt(2/tau)) 
 loglik <- log(pnorm(w)) - (abs(z)^(tau))/tau - log(sigma) - lgamma(1/tau) - ((1/tau)-1)*log(tau)
       if(log==FALSE) ft  <- exp(loglik) else ft <- loglik 
       ft
  }    
#------------------------------------------------------------------------------------------
pSEP2 <- function(q, mu = 0, sigma = 1, nu = 0, tau = 2, lower.tail = TRUE, log.p = FALSE)
 {   if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", "")) 
     if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))    
         lp <- pmax.int(length(q), length(mu), length(sigma), length(nu), length(tau))                                                                  
          q <- rep(q, length = lp)
      sigma <- rep(sigma, length = lp)
         mu <- rep(mu, length = lp)
         nu <- rep(nu, length = lp)
        tau <- rep(tau, length = lp)
        cdf <- rep(0, length = lp)
       for (i in 1:lp)
          {
         cdf[i] <- integrate(function(x) 
                 dSEP2(x, mu = 0, sigma = 1, nu = nu[i], tau = tau[i]), -Inf, (q[i]-mu[i])/sigma[i] )$value # ds rb 7-10-11
          }    
    if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf 
    if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf) 
    cdf
 }
#------------------------------------------------------------------------------------------
qSEP2 <-  function(p, mu = 0, sigma = 1, nu = 0, tau = 2, lower.tail = TRUE, log.p = FALSE)
  { 
    #---functions--------------------------------------------   
       h1 <- function(q)
       { 
     pSEP2(q , mu = mu[i], sigma = sigma[i], nu = nu[i], tau = tau[i]) - p[i]  
       }
       h <- function(q)
       { 
     pSEP2(q , mu = mu[i], sigma = sigma[i], nu = nu[i], tau = tau[i])  
       }
     #-----------------------------------------------------------------
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))      
    if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (lower.tail==TRUE) p <- p else p <- 1-p
    if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))     
         lp <-  max(length(p),length(mu),length(sigma),length(nu), length(tau))
          p <- rep(p, length = lp)                                                                     
      sigma <- rep(sigma, length = lp)
         mu <- rep(mu, length = lp)
         nu <- rep(nu, length = lp)
         tau <- rep(tau, length = lp)
          q <- rep(0,lp)  
         for (i in  seq(along=p)) 
         {
         if (h(mu[i])<p[i]) 
          { 
           interval <- c(mu[i], mu[i]+sigma[i])
           j <-2
           while (h(interval[2]) < p[i]) 
              {interval[2]<- mu[i]+j*sigma[i]
              j<-j+1 
              }
           } 
          else  
           {
           interval <-  c(mu[i]-sigma[i], mu[i])
           j <-2
           while (h(interval[1]) > p[i]) 
              {interval[1]<- mu[i]-j*sigma[i]
              j<-j+1 
              }
           }
        q[i] <- uniroot(h1, interval)$root
        #interval <- c(.Machine$double.xmin, 20)
         }
    q
   }
#----------------------------------------------------------------------------------------
rSEP2 <- function(n, mu=0, sigma=1, nu=0, tau=2)
  {
 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
    n <- ceiling(n)
    p <- runif(n)
    r <- qSEP2(p, mu = mu,sigma = sigma, nu = nu,tau = tau)
    r
  }
#-----------------------------------------------------------------  
