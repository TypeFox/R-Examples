# Monday, March 29, 2004 at 16:50 
# the first derivatives squares have been used here
SEP <- function (mu.link="identity", sigma.link="log", nu.link ="identity", tau.link="log")
{
    mstats <- checklink(   "mu.link", "Skew Exponential Power", substitute(mu.link), 
                           c("1/mu^2", "log", "identity"))
    dstats <- checklink("sigma.link", "Skew Exponential Power", substitute(sigma.link), 
                           c("inverse", "log", "identity"))
    vstats <- checklink(   "nu.link", "Skew Exponential Power", substitute(nu.link),    
                           c("1/nu^2", "log", "identity"))
    tstats <- checklink(  "tau.link", "Skew Exponential Power", substitute(tau.link),   
                           c("1/tau^2", "log", "identity")) 
    structure(
          list(family = c("SEP", "Skew Exponential Power"),
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
    -2*dSEP(y,mu,sigma,nu,tau,log=TRUE)
                        } ,                     
         rqres = expression(   
               rqres(pfun="pSEP", type="Continuous", y=y, mu=mu, 
                                sigma=sigma, nu=nu, tau=tau)) ,
    mu.initial = expression(mu <- (y+mean(y))/2),    #(y+mean(y))/2),# rep(mean(y),length(y)) 
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
dSEP <- function(x, mu = 0, sigma = 1, nu = 0, tau = 2, log = FALSE)
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
pSEP <- function(q, mu = 0, sigma = 1, nu = 0, tau = 2, lower.tail = TRUE, log.p = FALSE)
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
                 dSEP(x, mu = mu[i], sigma = sigma[i], nu = nu[i], tau = tau[i]), -Inf, q[i] )$value #ds br 7-10-11
          }    
    if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf 
    if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf) 
    cdf
 }
#------------------------------------------------------------------------------------------
qSEP <-  function(p, mu = 0, sigma = 1, nu = 0, tau = 2, lower.tail = TRUE, log.p = FALSE,
                  lower.limit = mu-5*sigma,
                  upper.limit = mu+5*sigma)
 {  
 
     #---Golded section search--------------------------------------------
       find.q.from.p <- function(p, mu, sigma, nu, tau,
                                 lower= lower.limit, 
                                 upper = upper.limit)
            {   
               usemode <- function(q,p)
                    { 
                      np <- pSEP(q , mu = mu, sigma = sigma, nu = nu, tau = tau )
                      fun <- (np-p)^2
                      fun
                    }      
            tol <-  0.000001
              r <- 0.61803399  
              b <- r*lower + (1-r)*upper 
             lo <- lower
             up <- upper 
             w1 <- TRUE 
           val1 <- if(up-b > b-lo) b   else b-(1-r)*(b-lo)
           val2 <- if(up-b > b-lo) b+(1-r)*(up-b) else b 
             f1 <- usemode(val1,p)
             f2 <- usemode(val2,p)
            while(w1) 
                 { if(f2 < f1) { lo <- val1 
                               val1 <- val2 
                               val2 <- r*val1+(1-r)*up
                                 f1 <- f2 
                                 f2 <- usemode(val2,p) 
                               } 
                   else        { up <- val2 
                               val2 <- val1
                               val1 <- r*val2+(1-r)*lo
                                 f2 <- f1 
                                 f1 <- usemode(val1,p)
                               }
                   w1 <- abs(up-lo) >  tol*(abs(val1)+abs(val2))
                 }             
            q <- if(f1<f2) val1 else val2
            q
            }
     #----------------------------------------------------------------- 
    if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
    if (lower.tail==TRUE) p <- p else p <- 1-p
         lp <-  pmax.int(length(p), length(mu), length(sigma), length(nu), length(tau))
          p <- rep(p, length = lp)
      sigma <- rep(sigma, length = lp)
         mu <- rep(mu, length = lp)
         nu <- rep(nu, length = lp)
        tau <- rep(tau, length = lp)
      upper <- rep(upper.limit, length = lp )
      lower <- rep(lower.limit, length = lp )
          q <- rep(0,lp)    
          for (i in 1:lp)
          {
          q[i] <- find.q.from.p(p[i], mu = mu[i], sigma = sigma[i], 
                                      nu = nu[i], tau = tau[i],
                                upper = upper[i], 
                                lower = lower[i])
          if (q[i]>=upper[i]) warning("q is at the upper limit, increase the upper.limit")
          if (q[i]<=lower[i]) warning("q is at the lower limit, decrease the lower.limit")
          }                                                                               
    q 
 }
#------------------------------------------------------------------------------------------
rSEP <- function(n, mu=0, sigma=1, nu=0, tau=2)
  {
 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    n <- ceiling(n)
    p <- runif(n)
    r <- qSEP(p, mu = mu,sigma = sigma, nu = nu,tau = tau)
    r
  }
#-----------------------------------------------------------------  
