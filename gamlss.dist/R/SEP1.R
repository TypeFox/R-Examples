# 28_11_2007
# the first derivatives squares have been used here
SEP1 <- function (mu.link="identity", sigma.link="log", nu.link ="identity", tau.link="log")
{
    mstats <- checklink(   "mu.link", "Skew exponential power (Azzalini type 1)", substitute(mu.link), 
                           c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Skew exponential power (Azzalini type 1)", substitute(sigma.link), 
                           c("inverse", "log", "identity", "own"))
    vstats <- checklink(   "nu.link", "Skew exponential power (Azzalini type 1)", substitute(nu.link),    
                           c("inverse", "log", "identity", "own"))
    tstats <- checklink(  "tau.link", "Skew exponential power (Azzalini type 1)", substitute(tau.link),   
                           c("inverse", "log", "identity", "own")) 
    structure(
          list(family = c("SEP1", "Skew exponential power (Azzalini type 1)"),
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
      w <- nu*z
      s <- ((abs(w))^tau)/tau 
   lpdf <- (1-(1/tau))*log(tau)-s-lgamma(1/tau)-log(2)
   lcdf <- log(0.5*(1+pgamma(s,shape=1/tau,scale=1)*sign(w)))
   dldm <- -(exp(lpdf-lcdf))*nu/sigma + sign(z)*(abs(z)^(tau-1))/sigma
   dldm
                    },
   d2ldm2 = function(y,mu,sigma,nu,tau){
      z <- (y-mu)/sigma
      w <- nu*z
      s <- ((abs(w))^tau)/tau 
   lpdf <- (1-(1/tau))*log(tau)-s-lgamma(1/tau)-log(2)
   lcdf <- log(0.5*(1+pgamma(s,shape=1/tau,scale=1)*sign(w)))
   dldm <- -(exp(lpdf-lcdf))*nu/sigma + sign(z)*(abs(z)^(tau-1))/sigma
   d2ldm2 <- -dldm*dldm
   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)
   d2ldm2  
                      }, 
   dldd = function(y,mu,sigma,nu,tau) {  
      z <- (y-mu)/sigma
      w <- nu*z
      s <- ((abs(w))^tau)/tau 
   lpdf <- (1-(1/tau))*log(tau)-s-lgamma(1/tau)-log(2)
   lcdf <- log(0.5*(1+pgamma(s,shape=1/tau,scale=1)*sign(w)))
   dldd <- -(exp(lpdf-lcdf))*nu*z/sigma + ((abs(z)^tau)-1)/sigma
   dldd
      } ,
   d2ldd2 = function(y,mu,sigma,nu,tau){
      z <- (y-mu)/sigma
      w <- nu*z
      s <- ((abs(w))^tau)/tau 
   lpdf <- (1-(1/tau))*log(tau)-s-lgamma(1/tau)-log(2)
   lcdf <- log(0.5*(1+pgamma(s,shape=1/tau,scale=1)*sign(w)))
   dldd <- -(exp(lpdf-lcdf))*nu*z/sigma + ((abs(z)^tau)-1)/sigma
      d2ldd2 <- -dldd*dldd
      d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
      d2ldd2
                      },   
     dldv = function(y,mu,sigma,nu,tau) { 
      z <- (y-mu)/sigma
      w <- nu*z
   dwdv <- w/nu 
      s <- ((abs(w))^tau)/tau 
   lpdf <- (1-(1/tau))*log(tau)-s-lgamma(1/tau)-log(2)
   lcdf <- log(0.5*(1+pgamma(s,shape=1/tau,scale=1)*sign(w)))
   dldv <- (exp(lpdf-lcdf))*dwdv
   dldv
                        } ,
    d2ldv2 = function(y,mu,sigma,nu,tau) { 
      z <- (y-mu)/sigma
      w <- nu*z
   dwdv <- w/nu 
      s <- ((abs(w))^tau)/tau 
   lpdf <- (1-(1/tau))*log(tau)-s-lgamma(1/tau)-log(2)
   lcdf <- log(0.5*(1+pgamma(s,shape=1/tau,scale=1)*sign(w)))
   dldv <- (exp(lpdf-lcdf))*dwdv
 d2ldv2 <-  -dldv*dldv             
 d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)  
 d2ldv2
                        },
      dldt = function(y,mu,sigma,nu,tau) {
      z <- (y-mu)/sigma
      w <- nu*z 
     t1 <- tau+0.00002
     t2 <- tau-0.00002
     s <-  ((abs(w))^tau)/tau
     s1 <- ((abs(w))^t1)/t1
     s2 <- ((abs(w))^t2)/t2
   dldf <- (log(tau)+tau-1+digamma(1/tau)-tau*((abs(z))^tau)*log(abs(z))+((abs(z))^tau))/(tau*tau)
dldcdf1 <- sign(w)*(log(1+pgamma(s1,shape=1/t1,scale=1))-log(1+pgamma(s2,shape=1/t2,scale=1)))/0.00004
  suppressWarnings(lcdf1 <- log(0.5 + w*exp((1-(1/t1))*log(t1)-lgamma(1/t1)-log(2))))
  suppressWarnings(lcdf2 <- log(0.5 + w*exp((1-(1/t2))*log(t2)-lgamma(1/t2)-log(2))))
dldcdf2 <- (lcdf1-lcdf2)/0.00004  
 dldcdf <- ifelse((s==0),dldcdf2,dldcdf1)
   dldt <- dldf + dldcdf
   dldt
                       } ,
      d2ldt2 = function(y,mu,sigma,nu,tau) 
                       { 
      z <- (y-mu)/sigma
      w <- nu*z 
     t1 <- tau+0.00002
     t2 <- tau-0.00002
      s <-  ((abs(w))^tau)/tau
     s1 <- ((abs(w))^t1)/t1
     s2 <- ((abs(w))^t2)/t2
   dldf <- (log(tau)+tau-1+digamma(1/tau)-tau*((abs(z))^tau)*log(abs(z))+((abs(z))^tau))/(tau*tau)
dldcdf1 <- sign(w)*(log(1+pgamma(s1,shape=1/t1,scale=1))-log(1+pgamma(s2,shape=1/t2,scale=1)))/0.00004
  suppressWarnings(lcdf1 <- log(0.5 + w*exp((1-(1/t1))*log(t1)-lgamma(1/t1)-log(2))))
  suppressWarnings(lcdf2 <- log(0.5 + w*exp((1-(1/t2))*log(t2)-lgamma(1/t2)-log(2))))
dldcdf2 <- (lcdf1-lcdf2)/0.00004  
 dldcdf <- ifelse((s==0),dldcdf2,dldcdf1)
   dldt <- dldf + dldcdf    
 d2ldt2 <-   -dldt*dldt   
 d2ldt2 <- ifelse(d2ldt2 < -1e-15, d2ldt2,-1e-15)  
 d2ldt2
                       },
  d2ldmdd = function(y,mu,sigma,nu,tau) {
      z <- (y-mu)/sigma
      w <- nu*z
      s <- ((abs(w))^tau)/tau 
   lpdf <- (1-(1/tau))*log(tau)-s-lgamma(1/tau)-log(2)
   lcdf <- log(0.5*(1+pgamma(s,shape=1/tau,scale=1)*sign(w)))
   dldm <- -(exp(lpdf-lcdf))*nu/sigma + sign(z)*(abs(z)^(tau-1))/sigma
   dldd <- -(exp(lpdf-lcdf))*nu*z/sigma + ((abs(z)^tau)-1)/sigma
    d2ldmdd <- -(dldm*dldd)
    d2ldmdd
    
                       },
  d2ldmdv = function(y,mu,sigma,nu,tau) {
     z <- (y-mu)/sigma
      w <- nu*z
      s <- ((abs(w))^tau)/tau 
   lpdf <- (1-(1/tau))*log(tau)-s-lgamma(1/tau)-log(2)
   lcdf <- log(0.5*(1+pgamma(s,shape=1/tau,scale=1)*sign(w)))
   dldm <- -(exp(lpdf-lcdf))*nu/sigma + sign(z)*(abs(z)^(tau-1))/sigma
   dwdv <- w/nu 
   dldv <- (exp(lpdf-lcdf))*dwdv
   d2ldmdv <- -(dldm*dldv)
   d2ldmdv
                       },    
  d2ldmdt = function(y,mu,sigma,nu,tau) {
     z <- (y-mu)/sigma
      w <- nu*z
      s <- ((abs(w))^tau)/tau 
   lpdf <- (1-(1/tau))*log(tau)-s-lgamma(1/tau)-log(2)
   lcdf <- log(0.5*(1+pgamma(s,shape=1/tau,scale=1)*sign(w)))
   dldm <- -(exp(lpdf-lcdf))*nu/sigma + sign(z)*(abs(z)^(tau-1))/sigma
     t1 <- tau+0.00002
     t2 <- tau-0.00002
     s1 <- ((abs(w))^t1)/t1
     s2 <- ((abs(w))^t2)/t2
   dldf <- (log(tau)+tau-1+digamma(1/tau)-tau*((abs(z))^tau)*log(abs(z))+((abs(z))^tau))/(tau*tau)
dldcdf1 <- sign(w)*(log(1+pgamma(s1,shape=1/t1,scale=1))-log(1+pgamma(s2,shape=1/t2,scale=1)))/0.00004
  lcdf1 <- log(0.5 + w*exp((1-(1/t1))*log(t1)-lgamma(1/t1)-log(2)))
  lcdf2 <- log(0.5 + w*exp((1-(1/t2))*log(t2)-lgamma(1/t2)-log(2)))
dldcdf2 <- (lcdf1-lcdf2)/0.00004  
 dldcdf <- ifelse((s==0),dldcdf2,dldcdf1)
   dldt <- dldf + dldcdf
   d2ldmdt <- -(dldm*dldt)
   d2ldmdt
                       },
  d2ldddv = function(y,mu,sigma,nu,tau) {
      z <- (y-mu)/sigma
      w <- nu*z
      s <- ((abs(w))^tau)/tau 
   lpdf <- (1-(1/tau))*log(tau)-s-lgamma(1/tau)-log(2)
   lcdf <- log(0.5*(1+pgamma(s,shape=1/tau,scale=1)*sign(w)))
   dldd <- -(exp(lpdf-lcdf))*nu*z/sigma + ((abs(z)^tau)-1)/sigma
   dwdv <- w/nu 
   dldv <- (exp(lpdf-lcdf))*dwdv
 d2ldddv <- -(dldd*dldv)
  d2ldddv
                       },
  d2ldddt = function(y,mu,sigma,nu,tau) {
      z <- (y-mu)/sigma
      w <- nu*z
      s <- ((abs(w))^tau)/tau 
   lpdf <- (1-(1/tau))*log(tau)-s-lgamma(1/tau)-log(2)
   lcdf <- log(0.5*(1+pgamma(s,shape=1/tau,scale=1)*sign(w)))
   dldd <- -(exp(lpdf-lcdf))*nu*z/sigma + ((abs(z)^tau)-1)/sigma
     t1 <- tau+0.00002
     t2 <- tau-0.00002
     s1 <- ((abs(w))^t1)/t1
     s2 <- ((abs(w))^t2)/t2
   dldf <- (log(tau)+tau-1+digamma(1/tau)-tau*((abs(z))^tau)*log(abs(z))+((abs(z))^tau))/(tau*tau)
dldcdf1 <- sign(w)*(log(1+pgamma(s1,shape=1/t1,scale=1))-log(1+pgamma(s2,shape=1/t2,scale=1)))/0.00004
  lcdf1 <- log(0.5 + w*exp((1-(1/t1))*log(t1)-lgamma(1/t1)-log(2)))
  lcdf2 <- log(0.5 + w*exp((1-(1/t2))*log(t2)-lgamma(1/t2)-log(2)))
dldcdf2 <- (lcdf1-lcdf2)/0.00004  
 dldcdf <- ifelse((s==0),dldcdf2,dldcdf1)
   dldt <- dldf + dldcdf
d2ldddt <- -(dldd*dldt)
d2ldddt   
                       },
  d2ldvdt = function(y,mu,sigma,nu,tau) {
      z <- (y-mu)/sigma
      w <- nu*z
   dwdv <- w/nu 
      s <- ((abs(w))^tau)/tau 
   lpdf <- (1-(1/tau))*log(tau)-s-lgamma(1/tau)-log(2)
   lcdf <- log(0.5*(1+pgamma(s,shape=1/tau,scale=1)*sign(w)))
   dldv <- (exp(lpdf-lcdf))*dwdv
     t1 <- tau+0.00002
     t2 <- tau-0.00002
     s1 <- ((abs(w))^t1)/t1
     s2 <- ((abs(w))^t2)/t2
   dldf <- (log(tau)+tau-1+digamma(1/tau)-tau*((abs(z))^tau)*log(abs(z))+((abs(z))^tau))/(tau*tau)
dldcdf1 <- sign(w)*(log(1+pgamma(s1,shape=1/t1,scale=1))-log(1+pgamma(s2,shape=1/t2,scale=1)))/0.00004
  lcdf1 <- log(0.5 + w*exp((1-(1/t1))*log(t1)-lgamma(1/t1)-log(2)))
  lcdf2 <- log(0.5 + w*exp((1-(1/t2))*log(t2)-lgamma(1/t2)-log(2)))
dldcdf2 <- (lcdf1-lcdf2)/0.00004  
 dldcdf <- ifelse((s==0),dldcdf2,dldcdf1)
   dldt <- dldf + dldcdf
d2ldvdt <- -(dldv*dldt)
d2ldvdt  
                       },
G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
                       { 
   -2*dSEP1(y,mu,sigma,nu,tau,log=TRUE)
                        } ,                     
         rqres = expression(   
               rqres(pfun="pSEP1", type="Continuous", y=y, mu=mu, 
                                sigma=sigma, nu=nu, tau=tau)) ,
    mu.initial = expression(mu <- (y+mean(y))/2),    # rep(mean(y),length(y)) 
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
dSEP1 <- function(x, mu = 0, sigma = 1, nu = 0, tau = 2, log = FALSE)
 {
          if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
      z <- (x-mu)/sigma
      w <- nu*z
     sz <- ((abs(z))^tau)/tau 
      s <- ((abs(w))^tau)/tau 
    lpdf <- (1-(1/tau))*log(tau)-sz-lgamma(1/tau)-log(2)
    lcdf1 <- log(0.5*(1+pgamma(s,shape=1/tau,scale=1)*sign(w)))
    cdf2 <- 0.5 + w*exp((1-(1/tau))*log(tau)-lgamma(1/tau)-log(2))
     suppressWarnings(lcdf2 <- log(cdf2))
    lcdf <- ifelse((s==0),lcdf2,lcdf1) # note this change, note ifelse used here 
 loglik <- lpdf + lcdf + log(2) - log(sigma)
       if(log==FALSE) ft  <- exp(loglik) else ft <- loglik 
       ft
  }    
#------------------------------------------------------------------------------------------
pSEP1 <- function(q, mu = 0, sigma = 1, nu = 0, tau = 2, lower.tail = TRUE, log.p = FALSE)
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
                 dSEP1(x, mu = 0, sigma = 1, nu = nu[i], tau = tau[i]), -Inf, (q[i]-mu[i])/sigma[i] )$value #DS BR 7-10-11
          }    
    if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf 
    if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf) 
    cdf
 }
#------------------------------------------------------------------------------------------
qSEP1 <-  function(p, mu = 0, sigma = 1, nu = 0, tau = 2, lower.tail = TRUE, log.p = FALSE)
  { 
    #---functions--------------------------------------------   
       h1 <- function(q)
       { 
     pSEP1(q , mu = mu[i], sigma = sigma[i], nu = nu[i], tau = tau[i]) - p[i] 
       }
       h <- function(q)
       { 
     pSEP1(q , mu = mu[i], sigma = sigma[i], nu = nu[i], tau = tau[i]) 
       }
     #-----------------------------------------------------------------
    #if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))      
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
#------------------------------------------------------------------------------------------
rSEP1 <- function(n, mu=0, sigma=1, nu=0, tau=2)
  {
 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    n <- ceiling(n)
    p <- runif(n)
    r <- qSEP1(p, mu = mu,sigma = sigma, nu = nu,tau = tau)
    r
  }
#-----------------------------------------------------------------  
