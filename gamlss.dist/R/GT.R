#  amended 27_11_2007
GT <- function (mu.link="identity", sigma.link="log", nu.link ="log", tau.link="log")
{
    mstats <- checklink("mu.link", "Generalized t", substitute(mu.link), 
                         c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Generalized t", substitute(sigma.link), 
                         c("inverse", "log", "identity", "own"))
    vstats <- checklink("nu.link", "Generalized t",substitute(nu.link), 
                         c("inverse", "log", "identity", "own"))
    tstats <- checklink("tau.link", "Generalized t ",substitute(tau.link), 
                         c("inverse", "log", "identity", "own")) 
    structure(
          list(family = c("GT",  "Generalized t"),
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
       zt <- (abs(z))^tau
       w  <- (nu*tau+1)/(nu+zt)
     dldm <- w*((abs(z))^(tau-1))*sign(z)/sigma
     dldm
                                    },
               d2ldm2 = function(y,mu,sigma,nu,tau){
        z <- (y-mu)/sigma
       zt <- (abs(z))^tau
       w  <- (nu*tau+1)/(nu+zt)
     dldm <- w*((abs(z))^(tau-1))*sign(z)/sigma
    d2ldm2 <- -dldm*dldm
    d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)                                    
    d2ldm2
                                   },
                 dldd = function(y,mu,sigma,nu,tau) {
        z <- (y-mu)/sigma
       zt <- (abs(z))^tau
       w  <- (nu*tau+1)/(nu+zt)
     dldd <- (w*zt-1)/sigma
     dldd
                                     } ,
               d2ldd2 = function(y,mu,sigma,nu,tau) {
        z <- (y-mu)/sigma
       zt <- (abs(z))^tau
       w  <- (nu*tau+1)/(nu+zt)
     dldd <- (w*zt-1)/sigma
    d2ldd2 <- -dldd*dldd
   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
   d2ldd2
                                     },
                 dldv = function(y,mu,sigma,nu,tau) {
        z <- (y-mu)/sigma
       zt <- (abs(z))^tau
        w <- (nu*tau+1)/(nu+zt)
     dldv <- (w*zt-1)/(nu*tau) - digamma(nu)+digamma(nu+(1/tau)) - log(1+(zt/nu))
     dldv
                                     } ,
               d2ldv2 = function(y,mu,sigma,nu,tau) { 
        z <- (y-mu)/sigma
       zt <- (abs(z))^tau
        w <- (nu*tau+1)/(nu+zt)
     dldv <- (w*zt-1)/(nu*tau) - digamma(nu) + digamma(nu+(1/tau)) - log(1+(zt/nu))
    d2ldv2 <- -dldv*dldv
   d2ldv2 <- ifelse(d2ldv2 < -1e-4, d2ldv2,-1e-4)  
   d2ldv2
                                     },
                 dldt = function(y,mu,sigma,nu,tau) { 
        z <- (y-mu)/sigma
       zt <- (abs(z))^tau
        w <- (nu*tau+1)/(nu+zt)
     dldt <- -(tau*w*zt*log(abs(z))) + log(1+(zt/nu)) 
     dldt <- dldt + digamma(1/tau)-digamma(nu+(1/tau))+log(nu)+tau
     dldt <- dldt/(tau^2)
     dldt
                                     } ,
               d2ldt2 = function(y,mu,sigma,nu,tau) {
        z <- (y-mu)/sigma
       zt <- (abs(z))^tau
        w <- (nu*tau+1)/(nu+zt)
     dldt <- -(tau*w*zt*log(abs(z))) + log(1+(zt/nu)) 
     dldt <- dldt + digamma(1/tau)-digamma(nu+(1/tau))+log(nu)+tau
     dldt <- dldt/(tau^2)
    d2ldt2 <- -dldt*dldt
    d2ldt2 <- ifelse(d2ldt2 < -1e-4, d2ldt2,-1e-4)                                    
    d2ldt2
                                     } ,
              d2ldmdd = function(y,mu,sigma,nu,tau)  
                     {
        z <- (y-mu)/sigma
       zt <- (abs(z))^tau
       w  <- (nu*tau+1)/(nu+zt)
     dldm <- w*((abs(z))^(tau-1))*sign(z)/sigma
     dldd <- (w*zt-1)/sigma
  d2ldmdd <- -(dldm*dldd)
  d2ldmdd                  
                      },
              d2ldmdv = function(y,mu,sigma,nu,tau)
                      {
        z <- (y-mu)/sigma
       zt <- (abs(z))^tau
       w  <- (nu*tau+1)/(nu+zt)
     dldm <- w*((abs(z))^(tau-1))*sign(z)/sigma
     dldv <- (w*zt-1)/(nu*tau) - digamma(nu)+digamma(nu+(1/tau)) - log(1+(zt/nu))
  d2ldmdv <- -(dldm*dldv)
  d2ldmdv
                      } ,
              d2ldmdt = function(y,mu,sigma,nu,tau)
                     {  
        z <- (y-mu)/sigma
       zt <- (abs(z))^tau
       w  <- (nu*tau+1)/(nu+zt)
     dldm <- w*((abs(z))^(tau-1))*sign(z)/sigma
     dldt <- -(tau*w*zt*log(abs(z))) + log(1+(zt/nu)) 
     dldt <- dldt + digamma(1/tau)-digamma(nu+(1/tau))+log(nu)+tau
     dldt <- dldt/(tau^2)
  d2ldmdt <- -(dldm*dldt)
  d2ldmdt
                      },
              d2ldddv = function(y,mu,sigma,nu,tau)  
                      {
        z <- (y-mu)/sigma
       zt <- (abs(z))^tau
       w  <- (nu*tau+1)/(nu+zt)
     dldd <- (w*zt-1)/sigma
     dldv <- (w*zt-1)/(nu*tau) - digamma(nu)+digamma(nu+(1/tau)) - log(1+(zt/nu))
  d2ldddv <- -(dldd*dldv)
  d2ldddv
                      },
              d2ldddt = function(y,mu,sigma,nu,tau)  
                      {
        z <- (y-mu)/sigma
       zt <- (abs(z))^tau
       w  <- (nu*tau+1)/(nu+zt)
     dldd <- (w*zt-1)/sigma
     dldt <- -(tau*w*zt*log(abs(z))) + log(1+(zt/nu)) 
     dldt <- dldt + digamma(1/tau)-digamma(nu+(1/tau))+log(nu)+tau
     dldt <- dldt/(tau^2)
  d2ldddt <- -(dldd*dldt)
  d2ldddt                   
                      },
              d2ldvdt = function(y,mu,sigma,nu,tau)  
                      {
        z <- (y-mu)/sigma
       zt <- (abs(z))^tau
        w <- (nu*tau+1)/(nu+zt)
     dldv <- (w*zt-1)/(nu*tau) - digamma(nu)+digamma(nu+(1/tau)) - log(1+(zt/nu))
     dldt <- -(tau*w*zt*log(abs(z))) + log(1+(zt/nu)) 
     dldt <- dldt + digamma(1/tau)-digamma(nu+(1/tau))+log(nu)+tau
     dldt <- dldt/(tau^2)
  d2ldvdt <- -(dldv*dldt)
  d2ldvdt    
                      }, 
          G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
                                  -2*dGT(y,mu,sigma,nu,tau,log=TRUE), 
                 rqres = expression(rqres(pfun="pGT", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)) ,
            mu.initial = expression(mu <- (y+mean(y))/2), 
         sigma.initial = expression(sigma<- rep(sd(y)/4, length(y))),
            nu.initial = expression(nu <- rep(5, length(y))), 
           tau.initial = expression(tau <-rep(2, length(y))), 
              mu.valid = function(mu) TRUE, 
           sigma.valid = function(sigma)  all(sigma > 0),
              nu.valid = function(nu) all(nu > 0), 
             tau.valid = function(tau) all(tau > 0), 
               y.valid = function(y)  TRUE
          ),
            class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------  
dGT <- function(x, mu=0, sigma=1, nu=3, tau=1.5, log=FALSE)
 {
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
          if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))  
           z <- (x-mu)/sigma
          zt <- (abs(z))^tau
          loglik <- log(tau)-log(2*sigma)-(1/tau)*log(nu)- lgamma(1/tau)-lgamma(nu)
          loglik <- loglik +lgamma(nu+(1/tau)) - (nu+(1/tau))*log(1+(zt/nu))
         loglik2 <- log(tau) - log(2*sigma) - lgamma(1/tau) - zt 
          if (length(nu)>1) loglik <- ifelse(nu<1000000, loglik, loglik2)
          else loglik <- if (nu<1000000) loglik else loglik2
          ft <- if(log==FALSE) exp(loglik) else loglik 
          ft
  }    
#-----------------------------------------------------------------  
pGT <- function(q, mu=0, sigma=1, nu=3, tau=1.5, lower.tail = TRUE, log.p = FALSE)
 {  
         if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
         if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
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
          endInt <- (q[i]-mu[i])/sigma[i]
         cdf[i] <- integrate(function(x) 
                 dGT(x, mu = 0, sigma = 1, nu = nu[i], tau = tau[i]), -Inf, endInt, rel.tol=.Machine$double.eps^0.4)$value
         if(endInt>0&&cdf<0.001) cdf[i] <- 1  # MS and BR 7-10-11    
          }    
         if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf 
         if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf) 
         cdf    
 }
#-----------------------------------------------------------------  
qGT <- function(p, mu=0, sigma=1, nu=3, tau=1.5, lower.tail = TRUE, log.p = FALSE)
  { 
    #---functions--------------------------------------------   
       h1 <- function(q)
       { 
     pGT(q , mu = mu[i], sigma = sigma[i], nu = nu[i], tau = tau[i]) - p[i]  
       }
       h <- function(q)
       { 
     pGT(q , mu = mu[i], sigma = sigma[i], nu = nu[i], tau = tau[i])  
       }
     #-----------------------------------------------------------------
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
        q[i] <- uniroot(h1, interval, tol=.Machine$double.eps^0.4)$root
        #interval <- c(.Machine$double.xmin, 20)
         }
    q
   }
#-----------------------------------------------------------------  
rGT <- function(n, mu=0, sigma=1, nu=3, tau=1.5)
  {
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
    if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))  
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))    
    n <- ceiling(n)
    p <- runif(n)
    r <- qGT(p,mu=mu,sigma=sigma,nu=nu,tau=tau)
    r
  }
