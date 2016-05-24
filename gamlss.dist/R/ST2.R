# 27_11_2007
ST2 <- function (mu.link="identity", sigma.link="log", nu.link ="identity", tau.link="log")
{
    mstats <- checklink(   "mu.link", "Skew t (Azzalini type 2)", substitute(mu.link), 
                           c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Skew t (Azzalini type 2)", substitute(sigma.link), 
                           c("inverse", "log", "identity", "own"))
    vstats <- checklink(   "nu.link", "Skew t (Azzalini type 2)", substitute(nu.link),    
                           c("inverse", "log", "identity", "own"))
    tstats <- checklink(  "tau.link", "Skew t (Azzalini type 2)", substitute(tau.link),   
                           c("inverse", "log", "identity", "own")) 
    structure(
          list(family = c("ST2", "Skew t (Azzalini type 2)"),
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
    lam <- ifelse(tau < 1000000, (tau+1)/(tau+(z^2)),1)
      w <- nu*(lam^0.5)*z
   dwdz <- nu*(lam^0.5)*(1-lam*(z^2)/(tau+1))
   dldm <- -(dt(w,df=tau+1)/pt(w,df=tau+1))*dwdz/sigma + lam*z/sigma
   dldm
                    },
   d2ldm2 = function(y,mu,sigma,nu,tau){
      z <- (y-mu)/sigma
    lam <- ifelse(tau < 1000000, (tau+1)/(tau+(z^2)),1)
      w <- nu*(lam^0.5)*z
   dwdz <- nu*(lam^0.5)*(1-lam*(z^2)/(tau+1))
   dldm <- -(dt(w,df=tau+1)/pt(w,df=tau+1))*dwdz/sigma + lam*z/sigma
   d2ldm2 <- -dldm*dldm
   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)
   d2ldm2  
                      }, 
   dldd = function(y,mu,sigma,nu,tau) {  
      z <- (y-mu)/sigma
lam <- ifelse(tau < 1000000, (tau+1)/(tau+(z^2)),1) 
      w <- nu*(lam^0.5)*z
   dwdz <- nu*(lam^0.5)*(1-lam*(z^2)/(tau+1))
   dldd <- -(dt(w,df=tau+1)/pt(w,df=tau+1))*z*dwdz/sigma + ((lam*(z^2))-1)/sigma
   dldd   
      } ,
   d2ldd2 = function(y,mu,sigma,nu,tau){
      z <- (y-mu)/sigma
lam <- ifelse(tau < 1000000, (tau+1)/(tau+(z^2)),1) 
      w <- nu*(lam^0.5)*z
   dwdz <- nu*(lam^0.5)*(1-lam*(z^2)/(tau+1))
   dldd <- -(dt(w,df=tau+1)/pt(w,df=tau+1))*z*dwdz/sigma + ((lam*(z^2))-1)/sigma
      d2ldd2 <- -dldd*dldd
      d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
      d2ldd2
                      },   
     dldv = function(y,mu,sigma,nu,tau) { 
      z <- (y-mu)/sigma
lam <- ifelse(tau < 1000000, (tau+1)/(tau+(z^2)),1) 
      w <- nu*(lam^0.5)*z
   dwdv <- w/nu 
   dldv <- (dt(w,df=tau+1)/pt(w,df=tau+1))*dwdv
   dldv
                        } ,
    d2ldv2 = function(y,mu,sigma,nu,tau) { 
      z <- (y-mu)/sigma
lam <- ifelse(tau < 1000000, (tau+1)/(tau+(z^2)),1) 
      w <- nu*(lam^0.5)*z
   dwdv <- w/nu 
   dldv <- (dt(w,df=tau+1)/pt(w,df=tau+1))*dwdv
      d2ldv2 <-  -dldv*dldv             
   d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)  
   d2ldv2
                        },
      dldt = function(y,mu,sigma,nu,tau) {
      z <- (y-mu)/sigma
lam <- ifelse(tau < 1000000, (tau+1)/(tau+(z^2)),1) 
      w <- nu*(lam^0.5)*z
      j <- (pt(w,df=tau+1.0005,log.p=TRUE)-pt(w,df=tau+0.9995,log.p=TRUE))/0.001
    dldt <- j + (digamma((tau+1)/2)-digamma(tau/2)-(1/tau)-log(1+(z^2)/tau)+lam*(z^2)/tau)/2
    dldt
                        } ,
      d2ldt2 = function(y,mu,sigma,nu,tau) 
                       { 
      z <- (y-mu)/sigma
lam <- ifelse(tau < 1000000, (tau+1)/(tau+(z^2)),1) 
      w <- nu*(lam^0.5)*z
      j <- (pt(w,df=tau+1.0005,log.p=TRUE)-pt(w,df=tau+0.9995,log.p=TRUE))/0.001
    dldt <- j + (digamma((tau+1)/2)-digamma(tau/2)-(1/tau)-log(1+(z^2)/tau)+lam*(z^2)/tau)/2
      d2ldt2 <-   -dldt*dldt   
# cat("mu", mu[1], "at iter", evalq(iter,envir=sys.frames()[[4]]), "\n")
#   if (evalq(iter,envir=sys.frames()[[4]])==40) browser()
   d2ldt2 <- ifelse(d2ldt2 < -1e-15, d2ldt2,-1e-15)  
   d2ldt2
                       },
  d2ldmdd = function(y,mu,sigma,nu,tau) {
      z <- (y-mu)/sigma
lam <- ifelse(tau < 1000000, (tau+1)/(tau+(z^2)),1) 
      w <- nu*(lam^0.5)*z
   dwdz <- nu*(lam^0.5)*(1-lam*(z^2)/(tau+1))
   dldm <- -(dt(w,df=tau+1)/pt(w,df=tau+1))*dwdz/sigma + lam*z/sigma
   dldd <- -(dt(w,df=tau+1)/pt(w,df=tau+1))*z*dwdz/sigma + ((lam*(z^2))-1)/sigma
   d2ldmdd <- -(dldm*dldd)
   d2ldmdd
                       },
  d2ldmdv = function(y,mu,sigma,nu,tau) {
      z <- (y-mu)/sigma
lam <- ifelse(tau < 1000000, (tau+1)/(tau+(z^2)),1) 
      w <- nu*(lam^0.5)*z
   dwdz <- nu*(lam^0.5)*(1-lam*(z^2)/(tau+1))
   dldm <- -(dt(w,df=tau+1)/pt(w,df=tau+1))*dwdz/sigma + lam*z/sigma
   dwdv <- w/nu 
   dldv <- (dt(w,df=tau+1)/pt(w,df=tau+1))*dwdv
  d2ldmdv <- -(dldm*dldv)
  d2ldmdv
                       },    
  d2ldmdt = function(y,mu,sigma,nu,tau) {
      z <- (y-mu)/sigma
lam <- ifelse(tau < 1000000, (tau+1)/(tau+(z^2)),1) 
      w <- nu*(lam^0.5)*z
   dwdz <- nu*(lam^0.5)*(1-lam*(z^2)/(tau+1))
   dldm <- -(dt(w,df=tau+1)/pt(w,df=tau+1))*dwdz/sigma + lam*z/sigma
      j <- (pt(w,df=tau+1.0005,log.p=TRUE)-pt(w,df=tau+0.9995,log.p=TRUE))/0.001
    dldt <- j + (digamma((tau+1)/2)-digamma(tau/2)-(1/tau)-log(1+(z^2)/tau)+lam*(z^2)/tau)/2
 d2ldmdt <- -(dldm*dldt)
 d2ldmdt
                       },
  d2ldddv = function(y,mu,sigma,nu,tau) {
      z <- (y-mu)/sigma
lam <- ifelse(tau < 1000000, (tau+1)/(tau+(z^2)),1) 
      w <- nu*(lam^0.5)*z
   dwdz <- nu*(lam^0.5)*(1-lam*(z^2)/(tau+1))
   dldd <- -(dt(w,df=tau+1)/pt(w,df=tau+1))*z*dwdz/sigma + ((lam*(z^2))-1)/sigma
   dwdv <- w/nu 
   dldv <- (dt(w,df=tau+1)/pt(w,df=tau+1))*dwdv
d2ldddv <- -(dldd*dldv)
d2ldddv
                       },
  d2ldddt = function(y,mu,sigma,nu,tau) {
      z <- (y-mu)/sigma
lam <- ifelse(tau < 1000000, (tau+1)/(tau+(z^2)),1) 
      w <- nu*(lam^0.5)*z
   dwdz <- nu*(lam^0.5)*(1-lam*(z^2)/(tau+1))
   dldd <- -(dt(w,df=tau+1)/pt(w,df=tau+1))*z*dwdz/sigma + ((lam*(z^2))-1)/sigma
      j <- (pt(w,df=tau+1.0005,log.p=TRUE)-pt(w,df=tau+0.9995,log.p=TRUE))/0.001
   dldt <- j + (digamma((tau+1)/2)-digamma(tau/2)-(1/tau)-log(1+(z^2)/tau)+lam*(z^2)/tau)/2
d2ldddt <- -(dldd*dldt)
d2ldddt  
                       },
  d2ldvdt = function(y,mu,sigma,nu,tau) {
      z <- (y-mu)/sigma
lam <- ifelse(tau < 1000000, (tau+1)/(tau+(z^2)),1) 
      w <- nu*(lam^0.5)*z
   dwdv <- w/nu 
   dldv <- (dt(w,df=tau+1)/pt(w,df=tau+1))*dwdv
      j <- (pt(w,df=tau+1.0005,log.p=TRUE)-pt(w,df=tau+0.9995,log.p=TRUE))/0.001
   dldt <- j + (digamma((tau+1)/2)-digamma(tau/2)-(1/tau)-log(1+(z^2)/tau)+lam*(z^2)/tau)/2
d2ldvdt <- -(dldv*dldt)
d2ldvdt  
                       },
 G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
                       { 
     -2*dST2(y,mu,sigma,nu,tau,log=TRUE)
                        } ,                     
         rqres = expression(   
               rqres(pfun="pST2", type="Continuous", y=y, mu=mu, 
                                sigma=sigma, nu=nu, tau=tau)) ,
    mu.initial = expression(mu <- (y+mean(y))/2),   
 sigma.initial = expression(sigma <- rep(sd(y)/4, length(y))),
    nu.initial = expression(nu <- rep(0.1, length(y))), 
   tau.initial = expression(tau <-rep(5, length(y))), 
      mu.valid = function(mu) TRUE, 
   sigma.valid = function(sigma)  all(sigma > 0),
      nu.valid = function(nu) TRUE , 
     tau.valid = function(tau) all(tau > 0), 
       y.valid = function(y)  TRUE
          ),
            class = c("gamlss.family","family"))
}
#------------------------------------------------------------------------------------------
dST2 <- function(x, mu = 0, sigma = 1, nu = 0, tau = 2, log = FALSE)
 {
          if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
      z <- (x-mu)/sigma
    if (length(tau)>1) lam <- ifelse(tau<1000000,(tau+1)/(tau+(z^2)),1) 
     else lam <- if (tau<1000000) (tau+1)/(tau+(z^2)) else 1
      w <- nu*(lam^0.5)*z
 loglik1 <- pt(w,df=tau+1,log.p=TRUE) + dt(z,df=tau,log =TRUE) + log(2) - log(sigma)
 loglik2 <- pNO(w,mu=0,sigma=1,log.p=TRUE) + dNO(z,mu=0,sigma=1,log =TRUE) + log(2) - log(sigma)
    if (length(tau)>1) loglik <- ifelse(tau<1000000, 
                      loglik1, 
                      loglik2) 
    else loglik <- if (tau<1000000) loglik1
                else  loglik2
       if(log==FALSE) ft  <- exp(loglik) else ft <- loglik 
       ft
  }    
#------------------------------------------------------------------------------------------
pST2 <- function(q, mu = 0, sigma = 1, nu = 0, tau = 2, lower.tail = TRUE, log.p = FALSE)
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
                 dST2(x, mu = 0, sigma = 1, nu = nu[i], tau = tau[i]), -Inf, (q[i]-mu[i])/sigma[i] )$value #ds br 7-10-11
          }    
    if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf 
    if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf) 
    cdf
 }
#------------------------------------------------------------------------------------------
#functions to find the inverse cdf numericaly
# for cdf in -Inf to +Inf
qST2 <- function(p, mu=1, sigma=1, nu = 0, tau = 2,  lower.tail = TRUE, log.p = FALSE)
  { 
    #---functions--------------------------------------------   
       h1 <- function(q)
       { 
     pST2(q , mu = mu[i], sigma = sigma[i], nu = nu[i], tau = tau[i]) - p[i] 
       }
       h <- function(q)
       { 
     pST2(q , mu = mu[i], sigma = sigma[i], nu = nu[i], tau = tau[i]) 
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
#------------------------------------------------------------------------------------------
rST2 <- function(n, mu=0, sigma=1, nu=0, tau=2)
  {
 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))    
    n <- ceiling(n)
    p <- runif(n)
    r <- qST2(p, mu = mu,sigma = sigma, nu = nu,tau = tau)
    r
  }
#-----------------------------------------------------------------  
