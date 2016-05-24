# 3/12/2004 RAR DS
# this is BEINF, the Beta distribution with probabilities for 0 and 1 (i.e. 4 parameters)
#------------------------------------------------------------------------------------------
BEINF <- function (mu.link = "logit", sigma.link = "logit", 
                   nu.link = "log", tau.link = "log")
{
    mstats <- checklink("mu.link", "BEINF", substitute(mu.link),
                            c("logit", "probit", "cloglog", "cauchit", "log", "own"))
    dstats <- checklink("sigma.link", "BEINF", substitute(sigma.link), 
                           c("logit", "probit", "cloglog", "cauchit", "log", "own"))   
    vstats <- checklink("nu.link", "BEINF", substitute(nu.link),    
                           c("inverse", "log", "identity", "own"))
    tstats <- checklink("tau.link", "BEINF", substitute(tau.link),   
                           c("inverse", "log", "identity", "own")) 
    structure(
          list(family = c("BEINF", "Beta Inflated"),
           parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE, tau=TRUE), 
                nopar = 4, 
                 type = "Mixed",
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
    dldm = function(y,mu,sigma) { a <- mu*(1-sigma^2)/(sigma^2)
                        b <- a*(1-mu)/mu
                   dldm <- ifelse(((y==0)|(y==1)),0,((1-sigma^2)/(sigma^2))*(-digamma(a)
                                    +digamma(b) +log(y) - log(1-y)))
                       dldm 
                      },
    d2ldm2 = function(y,mu,sigma) { a <- mu*(1-sigma^2)/(sigma^2)
                          b <- a*(1-mu)/mu
                     d2ldm2 <- ifelse( ((y==0)|(y==1)), 0 , 
                           -(((1-sigma^2)^2)/(sigma^4))*(trigamma(a) +trigamma(b)))
                      d2ldm2 
                        },
    dldd = function(y,mu,sigma) {  a <- mu*(1-sigma^2)/(sigma^2)
                         b <- a*(1-mu)/mu
                      dldd <- ifelse( ((y==0)|(y==1)), 0 , 
                            -(2/(sigma^3))*( mu*(-digamma(a)+digamma(a+b)+log(y))
                                  +(1-mu)*(-digamma(b)+digamma(a+b)+log(1-y)) ))
                      dldd 
                      }, 
    d2ldd2 = function(y,mu,sigma,nu,tau) { 
                         a <- mu*(1-sigma^2)/(sigma^2)
                         b <- a*(1-mu)/mu
                              d2ldd2 <- ifelse( ((y==0)|(y==1)), 0 , 
                            -(4/(sigma^6))*((mu^2)*trigamma(a) +((1-mu)^2)*trigamma(b)
                                    -trigamma(a+b)))
                   d2ldd2
                       }, 
     dldv = function(y,nu,tau)  {dldv <- ifelse(y==0,(1/nu),0) -(1/(1+nu+tau))
                         dldv
                        }, 
    d2ldv2 = function(nu,tau) {d2ldv2 <- -(1+tau)/(nu*((1+nu+tau)^2))
                         d2ldv2
                         },
     dldt = function(y,nu,tau) { dldt <- ifelse(y==1,(1/tau),0) -(1/(1+nu+tau))
                         dldt
                       }, 
                        
    d2ldt2 = function(nu,tau) {
                         d2ldt2 <- -(1+nu)/(tau*((1+nu+tau)^2))
                         d2ldt2
                        },
    d2ldmdd = function(y,mu,sigma) { a <-  mu*(1-sigma^2)/(sigma^2)
                           b <- a*(1-mu)/mu
                     d2ldmdd <- ifelse( ((y==0)|(y==1)), 0 , 
                       (2*(1-sigma^2)/(sigma^5))*(mu*trigamma(a)-(1-mu)*trigamma(b)))
                     d2ldmdd 
                         },
  d2ldmdv = function(y) {
                        d2ldmdv <- rep(0,length(y))
                        d2ldmdv
                       },

  d2ldmdt = function(y) {
                        d2ldmdt <- rep(0,length(y))
                        d2ldmdt
                       },

  d2ldddv = function(y) {
                        d2ldddv <- rep(0,length(y))
                        d2ldddv
                       },

  d2ldddt = function(y) {
                        d2ldddt <- rep(0,length(y)) 
                        d2ldddt 
                       },

  d2ldvdt = function(nu,tau) {
                        d2ldvdt <- 1/((1+nu+tau)^2)
                        d2ldvdt  
                       },

 G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
                        -2*dBEINF(y,mu,sigma,nu,tau,log=TRUE),                     
      rqres = expression(rqres(pfun="pBEINF", type="Mixed",  mass.p=c(0,1),  
                             prob.mp=cbind(nu/(1+nu+tau),tau/(1+nu+tau)), y=y, mu=mu, 
                             sigma=sigma, nu=nu, tau=tau)),
    mu.initial = expression(mu <- (y+mean(y))/2),    #(y+mean(y))/2),# rep(mean(y),length(y)) 
 sigma.initial = expression(sigma <- rep(0.5, length(y))),
    nu.initial = expression(nu <- rep(0.3, length(y))), 
   tau.initial = expression(tau <-rep(0.3, length(y))), 
      mu.valid = function(mu) all(mu > 0 & mu < 1) , 
   sigma.valid = function(sigma)  all(sigma > 0 & sigma < 1), 
      nu.valid = function(nu)  all(nu > 0) , 
     tau.valid = function(tau) all(tau > 0), 
       y.valid = function(y)  all(y >= 0 & y <= 1)
          ),
            class = c("gamlss.family","family"))
}
#------------------------------------------------------------------------------------------
dBEINF<-function(x, mu = 0.5, sigma = 0.1, 
                    nu = 0.1, tau = 0.1, log = FALSE)
 { 
          if (any(mu <= 0) | any(mu >= 1) )  
              stop(paste("mu must be between 0 and 1", "\n", "")) 
          if (any(sigma <= 0) | any(sigma >= 1))  
              stop(paste("sigma must be between 0 and 1", "\n", "")) 
          if (any(nu <= 0) )  
             stop(paste("nu must greated than 0", "\n", ""))           
          if (any(tau <= 0) )  
             stop(paste("tau must greated than 0", "\n", "")) 
          if (any(x < 0) | any(x > 1))  
             stop(paste("x must be 0<=x<=1, i.e. 0 to 1 inclusively", "\n", ""))  
              a <- mu*(1-sigma^2)/(sigma^2)
              b <- a*(1-mu)/mu
          logfy <- rep(0, length(x))
          logfy <- ifelse((x>0 & x<1), dbeta(x, shape1=a, shape2=b, ncp=0, log=TRUE), 0)
          logfy <- ifelse((x==0), log(nu), logfy)          
          logfy <- ifelse((x==1), log(tau) , logfy)
          logfy <- logfy - log(1+nu+tau)          
          if(log==FALSE) fy <- exp(logfy) else fy <- logfy
          fy
  }
#------------------------------------------------------------------------------------------
pBEINF <- function(q, mu = 0.5, sigma = 0.1, nu = 0.1, tau = 0.1, 
                      lower.tail = TRUE, log.p = FALSE)
  {     
         if (any(mu <= 0) | any(mu >= 1) )  
             stop(paste("mu must be between 0 and 1", "\n", "")) 
         if (any(sigma <= 0) | any(sigma >= 1))  
             stop(paste("sigma must be between 0 and 1", "\n", "")) 
         if (any(nu <= 0) )  
             stop(paste("nu must greated than 0", "\n", ""))           
         if (any(tau <= 0) )  
             stop(paste("tau must greated than 0", "\n", "")) 
         if (any(q < 0) | any(q > 1))  
             stop(paste("y must be 0<=y<=1, i.e. 0 to 1 inclusively", "\n", ""))  
            a <- mu*(1-sigma^2)/(sigma^2)
            b <- a*(1-mu)/mu
          cdf <- ifelse((q>0 & q<1), nu + pbeta(q, shape1=a, shape2=b, ncp=0, 
                                               lower.tail=TRUE,log.p=FALSE), 0)
          cdf <- ifelse((q==0), nu, cdf)          
          cdf <- ifelse((q==1), 1+nu+tau , cdf)
          cdf <- cdf/(1+nu+tau)          
          if(lower.tail==TRUE) cdf <- cdf else cdf=1-cdf
          if(log.p==FALSE) cdf <- cdf else cdf <- log(cdf)    
          cdf
   }
#------------------------------------------------------------------------------------------
qBEINF <- function(p, mu = 0.5, sigma = 0.1, nu = 0.1, tau = 0.1, 
                    lower.tail = TRUE, log.p = FALSE)
  {      if (any(mu <= 0) | any(mu >= 1) )  
            stop(paste("mu must be between 0 and 1", "\n", "")) 
         if (any(sigma <= 0) | any(sigma >= 1))  
            stop(paste("sigma must be between 0 and 1", "\n", ""))   
         if (any(nu <= 0) )  
            stop(paste("nu must greated than 0", "\n", ""))           
         if (any(tau <= 0) )  
            stop(paste("tau must greated than 0", "\n", "")) 
         if (any(p < 0) | any(p > 1))  
            stop(paste("p must be between 0 and 1", "\n", ""))    
         if (log.p==TRUE) p <- exp(p) else p <- p
         if (lower.tail==TRUE) p <- p else p <- 1-p
          a <- mu*(1-sigma^2)/(sigma^2)
          b <- a*(1-mu)/mu
          suppressWarnings(
          q <- ifelse((p<=(nu/(1+nu+tau))),0, qbeta((p-(nu/(1+nu+tau)))/(1/(1+nu+tau)), shape1=a, 
                                            shape2=b, lower.tail=TRUE, log.p=FALSE)))
          q <- ifelse((p>=((1+nu)/(1+nu+tau))), 1, q)        
          q
   }
#------------------------------------------------------------------------------------------
rBEINF <- function(n, mu = 0.5, sigma = 0.1, nu = 0.1, tau = 0.1)
  { if (any(mu <= 0) | any(mu >= 1) )  
        stop(paste("mu must be between 0 and 1", "\n", "")) 
    if (any(sigma <= 0) | any(sigma >= 1))  
        stop(paste("sigma must be between 0 and 1", "\n", ""))   
    if (any(nu <= 0) )  
        stop(paste("nu must greated than 0", "\n", ""))           
    if (any(tau <= 0) )  
        stop(paste("tau must greated than 0", "\n", "")) 
    if (any(n <= 0))  
        stop(paste("n must be a positive integer", "\n", ""))  
          n <- ceiling(n)
          p <- runif(n)
          r <- qBEINF(p, mu=mu, sigma=sigma, nu=nu, tau=tau)
          r
  }
#------------------------------------------------------------------------------------------
plotBEINF <- function( mu =.5 , sigma=.5, nu = 0.5, tau = 0.5, from = 0.001, to=0.999, n = 101, ...)
 {
  fy <- dBEINF( seq(from,to,length=n), mu = mu ,sigma = sigma, nu = nu, tau=tau)
  maxfy <- max(fy)
  pr<-c(dBEINF(0, mu=mu ,sigma=sigma,  nu=nu, tau=tau),dBEINF(1, mu=mu ,sigma=sigma,  nu=nu, tau=tau))
  allmax<- max(c(pr, maxfy))
  plot(function(y) dBEINF(y, mu = mu ,sigma = sigma, nu = nu, tau=tau), from = from, to = to, n = n, ylim=c(0, allmax), ... )  
  po<-c(0,1)
  points(po,pr,type="h")
  points(po,pr,type="p", col="blue")
 }
#-----------------------------------------------------------------------------------------
meanBEINF <- function(obj)
  {
  if ( obj$family[1]!="BEINF") stop("the object do not have a BEINF distribution")
  meanofY<-(fitted(obj,"tau")+fitted(obj,"mu"))/(1+fitted(obj,"nu")+fitted(obj,"tau"))
  meanofY
  }
#---------------------------------------------------------------------------------------- 
