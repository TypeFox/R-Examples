# 3/12/2004 RAR DS
# this is BEINF, the Beta distribution with probabilities for 0 and 1 (i.e. 4 parameters)
#------------------------------------------------------------------------------------------
BEINF1 <- function (mu.link = "logit", sigma.link = "logit", 
                   nu.link = "log")
{
    mstats <- checklink("mu.link", "BEINF1", substitute(mu.link),
                            c("logit", "probit", "cloglog", "cauchit", "log", "own"))
    dstats <- checklink("sigma.link", "BEINF1", substitute(sigma.link), 
                           c("logit", "probit", "cloglog", "cauchit", "log", "own"))   
    vstats <- checklink("nu.link", "BEINF1", substitute(nu.link),    
                           c("inverse", "log", "identity", "own"))
    structure(
          list(family = c("BEINF1", "Beta Inflated one"),
           parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE), 
                nopar = 3, 
                 type = "Mixed",
              mu.link = as.character(substitute(mu.link)),  
           sigma.link = as.character(substitute(sigma.link)), 
              nu.link = as.character(substitute(nu.link)),  
           mu.linkfun = mstats$linkfun, 
        sigma.linkfun = dstats$linkfun, 
           nu.linkfun = vstats$linkfun,
           mu.linkinv = mstats$linkinv, 
        sigma.linkinv = dstats$linkinv,
           nu.linkinv = vstats$linkinv, 
                mu.dr = mstats$mu.eta, 
             sigma.dr = dstats$mu.eta, 
                nu.dr = vstats$mu.eta,
    dldm = function(y,mu,sigma) { a <- mu*(1-sigma^2)/(sigma^2)
                        b <- a*(1-mu)/mu
                   dldm <- ifelse((y==1),0,((1-sigma^2)/(sigma^2))*(-digamma(a)
                                    +digamma(b) +log(y) - log(1-y)))
                       dldm 
                      },
    d2ldm2 = function(y,mu,sigma) { a <- mu*(1-sigma^2)/(sigma^2)
                          b <- a*(1-mu)/mu
                     d2ldm2 <- ifelse( (y==1), 0 , 
                           -(((1-sigma^2)^2)/(sigma^4))*(trigamma(a) +trigamma(b)))
                      d2ldm2 
                        },
    dldd = function(y,mu,sigma) {  a <- mu*(1-sigma^2)/(sigma^2)
                         b <- a*(1-mu)/mu
                      dldd <- ifelse( (y==1), 0 , 
                            -(2/(sigma^3))*( mu*(-digamma(a)+digamma(a+b)+log(y))
                                  +(1-mu)*(-digamma(b)+digamma(a+b)+log(1-y)) ))
                      dldd 
                      }, 
    d2ldd2 = function(y,mu,sigma,nu) { 
                         a <- mu*(1-sigma^2)/(sigma^2)
                         b <- a*(1-mu)/mu
                              d2ldd2 <- ifelse( (y==1), 0 , 
                            -(4/(sigma^6))*((mu^2)*trigamma(a) +((1-mu)^2)*trigamma(b)
                                    -trigamma(a+b)))
                   d2ldd2
                       }, 
     dldv = function(y,nu)  {
                         dldv <- ifelse(y==1,(1/nu),0) -(1/(1+nu))
                         dldv
                        }, 
    d2ldv2 = function(nu) {d2ldv2 <- -(1)/(nu*((1+nu)^2))
                         d2ldv2
                         },                   
    d2ldmdd = function(y,mu,sigma) { a <-  mu*(1-sigma^2)/(sigma^2)
                           b <- a*(1-mu)/mu
                     d2ldmdd <- ifelse( (y==1), 0 , 
                       (2*(1-sigma^2)/(sigma^5))*(mu*trigamma(a)-(1-mu)*trigamma(b)))
                     d2ldmdd 
                         },
  d2ldmdv = function(y) {
                        d2ldmdv <- rep(0,length(y))
                        d2ldmdv
                       },

  d2ldddv = function(y) {
                        d2ldddv <- rep(0,length(y))
                        d2ldddv
                       },

 G.dev.incr  = function(y,mu,sigma,nu,...) 
                        -2*dBEINF1(y,mu,sigma,nu,log=TRUE),                     
      rqres = expression(rqres(pfun="pBEINF1", type="Mixed",  mass.p=c(1),  
                             prob.mp=cbind(nu/(1+nu)), y=y, mu=mu, 
                             sigma=sigma, nu=nu)),
    mu.initial = expression(mu <- (y+mean(y))/2),    #(y+mean(y))/2),# rep(mean(y),length(y)) 
 sigma.initial = expression(sigma <- rep(0.5, length(y))),
    nu.initial = expression(nu <- rep(0.1, length(y))),  
      mu.valid = function(mu) all(mu > 0 & mu < 1) , 
   sigma.valid = function(sigma)  all(sigma > 0 & sigma < 1), 
      nu.valid = function(nu)  all(nu > 0) , 
       y.valid = function(y)  all(y > 0 & y <= 1)
          ),
            class = c("gamlss.family","family"))
}
#------------------------------------------------------------------------------------------
dBEINF1<-function(x, mu = 0.5, sigma = 0.1, nu = 0.1,  log = FALSE)
 { 
          if (any(mu <= 0) | any(mu >= 1) )  
              stop(paste("mu must be between 0 and 1", "\n", "")) 
          if (any(sigma <= 0) | any(sigma >= 1))  
              stop(paste("sigma must be between 0 and 1", "\n", "")) 
          if (any(nu <= 0) )  
             stop(paste("nu must greated than 0", "\n", ""))            
          if (any(x <= 0) | any(x > 1))
          {
             stop(paste("x must be 0< x<=1, i.e. (0 to 1] inclusively", "\n", ""))  
          } 
              a <- mu*(1-sigma^2)/(sigma^2)
              b <- a*(1-mu)/mu
          logfy <- rep(0, length(x))
          logfy <- ifelse((x>0), dbeta(x, shape1=a, shape2=b, ncp=0, log=TRUE), 0)
          logfy <- ifelse((x==1), log(nu), logfy)          
          logfy <- logfy - log(1+nu)          
          if(log==FALSE) fy <- exp(logfy) else fy <- logfy
          fy
  }
#------------------------------------------------------------------------------------------
pBEINF1 <- function(q, mu = 0.5, sigma = 0.1, nu = 0.1, 
                      lower.tail = TRUE, log.p = FALSE)
  {     
         if (any(mu <= 0) | any(mu >= 1) )  
             stop(paste("mu must be between 0 and 1", "\n", "")) 
         if (any(sigma <= 0) | any(sigma >= 1))  
             stop(paste("sigma must be between 0 and 1", "\n", "")) 
         if (any(nu <= 0) )  
             stop(paste("nu must greated than 0", "\n", ""))           
         if (any(q < 0) | any(q > 1))  
             stop(paste("y must be 0<=y<=1, i.e. 0 to 1 inclusively", "\n", ""))  
            a <- mu*(1-sigma^2)/(sigma^2)
            b <- a*(1-mu)/mu
           cdf <- ifelse((q>0 & q<1),  pbeta(q, shape1=a, shape2=b, ncp=0, 
                                               lower.tail=TRUE,log.p=FALSE), 0)          
           cdf <- ifelse((q==1), 1+nu , cdf)
           cdf <- cdf/(1+nu)               
          if(lower.tail==TRUE) cdf <- cdf else cdf=1-cdf
          if(log.p==FALSE) cdf <- cdf else cdf <- log(cdf)    
          cdf
   }
#------------------------------------------------------------------------------------------
qBEINF1 <- function(p, mu = 0.5, sigma = 0.1, nu = 0.1, 
                    lower.tail = TRUE, log.p = FALSE)
  {      if (any(mu <= 0) | any(mu >= 1) )  
            stop(paste("mu must be between 0 and 1", "\n", "")) 
         if (any(sigma <= 0) | any(sigma >= 1))  
            stop(paste("sigma must be between 0 and 1", "\n", ""))   
         if (any(nu <= 0) )  
            stop(paste("nu must greated than 0", "\n", ""))           
         if (any(p <= 0) | any(p >= 1))  
            stop(paste("p must be between (0 and 1]", "\n", ""))    
         if (log.p==TRUE) p <- exp(p) else p <- p
         if (lower.tail==TRUE) p <- p else p <- 1-p
          a <- mu*(1-sigma^2)/(sigma^2)
          b <- a*(1-mu)/mu 
          
         # suppressWarnings(
         # q <- ifelse((p<=(nu/(1+nu+tau))),0, qbeta((p-(nu/(1+nu+tau)))/(1/(1+nu+tau)), shape1=a, 
         #                                   shape2=b, lower.tail=TRUE, log.p=FALSE)))
         # q <- ifelse((p>=((1+nu)/(1+nu+tau))), 1, q)        
          
          
          suppressWarnings(
          q <- ifelse((p>=(1/(1+nu))),1, qbeta((p*(1+nu)), shape1=a, 
                                            shape2=b, lower.tail=TRUE, log.p=FALSE)))   
          q
   }
#------------------------------------------------------------------------------------------
rBEINF1 <- function(n, mu = 0.5, sigma = 0.1, nu = 0.1)
  { if (any(mu <= 0) | any(mu >= 1) )  
        stop(paste("mu must be between 0 and 1", "\n", "")) 
    if (any(sigma <= 0) | any(sigma >= 1))  
        stop(paste("sigma must be between 0 and 1", "\n", ""))   
    if (any(nu <= 0) )  
        stop(paste("nu must greated than 0", "\n", ""))           
    if (any(n <= 0))  
        stop(paste("n must be a positive integer", "\n", ""))  
          n <- ceiling(n)
          p <- runif(n)
          r <- qBEINF1(p, mu=mu, sigma=sigma, nu=nu)
          r
  }
#------------------------------------------------------------------------------------------
plotBEINF1 <- function( mu =.5 , sigma=.5, nu = 0.5, from = 0.0001, to=.9999, n = 101, ...)
 { 
  fy <- dBEINF1( seq(from,to,length=n), mu = mu ,sigma = sigma, nu = nu)
  maxfy <- max(fy)
  pr<-c(dBEINF1(1, mu=mu ,sigma=sigma,  nu=nu))
  allmax<- max(c(pr, maxfy))
  plot(function(y) dBEINF1(y, mu = mu ,sigma = sigma, nu = nu), from = from, to = to, n = n, ylim=c(0, allmax), ... )  
  po<-c(1)
  points(po,pr,type="h")
  points(po,pr,type="p", col="blue")
 }
#-----------------------------------------------------------------------------------------
meanBEINF1 <- function(obj)
  {
  if ( obj$family[1]!="BEINF1") stop("the object do not have a BEINF1 distribution")
  meanofY<-(fitted(obj,"nu")+fitted(obj,"mu"))/(1+fitted(obj,"nu"))
  meanofY
  }
#---------------------------------------------------------------------------------------- 
