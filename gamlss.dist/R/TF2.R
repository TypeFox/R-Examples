#----------------------------------------------------------------------------------------
# MS + BR on 10-07-2012 at 5pm
# this t family has mean mu and variance sigma^2
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
# for testing
#library(gamlss)
#source('/Volumes/Data/Users/stasinom/Documents/gamlss/projects/DISTRIBUTIONS/testingDistributions/testContDist.R')
#testContDist("TF2", y.range=c(-Inf,Inf), mu.range = c(-Inf,Inf), sigma.range=c(0, Inf),  nu.range = c(0, Inf), mu.val = c(0), sigma.val=c(1,2,10), nu=c(2.5, 4, 10))
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
TF2 <- function(mu.link="identity", sigma.link="log", nu.link ="logshiftto2")
{
    mstats <- checklink("mu.link", "t Family", substitute(mu.link), c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "t Family", substitute(sigma.link), c("inverse", "log", "identity", "own"))
    vstats <- checklink("nu.link", "t Family", substitute(nu.link), c("logshiftto2","inverse", "log", "identity", "own"))
    structure(
          list(family = c("TF2", "t Family 2"),
           parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE), 
                nopar = 3, 
                 type = "Continuous",
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
                 dldm = function(y,mu,sigma,nu) {
                                  sigma1 <- (sqrt((nu-2)/nu))*sigma                
                                    dldm <- TF()$dldm(y, mu, sigma1, nu)
                                    dldm
                                    },
               d2ldm2 = function(sigma,nu) {
                                    sigma1 <- (sqrt((nu-2)/nu))*sigma                
                                    d2ldm2 <- -(nu+1)/((nu+3)*(sigma1^2))
                                    d2ldm2
                                    },    
                 dldd = function(y,mu,sigma,nu) {
                                  sigma1 <- (sqrt((nu-2)/nu))*sigma   
                                   ds1dd <-  (sqrt((nu-2)/nu))        
                                    dldd <- ds1dd*TF()$dldd(y, mu, sigma1, nu)
                                    dldd
                                    },
               d2ldd2 = function(sigma,nu) {
                                 sigma1 <- (sqrt((nu-2)/nu))*sigma   
                                   ds1dd <-  (sqrt((nu-2)/nu)) 
                                       s2 <- sigma1^2
                                  d2ldd2 <- -(ds1dd^2)*((2*nu)/((nu+3)*s2))
                                  d2ldd2 
                                    },
                 dldv = function(y,mu,sigma,nu) {
                                  sigma1 <- (sqrt((nu-2)/nu))*sigma   
                                   ds1dv <-  (sqrt(nu/(nu-2)))*sigma/(nu^2)        
                                    dldv <- TF()$dldv(y, mu, sigma1, nu) +
                                            ds1dv*TF()$dldd(y, mu, sigma1, nu) 
                                    dldv
                                    },
               d2ldv2 = function(y,mu,sigma,nu)  {
                                  sigma1 <- (sqrt((nu-2)/nu))*sigma                 
                                   ds1dv <-  (sqrt(nu/(nu-2)))*sigma/(nu^2)        
                                      v2 <- nu/2
                                      v3 <-(nu+1)/2
                                  d2ldv2 <- trigamma(v3)-trigamma(v2)+(2*(nu+5))/(nu*(nu+1)*(nu+3))
                                  d2ldv2 <- d2ldv2/4
                                  d2ldv2 <- d2ldv2 + (ds1dv^2)*TF()$d2ldd2(sigma1, nu) 
                                  d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)
                                  d2ldv2  
                                    },
              d2ldmdd = function(y)  rep(0,length(y)),
              d2ldmdv = function(y)  rep(0,length(y)),
              d2ldddv = function(sigma,nu) {
                                   sigma1 <- (sqrt((nu-2)/nu))*sigma                 
                                    ds1dd <-  (sqrt((nu-2)/nu)) 
                                    ds1dv <-  (sqrt(nu/(nu-2)))*sigma/(nu^2)        
                                  d2ldddv <-  ds1dd*2/(sigma1*(nu+3)*(nu+1)) 
                                  d2ldddv <- d2ldddv + ds1dd*ds1dv*TF()$d2ldd2(sigma1, nu)    
                                    d2ldddv
                                    }, 
          G.dev.incr  = function(y,mu,sigma,nu,...) -2*dTF2(y,mu,sigma,nu,log=TRUE),                           
                rqres = expression(
                rqres(pfun="pTF2", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu) 
                                   ),
           mu.initial =  expression(mu <- (y+mean(y))/2), 
         sigma.initial =  expression(sigma <- rep(sd(y),length(y))),
            nu.initial =  expression(   nu <- rep(4, length(y))), 
              mu.valid = function(mu) TRUE , 
           sigma.valid = function(sigma)  all(sigma > 0),
              nu.valid = function(nu) all(nu > 0), 
               y.valid = function(y)  TRUE
          ),
            class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
dTF2<-function(x, mu=0, sigma=1, nu=10, log=FALSE)
 {  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))  
    sigma1 <- (sqrt((nu-2)/nu))*sigma
    if (length(nu)>1) lik <- ifelse(nu>1000000, 
                       dNO(x,mu=mu,sigma=sigma1,log=FALSE), 
                       (1/sigma1)*dt((x-mu)/sigma1, df=nu, log =FALSE)) 
    else lik <- if (nu>1000000) dNO(x,mu=mu,sigma=sigma1,log=FALSE)
                else  (1/sigma1)*dt((x-mu)/sigma1, df=nu, log =FALSE)
    fy <- if(log==FALSE) lik else log(lik)
    fy 
  }
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------- 
pTF2 <- function(q, mu=0, sigma=1, nu=10, lower.tail = TRUE, log.p = FALSE)
  { if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))  
    sigma1 <- (sqrt((nu-2)/nu))*sigma
    if (length(nu)>1) cdf <- ifelse(nu>1000000, 
                                   pNO(q, mu=mu, sigma=sigma1, lower.tail=lower.tail, log.p=log.p),
                                    pt((q-mu)/sigma1, df=nu,   lower.tail=lower.tail, log.p=log.p))
    else cdf <- if (nu>1000000)    pNO(q, mu=mu, sigma=sigma1, lower.tail=lower.tail, log.p=log.p)
                   else pt((q-mu)/sigma1, df=nu, lower.tail = lower.tail, log.p = log.p)
    cdf
   }
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
qTF2 <- function(p, mu=0, sigma=1, nu=10, lower.tail = TRUE, log.p = FALSE)
  { if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))  
    if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
    sigma1 <- (sqrt((nu-2)/nu))*sigma
    if (length(nu)>1) q <- ifelse(nu>1000000, 
                                 qNO(p, mu=mu, sigma=sigma1, lower.tail=lower.tail, log.p=log.p),
                                  mu+sigma1*qt(p,df=nu, lower.tail = lower.tail))
    else q <- if (nu>1000000)    qNO(p, mu=mu, sigma=sigma1, lower.tail=lower.tail, log.p=log.p)
              else               mu+sigma1*qt(p,df=nu, lower.tail = lower.tail)
    q
   }

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
rTF2 <- function(n, mu=0, sigma=1, nu=10)
  { if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))  
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
    n <- ceiling(n)
    p <- runif(n)
    r <- qTF2(p, mu=mu, sigma=sigma, nu=nu)
    r
  }
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------