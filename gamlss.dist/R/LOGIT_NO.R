# Checked 1-7-12
#source("/Volumes/Data/Users/stasinom/Documents/gamlss/projects/DISTRIBUTIONS/testingDistributions/testContDist.R")
#testContDist("LOGITNO", y.range=c(0,1), mu.range = c(0,1),sigma.range=c(0, Inf), mu.val = c(0.01,.1,0.6,.99), sigma.val=c(1, 2, 5))

LOGITNO<-function (mu.link ="logit", sigma.link="log") 
{
    mstats <- checklink("mu.link", "logitNO", substitute(mu.link),
                         c("logit", "probit", "cloglog", "cauchit", "log", "own"))
    dstats <- checklink("sigma.link", "logitNO", substitute(sigma.link), 
                                              c("inverse", "log", "identity"))        
    structure(
          list(family = c("LOGITNO", "Logit Normal"),
           parameters = list(mu=TRUE, sigma=TRUE),
                nopar = 2, 
                 type = "Continuous",
              mu.link = as.character(substitute(mu.link)), 
           sigma.link = as.character(substitute(sigma.link)), 
           mu.linkfun = mstats$linkfun, 
        sigma.linkfun = dstats$linkfun, 
           mu.linkinv = mstats$linkinv, 
        sigma.linkinv = dstats$linkinv,
                mu.dr = mstats$mu.eta, 
             sigma.dr = dstats$mu.eta,
                 dldm = function(y,mu,sigma){
                 	      z <- log(y/(1-y))
                 	   mu_z <- log(mu/(1-mu))
                 	   dldm <- (1/(sigma^2)*(z-mu_z))*(1/(mu*(1-mu)))
                        #(log(y/(1-y))-log(mu/(1-mu)))/(sigma^2*mu*(1-mu))
                        dldm
                           },
               d2ldm2 = function(sigma,mu) {(-1/((sigma^2)*(mu^2)*((1-mu)^2)))},
                 dldd = function(y,mu,sigma)
                                {
                                 z <- log(y/(1-y))
                 	          mu_z <- log(mu/(1-mu))	
                                dldd <-  ((z-mu_z)^2-sigma^2)/(sigma^3)                                
                                #(-1/sigma)+((log(y/(1-y))-log(mu/(1-mu))))^2/(sigma^3)
                                dldd	
                                },
               d2ldd2 = function(sigma) -(2/(sigma^2)),
              d2ldmdd = function(y) rep(0,length(y)),
          G.dev.incr  = function(y,mu,sigma,...) -2*dLOGITNO(y,mu,sigma,log=TRUE),                         
                rqres = expression(rqres(pfun="pLOGITNO", type="Continuous", y=y, mu=mu, sigma=sigma)),
           mu.initial = expression({   mu <- rep(.5, length(y)) }),
        sigma.initial = expression({sigma <- rep(sd(log(y/(1-y))),length(y)) }), # 
             mu.valid = function(mu) all(mu > 0 & mu < 1) ,
          sigma.valid = function(sigma) all(sigma > 0), 
              y.valid = function(y)  all(y > 0 & y < 1)
          ),
            class = c("gamlss.family","family")
          )
}

#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
dLOGITNO<-function(x, mu=0.5, sigma=1, log=FALSE)
 { 
        if (any(mu < 0) | any(mu > 1))  stop(paste("mu must be between 0 and 1", "\n", "")) 
        if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
     z <- log(x/(1-x))
    loglik <- dnorm(z, mean=log(mu/(1-mu)), sd=sigma, log=TRUE)-log(x)-log(1-x)
    if(log==FALSE) fy  <- exp(loglik) else fy <- loglik 
    fy
  }
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
pLOGITNO<- function(q, mu=0.5, sigma=1, lower.tail = TRUE, log.p = FALSE)
  {     
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    qz <- log(q/(1-q))
    cdf <- pnorm(qz, mean=log(mu/(1-mu)), sd=sigma, lower.tail = lower.tail, log.p = log.p)
    cdf
   }
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
qLOGITNO<- function(p, mu=0.5, sigma=1, lower.tail = TRUE, log.p = FALSE)
  { if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", "")) 
    qz <- qnorm(p, mean=log(mu/(1-mu)), sd=sigma, lower.tail = lower.tail )
    qy <- 1/(1+exp(-qz))
    qy
   }
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
rLOGITNO<- function(n, mu=0.5, sigma=1)
  { 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    rz <- rnorm(n, mean=log(mu/(1-mu)), sd=sigma)
    ry <- 1/(1+exp(-rz))
    ry
  }
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------