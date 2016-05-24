# MS+KA+RS Tuesday, October 22, 2002 at 11:29 
# last change Thursday, January 16, 2003 at 12:06 MS
# the sigma is the standard deviation
# pdf is out now 
NO <-function (mu.link ="identity", sigma.link="log") 
{
    mstats <- checklink(   "mu.link", "Normal", substitute(mu.link),    
                                              c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Normal", substitute(sigma.link), 
                                              c("inverse", "log", "identity", "own"))        
    structure(
          list(family = c("NO", "Normal"),
           parameters = list(mu=TRUE,sigma=TRUE),
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
                 dldm = function(y,mu,sigma) (1/sigma^2)*(y-mu),
               d2ldm2 = function(sigma) -(1/sigma^2),
                 dldd = function(y,mu,sigma)  ((y-mu)^2-sigma^2)/(sigma^3),
               d2ldd2 = function(sigma) -(2/(sigma^2)),
              d2ldmdd = function(y) rep(0,length(y)),
          G.dev.incr  = function(y,mu,sigma,...) -2*dNO(y,mu,sigma,log=TRUE),                         
                rqres = expression(rqres(pfun="pNO", type="Continuous", y=y, mu=mu, sigma=sigma)),
           mu.initial = expression({ mu <- (y+mean(y))/2 }),
        sigma.initial = expression({sigma <- rep(sd(y),length(y))}), 
             mu.valid = function(mu) TRUE , 
          sigma.valid = function(sigma) all(sigma > 0), 
              y.valid = function(y)  TRUE
          ),
            class = c("gamlss.family","family"))
}

#----------------------------------------------------------------------------------------
dNO<-function(x, mu=0, sigma=1, log=FALSE)
 { 
         # if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    fy <- dnorm(x, mean=mu, sd=sigma, log=log)
    fy
  }
#----------------------------------------------------------------------------------------
pNO <- function(q, mu=0, sigma=1, lower.tail = TRUE, log.p = FALSE)
  {     
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    cdf <- pnorm(q, mean=mu, sd=sigma, lower.tail = lower.tail, log.p = log.p)
    cdf
   }
#----------------------------------------------------------------------------------------
qNO <- function(p, mu=0, sigma=1, lower.tail = TRUE, log.p = FALSE)
  { if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", "")) 
    q <- qnorm(p, mean=mu, sd=sigma, lower.tail = lower.tail )
    q
   }

rNO <- function(n, mu=0, sigma=1)
  { 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    r <- rnorm(n, mean=mu, sd=sigma)
    r
  }
#----------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------
# MS Saturday, April 6, 2002 at 14:40
# Thursday, January 16, 2003 at 12:20
# the version of the Normal with sigma as the variance 
#----------------------------------------------------------------------------------------
 NO2 <-function (mu.link ="identity", sigma.link="log") 
{
    mstats <- checklink(   "mu.link", "Normal", substitute(mu.link),    
                                              c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Normal", substitute(sigma.link), 
                                              c("inverse", "log", "identity", "own"))        
    structure(
          list(family = c("NO2","Normal with variance"),
           parameters = list(mu=TRUE,sigma=TRUE), 
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
                 dldm = function(y,mu,sigma) (1/sigma)*(y-mu),
               d2ldm2 = function(sigma) -(1/sigma),
                 dldd = function(y,mu,sigma)  0.5*((y-mu)^2-sigma)/(sigma^2),
               d2ldd2 = function(sigma)  -(1/(2*sigma^2)),
              d2ldmdd = function(y) rep(0,length(y)),
          G.dev.incr  = function(y,mu,sigma,...) -2*dNO2(y,mu,sigma,log=TRUE),                                                          
                rqres = expression(rqres(pfun="pNO2", type="Continuous", y=y, mu=mu, sigma=sigma)),
           mu.initial = expression({  mu <- y+0.00001  }),
        sigma.initial = expression({sigma <- rep(var(y),length(y)) }), 
             mu.valid = function(mu) TRUE , 
          sigma.valid = function(sigma)  all(sigma > 0), 
              y.valid = function(y)  TRUE
          ),
            class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------
dNO2<-function(x, mu=0, sigma=1, log=FALSE)
 { 
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    fy <- dnorm(x, mean=mu, sd=sqrt(sigma), log=log)
    fy
  }
#----------------------------------------------------------------------------------------
pNO2 <- function(q, mu=0, sigma=1, lower.tail = TRUE, log.p = FALSE)
  {     
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    cdf <- pnorm(q, mean=mu, sd=sqrt(sigma), lower.tail = lower.tail, log.p = log.p)
    cdf
   }
#----------------------------------------------------------------------------------------
qNO2 <- function(p, mu=0, sigma=1, lower.tail = TRUE, log.p = FALSE)
  { if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", "")) 
    q <- qnorm(p, mean=mu, sd=sqrt(sigma), lower.tail = lower.tail )
    q
   }
#----------------------------------------------------------------------------------------
rNO2 <- function(n, mu=0, sigma=1)
  { 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    r <- rnorm(n, mean=mu, sd=sqrt(sigma))
    r
  }
#----------------------------------------------------------------------------------------
