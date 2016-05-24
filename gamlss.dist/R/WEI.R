# last change MS Wednesday, September 17, 2003 at 08:51
WEI <- function (mu.link="log", sigma.link="log") 
{
    mstats <- checklink("mu.link", "Weibull", substitute(mu.link), c("inverse", "log", "identity", "own" ))# dummy
    dstats <- checklink("sigma.link", "Weibull", substitute(sigma.link), c("inverse", "log", "identity", "own"))
    structure(
          list(family = c("WEI", "Weibull"),
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
                 dldm = function(y,mu,sigma) ((y/mu)^sigma -1)*(sigma/mu),
               d2ldm2 = function(mu,sigma) - sigma^2/mu^2,
                 dldd = function(y,mu,sigma) 1/sigma - log(y/mu)*((y/mu)^sigma-1) ,
               d2ldd2 = function(sigma) -1.82368/sigma^2 ,
              d2ldmdd = function(mu) 0.422784/mu,
          G.dev.incr  = function(y,mu,sigma,...) -2*dWEI(y, mu ,sigma, log=TRUE), 
                rqres = expression(rqres(pfun="pWEI", type="Continuous", y=y, mu=mu, sigma=sigma)),
           mu.initial = expression( { log.Y.m <- log(y) 
                                      var.Y.v <- var(log(y))
                                      sd.Y.s <- 1.283/sqrt(var.Y.v)
                                      mu <- exp(log.Y.m + 0.5772/sd.Y.s)
                                      }),
        sigma.initial = expression({  var.logY <- var(log(y))
                                      s.Y.s <- 1.283/sqrt(var.logY) 
                                    sigma <-  rep(s.Y.s,length(y))}),
              mu.valid = function(mu) all(mu > 0) , 
           sigma.valid = function(sigma)  all(sigma > 0), 
               y.valid = function(y)  all(y > 0)
          ),
            class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------
dWEI<-function(x, mu=1, sigma=1, log=FALSE)
 { 
          if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(x < 0))  stop(paste("x must be positive", "\n", ""))  
      fy <- dweibull(x, scale=mu, shape=sigma, log =log)
      fy 
  }
#----------------------------------------------------------------------------------------
pWEI <- function(q, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE)
  {     
          if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(q < 0))  stop(paste("y must be positive", "\n", ""))  
    cdf <- pweibull(q, scale=mu, shape=sigma, lower.tail = lower.tail, log.p = log.p)
    cdf
   }
#----------------------------------------------------------------------------------------
qWEI <- function(p, mu=1, sigma=1,  lower.tail = TRUE, log.p = FALSE)
  { if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
    q <- qweibull(p, scale=mu, shape=sigma, lower.tail = lower.tail)
    q
   }
#----------------------------------------------------------------------------------------
rWEI <- function(n, mu=1, sigma=1)
  { if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
    r <-  rweibull(n, scale=mu, shape=sigma)
    r
  }
#----------------------------------------------------------------------------------------
