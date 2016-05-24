#-------------------------------------------------------------
WEI2 <- function (mu.link ="log", sigma.link="log") 
{
    mstats <- checklink("mu.link", "Weibull.2", substitute(mu.link), c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Weibull.2", substitute(sigma.link), c("inverse", "log", "identity", "own"))
    structure(
          list(family = c("WEI2", "Weibull type 2"),
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
                 dldm = function(y,mu,sigma) (1/mu)-(y^sigma),
               d2ldm2 = function(mu)  -1/mu^2,
                 dldd = function(y,mu,sigma) (1/sigma) + log(y) -log(y)*mu*(y^sigma),
               d2ldd2 = function(mu,sigma) -(1/sigma^2)*(1.64493+(0.42278-log(mu))^2) ,
              d2ldmdd = function(mu,sigma) -(1/(sigma*mu))*(0.42278-log(mu)) ,
          G.dev.incr  = function(y,mu,sigma,...) -2*dWEI2(y,mu,sigma,log=TRUE) ,  
                rqres = expression(rqres(pfun="pWEI2", type="Continuous", y=y, mu=mu, sigma=sigma)),
           mu.initial = expression( { log.Y.m <- log(y) 
                                      var.log.Y <- var(log(y))
                                      s.Y.s <- 1.283/sqrt(var.log.Y)
                                      mu <- exp(-s.Y.s*(log.Y.m + 0.5772/s.Y.s))
                                      }),
        sigma.initial = expression({  var.log.Y <- var(log(y))
                                      s.Y.s <- 1.283/sqrt(var.log.Y) 
                                    sigma <-  rep(s.Y.s,length(y))}),
              mu.valid = function(mu) all(mu > 0), 
           sigma.valid = function(sigma)  all(sigma > 0), 
               y.valid = function(y)  all(y > 0)
          ),
            class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------
dWEI2<-function(x, mu=1, sigma=1, log=FALSE)
 { 
          muinv <- function(mubar,sigma) {  mu <- mubar^(-1/sigma); mu}
          if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(x < 0))  stop(paste("x must be positive", "\n", ""))  
      fy <- dweibull(x, scale=muinv(mu,sigma), shape=sigma, log=log)
      fy 
  }
#----------------------------------------------------------------------------------------  
pWEI2 <- function(q, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE)
  {     
         muinv <- function(mubar,sigma) {  mu <- mubar^(-1/sigma); mu}      
          if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(q < 0))  stop(paste("y must be positive", "\n", ""))  
    cdf <- pweibull(q, scale=muinv(mu,sigma), shape=sigma, lower.tail = lower.tail, log.p = log.p)
    cdf
   }
#----------------------------------------------------------------------------------------
qWEI2 <- function(p, mu=1, sigma=1,  lower.tail = TRUE, log.p = FALSE)
  {  muinv <- function(mubar,sigma) {  mu <- mubar^(-1/sigma); mu}   
    if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
    q <- qweibull(p, scale=muinv(mu,sigma), shape=sigma, lower.tail = lower.tail)
    q
   }
#----------------------------------------------------------------------------------------
rWEI2 <- function(n, mu=1, sigma=1)
  {  muinv <- function(mubar,sigma) {  mu <- mubar^(-1/sigma); mu}   
    if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
    r <-  rweibull(n, scale=muinv(mu,sigma), shape=sigma)
    r
  }
#----------------------------------------------------------------------------------------
