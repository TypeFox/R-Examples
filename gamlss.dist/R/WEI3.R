# RR
# the weibull distribution parameterised with mu as the mean
#----------------------------------------------------------------------------------------
WEI3 <- function (mu.link="log", sigma.link="log") 
{
    mstats <- checklink("mu.link", "WEI3bull", substitute(mu.link), c("inverse", "log", "identity"))# dummy
    dstats <- checklink("sigma.link", "WEI3bull", substitute(sigma.link), c("inverse", "log", "identity"))
    structure(
          list(family = c("WEI3", "Weibull type 3"),
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
                 dldm = function(y,mu,sigma) ((y*gamma((1/sigma)+1)/mu)^sigma -1)*(sigma/mu),
               d2ldm2 = function(sigma,mu) - sigma^2/mu^2,
                 dldd = function(y,mu,sigma) 1/sigma - log(y*gamma((1/sigma)+1)/mu)*
                         ((y*gamma((1/sigma)+1)/mu)^sigma-1)+
                        (digamma((1/sigma)+1))*((y*gamma((1/sigma)+1)/mu)^sigma-1)/sigma,
               d2ldd2 = function(sigma) -(1.644934+(0.422784-digamma((1/sigma)+1))^2)/(sigma*sigma) ,
              d2ldmdd = function(mu,sigma) (0.422784-digamma((1/sigma)+1))/mu,
          G.dev.incr  = function(y,mu,sigma,...) -2*dWEI3(y, mu ,sigma, log=TRUE), 
                rqres = expression(rqres(pfun="pWEI3", type="Continuous", y=y, mu=mu, sigma=sigma)),
           mu.initial = expression( {  mu <- y+0.01
                                      }),
        sigma.initial = expression({  var.log.Y <- var(log(y))
                                      s.Y.s <- 1.283/sqrt(var.log.Y) 
                                    sigma <-  rep(s.Y.s,length(y))}),
              mu.valid = function(mu) all(mu > 0) , 
           sigma.valid = function(sigma)  all(sigma > 0), 
               y.valid = function(y)  all(y > 0)
          ),
            class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------
dWEI3<-function(x, mu=1, sigma=1, log=FALSE)
 { 
          if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(x < 0))  stop(paste("x must be positive", "\n", ""))  
      mu2 <- mu/gamma((1/sigma)+1)
      fy <- dweibull(x, scale=mu2, shape=sigma, log =log)
      fy 
  }
#----------------------------------------------------------------------------------------  
pWEI3 <- function(q, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE)
  {     
          if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(q < 0))  stop(paste("y must be positive", "\n", ""))  
    mu2 <- mu/gamma((1/sigma)+1)    
    cdf <- pweibull(q, scale=mu2, shape=sigma, lower.tail = lower.tail, log.p = log.p)
    cdf
   }
#----------------------------------------------------------------------------------------
qWEI3 <- function(p, mu=1, sigma=1,  lower.tail = TRUE, log.p = FALSE)
  { if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
    mu2 <- mu/gamma((1/sigma)+1)
    q <- qweibull(p, scale=mu2, shape=sigma, lower.tail = lower.tail)
    q
   }
#----------------------------------------------------------------------------------------
rWEI3 <- function(n, mu=1, sigma=1)
  { if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
    mu2 <- mu/gamma((1/sigma)+1)
    r <-  rweibull(n, scale=mu2, shape=sigma)
    r
  }
#----------------------------------------------------------------------------------------
