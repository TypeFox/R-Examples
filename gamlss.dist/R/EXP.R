# MS+PA Wednesday, April 3, 2002 at 09:00
EXP <- function (mu.link ="log") 
{
    mstats <- checklink("mu.link", "Exponential", substitute(mu.link), c("inverse", "log", "sqrt","identity")) 
    structure(
          list(family = c("EXP","Exponential"),
           parameters = list(mu=TRUE), 
                nopar = 1,
                type = "Continuous", 
              mu.link = as.character(substitute(mu.link)),  
           mu.linkfun = mstats$linkfun, 
           mu.linkinv = mstats$linkinv,
                mu.dr = mstats$mu.eta, 
                 dldm = function(y,mu) ((y-mu)/mu^2),
               d2ldm2 = function(mu) (-1/mu^2), 
          G.dev.incr  = function(y,mu,...)  -2*dEXP(x = y, mu = mu, log = TRUE) , 
                rqres = expression(rqres(pfun="pEXP", type="Continuous", y=y, mu=mu)), 
           mu.initial = expression(mu <- (y+mean(y))/2),
             mu.valid = function(mu) all(mu > 0) ,  
              y.valid = function(y) all(y > 0)
          ),
            class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------
dEXP<-function(x, mu = 1, log = FALSE)
 { 
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(x < 0) )  stop(paste("x must be >0", "\n", ""))  
          fy <- dexp(x = x, rate =1/mu, log = log)
          fy
  }
  
#----------------------------------------------------------------------------------------
pEXP <- function(q, mu = 1, lower.tail = TRUE, log.p = FALSE)
  {     
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(q < 0) )  stop(paste("y must be > 0", "\n", ""))  
          cdf <- pexp(q, rate =1/mu, lower.tail = lower.tail, log.p = log.p)
          cdf
   }
#----------------------------------------------------------------------------------------
qEXP <- function(p, mu = 1, lower.tail = TRUE, log.p = FALSE)
  {      
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(p < 0) | any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
          q <- qexp(p, rate = 1/mu, lower.tail = lower.tail, log.p = log.p)
          q
   }
#----------------------------------------------------------------------------------------
rEXP <- function(n, mu = 1)
  { 
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
          r <- rexp(n, rate =1/mu)
          r
  }
#----------------------------------------------------------------------------------------
