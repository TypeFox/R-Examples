# MS+KA Wednesday, April 3, 2002 at 09:27
# BR December 2004
# last change  BR December 2004
# last modification Tuesday, March 28, 2006 at 11:14 DS
PO <- function (mu.link = "log") 
{
    mstats <- checklink("mu.link", "Poisson", substitute(mu.link),c("inverse", "log", "sqrt", "identity")) 
    structure(
          list(family = c("PO", "Poisson"),
           parameters = list(mu = TRUE), # the mean
                nopar = 1, 
                 type = "Discrete", 
              mu.link = as.character(substitute(mu.link)), 
           mu.linkfun = mstats$linkfun, 
           mu.linkinv = mstats$linkinv, 
                mu.dr = mstats$mu.eta, 
                 dldm = function(y,mu) (y-mu)/mu,
               d2ldm2 = function(mu) -1/mu,
          G.dev.incr  = function(y,mu,...) -2*dPO(x = y, mu = mu, log = TRUE),
                rqres = expression(rqres(pfun="pPO", type="Discrete", ymin=0, y=y, mu=mu)), 
            mu.initial =expression({mu <- (y +mean(y))/2 } ),
              mu.valid = function(mu) all(mu > 0), 
               y.valid = function(y)  all(y >= 0) 
          ),
            class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------------------------------
dPO<-function(x, mu = 1, log = FALSE)
 { 
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(x < 0) )  stop(paste("x must be >=0", "\n", ""))  
          fy <- dpois(x = x, lambda = mu, log = log)
          fy
  }
  
#----------------------------------------------------------------------------------------
pPO <- function(q, mu = 1, lower.tail = TRUE, log.p = FALSE)
  {     
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(q < 0) )  stop(paste("y must be >=0", "\n", ""))  
          cdf <- ppois(q, lambda = mu, lower.tail = lower.tail, log.p = log.p)
          cdf
   }
#----------------------------------------------------------------------------------------
qPO <- function(p, mu = 1, lower.tail = TRUE, log.p = FALSE)
  {      
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(p < 0) | any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
          q <- qpois(p, lambda = mu, lower.tail = lower.tail, log.p = log.p)
          q
   }
#----------------------------------------------------------------------------------------
rPO <- function(n, mu = 1)
  { 
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
          r <- rpois(n, lambda = mu)
          r
  }
#----------------------------------------------------------------------------------------
