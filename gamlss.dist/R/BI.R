# RAR+KA+MS Friday, April 5, 2002 at 15:10
# last change Tuesday, December 14, 2004 Saturday, October 8, 2005 
BI <- function (mu.link = "logit") 
{
    mstats <- checklink("mu.link", "Binomial", substitute(mu.link),
                         c("logit", "probit", "cloglog", "cauchit", "log", "own"))# ms 8-10-05
    structure(
          list(family = c("BI", "Binomial"),
           parameters = list(mu=TRUE), 
                nopar = 1,
                 type = "Discrete",
              mu.link = as.character(substitute(mu.link)),  
           mu.linkfun = mstats$linkfun, 
           mu.linkinv = mstats$linkinv, 
                mu.dr = mstats$mu.eta,
                 dldm = function(y, mu, bd) (y-bd*mu)/(mu*(1-mu)),
               d2ldm2 = function(mu,bd) -(bd/(mu*(1-mu))),
          G.dev.incr  = function(y,mu,bd,...)  -2*dBI(y,bd,mu,log=TRUE),
                rqres = expression(
                  rqres(pfun="pBI", type="Discrete", ymin=0, y=y, mu=mu, bd=bd)
                                   ), #
            mu.initial = expression({mu <- (y + 0.5)/(bd + 1)}),
              mu.valid = function(mu) all(mu > 0) && all(mu < 1),  
               y.valid = function(y)  all(y >= 0)
          ),
            class = c("gamlss.family","family") )
}
#------------------------------------------------------------------------------------------
dBI<-function(x, bd = 1, mu = 0.5, log = FALSE)
 { 
    if (any(mu < 0) | any(mu > 1))  stop(paste("mu must be between 0 and 1", "\n", "")) 
    if (any(x < 0) )  stop(paste("x must be >=0", "\n", ""))
    if (any(bd < x))  stop(paste("x  must be <=  than the binomial denominator", bd, "\n"))   
    fy <- dbinom(x, size = bd, prob = mu, log = log)
    fy
  }
#------------------------------------------------------------------------------------------
pBI <- function(q, bd=1, mu=0.5, lower.tail = TRUE, log.p = FALSE)
  {     
    if (any(mu < 0) | any(mu > 1))  stop(paste("mu must be between 0 and 1", "\n", "")) 
    if (any(q < 0) )  stop(paste("q must be >=0", "\n", ""))
    if (any(bd < q))  stop(paste("q  must be <=  than the binomial denominator", bd, "\n"))   
    cdf <- pbinom(q, size = bd, prob = mu, lower.tail=lower.tail, log.p=log.p)
    cdf
   }
#------------------------------------------------------------------------------------------
qBI <- function(p, bd = 1, mu = 0.5,  lower.tail = TRUE, log.p = FALSE)
  {      
    if (any(mu < 0) | any(mu > 1))  stop(paste("mu must be between 0 and 1", "\n", "")) 
    if (any( p < 0) | any(p >  1))  stop(paste(" p must be between 0 and 1", "\n", ""))    
     q <- qbinom(p, size = bd, prob = mu, lower.tail = lower.tail, log.p = log.p)
     q
   }
#------------------------------------------------------------------------------------------
rBI <- function(n, bd = 1, mu = 0.5)
  { 
     if (any(mu < 0) | any(mu > 1))  stop(paste("mu must be between 0 and 1", "\n", "")) 
     if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
      n <- ceiling(n)
      p <- runif(n)
      r <- qbinom(p, size = bd, prob = mu)
      r
  }
#-------------------------------------------------------------------------------------------
