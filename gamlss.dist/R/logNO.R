# MS + BR Thursday, April 13, 2006 
LOGNO <- function (mu.link="identity", sigma.link="log") 
{
    mstats <- checklink("mu.link", "Log Normal", substitute(mu.link), c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Log Normal", substitute(sigma.link), c("inverse", "log", "identity", "own"))
    
    structure(
          list(family = c("LOGNO", "Log Normal"),
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
                  dldm = function(y,mu,sigma) {
                                        dldm <- (log(y)-mu)/sigma^2
                                        dldm
                                    },
               d2ldm2 = function(sigma) -1/sigma^2,
                 dldd = function(y,mu,sigma)  {
                                        dldd <- (1/(sigma^3))*((log(y)-mu)^2-sigma^2)
                                        dldd
                                    },
               d2ldd2 = function(sigma) -2/sigma^2,
              d2ldmdd = function(y)  rep(0,length(y)),
          G.dev.incr  = function(y,mu,sigma,...) -2*dLOGNO(x=y, mu=mu,  sigma = sigma, log = TRUE), #
                            rqres = expression(  rqres(pfun="pLOGNO", type="Continuous", y=y, mu=mu, sigma=sigma)),
           mu.initial = expression({    mu <- (log(y)+mean(log(y)))/2  }),
         sigma.initial = expression({sigma <- rep(sd(log(y)),length(y)) }),  
              mu.valid = function(mu) all(mu > 0), 
           sigma.valid = function(sigma)  all(sigma > 0),
               y.valid = function(y)  all(y > 0)
          ),
            class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------

dLOGNO <- function(x, mu=0, sigma=1,  log = FALSE)
 {
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
          if (any(x < 0) )  stop(paste("x must be >=0", "\n", ""))  
          fy <- dlnorm(x=x, meanlog = mu, sdlog = sigma, log = log) 
          fy
  }    
#----------------------------------------------------------------------------------------
pLOGNO <- function(q, mu=0, sigma=1,  lower.tail = TRUE, log.p = FALSE)
 {  
    #      if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
          if (any(q < 0) )  stop(paste("y must be >=0", "\n", ""))  
          cdf <- plnorm(q=q, meanlog = mu, sdlog = sigma, lower.tail = lower.tail, log.p = log.p)
          cdf
 }
#----------------------------------------------------------------------------------------
qLOGNO <- function(p, mu=0, sigma=1,  lower.tail = TRUE, log.p = FALSE )
 { 
     #if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
          if (any(p < 0) | any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
          q <-qlnorm(p=p, meanlog = mu, sdlog = sigma, lower.tail = lower.tail, log.p = log.p)
          q
 }
#----------------------------------------------------------------------------------------
rLOGNO <- function(n, mu=0, sigma=1)
  {
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
          if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
          r <- rlnorm(n=n, meanlog = mu, sdlog = sigma)
  }
#----------------------------------------------------------------------------------------
