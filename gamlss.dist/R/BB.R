# MS , RAR, KA amended 21_04_2010
BB <- function (mu.link = "logit", sigma.link = "log") 
{
    mstats <- checklink("mu.link", "Beta Binomial", substitute(mu.link),
                       c("logit", "probit", "cloglog", "cauchit", "log", "own"))
    dstats <- checklink("sigma.link", "Beta Binomial", substitute(sigma.link), 
                        c("inverse", "log", "identity", "sqrt", "own"))   
    structure(
          list(family = c("BB", "Beta Binomial"),
           parameters = list(mu=TRUE,sigma=TRUE), 
                nopar = 2,  
                 type = "Discrete", 
              mu.link = as.character(substitute(mu.link)),
           sigma.link = as.character(substitute(sigma.link)), 
           mu.linkfun = mstats$linkfun, 
        sigma.linkfun = dstats$linkfun, 
           mu.linkinv = mstats$linkinv, 
        sigma.linkinv = dstats$linkinv,
                mu.dr = mstats$mu.eta, 
             sigma.dr = dstats$mu.eta, 
                 dldm = function(y,mu,sigma,bd) (1/sigma)*(digamma(y+(mu/sigma))
                                    -digamma(bd+((1-mu)/sigma)-y)
                                    -digamma(mu/sigma)
                                    +digamma((1-mu)/sigma)),
               d2ldm2 = function(y,mu,sigma,bd) (1/(sigma)^2)*(trigamma(y+(mu/sigma))
                                                  +trigamma(bd+((1-mu)/sigma)-y)
                                                  -trigamma(mu/sigma)
                                                  -trigamma((1-mu)/sigma)),
                 dldd = function(y,mu,sigma,bd) {k <- 1/sigma
                                 dldd <- -(k^2)*(digamma(k)+mu*digamma(y+mu*k)+(1-mu)*
                                         digamma(bd+(1-mu)*k-y)-mu*digamma(mu*k)-(1-mu)*
                                         digamma((1-mu)*k)-digamma(bd+k))
                                 dldd 
                                   }, 
               d2ldd2 = function(y,mu,sigma,bd) {k <- 1/sigma
                                dldd <- -(k^2)*(digamma(k)+mu*digamma(y+mu*k)+(1-mu)*
                                         digamma(bd+(1-mu)*k-y)-mu*digamma(mu*k)-(1-mu)*
                                         digamma((1-mu)*k)-digamma(bd+k))
                              d2ldd2 <- -dldd^2
                              d2ldd2 
                                   }, 
              d2ldmdd = function(y) rep(0,length(y)),
           G.dev.incr = function(y,mu,sigma,bd,...) -2*dBB(y,mu,sigma,bd, log = TRUE), 
                rqres = expression(
                 rqres(pfun="pBB", type="Discrete", ymin=0, y=y, mu=mu, sigma=sigma, bd=bd)
                                   ), 
            mu.initial = expression (   mu <- (y+0.5)/(bd+1)),      
         sigma.initial = expression (sigma <- rep(1,length(y))),
              mu.valid = function(mu) all(mu > 0) && all(mu < 1), 
           sigma.valid = function(sigma)  all(sigma > 0), 
               y.valid = function(y)  all(y >= 0)
          ),
            class = c("gamlss.family","family"))
}
#------------------------------------------------------------------------------------------
dBB <- function(x, mu = 0.5, sigma = 1, bd = 10, log = FALSE)
 { 
    if (any(mu < 0) | any(mu > 1))   stop(paste("mu must be between 0 and 1 ", "\n", "")) 
    if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", ""))
    if (any(sigma < 1e-10)) warning(" values of sigma in BB less that 1e-10 are set to 1e-10" )
    sigma <- ifelse((sigma < 1e-10),1e-10,sigma)
    if (any(x < 0) )  stop(paste("x must be >=0", "\n", ""))  
    if (any(bd < x))  stop(paste("x  must be <=  than the binomial denominator", bd, "\n"))
      logfy <-   (lgamma(bd+1)-lgamma(x+1)-lgamma(bd-x+1)
                  +lgamma((1/sigma))+lgamma(x+mu*(1/sigma))
                  +lgamma(bd+((1-mu)/sigma)-x)-lgamma(mu*(1/sigma))
                  -lgamma((1-mu)/sigma)-lgamma(bd+(1/sigma)))
        if (length(sigma)>1) logfy2 <- ifelse(sigma>0.0001, logfy, 
                                          dBI(x, mu = mu, bd=bd, log = TRUE) )
        else logfy2 <- if (sigma<0.0001)  dBI(x, mu = mu, bd=bd, log = TRUE) 
                   else logfy
        fy <- if(log == FALSE) exp(logfy2) else logfy2
        fy
  }
#------------------------------------------------------------------------------------------
pBB <- function(q, mu = 0.5, sigma = 1, bd = 10, lower.tail = TRUE, log.p = FALSE)
  {     
    if (any(mu <= 0) | any(mu >= 1))   stop(paste("mu must be between 0 and 1 ", "\n", "")) 
    if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
    if (any(q < 0) )  stop(paste("y must be >=0", "\n", ""))
    if (any(bd < q))  stop(paste("y  must be <=  than the binomial denominator", bd, "\n"))    
        ly <- length(q)                                                       
       FFF <- rep(0,ly)                         
    nsigma <- rep(sigma, length = ly)
       nmu <- rep(mu, length = ly) 
       nbd <- rep(bd, length = ly)                                                         
         j <- seq(along=q) 
   for (i in j)                                                          
      {                                                                 
        y.y <- q[i]                                                   
         nn <- nbd[i]                                                  
         mm <- nmu[i]
       nsig <- nsigma[i]                                                     
     allval <- seq(0,y.y)
     pdfall <- dBB(allval, mu = mm, sigma = nsig, bd = nn, log = FALSE)
     FFF[i] <- sum(pdfall)                                             
      }  
      cdf <- FFF
      cdf <- if(lower.tail==TRUE) cdf else 1-cdf
      cdf <- if(log.p==FALSE) cdf else log(cdf)                                                                    
      if (length(sigma)>1) cdf2 <- ifelse(sigma>0.0001, cdf, 
                                          pBI(q, mu = mu, bd=bd, lower.tail=lower.tail, log.p = log.p) )
      else cdf2 <- if (sigma<0.0001) pBI(q, mu = mu, bd=bd, lower.tail=lower.tail, log.p = log.p)
                   else cdf
      cdf2
  }
#------------------------------------------------------------------------------------------
qBB <- function(p, mu=0.5, sigma=1, bd=10, lower.tail = TRUE, log.p = FALSE, fast = FALSE)
  {      
          if (any(mu <= 0) | any(mu >= 1))   stop(paste("mu must be between 0 and 1 ", "\n", "")) 
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
          if (any(p < 0) | any(p > 1.0001))  stop(paste("p must be between 0 and 1", "\n", "")) 
          if (log.p==TRUE) p <- exp(p) else p <- p
          if (lower.tail==TRUE) p <- p else p <- 1-p  
        ly <- length(p)                                                       
       QQQ <- rep(0,ly)                         
    nsigma <- rep(sigma, length = ly)
       nmu <- rep(mu, length = ly) 
       nbd <- rep(bd, length = ly)  
    for (i in seq(along=p))                                                          
      { 
    cumpro <- 0                                                                      
          for (j in seq(from = 0, to = nbd[i]))
            {
       cumpro <-  if (fast == FALSE) pBB(j, mu = nmu[i], sigma = nsigma[i] , bd = nbd[i] , log.p = FALSE)
                 else  cumpro+dBB(j, mu = nmu[i], sigma = nsigma[i] , bd = nbd[i] , log = FALSE)# the above is faster 
       QQQ[i] <- j 
       if  (p[i] <= cumpro ) break 
            } 
      }           
          invcdf <- QQQ
         if (length(sigma)>1) invcdf2 <- ifelse(sigma>0.0001, invcdf, 
                                          qBI(p, mu = mu, bd=bd, lower.tail=TRUE))
        else invcdf2 <- if (sigma<0.0001) qBI(p, mu = mu, bd=bd, lower.tail=TRUE)
                   else invcdf
     invcdf2    
   }
#------------------------------------------------------------------------------------------
rBB <- function(n, mu = 0.5, sigma = 1, bd = 10, fast = FALSE )
  { 
          if (any(mu <= 0) | any(mu >= 1))   stop(paste("mu must be between 0 and 1 ", "\n", "")) 
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
          if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
          n <- ceiling(n)
          p <- runif(n)
          r <- qBB(p, mu=mu, sigma=sigma, bd=bd, fast=fast)
          r
  }
#-------------------------------------------------------------------------------------------
