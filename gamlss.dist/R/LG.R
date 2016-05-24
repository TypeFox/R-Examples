# I think this is working correctly 01/03/10
LG <- function (mu.link = "logit") 
{
    mstats <- checklink("mu.link", "LG", substitute(mu.link),c("logit", "probit", "cloglog", "cauchit", "log", "own"))
    structure(
          list(family = c("LG", "Logarithmic"),
           parameters = list(mu = TRUE), # the mean
                nopar = 1, 
                 type = "Discrete", 
              mu.link = as.character(substitute(mu.link)), 
           mu.linkfun = mstats$linkfun, 
           mu.linkinv = mstats$linkinv, 
                mu.dr = mstats$mu.eta, 
                 dldm = function(y,mu) (y/mu)+1/((1-mu)*log(1-mu)),
               d2ldm2 = function(y,mu) 
                        {
                         dldm <- (y/mu)+1/((1-mu)*log(1-mu))
                       d2ldm2 <- -dldm^2
                       d2ldm2
                        },
          G.dev.incr  = function(y,mu,...) -2*dLG(x = y, mu = mu, log = TRUE),
                rqres = expression(rqres(pfun="pLG", type="Discrete", ymin=1, y=y, mu=mu)), 
            mu.initial =expression({mu <- 0.9 } ),
              mu.valid = function(mu) all(mu > 0  & mu < 1), 
               y.valid = function(y)  all(y > 0) 
          ),
            class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------------------------------
dLG<-function(x, mu = 0.5, log = FALSE)
 { 
          if (any(mu <= 0) | any(mu >= 1) )  stop(paste("mu must be greater than 0 and less than 1", "\n", ""))
          if (any(x <= 0) )  stop(paste("x must be >0", "\n", ""))
       logfy <- x*log(mu)-log(x)-log(-log(1-mu))
      if(log == FALSE) fy <- exp(logfy) else fy <- logfy
          fy
  }
#----------------------------------------------------------------------------------------
pLG <- function(q, mu = 0.5, lower.tail = TRUE, log.p = FALSE)
  {     
          if (any(mu <= 0) | any(mu >= 1) )  stop(paste("mu must be greater than 0 and less than 1", "\n", ""))
   if (any(q <= 0) )  stop(paste("q must be >0", "\n", ""))    
        ly <- length(q)                                                       
       FFF <- rep(0,ly)                         
       nmu <- rep(mu, length = ly)                                                       
         j <- seq(along=q) 
   for (i in j)                                                          
      {                                                                 
        y.y <- q[i]                                                                                      
         mm <- nmu[i]                                      
     allval <- seq(1,y.y)
     pdfall <- dLG(allval, mu = mm, log = FALSE)
     FFF[i] <- sum(pdfall)                                             
      }  
      cdf <- FFF
      cdf <- if(lower.tail==TRUE) cdf else 1-cdf
      cdf <- if(log.p==FALSE) cdf else log(cdf)                                                                    
      cdf
  }
#----------------------------------------------------------------------------------------
qLG <- function(p, mu=0.5,  lower.tail = TRUE, log.p = FALSE,  
                 max.value = 10000)
  {      
          if (any(mu <= 0) | any(mu >= 1) )  stop(paste("mu must be greater than 0 and less than 1", "\n", "")) 
          if (any(p < 0) | any(p > 1.0001))  stop(paste("p must be between 0 and 1", "\n", "")) 
          if (log.p==TRUE) p <- exp(p) else p <- p
          if (lower.tail==TRUE) p <- p else p <- 1-p    
           ly <- length(p)                                                       
          QQQ <- rep(0,ly)                         
          nmu <- rep(mu, length = ly)                   
       for (i in seq(along=p))                                                          
      {
       cumpro <- 0                                                                         
     if (p[i]+0.000000001 >= 1) QQQ[i] <- Inf
     else  
        {  
            for (j in seq(from = 1, to = max.value))
            {
            cumpro <-  pLG(j, mu = nmu[i], log.p = FALSE) 
           QQQ[i] <- j 
       if  (p[i] <= cumpro ) break 
            } 
        }
      }          
          QQQ   
   } 
#----------------------------------------------------------------------------------------
rLG <- function(n, mu = 0.5)
  { 
          if (any(mu <= 0) | any(mu >= 1) )  stop(paste("mu must be greater than 0 and less than 1", "\n", "")) 
          if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
          n <- ceiling(n)
          p <- runif(n)
          r <- qLG(p, mu=mu)
          r
  }
#----------------------------------------------------------------------------------------
