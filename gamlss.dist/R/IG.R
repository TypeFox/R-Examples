# amended 27_11_2007
IG <-function (mu.link = "log", sigma.link = "log") 
{   
    mstats <- checklink("mu.link", "Inverse Gaussian", substitute(mu.link), c("1/mu^2", "inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Inverse Gaussian", substitute(sigma.link),  c("inverse", "log", "identity", "own"))
    structure(
          list(family = c("IG", "Inverse Gaussian"),
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
                 dldm = function(y, mu, sigma) (y-mu)/((sigma^2)*(mu^3)),
               d2ldm2 = function(mu,sigma) -1/((mu^3)*(sigma^2)),
                 dldd = function(y,mu,sigma) (-1/sigma) +((y-mu)^2)/(y*(sigma^3)*(mu^2)),
               d2ldd2 = function(sigma) -2/(sigma^2),
              d2ldmdd = function(y)  rep(0,length(y)),
          G.dev.incr  = function(y,mu,sigma,...) 
                               {  -2*dIG(y,mu,sigma,log=TRUE)},
                rqres = expression(rqres(pfun="pIG", type="Continuous", y=y, mu=mu, sigma=sigma)),
            mu.initial = expression( mu <- (y+mean(y))/2) ,
         sigma.initial = expression(sigma <- sd(y)/(mean(y))^1.5 ), 
              mu.valid = function(mu) all(mu > 0), 
           sigma.valid = function(sigma)  all(sigma > 0), 
               y.valid = function(y)  all(y > 0)
            ),
            class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------  
dIG<-function(x, mu = 1, sigma = 1, log=FALSE)
 {        if (any(mu < 0))  stop(paste("mu must be positive", "\n", "")) 
          if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(x < 0))  stop(paste("x must be positive", "\n", ""))  
 log.lik <- (-0.5*log(2*pi)-log(sigma)-(3/2)*log(x)-((x-mu)^2)/(2*sigma^2*(mu^2)*x) )
     if(log==FALSE) fy  <- exp(log.lik) else fy <- log.lik
      fy 
  }
#----------------------------------------------------------------------------------------  
pIG <- function(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
  {    #  browser() 
    if (any(mu < 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(q < 0))  stop(paste("y must be positive", "\n", ""))  
      lq <- length(q)                                                                    
   sigma <- rep(sigma, length = lq)
      mu <- rep(mu, length = lq)           

    cdf1 <- pnorm(((q/mu)-1)/(sigma*sqrt(q))) 
    lcdf2 <- (2/(mu*sigma^2))+pnorm((-((q/mu)+1))/(sigma*sqrt(q)),log.p=TRUE)
     cdf <- cdf1+ exp(lcdf2)
    if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf 
    if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf) 
    cdf
   }
#----------------------------------------------------------------------------------------  
qIG <- function(p, mu=1, sigma=1,  lower.tail = TRUE, log.p = FALSE)
 {
    #---functions--------------------------------------------   
       h1 <- function(q)
       { 
     pIG(q , mu = mu[i], sigma = sigma[i])-p[i]   
       }
       h <- function(q)
       { 
     pIG(q , mu = mu[i], sigma = sigma[i])   
       }
     #-------------------------------------------------------
    if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))      
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (lower.tail==TRUE) p <- p else p <- 1-p
    if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))     
        lp <-  max(length(p),length(mu),length(sigma))
          p <- rep(p, length = lp)                                                                      
      sigma <- rep(sigma, length = lp)
         mu <- rep(mu, length = lp)
          q <- rep(0,lp)      
         for (i in  seq(along=p)) 
         {
         if (h(mu[i])<p[i]) 
          { 
           interval <- c(mu[i], mu[i]+sigma[i])
           j <-2
           while (h(interval[2]) < p[i]) 
              {interval[2]<- mu[i]+j*sigma[i]
              j<-j+1 
              }
           } 
          else  
           {
           interval <-  interval <- c(.Machine$double.xmin, mu[i])
           }
        q[i] <- uniroot(h1, interval)$root
         }
    q
   }

#----------------------------------------------------------------------------------------  
rIG <- function(n, mu=1, sigma=1, ...)
  { 
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
    n <- ceiling(n)
    p <- runif(n)
    r <- qIG(p,mu=mu,sigma=sigma, ...)
    r
  }
