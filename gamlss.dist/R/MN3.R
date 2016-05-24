# 3/12/2004 BR
# last change 2012 23-09-12
# this is MN3, the Multinomial distribution with probabilities for 1, 2 and 3
#------------------------------------------------------------------------------------------
MN3 <- function (mu.link = "log", sigma.link = "log")
{
    mstats <- checklink("mu.link", "MN3", substitute(mu.link),    
                           c("log"))
    dstats <- checklink("sigma.link", "MN3", substitute(sigma.link),   
                           c("log")) 
    structure(
          list(family = c("MN3", "Multinomial 3 levels"),
           parameters = list(mu = TRUE, sigma = TRUE), 
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
     dldm = function(y,mu,sigma) {dldm <- ifelse(y==1,(1/mu),0) -(1/(1+mu+sigma))
                                     dldm}, 
   d2ldm2 = function(mu,sigma) {d2ldm2 <- -(1+sigma)/(mu*((1+mu+sigma)^2))
                                         d2ldm2},
     dldd = function(y,mu,sigma) {dldd <- ifelse(y==2,(1/sigma),0) -(1/(1+mu+sigma))
                                     dldd}, 
   d2ldd2 = function(mu,sigma) {d2ldd2 <- -(1+mu)/(sigma*((1+mu+sigma)^2))
                                     d2ldd2},
  d2ldmdd = function(y,mu,sigma) {
                       d2ldmdd <- 1/((1+mu+sigma)^2)  
                       },
G.dev.incr  = function(y,mu,sigma,...) 
                       {G.dev.incr <- -2*dMN3(y, mu, sigma, log = TRUE)} ,                     
         rqres = expression(   
               {                                                               
               uval <- ifelse(y==1,(mu/(1+mu+sigma))*runif(length(y),0,1),0)
               uval <- ifelse(y==2, (mu/(1+mu+sigma))+(sigma/(1+mu+sigma))*
                                                 runif(length(y),0,1) ,uval)
               uval <- ifelse(y==3,((mu+sigma)/(1+mu+sigma))+(1/(1+mu+sigma))*
                                                 runif(length(y),0,1),uval)
              rqres <- qnorm(uval)  
                }) ,
      mu.initial = expression(mu <- rep(0.3, length(y))), 
   sigma.initial = expression(sigma <-rep(0.3, length(y))), 
        mu.valid = function(mu) all(mu > 0) , 
     sigma.valid = function(sigma)  all(sigma > 0), 
         y.valid = function(y)   if (is.factor(y))  nlevels(y)==3
                                 else all(y >= 1 & y <= 3)
          ),
           class = c("gamlss.family","family"))
}
#------------------------------------------------------------------------------------------
dMN3<-function(x, mu=1, sigma=1, log=FALSE)
 { 
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))           
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0", "\n", "")) 
          if (is.factor(x)) 
               {
                if (nlevels(x)!=3) stop(paste("x levels must be 3", "\n", "")) 
               }
          else {  
               if (any(x < 1) | any(x > 3))  stop(paste("x must be 1, 2 or 3", "\n", ""))
               }   
          logfy <- rep(0, length(x))
          logfy <- ifelse((x==1), log(mu), logfy)          
          logfy <- ifelse((x==2), log(sigma) , logfy)
          logfy <- logfy - log(1+mu+sigma)          
          if(log == FALSE) fy <- exp(logfy) else fy <- logfy
          fy
  }
#------------------------------------------------------------------------------------------
pMN3 <- function(q, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE)
  {     
         if (any(mu <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))           
         if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0", "\n", "")) 
         if (any(q < 1) | any(q > 3))  stop(paste("y must be 1, 2 or 3", "\n", ""))  
          cdf <- ifelse((q==1), mu, 0)          
          cdf <- ifelse ((q==2),mu+sigma,cdf)
          cdf <- cdf/(1+mu+sigma)          
          cdf <- ifelse((q==3), 1 , cdf)
          if(lower.tail==TRUE) cdf <- cdf else cdf=1-cdf
          if(log.p==FALSE) cdf <- cdf else cdf <- log(cdf)    
          cdf
   }
#------------------------------------------------------------------------------------------
qMN3 <- function(p, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE)
  {      
         if (any(mu <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))           
         if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0", "\n", "")) 
         if (any(p <= 0) | any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
          q <- 1
          q <- ifelse((p>=(mu/(1+mu+sigma))), 2, q) 
          q <- ifelse((p>=((mu+sigma)/(1+mu+sigma))), 3, q)          
          q
   }
#------------------------------------------------------------------------------------------
rMN3 <- function(n, mu=1, sigma=1)
  { 
    if (any(mu <= 0) )  stop(paste("mu must greated than 0", "\n", ""))           
    if (any(sigma <= 0) )  stop(paste("sigma must greated than 0", "\n", "")) 
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
          n <- ceiling(n)
          p <- runif(n)
          r <- qMN3(p, mu=mu, sigma=sigma)
          r
  }
#------------------------------------------------------------------------------------------
