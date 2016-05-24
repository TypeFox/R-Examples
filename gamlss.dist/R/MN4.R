# 6/12/2004
# this is MN4, the Multinomial distribution with probabilities for y values 1, 2, 3 and 4.
MN4 <- function (mu.link="log", sigma.link="log", nu.link="log")
{
    mstats <- checklink("mu.link", "MN4", substitute(mu.link),    
                           c("1/mu^2", "log", "identity"))
    dstats <- checklink("sigma.link", "MN4", substitute(sigma.link),   
                           c("1/sigma^2", "log", "identity")) 
    vstats <- checklink("nu.link", "MN4", substitute(nu.link),   
                           c("1/nu^2", "log", "identity")) 


    structure(
          list(family = c("MN4", "Multinomial 4 levels"),
           parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE), 
                nopar = 3, 
                 type = "Discrete",
              mu.link = as.character(substitute(mu.link)),  
           sigma.link = as.character(substitute(sigma.link)), 
              nu.link = as.character(substitute(nu.link)), 
           mu.linkfun = mstats$linkfun, 
        sigma.linkfun = dstats$linkfun, 
           nu.linkfun = vstats$linkfun,
           mu.linkinv = mstats$linkinv, 
        sigma.linkinv = dstats$linkinv,
           nu.linkinv = vstats$linkinv,
                mu.dr = mstats$mu.eta, 
             sigma.dr = dstats$mu.eta, 
                nu.dr = vstats$mu.eta,
     dldm = function(y,mu,sigma,nu) {dldm <- ifelse(y==1,(1/mu),0) -(1/(1+mu+sigma+nu))
                                     dldm}, 
     dldd = function(y,mu,sigma,nu) {dldd <- ifelse(y==2,(1/sigma),0) -(1/(1+mu+sigma+nu))
                                     dldd}, 
     dldv = function(y,mu,sigma,nu) {dldv <- ifelse(y==3,(1/nu),0) -(1/(1+mu+sigma+nu))
                                     dldv}, 
    d2ldm2 = function(mu,sigma,nu) {d2ldm2 <- -(1+sigma+nu)/(mu*((1+mu+sigma+nu)^2))
                                         d2ldm2},
    d2ldd2 = function(mu,sigma,nu) {d2ldd2 <- -(1+mu+nu)/(sigma*((1+mu+sigma+nu)^2))
                                     d2ldd2},
    d2ldv2 = function(mu,sigma,nu) {d2ldv2 <- -(1+mu+sigma)/(nu*((1+mu+sigma+nu)^2))
                                     d2ldv2},
  d2ldmdd = function(mu,sigma,nu) {d2ldmdd <- 1/((1+mu+sigma+nu)^2)  
                                    d2ldmdd},
  d2ldmdv = function(mu,sigma,nu) {d2ldmdv <- 1/((1+mu+sigma+nu)^2)  
                                    d2ldmdv},
  d2ldmdt = function(mu,sigma,nu) {d2ldmdt <- 1/((1+mu+sigma+nu)^2)  
                                    d2ldmdt},
  d2ldddv = function(mu,sigma,nu) {d2ldddv <- 1/((1+mu+sigma+nu)^2)  
                                    d2ldddv},

G.dev.incr  = function(y,mu,sigma,nu,...) 
                       { 
G.dev.incr <- -2*dMN4(y,mu,sigma,nu,log=TRUE)
                        } ,                     
         rqres = expression(   
               { stot <- 1+mu+sigma+nu                                                              
               uval <- ifelse(y==1,(mu/(stot))*runif(length(y),0,1),0)
               uval <- ifelse(y==2, (mu/(stot))+(sigma/(stot))*runif(length(y),0,1) ,uval)
               uval <- ifelse(y==3, ((mu+sigma)/(stot))+(nu/(stot))*runif(length(y),0,1) ,uval)
               uval <- ifelse(y==4, ((mu+sigma+nu)/(stot))+(1/(stot))*runif(length(y),0,1) ,uval)
               rqres <- qnorm(uval)  
                }) ,
    mu.initial = expression(mu <- rep(0.25, length(y))), 
   sigma.initial = expression(sigma <-rep(0.25, length(y))), 
    nu.initial = expression(nu <- rep(0.25, length(y))), 

      mu.valid = function(mu) all(mu > 0) , 
   sigma.valid = function(sigma)  all(sigma > 0), 
      nu.valid = function(nu) all(nu > 0) , 
       y.valid = function(y)  all(y >= 1 & y <= 4)
          ),
            class = c("gamlss.family","family"))
}

dMN4<-function(x, mu=1, sigma=1, nu=1, log=FALSE)
 { 
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))           
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0", "\n", "")) 
          if (any(nu <= 0) )  stop(paste("nu must be greater than 0", "\n", ""))           
          if (any(x < 1) | any(x > 5))  stop(paste("x must be 1, 2, 3 or 4", "\n", ""))  
          logfy <- rep(0, length(x))
          logfy <- ifelse((x==1), log(mu), logfy)          
          logfy <- ifelse((x==2), log(sigma) , logfy)
          logfy <- ifelse((x==3), log(nu) , logfy)
          logfy <- logfy - log(1+mu+sigma+nu)          
          if(log==FALSE) fy <- exp(logfy) else fy <- logfy
          fy

  }
  

pMN4 <- function(q, mu=1, sigma=1, nu=1, lower.tail = TRUE, log.p = FALSE)
  {     
         if (any(mu <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))           
         if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0", "\n", "")) 
         if (any(nu <= 0) )  stop(paste("nu must be greater than 0", "\n", ""))           
         if (any(q < 1) | any(q > 4))  stop(paste("y must be 1, 2, 3 or 4", "\n", ""))  
          cdf <- ifelse((q==1), mu, 0)          
          cdf <- ifelse((q==2), mu+sigma,cdf)
          cdf <- ifelse((q==3), mu+sigma+nu , cdf)
          cdf <- cdf/(1+mu+sigma+nu)          
          cdf <- ifelse((q==4), 1 , cdf)
          if(lower.tail==TRUE) cdf <- cdf else cdf=1-cdf
          if(log.p==FALSE) cdf <- cdf else cdf <- log(cdf)    
          cdf
   }

qMN4 <- function(p, mu=1, sigma=1, nu=1, lower.tail = TRUE, log.p = FALSE)
  {      
         if (any(mu <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))           
         if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0", "\n", "")) 
         if (any(nu <= 0) )  stop(paste("nu must be greater than 0", "\n", ""))           
         if (any(p <= 0) | any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
          q <- 1
          q <- ifelse((p>=(mu/(1+mu+sigma+nu))), 2, q) 
          q <- ifelse((p>=((mu+sigma)/(1+mu+sigma+nu))), 3, q)          
          q <- ifelse((p>=((mu+sigma+nu)/(1+mu+sigma+nu))), 4, q)          
          q
   }

rMN4 <- function(n, mu=1, sigma=1, nu=1)
  { 
    if (any(mu <= 0) )  stop(paste("mu must greated than 0", "\n", ""))           
    if (any(sigma <= 0) )  stop(paste("sigma must greated than 0", "\n", "")) 
    if (any(nu <= 0) )  stop(paste("nu must be greater than 0", "\n", ""))           
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
          n <- ceiling(n)
          p <- runif(n)
          r <- qMN4(p, mu=mu, sigma=sigma, nu=nu)
          r
  }
