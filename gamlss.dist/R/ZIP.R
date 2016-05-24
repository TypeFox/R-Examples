# 3/12/2004
#------------------------------------------------------------------------------------------
# this is ZIP, the Poisson Zero Inflated distribution with extra probability for 0
ZIP <- function (mu.link = "log", sigma.link = "logit")
{
    mstats <- checklink("mu.link", "ZIP", substitute(mu.link),    
                           c("1/mu^2", "log", "identity"))
    dstats <- checklink("sigma.link", "ZIP", substitute(sigma.link), 
                           c("logit", "probit", "cloglog", "cauchit", "log", "own"))   
    structure(
          list(family = c("ZIP", "Poisson Zero Inflated"),
           parameters = list(mu=TRUE, sigma=TRUE), 
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
     dldm = function(y,mu,sigma) {dldm0 <- -(1-sigma)*(((1-sigma)+sigma*exp(mu))^(-1))
                         dldm <- ifelse(y==0, dldm0, (y/mu)-1)
                         dldm}, 
    d2ldm2 = function(y,mu,sigma) {dldm0 <- -(1-sigma)*(((1-sigma)+sigma*exp(mu))^(-1))
                          dldm <- ifelse(y==0, dldm0, (y/mu)-1)
                        d2ldm2 <- -dldm*dldm        
 # d2ldm2 <- (sigma*(1-sigma)*(((1-sigma)+sigma*exp(mu))^(-1)))-(1-sigma)/mu
                        d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)  
                         d2ldm2},
     dldd = function(y,mu,sigma) {dldd0 <- (1-exp(-mu))*((sigma+(1-sigma)*exp(-mu))^(-1)) 
                         dldd <- ifelse(y==0, dldd0, -1/(1-sigma))
                         dldd}, 
    d2ldd2 = function(y,mu,sigma) {dldd0 <- (1-exp(-mu))*((sigma+(1-sigma)*exp(-mu))^(-1)) 
                          dldd <- ifelse(y==0, dldd0, -1/(1-sigma))
                        d2ldd2 <- -dldd*dldd   
                       # py0 <- sigma+(1-sigma)*exp(-mu)
                       # d2ldd2 <- -((1-exp(-mu))^2)*py0 - (1-py0)/((1-sigma)^2)
                        d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
                         d2ldd2},
  d2ldmdd = function(y,mu,sigma) {dldm0 <- -(1-sigma)*(((1-sigma)+sigma*exp(mu))^(-1))
                         dldm <- ifelse(y==0, dldm0, (y/mu)-1)
                        dldd0 <- (1-exp(-mu))*((sigma+(1-sigma)*exp(-mu))^(-1)) 
                         dldd <- ifelse(y==0, dldd0, -1/(1-sigma)) 
                      d2ldmdd <- -dldm*dldd                     
#d2ldmdd <- exp(-mu) + (1-sigma)*exp(-mu)*(1-exp(-mu))*((sigma+(1-sigma)*exp(-mu))^(-1))  
                      d2ldmdd   
                       },
 G.dev.incr  = function(y,mu,sigma,...) -2*dZIP(y,mu,sigma,log=TRUE),                      
         rqres = expression(rqres(pfun="pZIP", type="Discrete", ymin=0, y=y, mu=mu, sigma=sigma)) ,
    mu.initial = expression(mu <- (y+mean(y))/2),     #rep(mean(y),length(y)) ), 
   sigma.initial = expression(sigma <-rep(0.1, length(y))), 
      mu.valid = function(mu) all(mu > 0) , 
   sigma.valid = function(sigma)  all(sigma > 0 & sigma < 1), 
       y.valid = function(y)  all(y >= 0)
          ),
            class = c("gamlss.family","family"))
}
#------------------------------------------------------------------------------------------
dZIP<-function(x, mu = 5, sigma = 0.1, log = FALSE)
 { 
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))           
          if (any(sigma <= 0) | any(sigma >= 1) )  stop(paste("sigma must be between 0 and 1", "\n", "")) 
          if (any(x < 0) )  stop(paste("x must be 0 or greater than 0", "\n", ""))   
          logfy <- rep(0, length(x))
          logfy <- ifelse((x==0), log(sigma+(1-sigma)*exp(-mu)), (log(1-sigma) -mu +x*log(mu) -lgamma(x+1)))          
          if(log == FALSE) fy <- exp(logfy) else fy <- logfy
          fy
  }
#------------------------------------------------------------------------------------------
pZIP <- function(q, mu = 5, sigma = 0.1, lower.tail = TRUE, log.p = FALSE)
  {     
         if (any(mu <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))           
         if (any(sigma <= 0) | any(sigma >= 1) )  stop(paste("sigma must be between 0 and 1", "\n", "")) 
         if (any(q < 0) )  stop(paste("y must be 0 or greater than 0", "\n", ""))  
         cdf <- rep(0,length(q))
         cdf <- ppois(q, lambda = mu, lower.tail = TRUE, log.p = FALSE)
         cdf <- sigma + (1-sigma)*cdf
         if(lower.tail == TRUE) cdf <- cdf else cdf <-1-cdf
         if(log.p==FALSE) cdf <- cdf else cdf <- log(cdf)    
         cdf
   }
#-----------------------------------------------------------------------------------------
qZIP <- function(p, mu = 5, sigma = 0.1, lower.tail = TRUE, log.p = FALSE)
  {      
         if (any(mu <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))           
         if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0", "\n", "")) 
         if (any(p <= 0) | any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", "")) 
         if (log.p == TRUE) p <- exp(p)   else p <- p
         if (lower.tail == TRUE)  p <- p  else p <- 1 - p
          pnew <- ((p-sigma)/(1-sigma))-(1e-7)# added Monday, March 15, 2010 
          pnew <- ifelse(pnew > 0, pnew,0 )
          q <- qpois(pnew, lambda = mu )
          q
   }
#-----------------------------------------------------------------------------------------
rZIP <- function(n, mu=5, sigma=0.1)
  { 
    if (any(mu <= 0) )  stop(paste("mu must greated than 0", "\n", ""))           
    if (any(sigma <= 0) )  stop(paste("sigma must greated than 0", "\n", "")) 
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
          n <- ceiling(n)
          p <- runif(n)
          r <- qZIP(p, mu = mu, sigma = sigma)
          r
  }
#-----------------------------------------------------------------------------------------
