# 6/12/2004
# this is MN5, the Multinomial distribution with probabilities for y values 1, 2, 3, 4 and 5.
MN5 <- function (mu.link="log", sigma.link="log", nu.link="log", tau.link="log")
{
    mstats <- checklink("mu.link", "MN5", substitute(mu.link),    
                           c("1/mu^2", "log", "identity"))
    dstats <- checklink("sigma.link", "MN5", substitute(sigma.link),   
                           c("1/sigma^2", "log", "identity")) 
    vstats <- checklink("nu.link", "MN5", substitute(nu.link),   
                           c("1/nu^2", "log", "identity")) 
    tstats <- checklink("tau.link", "MN5", substitute(tau.link),   
                           c("1/tau^2", "log", "identity")) 


    structure(
          list(family = c("MN5", "Multinomial 5 levels"),
           parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE, tau=TRUE), 
                nopar = 4, 
                 type = "Discrete",
              mu.link = as.character(substitute(mu.link)),  
           sigma.link = as.character(substitute(sigma.link)), 
              nu.link = as.character(substitute(nu.link)), 
             tau.link = as.character(substitute(tau.link)), 
           mu.linkfun = mstats$linkfun, 
        sigma.linkfun = dstats$linkfun, 
           nu.linkfun = vstats$linkfun,
           tau.linkfun = tstats$linkfun,  
           mu.linkinv = mstats$linkinv, 
        sigma.linkinv = dstats$linkinv,
           nu.linkinv = vstats$linkinv,
           tau.linkinv = tstats$linkinv, 
                mu.dr = mstats$mu.eta, 
             sigma.dr = dstats$mu.eta, 
                nu.dr = vstats$mu.eta,
               tau.dr = tstats$mu.eta, 


     dldm = function(y,mu,sigma,nu,tau) {dldm <- ifelse(y==1,(1/mu),0) -(1/(1+mu+sigma+nu+tau))
                                     dldm}, 
                        
     dldd = function(y,mu,sigma,nu,tau) {dldd <- ifelse(y==2,(1/sigma),0) -(1/(1+mu+sigma+nu+tau))
                                     dldd}, 

     dldv = function(y,mu,sigma,nu,tau) {dldv <- ifelse(y==3,(1/nu),0) -(1/(1+mu+sigma+nu+tau))
                                     dldv}, 

     dldt= function(y,mu,sigma,nu,tau) {dldt <- ifelse(y==4,(1/tau),0) -(1/(1+mu+sigma+nu+tau))
                                     dldt}, 
                        

    d2ldm2 = function(mu,sigma,nu,tau) {d2ldm2 <- -(1+sigma+nu+tau)/(mu*((1+mu+sigma+nu+tau)^2))
                                         d2ldm2},

    d2ldd2 = function(mu,sigma,nu,tau) {d2ldd2 <- -(1+mu+nu+tau)/(sigma*((1+mu+sigma+nu+tau)^2))
                                     d2ldd2},

    d2ldv2 = function(mu,sigma,nu,tau) {d2ldv2 <- -(1+mu+sigma+tau)/(nu*((1+mu+sigma+nu+tau)^2))
                                     d2ldv2},

    d2ldt2 = function(y,mu,sigma,nu,tau) {d2ldt2 <- -(1+mu+sigma+nu)/(tau*((1+mu+sigma+nu+tau)^2))
                                     d2ldt2},


  d2ldmdd = function(mu,sigma,nu,tau) {d2ldmdd <- 1/((1+mu+sigma+nu+tau)^2)  
                                    d2ldmdd},

  d2ldmdv = function(mu,sigma,nu,tau) {d2ldmdv <- 1/((1+mu+sigma+nu+tau)^2)  
                                    d2ldmdv},

  d2ldmdt = function(mu,sigma,nu,tau) {d2ldmdt <- 1/((1+mu+sigma+nu+tau)^2)  
                                    d2ldmdt},

  d2ldddv = function(mu,sigma,nu,tau) {d2ldddv <- 1/((1+mu+sigma+nu+tau)^2)  
                                    d2ldddv},

  d2ldddt = function(mu,sigma,nu,tau) {d2ldddt <- 1/((1+mu+sigma+nu+tau)^2)  
                                    d2ldddt},

  d2ldvdt = function(mu,sigma,nu,tau) {d2ldvdt <- 1/((1+mu+sigma+nu+tau)^2)  
                                    d2ldvdt},


G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
                       { 
G.dev.incr <- -2*dMN5(y,mu,sigma,nu,tau,log=TRUE)
                        } ,                     
         rqres = expression(   
               { stot <- 1+mu+sigma+nu+tau                                                              
               uval <- ifelse(y==1,(mu/(stot))*runif(length(y),0,1),0)
               uval <- ifelse(y==2, (mu/(stot))+(sigma/(stot))*runif(length(y),0,1) ,uval)
               uval <- ifelse(y==3, ((mu+sigma)/(stot))+(nu/(stot))*runif(length(y),0,1) ,uval)
               uval <- ifelse(y==4, ((mu+sigma+nu)/(stot))+(tau/(stot))*runif(length(y),0,1) ,uval)
               uval <- ifelse(y==5, ((mu+sigma+nu+tau)/(stot))+(1/(stot))*runif(length(y),0,1) ,uval)
               rqres <- qnorm(uval)  
                }) ,
    mu.initial = expression(mu <- rep(0.2, length(y))), 
   sigma.initial = expression(sigma <-rep(0.2, length(y))), 
    nu.initial = expression(nu <- rep(0.2, length(y))), 
   tau.initial = expression(tau <-rep(0.2, length(y))), 

      mu.valid = function(mu) all(mu > 0) , 
   sigma.valid = function(sigma)  all(sigma > 0), 
      nu.valid = function(nu) all(nu > 0) , 
     tau.valid = function(tau)  all(tau > 0), 
       y.valid = function(y)  all(y >= 1 & y <= 5)
          ),
            class = c("gamlss.family","family"))
}

dMN5<-function(x, mu=1, sigma=1, nu=1, tau=1, log=FALSE)
 { 
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))           
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0", "\n", "")) 
          if (any(nu <= 0) )  stop(paste("nu must be greater than 0", "\n", ""))           
          if (any(tau <= 0) )  stop(paste("tau must be greater than 0", "\n", "")) 
          if (any(x < 1) | any(x > 5))  stop(paste("x must be 1, 2, 3, 4 or 5", "\n", ""))  
          logfy <- rep(0, length(x))
          logfy <- ifelse((x==1), log(mu), logfy)          
          logfy <- ifelse((x==2), log(sigma) , logfy)
          logfy <- ifelse((x==3), log(nu) , logfy)
          logfy <- ifelse((x==4), log(tau) , logfy)
          logfy <- logfy - log(1+mu+sigma+nu+tau)          
          if(log==FALSE) fy <- exp(logfy) else fy <- logfy
          fy

  }
  

pMN5 <- function(q, mu=1, sigma=1, nu=1, tau=1, lower.tail = TRUE, log.p = FALSE)
  {     
         if (any(mu <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))           
         if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0", "\n", "")) 
         if (any(nu <= 0) )  stop(paste("nu must be greater than 0", "\n", ""))           
         if (any(tau <= 0) )  stop(paste("tau must be greater than 0", "\n", "")) 
         if (any(q < 1) | any(q > 5))  stop(paste("y must be 1, 2, 3, 4 or 5", "\n", ""))  
          cdf <- ifelse((q==1), mu, 0)          
          cdf <- ifelse((q==2), mu+sigma,cdf)
          cdf <- ifelse((q==3), mu+sigma+nu , cdf)
          cdf <- ifelse((q==4), mu+sigma+nu+tau , cdf)
          cdf <- cdf/(1+mu+sigma+nu+tau)          
          cdf <- ifelse((q==5), 1 , cdf)
          if(lower.tail==TRUE) cdf <- cdf else cdf=1-cdf
          if(log.p==FALSE) cdf <- cdf else cdf <- log(cdf)    
          cdf
   }

qMN5 <- function(p, mu=1, sigma=1, nu=1, tau=1, lower.tail = TRUE, log.p = FALSE)
  {      
         if (any(mu <= 0) )  stop(paste("mu must be greater than 0", "\n", ""))           
         if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0", "\n", "")) 
         if (any(nu <= 0) )  stop(paste("nu must be greater than 0", "\n", ""))           
         if (any(tau <= 0) )  stop(paste("tau must be greater than 0", "\n", "")) 
         if (any(p <= 0) | any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
          q <- 1
          q <- ifelse((p>=(mu/(1+mu+sigma+nu+tau))), 2, q) 
          q <- ifelse((p>=((mu+sigma)/(1+mu+sigma+nu+tau))), 3, q)          
          q <- ifelse((p>=((mu+sigma+nu)/(1+mu+sigma+nu+tau))), 4, q)          
          q <- ifelse((p>=((mu+sigma+nu+tau)/(1+mu+sigma+nu+tau))), 5, q)          
          q
   }

rMN5 <- function(n, mu=1, sigma=1, nu=1, tau=1)
  { 
    if (any(mu <= 0) )  stop(paste("mu must greated than 0", "\n", ""))           
    if (any(sigma <= 0) )  stop(paste("sigma must greated than 0", "\n", "")) 
    if (any(nu <= 0) )  stop(paste("nu must be greater than 0", "\n", ""))           
    if (any(tau <= 0) )  stop(paste("tau must be greater than 0", "\n", "")) 
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
          n <- ceiling(n)
          p <- runif(n)
          r <- qMN5(p, mu=mu, sigma=sigma, nu=nu, tau=tau)
          r
  }
