#  27_11_2007
# the first derivatives squares have been used here
SHASH <- function (mu.link="identity", sigma.link="log", nu.link ="log", tau.link="log")
{
    mstats <- checklink(   "mu.link", "Sinh-Arcsinh", substitute(mu.link), 
                           c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Sinh-Arcsinh", substitute(sigma.link), 
                           c("inverse", "log", "identity", "own"))
    vstats <- checklink(   "nu.link", "Sinh-Arcsinh", substitute(nu.link),    
                           c("inverse", "log", "identity", "own"))
    tstats <- checklink(  "tau.link", "Sinh-Arcsinh", substitute(tau.link),   
                           c("inverse", "log", "identity", "own")) 
    structure(
          list(family = c("SHASH", "Sinh-Arcsinh"),
           parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE, tau=TRUE), 
                nopar = 4, 
                 type = "Continuous",
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
    dldm = function(y,mu,sigma,nu,tau) {
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
      h <- (1/2)*((tau^2)*exp(tau*asinh(z))-(nu^2)*exp(-nu*asinh(z)))
   dldz <- -z/(1+(z^2))
   dldr <- -r
   dldc <- 1/c
   dcdz <- h*((1+(z^2))^(-1/2))
   drdz <- c*((1+(z^2))^(-1/2))
   dzdm <- -1/sigma
   dldm <- (1/sigma)*((1+(z^2))^(-1/2))*((-h/c)+(r*c)+z*((1+(z^2))^(-1/2)))
   dldm <- (dldr*drdz+dldc*dcdz+dldz)*dzdm
   dldm 
                      },
   d2ldm2 = function(y,mu,sigma,nu,tau){
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
      h <- (1/2)*((tau^2)*exp(tau*asinh(z))-(nu^2)*exp(-nu*asinh(z)))
   dldz <- -z/(1+(z^2))
   dldr <- -r
   dldc <- 1/c
   dcdz <- h*((1+(z^2))^(-1/2))
   drdz <- c*((1+(z^2))^(-1/2))
   dzdm <- -1/sigma
   dldm <- (1/sigma)*((1+(z^2))^(-1/2))*((-h/c)+(r*c)+z*((1+(z^2))^(-1/2)))
   dldm <- (dldr*drdz+dldc*dcdz+dldz)*dzdm
 d2ldm2 <- -dldm*dldm
 d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15) 
 d2ldm2
                      },     
   dldd = function(y,mu,sigma,nu,tau) {  
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
      h <- (1/2)*((tau^2)*exp(tau*asinh(z))-(nu^2)*exp(-nu*asinh(z)))
   dldz <- -z/(1+(z^2))
   dldr <- -r
   dldc <- 1/c
   dcdz <- h*((1+(z^2))^(-1/2))
   drdz <- c*((1+(z^2))^(-1/2))
   dzdd <- -z/sigma
   dldd <- (dldr*drdz+dldc*dcdz+dldz)*dzdd-1/sigma 
   dldd                    
                      } ,
   d2ldd2 = function(y,mu,sigma,nu,tau){
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
      h <- (1/2)*((tau^2)*exp(tau*asinh(z))-(nu^2)*exp(-nu*asinh(z)))
   dldz <- -z/(1+(z^2))
   dldr <- -r
   dldc <- 1/c
   dcdz <- h*((1+(z^2))^(-1/2))
   drdz <- c*((1+(z^2))^(-1/2))
   dzdd <- -z/sigma
   dldd <- (dldr*drdz+dldc*dcdz+dldz)*dzdd-1/sigma    
 d2ldd2 <- -dldd*dldd
 d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
 d2ldd2
                      },   
     dldv = function(y,mu,sigma,nu,tau) { 
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
   dldr <- -r
   dldc <- 1/c
   drdv <- (1/2)*asinh(z)*exp(-nu*asinh(z))
   dcdv <- (1/2)*(1-nu*asinh(z))*exp(-nu*asinh(z))
      dldv <- dldr*drdv+dldc*dcdv
      dldv
                        } ,
    d2ldv2 = function(y,mu,sigma,nu,tau) { 
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
   dldr <- -r
   dldc <- 1/c
   drdv <- (1/2)*asinh(z)*exp(-nu*asinh(z))
   dcdv <- (1/2)*(1-nu*asinh(z))*exp(-nu*asinh(z))
   dldv <- dldr*drdv+dldc*dcdv
   d2ldv2 <-  -dldv*dldv             
   d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)  
   d2ldv2
                        },
      dldt = function(y,mu,sigma,nu,tau) {
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
   dldr <- -r
   dldc <- 1/c
   drdt <- (1/2)*asinh(z)*exp(tau*asinh(z))
   dcdt <- (1/2)*(1+tau*asinh(z))*exp(tau*asinh(z))
   dldt <- dldr*drdt+dldc*dcdt
   dldt   
                        } ,
      d2ldt2 = function(y,mu,sigma,nu,tau) { 
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
   dldr <- -r
   dldc <- 1/c
   drdt <- (1/2)*asinh(z)*exp(tau*asinh(z))
   dcdt <- (1/2)*(1+tau*asinh(z))*exp(tau*asinh(z))
   dldt <- dldr*drdt+dldc*dcdt
   d2ldt2 <-   -dldt*dldt   
   d2ldt2 <- ifelse(d2ldt2 < -1e-15, d2ldt2,-1e-15)                                    
   d2ldt2
                            } ,
       d2ldmdd = function(y,mu,sigma,nu,tau)## ok
               {
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
      h <- (1/2)*((tau^2)*exp(tau*asinh(z))-(nu^2)*exp(-nu*asinh(z)))
   dldz <- -z/(1+(z^2))
   dldr <- -r
   dldc <- 1/c
   dcdz <- h*((1+(z^2))^(-1/2))
   drdz <- c*((1+(z^2))^(-1/2))
   dzdm <- -1/sigma
   dldm <- (1/sigma)*((1+(z^2))^(-1/2))*((-h/c)+(r*c)+z*((1+(z^2))^(-1/2)))
   dldm <- (dldr*drdz+dldc*dcdz+dldz)*dzdm     
   dzdd <- -z/sigma
   dldd <- (dldr*drdz+dldc*dcdz+dldz)*dzdd-1/sigma           
d2ldmdd <- -(dldm*dldd)
d2ldmdd
               },
       d2ldmdv = function(y,mu,sigma,nu,tau)# OK
               { 
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
      h <- (1/2)*((tau^2)*exp(tau*asinh(z))-(nu^2)*exp(-nu*asinh(z)))
   dldz <- -z/(1+(z^2))
   dldr <- -r
   dldc <- 1/c
   dcdz <- h*((1+(z^2))^(-1/2))
   drdz <- c*((1+(z^2))^(-1/2))
   dzdm <- -1/sigma
   dldm <- (1/sigma)*((1+(z^2))^(-1/2))*((-h/c)+(r*c)+z*((1+(z^2))^(-1/2)))
   dldm <- (dldr*drdz+dldc*dcdz+dldz)*dzdm                   
   drdv <- (1/2)*asinh(z)*exp(-nu*asinh(z))
   dcdv <- (1/2)*(1-nu*asinh(z))*exp(-nu*asinh(z))
   dldv <- dldr*drdv+dldc*dcdv
d2ldmdv <- -(dldm*dldv)
d2ldmdv
               },
       d2ldmdt = function(y,mu,sigma,nu,tau) #ok
               {
          z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
      h <- (1/2)*((tau^2)*exp(tau*asinh(z))-(nu^2)*exp(-nu*asinh(z)))
   dldz <- -z/(1+(z^2))
   dldr <- -r
   dldc <- 1/c
   dcdz <- h*((1+(z^2))^(-1/2))
   drdz <- c*((1+(z^2))^(-1/2))
   dzdm <- -1/sigma
   dldm <- (1/sigma)*((1+(z^2))^(-1/2))*((-h/c)+(r*c)+z*((1+(z^2))^(-1/2)))
   dldm <- (dldr*drdz+dldc*dcdz+dldz)*dzdm     
   drdt <- (1/2)*asinh(z)*exp(tau*asinh(z))
   dcdt <- (1/2)*(1+tau*asinh(z))*exp(tau*asinh(z))
   dldt <- dldr*drdt+dldc*dcdt                          
d2ldmdt <- -(dldm*dldt)
d2ldmdt
               },
       d2ldddv = function(y,mu,sigma,nu,tau) #ok
               {               
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
      h <- (1/2)*((tau^2)*exp(tau*asinh(z))-(nu^2)*exp(-nu*asinh(z)))
   dldz <- -z/(1+(z^2))
   dldr <- -r
   dldc <- 1/c
   dcdz <- h*((1+(z^2))^(-1/2))
   drdz <- c*((1+(z^2))^(-1/2))
   dzdd <- -z/sigma
   dldd <- (dldr*drdz+dldc*dcdz+dldz)*dzdd-1/sigma       
   drdv <- (1/2)*asinh(z)*exp(-nu*asinh(z))
   dcdv <- (1/2)*(1-nu*asinh(z))*exp(-nu*asinh(z))
   dldv <- dldr*drdv+dldc*dcdv
d2ldddv <- -(dldd*dldv)
d2ldddv 
               },
       d2ldddt = function(y,mu,sigma,nu,tau) #ok
               {
     z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
      h <- (1/2)*((tau^2)*exp(tau*asinh(z))-(nu^2)*exp(-nu*asinh(z)))
   dldz <- -z/(1+(z^2))
   dldr <- -r
   dldc <- 1/c
   dcdz <- h*((1+(z^2))^(-1/2))
   drdz <- c*((1+(z^2))^(-1/2))
   dzdd <- -z/sigma
   dldd <- (dldr*drdz+dldc*dcdz+dldz)*dzdd-1/sigma   
   drdt <- (1/2)*asinh(z)*exp(tau*asinh(z))
   dcdt <- (1/2)*(1+tau*asinh(z))*exp(tau*asinh(z))
   dldt <- dldr*drdt+dldc*dcdt   
d2ldddt <- -(dldd*dldt) 
d2ldddt 
               },
       d2ldvdt = function(y,mu,sigma,nu,tau) #ok
               { 
      z <- (y-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
   dldr <- -r
   dldc <- 1/c
   drdv <- (1/2)*asinh(z)*exp(-nu*asinh(z))
   dcdv <- (1/2)*(1-nu*asinh(z))*exp(-nu*asinh(z))
   dldv <- dldr*drdv+dldc*dcdv                         
   drdt <- (1/2)*asinh(z)*exp(tau*asinh(z))
   dcdt <- (1/2)*(1+tau*asinh(z))*exp(tau*asinh(z))
   dldt <- dldr*drdt+dldc*dcdt             
d2ldvdt <- -(dldv*dldt) 
d2ldvdt 
               },
 G.dev.incr  = function(y,mu,sigma,nu,tau,...) -2*dSHASH(y,mu,sigma,nu,tau,log=TRUE),                 
         rqres = expression(   
                   rqres(pfun="pSHASH", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)
                           ),
    mu.initial = expression(mu <- (y+mean(y))/2),   
 sigma.initial = expression(sigma<- rep(sd(y)/5, length(y))),
    nu.initial = expression(nu <- rep(.5, length(y))), 
   tau.initial = expression(tau <-rep(.5, length(y))), 
      mu.valid = function(mu) TRUE, 
   sigma.valid = function(sigma)  all(sigma > 0),
      nu.valid = function(nu) TRUE , 
     tau.valid = function(tau) all(tau > 0), 
       y.valid = function(y)  TRUE
          ),
            class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------
dSHASH <- function(x, mu = 0, sigma = 1, nu = .5, tau = .5, log = FALSE)
 {
          if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
          if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))
      z <- (x-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      c <- (1/2)*(tau*exp(tau*asinh(z))+nu*exp(-nu*asinh(z)))
 loglik <- -log(sigma)-log(2*pi)/2-log(1+(z^2))/2+log(c)-(r^2)/2
       if(log==FALSE) ft  <- exp(loglik) else ft <- loglik 
       ft
  }    
#-----------------------------------------------------------------  
pSHASH <- function(q, mu = 0, sigma = 1, nu = .5, tau = .5, lower.tail = TRUE, log.p = FALSE)
 {  
         if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
         if (any(tau < 0))  stop(paste("tau must be positive", "\n", "")) 
         if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))          
      z <- (q-mu)/sigma 
      r <- (1/2)*(exp(tau*asinh(z))-exp(-nu*asinh(z)))
      p <- pNO(r)
      if(lower.tail==TRUE) p  <- p else  p <- 1-p 
      if(log.p==FALSE) p  <- p else  p <- log(p) 
      p
 }
#-----------------------------------------------------------------  
qSHASH <-  function(p, mu=0, sigma=1, nu=.5, tau=.5, lower.tail = TRUE, log.p = FALSE)
  { 
    #---functions--------------------------------------------   
       h1 <- function(q)
       { 
     pSHASH(q , mu = mu[i], sigma = sigma[i], nu = nu[i], tau = tau[i]) - p[i] 
       }
       h <- function(q)
       { 
     pSHASH(q , mu = mu[i], sigma = sigma[i], nu = nu[i], tau = tau[i]) 
       }
     #-----------------------------------------------------------------
    #if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))      
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (lower.tail==TRUE) p <- p else p <- 1-p
    if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))   
         lp <-  max(length(p),length(mu),length(sigma),length(nu), length(tau))
          p <- rep(p, length = lp)                                                                   
      sigma <- rep(sigma, length = lp)
         mu <- rep(mu, length = lp)
         nu <- rep(nu, length = lp)
         tau <- rep(tau, length = lp)
           q <- rep(0,lp)
         for (i in 1:lp) 
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
           interval <-  c(mu[i]-sigma[i], mu[i])
           j <-2
           while (h(interval[1]) > p[i]) 
              {interval[1]<- mu[i]-j*sigma[i]
              j<-j+1 
              }
           }
        q[i] <- uniroot(h1, interval)$root
        #interval <- c(.Machine$double.xmin, 20)
         }
    q
   }
#-----------------------------------------------------------------  
rSHASH <- function(n, mu=0, sigma=1, nu=.5, tau=.5)
  {
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    n <- ceiling(n)
    p <- runif(n)
    r <- qSHASH(p,mu=mu,sigma=sigma,nu=nu,tau=tau)
   r
  }
#-----------------------------------------------------------------  
