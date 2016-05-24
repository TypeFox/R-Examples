# 27_11_2007
# the first derivatives squares have been used here
JSU <- function (mu.link="identity", sigma.link="log", nu.link ="identity", tau.link="log")
{
    mstats <- checklink(   "mu.link", "Johnson SU", substitute(mu.link), 
                           c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Johnson SU", substitute(sigma.link), 
                           c("inverse", "log", "identity", "own"))
    vstats <- checklink(   "nu.link", "Johnson SU", substitute(nu.link),    
                           c("inverse", "log", "identity", "own"))
    tstats <- checklink(  "tau.link", "Johnson SU", substitute(tau.link),   
                           c("inverse", "log", "identity", "own")) 
    structure(
          list(family = c("JSU", "Johnson SU"),
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
   rtau <- 1/tau
     w <- ifelse(rtau<0.0000001,1,exp(rtau^2))
  omega <- -nu*rtau 
      c <- (.5*(w-1)*(w*cosh(2*omega)+1))^(-0.5)  
      z <- (y-(mu+c*sigma*w^(.5)*sinh(omega)))/(c*sigma)
      r <- -nu + asinh(z)/rtau
   dldm <- (z/(z*z+1)+(r/(rtau*(z*z+1)^(.5))))/(c*sigma)
   dldm
                      },
   d2ldm2 = function(y,mu,sigma,nu,tau){
   rtau <- 1/tau
     w <- ifelse(rtau<0.0000001,1,exp(rtau^2))
  omega <- -nu*rtau 
      c <- (.5*(w-1)*(w*cosh(2*omega)+1))^(-0.5)  
      z <- (y-(mu+c*sigma*w^(.5)*sinh(omega)))/(c*sigma)
      r <- -nu + asinh(z)/rtau
   dldm <- (z/(z*z+1)+(r/(rtau*(z*z+1)^(.5))))/(c*sigma)
   d2ldm2 <- -dldm*dldm
    d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)                                    
    d2ldm2
                      },     
   dldd = function(y,mu,sigma,nu,tau) {  
   rtau <- 1/tau
     w <- ifelse(rtau<0.0000001,1,exp(rtau^2))
  omega <- -nu*rtau 
      c <- (.5*(w-1)*(w*cosh(2*omega)+1))^(-0.5)  
      z <- (y-(mu+c*sigma*w^(.5)*sinh(omega)))/(c*sigma)
      r <- -nu + asinh(z)/rtau
      dldd <- (z+w^(0.5)*sinh(omega))*(z/(z*z+1)+(r/(rtau*(z*z+1)^(.5))))
      dldd <- (dldd-1)/sigma  
      dldd                    
      } ,
   d2ldd2 = function(y,mu,sigma,nu,tau){
   rtau <- 1/tau
     w <- ifelse(rtau<0.0000001,1,exp(rtau^2))
  omega <- -nu*rtau 
      c <- (.5*(w-1)*(w*cosh(2*omega)+1))^(-0.5)  
      z <- (y-(mu+c*sigma*w^(.5)*sinh(omega)))/(c*sigma)
      r <- -nu + asinh(z)/rtau
      dldd <- (z+w^(0.5)*sinh(omega))*(z/(z*z+1)+(r/(rtau*(z*z+1)^(.5))))
      dldd <- (dldd-1)/sigma     
      d2ldd2 <- -dldd*dldd
   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
   d2ldd2
                      },   
     dldv = function(y,mu,sigma,nu,tau) { 
   rtau <- 1/tau
     w <- ifelse(rtau<0.0000001,1,exp(rtau^2))
  omega <- -nu*rtau 
      c <- (.5*(w-1)*(w*cosh(2*omega)+1))^(-0.5)  
      z <- (y-(mu+c*sigma*w^(.5)*sinh(omega)))/(c*sigma)
      r <- -nu + asinh(z)/rtau
      dlogcdv <- (rtau*w*sinh(2*omega))/(w*cosh(2*omega)+1)
      dzdv <- -(z+w^(.5)*sinh(omega))*dlogcdv+(rtau*w^(.5)*cosh(omega))
      dldv <- -dlogcdv-(z/(z*z+1)+(r/(rtau*(z*z+1)^(.5))))*dzdv+r
      dldv
                        } ,
    d2ldv2 = function(y,mu,sigma,nu,tau) { 
   rtau <- 1/tau
    w <- ifelse(rtau<0.0000001,1,exp(rtau^2))
  omega <- -nu*rtau 
      c <- (.5*(w-1)*(w*cosh(2*omega)+1))^(-0.5)  
      z <- (y-(mu+c*sigma*w^(.5)*sinh(omega)))/(c*sigma)
      r <- -nu + asinh(z)/rtau
      dlogcdv <- (rtau*w*sinh(2*omega))/(w*cosh(2*omega)+1)
      dzdv <- -(z+w^(.5)*sinh(omega))*dlogcdv+(rtau*w^(.5)*cosh(omega))
      dldv <- -dlogcdv-(z/(z*z+1)+(r/(rtau*(z*z+1)^(.5))))*dzdv+r
      d2ldv2 <-  -dldv*dldv             
   d2ldv2 <- ifelse(d2ldv2 < -1e-4, d2ldv2,-1e-4)  
   d2ldv2
                        },
      dldt = function(y,mu,sigma,nu,tau) {
   rtau <- 1/tau
     w <- ifelse(rtau<0.0000001,1,exp(rtau^2))
  omega <- -nu*rtau 
      c <- (.5*(w-1)*(w*cosh(2*omega)+1))^(-0.5)  
      z <- (y-(mu+c*sigma*w^(.5)*sinh(omega)))/(c*sigma)
      r <- -nu + asinh(z)/rtau
      dlogcdt <- -rtau*w*((1/(w-1))+((cosh(2*omega))/(w*cosh(2*omega)+1)))
      dlogcdt <- dlogcdt+((nu*w*sinh(2*omega))/(w*cosh(2*omega)+1))
      dzdt <- -(z+w^(.5)*sinh(omega))*dlogcdt-rtau*w^(.5)*sinh(omega)+nu*w^(.5)*cosh(omega)  
      dldt <- -dlogcdt-(1/rtau)-(z/(z*z+1)+(r/(rtau*(z*z+1)^(.5))))*dzdt+(r*(r+nu))/rtau
      dldt <- -dldt*rtau*rtau
      dldt
                        } ,
      d2ldt2 = function(y,mu,sigma,nu,tau) { 
   rtau <- 1/tau
     w <- ifelse(rtau<0.0000001,1,exp(rtau^2))
  omega <- -nu*rtau 
      c <- (.5*(w-1)*(w*cosh(2*omega)+1))^(-0.5)  
      z <- (y-(mu+c*sigma*w^(.5)*sinh(omega)))/(c*sigma)
      r <- -nu + asinh(z)/rtau
      dlogcdt <- -rtau*w*((1/(w-1))+((cosh(2*omega))/(w*cosh(2*omega)+1)))
      dlogcdt <- dlogcdt+((nu*w*sinh(2*omega))/(w*cosh(2*omega)+1))
      dzdt <- -(z+w^(.5)*sinh(omega))*dlogcdt-rtau*w^(.5)*sinh(omega)+nu*w^(.5)*cosh(omega)  
      dldt <- -dlogcdt-(1/rtau)-(z/(z*z+1)+(r/(rtau*(z*z+1)^(.5))))*dzdt+(r*(r+nu))/rtau
      dldt <- -dldt*rtau*rtau
      d2ldt2 <-   -dldt*dldt   
    d2ldt2 <- ifelse(d2ldt2 < -1e-4, d2ldt2,-1e-4)                                    
    d2ldt2
                            } ,
  d2ldmdd = function(y,mu,sigma,nu,tau) {
   rtau <- 1/tau
      w <- ifelse(rtau<0.0000001,1,exp(rtau^2))
   omega <- -nu*rtau 
       c <- (.5*(w-1)*(w*cosh(2*omega)+1))^(-0.5)  
       z <- (y-(mu+c*sigma*w^(.5)*sinh(omega)))/(c*sigma)
       r <- -nu + asinh(z)/rtau
    dldm <- (z/(z*z+1)+(r/(rtau*(z*z+1)^(.5))))/(c*sigma)
    dldd <- (z+w^(0.5)*sinh(omega))*(z/(z*z+1)+(r/(rtau*(z*z+1)^(.5))))
    dldd <- (dldd-1)/sigma  
 d2ldmdd <- -(dldm*dldd)
 d2ldmdd
                       },
  d2ldmdv = function(y,mu,sigma,nu,tau) {
   rtau <- 1/tau
      w <- ifelse(rtau<0.0000001,1,exp(rtau^2))
   omega <- -nu*rtau 
       c <- (.5*(w-1)*(w*cosh(2*omega)+1))^(-0.5)  
       z <- (y-(mu+c*sigma*w^(.5)*sinh(omega)))/(c*sigma)
       r <- -nu + asinh(z)/rtau
    dldm <- (z/(z*z+1)+(r/(rtau*(z*z+1)^(.5))))/(c*sigma)
 dlogcdv <- (rtau*w*sinh(2*omega))/(w*cosh(2*omega)+1)
    dzdv <- -(z+w^(.5)*sinh(omega))*dlogcdv+(rtau*w^(.5)*cosh(omega))
    dldv <- -dlogcdv-(z/(z*z+1)+(r/(rtau*(z*z+1)^(.5))))*dzdv+r   
 d2ldmdv <- -(dldm*dldv)
 d2ldmdv
                       },
  d2ldmdt = function(y,mu,sigma,nu,tau) {
   rtau <- 1/tau
      w <- ifelse(rtau<0.0000001,1,exp(rtau^2))
   omega <- -nu*rtau 
       c <- (.5*(w-1)*(w*cosh(2*omega)+1))^(-0.5)  
       z <- (y-(mu+c*sigma*w^(.5)*sinh(omega)))/(c*sigma)
       r <- -nu + asinh(z)/rtau
    dldm <- (z/(z*z+1)+(r/(rtau*(z*z+1)^(.5))))/(c*sigma)
 dlogcdt <- -rtau*w*((1/(w-1))+((cosh(2*omega))/(w*cosh(2*omega)+1)))
 dlogcdt <- dlogcdt+((nu*w*sinh(2*omega))/(w*cosh(2*omega)+1))
    dzdt <- -(z+w^(.5)*sinh(omega))*dlogcdt-rtau*w^(.5)*sinh(omega)+nu*w^(.5)*cosh(omega)  
    dldt <- -dlogcdt-(1/rtau)-(z/(z*z+1)+(r/(rtau*(z*z+1)^(.5))))*dzdt+(r*(r+nu))/rtau    
      dldt <- -dldt*rtau*rtau
 d2ldmdt <- -(dldm*dldt)
 d2ldmdt
                       },
  d2ldddv = function(y,mu,sigma,nu,tau) {
   rtau <- 1/tau
     w <- ifelse(rtau<0.0000001,1,exp(rtau^2))
  omega <- -nu*rtau 
      c <- (.5*(w-1)*(w*cosh(2*omega)+1))^(-0.5)  
      z <- (y-(mu+c*sigma*w^(.5)*sinh(omega)))/(c*sigma)
      r <- -nu + asinh(z)/rtau
   dldd <- (z+w^(0.5)*sinh(omega))*(z/(z*z+1)+(r/(rtau*(z*z+1)^(.5))))
   dldd <- (dldd-1)/sigma       
dlogcdv <- (rtau*w*sinh(2*omega))/(w*cosh(2*omega)+1)
   dzdv <- -(z+w^(.5)*sinh(omega))*dlogcdv+(rtau*w^(.5)*cosh(omega))
   dldv <- -dlogcdv-(z/(z*z+1)+(r/(rtau*(z*z+1)^(.5))))*dzdv+r         
d2ldddv <- -(dldd*dldv)
d2ldddv
                       },
  d2ldddt = function(y,mu,sigma,nu,tau) {
   rtau <- 1/tau
     w <- ifelse(rtau<0.0000001,1,exp(rtau^2))
  omega <- -nu*rtau 
      c <- (.5*(w-1)*(w*cosh(2*omega)+1))^(-0.5)  
      z <- (y-(mu+c*sigma*w^(.5)*sinh(omega)))/(c*sigma)
      r <- -nu + asinh(z)/rtau
   dldd <- (z+w^(0.5)*sinh(omega))*(z/(z*z+1)+(r/(rtau*(z*z+1)^(.5))))
   dldd <- (dldd-1)/sigma     
dlogcdt <- -rtau*w*((1/(w-1))+((cosh(2*omega))/(w*cosh(2*omega)+1)))
dlogcdt <- dlogcdt+((nu*w*sinh(2*omega))/(w*cosh(2*omega)+1))
   dzdt <- -(z+w^(.5)*sinh(omega))*dlogcdt-rtau*w^(.5)*sinh(omega)+nu*w^(.5)*cosh(omega)  
   dldt <- -dlogcdt-(1/rtau)-(z/(z*z+1)+(r/(rtau*(z*z+1)^(.5))))*dzdt+(r*(r+nu))/rtau    
      dldt <- -dldt*rtau*rtau
d2ldddt <- -(dldd*dldt)
d2ldddt  
                       },
  d2ldvdt = function(y,mu,sigma,nu,tau) {
   rtau <- 1/tau
     w <- ifelse(rtau<0.0000001,1,exp(rtau^2))
  omega <- -nu*rtau 
      c <- (.5*(w-1)*(w*cosh(2*omega)+1))^(-0.5)  
      z <- (y-(mu+c*sigma*w^(.5)*sinh(omega)))/(c*sigma)
      r <- -nu + asinh(z)/rtau
dlogcdv <- (rtau*w*sinh(2*omega))/(w*cosh(2*omega)+1)
   dzdv <- -(z+w^(.5)*sinh(omega))*dlogcdv+(rtau*w^(.5)*cosh(omega))
   dldv <- -dlogcdv-(z/(z*z+1)+(r/(rtau*(z*z+1)^(.5))))*dzdv+r
dlogcdt <- -rtau*w*((1/(w-1))+((cosh(2*omega))/(w*cosh(2*omega)+1)))
dlogcdt <- dlogcdt+((nu*w*sinh(2*omega))/(w*cosh(2*omega)+1))
   dzdt <- -(z+w^(.5)*sinh(omega))*dlogcdt-rtau*w^(.5)*sinh(omega)+nu*w^(.5)*cosh(omega)  
   dldt <- -dlogcdt-(1/rtau)-(z/(z*z+1)+(r/(rtau*(z*z+1)^(.5))))*dzdt+(r*(r+nu))/rtau          
      dldt <- -dldt*rtau*rtau
d2ldvdt <- -(dldv*dldt)
d2ldvdt   
                       },
 G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
                       { 
    -2*dJSU(y,mu,sigma,nu,tau,log=TRUE)
                        } ,                     
         rqres = expression(   
               rqres(pfun="pJSU", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)) ,
    mu.initial = expression(mu <- (y+mean(y))/2),  
 sigma.initial = expression(sigma<- rep(sd(y)/4, length(y))),
    nu.initial = expression(nu <- rep(0, length(y))), 
   tau.initial = expression(tau <-rep(1, length(y))), 
      mu.valid = function(mu) TRUE, 
   sigma.valid = function(sigma)  all(sigma > 0),
      nu.valid = function(nu) TRUE , 
     tau.valid = function(tau) all(tau > 0), 
       y.valid = function(y)  TRUE
          ),
            class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------
dJSU <- function(x, mu = 0, sigma = 1, nu = 1, tau = .5, log = FALSE)
 {
          if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
   rtau <- 1/tau
    if (length(tau)>1) 
           w <- ifelse(rtau<0.0000001,1,exp(rtau^2))       
    else w <- if (rtau<0.0000001) 1  else  exp(rtau^2)
  omega <- -nu*rtau 
      c <- (.5*(w-1)*(w*cosh(2*omega)+1))^(-0.5)  
      z <- (x-(mu+c*sigma*w^(.5)*sinh(omega)))/(c*sigma)
      r <- -nu + asinh(z)/rtau
 loglik <- -log(sigma)-log(c)-log(rtau)-.5*log(z*z+1)-.5*log(2*pi)-.5*r*r
       if(log==FALSE) ft  <- exp(loglik) else ft <- loglik 
       ft
  }    
#----------------------------------------------------------------------------------------  
pJSU <- function(q, mu = 0, sigma = 1, nu = 1, tau = .5, lower.tail = TRUE, log.p = FALSE)
 {  
         if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
         if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))           
   rtau <- 1/tau
      w <- ifelse(rtau<0.0000001,1,exp(rtau^2))
  omega <- -nu*rtau 
      c <- (.5*(w-1)*(w*cosh(2*omega)+1))^(-0.5)  
      z <- (q-(mu+c*sigma*w^(.5)*sinh(omega)))/(c*sigma)
      r <- -nu + asinh(z)/rtau    
      p <- pNO(r,0,1)
      if(lower.tail==TRUE) p  <- p else  p <- 1-p 
      if(log.p==FALSE) p  <- p else  p <- log(p) 
      p
 }
#----------------------------------------------------------------------------------------  
qJSU <-  function(p, mu=0, sigma=1, nu=0, tau=.5, lower.tail = TRUE, log.p = FALSE)
 {   
    if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
    if (lower.tail==TRUE) p <- p else p <- 1-p
   rtau <- 1/tau
    r <- qNO(p,0,1)
    z <- sinh(rtau*(r+nu))
    w <- ifelse(rtau<0.0000001,1,exp(rtau^2))
omega <- -nu*rtau 
    c <- (.5*(w-1)*(w*cosh(2*omega)+1))^(-0.5)     
    q <- (mu+c*sigma*w^(.5)*sinh(omega))+c*sigma*z   
    q
 }
#-----------------------------------------------------------------  
rJSU <- function(n, mu=0, sigma=1, nu=0, tau=.5)
  {
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    n <- ceiling(n)
    p <- runif(n)
    r <- qJSU(p,mu=mu,sigma=sigma,nu=nu,tau=tau)
    r
  }
#-----------------------------------------------------------------  
