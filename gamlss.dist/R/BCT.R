# last change : MS Wednesday, September 17, 2003 at 08:45
BCT <- function (mu.link="identity", sigma.link="log", nu.link ="identity", tau.link="log")
{
    mstats <- checklink("mu.link", "Box Cox t ", substitute(mu.link), 
                         c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Box Cox t ", substitute(sigma.link), 
                         c("inverse", "log", "identity", "own"))
    vstats <- checklink("nu.link", "Box Cox t ",substitute(nu.link), 
                         c("inverse", "log", "identity", "own"))
    tstats <- checklink("tau.link", "Box Cox t ",substitute(tau.link), 
                         c("inverse", "log", "identity", "own")) 
    structure(
          list(family = c("BCT",  "Box-Cox t"),
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
       z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
       w <- (tau+1)/(tau+z^2)
    dldm <- (w*z)/(mu*sigma)+(nu/mu)*(w*(z^2)-1)
    dldm
                                    },
               d2ldm2 = function(mu,sigma,nu,tau){
   d2ldm2 <- -(tau+2*nu*nu*sigma*sigma*tau+1)/(tau+3)
   d2ldm2 <- d2ldm2/(mu*mu*sigma*sigma)
                                   },
                 dldd = function(y,mu,sigma,nu,tau) {
        z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
        w <- (tau+1)/(tau+z^2)
        h <- dt(1/(sigma*abs(nu)),df=tau)/pt(1/(sigma*abs(nu)),df=tau)
     dldd <- (w*(z^2)-1)/sigma + h/(sigma^2*abs(nu))
     dldd
                                     } ,
               d2ldd2 = function(sigma, tau) -2*tau/(sigma^2*(tau+3)),
                 dldv = function(y,mu,sigma,nu,tau) {
        z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
        w <- (tau+1)/(tau+z^2)
        h <- dt(1/(sigma*abs(nu)),df=tau)/pt(1/(sigma*abs(nu)),df=tau)               
     dldv <- ((w*z^2)/nu)-log(y/mu) *(w*z^2+((w*z)/(sigma*nu))-1) 
     dldv <- dldv+sign(nu)*h/(sigma*nu^2) 
                                     } ,
               d2ldv2 = function(sigma) { 
                                    -7*(sigma^2)/4 
                                     },
                 dldt = function(y,mu,sigma,nu,tau) { 
         z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
         w <- (tau+1)/(tau+z^2)
         j <- (log(pt(1/(sigma*abs(nu)),df=tau+0.01))
               -log(pt(1/(sigma*abs(nu)),df=tau)))/0.01
      dldt <- -0.5*log(1+(z^2)/tau)+(w*(z^2))/(2*tau)
      dldt <- dldt+0.5*digamma((tau+1)/2)-0.5*digamma(tau/2)-1/(2*tau)-j
                                     } ,
               d2ldt2 = function(tau) {
    d2ldt2 <- trigamma((tau+1)/2) -trigamma(tau/2) +2*(tau+5)/(tau*(tau+1)*(tau+3))
    d2ldt2 <- d2ldt2/4
    d2ldt2 <- ifelse(d2ldt2 < -1e-15, d2ldt2,-1e-15)                                    
    d2ldt2
                                     } ,
              d2ldmdd = function(mu,sigma,nu,tau) -(2*nu*tau)/(mu*sigma*(tau+3)),
              d2ldmdv = function(mu,tau) (tau-3)/(2*mu*(tau+3)),
              d2ldmdt = function(mu,nu,tau) (2*nu)/(mu*(tau+1)*(tau+3)),
              d2ldddv = function(sigma,nu,tau) -(sigma*nu*tau)/(tau+3),
              d2ldddt = function(sigma,tau) 2/(sigma*(tau+1)*(tau+3)),
              d2ldvdt = function(sigma,nu,tau) (2*sigma^2*nu)/(tau^2), 
          G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
                                  -2*dBCT(y,mu,sigma,nu,tau,log=TRUE), 
                 rqres = expression(rqres(pfun="pBCT", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)) ,
            mu.initial = expression(mu <- (y+mean(y))/2), 
         sigma.initial = expression(sigma<- rep(0.1, length(y))),
            nu.initial = expression(nu <- rep(0.5, length(y))), 
           tau.initial = expression(tau <-rep(10, length(y))), 
              mu.valid = function(mu) all(mu > 0), 
           sigma.valid = function(sigma)  all(sigma > 0),
              nu.valid = function(nu) TRUE , 
             tau.valid = function(tau) all(tau > 0), 
               y.valid = function(y)  all(y > 0)
          ),
            class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------  
dBCT <- dBCTo <- function(x, mu=5, sigma=0.1, nu=1, tau=2, log=FALSE)
 {
          if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))  
          if (any(x < 0))  stop(paste("x must be positive", "\n", "")) 
          if(length(nu)>1)  z <- ifelse(nu != 0,(((x/mu)^nu-1)/(nu*sigma)),log(x/mu)/sigma)
          else   if (nu != 0) z <- (((x/mu)^nu-1)/(nu*sigma)) else z <- log(x/mu)/sigma
          loglik <- (nu-1)*log(x)-nu*log(mu)-log(sigma)
             fTz <-  lgamma((tau+1)/2)-lgamma(tau/2)-0.5*log(tau)-lgamma(0.5)
             fTz <- fTz-((tau+1)/2)* log(1+(z*z)/tau)
          loglik <- loglik+fTz-log(pt(1/(sigma*abs(nu)),df=tau))
          if (length(tau)>1) loglik <- ifelse(tau>1000000, dBCCG(x,mu,sigma,nu,log=TRUE), loglik) # MS Wednesday, April 12, 2006
          else if (tau>1000000) loglik <- dBCCG(x,mu,sigma,nu,log=TRUE)
             ft <- if(log==FALSE) exp(loglik) else loglik 
       ft
  }    
#-----------------------------------------------------------------  
pBCT <- pBCTo <- function(q, mu=5, sigma=0.1, nu=1, tau=2, lower.tail = TRUE, log.p = FALSE)
 {  
          if (any(mu < 0))  stop(paste("mu must be positive", "\n", "")) 
          if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
          if (any(q < 0))  stop(paste("y must be positive", "\n", ""))  
         if(length(nu)>1)  z <- ifelse(nu != 0,(((q/mu)^nu-1)/(nu*sigma)),log(q/mu)/sigma)
         else   if (nu != 0) z <- (((q/mu)^nu-1)/(nu*sigma)) else z <- log(q/mu)/sigma
         FYy1 <- pt(z,tau)
         if(length(nu)>1)  FYy2 <- ifelse(nu > 0, pt(-1/(sigma*abs(nu)),df=tau) ,0)
         else   if (nu>0)  FYy2 <-  pt(-1/(sigma*abs(nu)),df=tau) else FYy2 <- 0
         FYy3 <- pt(1/(sigma*abs(nu)),df=tau)
         FYy  <- (FYy1-FYy2)/FYy3
         if(lower.tail==TRUE) FYy  <- FYy else  FYy <- 1-FYy 
         if(log.p==FALSE) FYy  <- FYy else  FYy<- log(FYy) 
         FYy     
 }
#-----------------------------------------------------------------  
qBCT <- qBCTo <- function(p, mu=5, sigma=0.1, nu=1, tau=2, lower.tail = TRUE, log.p = FALSE)
 { 
    if (any(mu < 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
    if (lower.tail==TRUE) p <- p else p <- 1-p
    if(length(nu)>1) 
         {
         z <- ifelse((nu<=0),qt(p*pt(1/(sigma*abs(nu)),tau),tau),qt(1-(1-p)*pt(1/(sigma*abs(nu)),tau),tau)) 
         }
    else {
    z <- if (nu<=0) qt(p*pt(1/(sigma*abs(nu)),tau),tau)
         else       qt(1-(1-p)*pt(1/(sigma*abs(nu)),tau),tau)                     
         }                 
       if(length(nu)>1)  ya <- ifelse(nu != 0,mu*(nu*sigma*z+1)^(1/nu),mu*exp(sigma*z))
       else   if (nu != 0) ya <- mu*(nu*sigma*z+1)^(1/nu) else ya <- mu*exp(sigma*z)
       ya
 }
#-----------------------------------------------------------------  
rBCT <- rBCTo <- function(n, mu=5, sigma=0.1, nu=1, tau=2)
  {
    if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))  
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))    
    n <- ceiling(n)
    p <- runif(n)
    r <- qBCT(p,mu=mu,sigma=sigma,nu=nu,tau=tau)
    r
  }
#-----------------------------------------------------------------  
BCTuntr <- function (mu.link="identity", sigma.link="log", nu.link ="identity", tau.link="log")
{
    mstats <- checklink("mu.link", "Box Cox t ", substitute(mu.link), c("1/mu^2", "log", "identity"))
    dstats <- checklink("sigma.link", "Box Cox t ", substitute(sigma.link), c("inverse", "log", "identity"))
    vstats <- checklink("nu.link", "Box Cox t ",substitute(nu.link), c("1/nu^2", "log", "identity"))
    tstats <- checklink("tau.link", "Box Cox t ",substitute(tau.link), c("1/tau^2", "log", "identity"))
       
    structure(
          list(family = c("BCTuntr", "Box-Cox t untrucated"),
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
                        k <- 1/sigma
                        z <- ifelse(nu != 0,1,0)*(((y/mu)^nu-1)*(k/nu))
                        z <- z+ifelse(nu == 0,1,0)*log(y/mu)*k
                        u <- (tau+1)/(tau+z^2)
                     dldm <- k*u*(z/mu)+(nu/mu)*(u*(z*z)-1)
                     dldm
                                    },
               d2ldm2 = function(mu,sigma,nu,tau){
                   d2ldm2 <- -(tau+2*nu*nu*sigma*sigma*tau+1)/(tau+3)
                   d2ldm2 <- d2ldm2/(mu*mu*sigma*sigma)
                                   },
                 dldd = function(y,mu,sigma,nu,tau) {
                        k <- 1/sigma
                        z <- ifelse(nu != 0,1,0)*(((y/mu)^nu-1)*(k/nu))
                        z <- z+ifelse(nu == 0,1,0)*log(y/mu)*k
                        u <- (tau+1)/(tau+z^2)
                     dldd <- k*(u*(z*z)-1)
                     dldd
                                     } ,
               d2ldd2 = function(sigma,tau) -2*tau/(sigma^2*(tau+3)),
                 dldv = function(y,mu,sigma,nu,tau) {
                        k <- 1/sigma
                        z <- ifelse(nu != 0,1,0)*(((y/mu)^nu-1)*(k/nu))
                        z <- z+ifelse(nu == 0,1,0)*log(y/mu)*k
                        u <- (tau+1)/(tau+z^2)
                   lnybym <- log(y/mu)
                     dldv <- ((u*z)/nu)*(z-(lnybym*k))
                     dldv <- dldv-lnybym*(u*z*z-1) 
                                     } ,
               d2ldv2 = function(sigma) {
                          -7*(sigma^2)/4 
                                     },
                 dldt = function(y,mu,sigma,nu,tau) {
                       k <- 1/sigma
                       z <- ifelse(nu != 0,1,0)*(((y/mu)^nu-1)*(k/nu))
                       z <- z+ifelse(nu == 0,1,0)*log(y/mu)*k
                       u <- (tau+1)/(tau+z^2)
                      #t2 <-  tau/2
                      #t3 <- (tau+1)/2
                    dldt <- -0.5*log(1+(z*z)/tau)+0.5*u*(z*z)/tau
                    dldt <- dldt+0.5*digamma((tau+1)/2)-0.5*digamma(tau/2)-1/(2*tau)
                                     } ,
               d2ldt2 = function(tau) {
                      t2 <-  tau/2
                      t3 <- (tau+1)/2
                  d2ldt2 <- trigamma(t3) -trigamma(t2) +2*(tau+5)/(tau*(tau+1)*(tau+3))
                  d2ldt2 <- d2ldt2/4
                  d2ldt2 <- ifelse(d2ldt2 < -1e-15, d2ldt2,-1e-15)                                    
                         d2ldt2
                                     } ,
              d2ldmdd = function(mu,sigma,nu,tau) -(2*nu*tau)/(mu*sigma*(tau+3)),
              d2ldmdv = function(mu,tau) (tau-3)/(2*mu*(tau+3)),
              d2ldmdt = function(mu,nu,tau) (2*nu)/(mu*(tau+1)*(tau+3)),
              d2ldddv = function(sigma,nu,tau) -(sigma*nu*tau)/(tau+3),
              d2ldddt = function(sigma,tau) 2/(sigma*(tau+1)*(tau+3)),
              d2ldvdt = function(sigma,nu,tau) (2*sigma^2*nu)/(tau^2), 
          G.dev.incr  = function(y,mu,sigma,nu,tau,...) { 
                     tf <- pt(1/(sigma*abs(nu)),tau)
                     if(any(tf<.99)) 
                      warning(paste("The truncation factor F(1/(sigma*|nu|)) \n", 
                      "of BCT is less than 0.99 for some observations",
                            "\n"))       
                       k <- 1/sigma
                       z <- ifelse(nu != 0,1,0)*(((y/mu)^nu-1)*(k/nu))
                       z <- z + ifelse(nu == 0,1,0)*log(y/mu)*k
                     lik <- (nu*log(y/mu)-log(y)+log(k))
                     lik <-  lik+lgamma((tau+1)/2)-lgamma(tau/2)-0.5*log(tau)-lgamma(0.5)
                     lik <- lik-((tau+1)/2)* log(1+(z*z)/tau)
              G.dev.incr <- -2*lik
              G.dev.incr
                                                        } , 
                rqres =expression(   {
                       k <- 1/sigma
                       z <- ifelse(nu != 0,1,0)*(((y/mu)^nu-1)*(k/nu))
                       z <- z + ifelse(nu == 0,1,0)*log(y/mu)*k
                   rqres <- qnorm(pt(z,tau))
                                                           }) ,
            mu.initial = expression(mu <- y), 
         sigma.initial = expression(sigma<- rep(0.1, length(y))),
            nu.initial = expression(nu <- rep(0.5, length(y))), 
           tau.initial = expression(tau <-rep(10, length(y))), 
              mu.valid = function(mu) all(mu > 0), 
           sigma.valid = function(sigma)  all(sigma > 0),
              nu.valid = function(nu) TRUE , 
             tau.valid = function(tau) all(tau > 0), 
               y.valid = function(y)  all(y > 0)
          ),
            class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------------------------------
BCTo <- function (mu.link="log", sigma.link="log", nu.link ="identity", tau.link="log")
{
    mstats <- checklink("mu.link", "Box Cox t ", substitute(mu.link), 
                         c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Box Cox t ", substitute(sigma.link), 
                         c("inverse", "log", "identity", "own"))
    vstats <- checklink("nu.link", "Box Cox t ",substitute(nu.link), 
                         c("inverse", "log", "identity", "own"))
    tstats <- checklink("tau.link", "Box Cox t ",substitute(tau.link), 
                         c("inverse", "log", "identity", "own")) 
    structure(
          list(family = c("BCTo",  "Box-Cox-t-orig."),
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
       z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
       w <- (tau+1)/(tau+z^2)
    dldm <- (w*z)/(mu*sigma)+(nu/mu)*(w*(z^2)-1)
    dldm
                                    },
               d2ldm2 = function(mu,sigma,nu,tau){
   d2ldm2 <- -(tau+2*nu*nu*sigma*sigma*tau+1)/(tau+3)
   d2ldm2 <- d2ldm2/(mu*mu*sigma*sigma)
                                   },
                 dldd = function(y,mu,sigma,nu,tau) {
        z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
        w <- (tau+1)/(tau+z^2)
        h <- dt(1/(sigma*abs(nu)),df=tau)/pt(1/(sigma*abs(nu)),df=tau)
     dldd <- (w*(z^2)-1)/sigma + h/(sigma^2*abs(nu))
     dldd
                                     } ,
               d2ldd2 = function(sigma, tau) -2*tau/(sigma^2*(tau+3)),
                 dldv = function(y,mu,sigma,nu,tau) {
        z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
        w <- (tau+1)/(tau+z^2)
        h <- dt(1/(sigma*abs(nu)),df=tau)/pt(1/(sigma*abs(nu)),df=tau)               
     dldv <- ((w*z^2)/nu)-log(y/mu) *(w*z^2+((w*z)/(sigma*nu))-1) 
     dldv <- dldv+sign(nu)*h/(sigma*nu^2) 
                                     } ,
               d2ldv2 = function(sigma) { 
                                    -7*(sigma^2)/4 
                                     },
                 dldt = function(y,mu,sigma,nu,tau) { 
         z <- ifelse(nu != 0,(((y/mu)^nu-1)/(nu*sigma)),log(y/mu)/sigma)
         w <- (tau+1)/(tau+z^2)
         j <- (log(pt(1/(sigma*abs(nu)),df=tau+0.01))
               -log(pt(1/(sigma*abs(nu)),df=tau)))/0.01
      dldt <- -0.5*log(1+(z^2)/tau)+(w*(z^2))/(2*tau)
      dldt <- dldt+0.5*digamma((tau+1)/2)-0.5*digamma(tau/2)-1/(2*tau)-j
                                     } ,
               d2ldt2 = function(tau) {
    d2ldt2 <- trigamma((tau+1)/2) -trigamma(tau/2) +2*(tau+5)/(tau*(tau+1)*(tau+3))
    d2ldt2 <- d2ldt2/4
    d2ldt2 <- ifelse(d2ldt2 < -1e-15, d2ldt2,-1e-15)                                    
    d2ldt2
                                     } ,
              d2ldmdd = function(mu,sigma,nu,tau) -(2*nu*tau)/(mu*sigma*(tau+3)),
              d2ldmdv = function(mu,tau) (tau-3)/(2*mu*(tau+3)),
              d2ldmdt = function(mu,nu,tau) (2*nu)/(mu*(tau+1)*(tau+3)),
              d2ldddv = function(sigma,nu,tau) -(sigma*nu*tau)/(tau+3),
              d2ldddt = function(sigma,tau) 2/(sigma*(tau+1)*(tau+3)),
              d2ldvdt = function(sigma,nu,tau) (2*sigma^2*nu)/(tau^2), 
          G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
                                  -2*dBCTo(y,mu,sigma,nu,tau,log=TRUE), 
                 rqres = expression(rqres(pfun="pBCTo", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)) ,
            mu.initial = expression(mu <- (y+mean(y))/2), 
         sigma.initial = expression(sigma<- rep(0.1, length(y))),
            nu.initial = expression(nu <- rep(0.5, length(y))), 
           tau.initial = expression(tau <-rep(10, length(y))), 
              mu.valid = function(mu) all(mu > 0), 
           sigma.valid = function(sigma)  all(sigma > 0),
              nu.valid = function(nu) TRUE , 
             tau.valid = function(tau) all(tau > 0), 
               y.valid = function(y)  all(y > 0)
          ),
            class = c("gamlss.family","family"))
}
#------------------------------------------------------------------------
