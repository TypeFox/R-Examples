# updated 8_11_2007 fits very well both RS and mixed 
ST4 <- function (mu.link="identity", sigma.link="log", nu.link ="log", tau.link="log")
{
    mstats <- checklink("mu.link", "skew t type 4", substitute(mu.link), 
                         c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "skew t type 4", substitute(sigma.link), 
                         c("inverse", "log", "identity", "own"))
    vstats <- checklink("nu.link", "skew t type 4",substitute(nu.link), 
                         c("inverse", "log", "identity", "own"))
    tstats <- checklink("tau.link", "skew t type 4 ",substitute(tau.link), 
                         c("inverse", "log", "identity", "own")) 
    structure(
          list(family = c("ST4",  "skew t type 4"),
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
        s <- sigma
     dsq  <- ((y-mu)/s)^2
    w1 <- ifelse(nu < 1000000, (nu+1)/(nu+dsq),1)
    w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq),1)
     dldm <- ifelse(y < mu, (w1*(y-mu))/(s^2) , (w2*(y-mu))/(s^2))
     dldm
                                    },
               d2ldm2 = function(y,mu,sigma,nu,tau){
        s <- sigma
     dsq  <- ((y-mu)/s)^2
    w1 <- ifelse(nu < 1000000, (nu+1)/(nu+dsq),1)
    w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq),1)
     dldm <- ifelse(y < mu, (w1*(y-mu))/(s^2) , (w2*(y-mu))/(s^2))
    d2ldm2 <- -dldm*dldm
   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)  
   d2ldm2
                                   },
                 dldd = function(y,mu,sigma,nu,tau) {
        s <- sigma
     dsq  <- ((y-mu)/s)^2
    w1 <- ifelse(nu < 1000000, (nu+1)/(nu+dsq),1)
    w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq),1)
     dldd <- ifelse(y < mu, (w1*dsq-1)/(sigma) , (w2*dsq-1)/(sigma) )
     dldd
                                     } ,
               d2ldd2 = function(y,mu,sigma,nu,tau) {
        s <- sigma
     dsq  <- ((y-mu)/s)^2
    w1 <- ifelse(nu < 1000000, (nu+1)/(nu+dsq),1)
    w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq),1)
     dldd <- ifelse(y < mu, (w1*dsq-1)/(sigma) , (w2*dsq-1)/(sigma) )
    d2ldd2 <- -dldd*dldd
   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
   d2ldd2
                                     },
                 dldv = function(y,mu,sigma,nu,tau) {
        s <- sigma
     dsq  <- ((y-mu)/s)^2
    w1 <- ifelse(nu < 1000000, (nu+1)/(nu+dsq),1)
      lk1 <- lgamma(0.5) + lgamma(nu/2) - lgamma((nu+1)/2) + 0.5*log(nu)
    lk1 <- ifelse(nu < 1000000, lk1, 0.5*log(2*pi))
      lk2 <- lgamma(0.5) + lgamma(tau/2) - lgamma((tau+1)/2) + 0.5*log(tau)
    lk2 <- ifelse(tau < 1000000, lk2, 0.5*log(2*pi))
       lk <- lk2 - lk1
        k <- exp(lk)       
    dldva <- -0.5*log(1+dsq/nu)+(w1*dsq)/(2*nu)
     dldv <- ifelse(y < mu, dldva , 0)
     dldv <- dldv+0.5*(digamma((nu+1)/2)-digamma(nu/2)-(1/nu))*(1/(1+k))
     dldv
                                     } ,
               d2ldv2 = function(y,mu,sigma,nu,tau) { 
        s <- sigma
     dsq  <- ((y-mu)/s)^2
    w1 <- ifelse(nu < 1000000, (nu+1)/(nu+dsq),1)
      lk1 <- lgamma(0.5) + lgamma(nu/2) - lgamma((nu+1)/2) + 0.5*log(nu)
    lk1 <- ifelse(nu < 1000000, lk1, 0.5*log(2*pi))
      lk2 <- lgamma(0.5) + lgamma(tau/2) - lgamma((tau+1)/2) + 0.5*log(tau)
    lk2 <- ifelse(tau < 1000000, lk2, 0.5*log(2*pi))
       lk <- lk2 - lk1
        k <- exp(lk)       
    dldva <- -0.5*log(1+dsq/nu)+(w1*dsq)/(2*nu)
     dldv <- ifelse(y < mu, dldva , 0)
     dldv <- dldv+0.5*(digamma((nu+1)/2)-digamma(nu/2)-(1/nu))*(1/(1+k))
    d2ldv2 <- -dldv*dldv
   d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)  
   d2ldv2
                                     },
                 dldt = function(y,mu,sigma,nu,tau) { 
        s <- sigma
     dsq  <- ((y-mu)/s)^2
    w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq),1)
      lk1 <- lgamma(0.5) + lgamma(nu/2) - lgamma((nu+1)/2) + 0.5*log(nu)
    lk1 <- ifelse(nu < 1000000, lk1, 0.5*log(2*pi))
      lk2 <- lgamma(0.5) + lgamma(tau/2) - lgamma((tau+1)/2) + 0.5*log(tau)
    lk2 <- ifelse(tau < 1000000, lk2, 0.5*log(2*pi))
       lk <- lk2 - lk1
        k <- exp(lk)       
    dldtb <- -0.5*log(1+dsq/tau)+(w2*dsq)/(2*tau)
     dldt <- ifelse(y < mu, 0 , dldtb)
     dldt <- dldt+0.5*(digamma((tau+1)/2)-digamma(tau/2)-(1/tau))*(k/(1+k))
     dldt
                                     } ,
               d2ldt2 = function(y,mu,sigma,nu,tau) {
        s <- sigma
     dsq  <- ((y-mu)/s)^2
    w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq),1)
      lk1 <- lgamma(0.5) + lgamma(nu/2) - lgamma((nu+1)/2) + 0.5*log(nu)
    lk1 <- ifelse(nu < 1000000, lk1, 0.5*log(2*pi))
      lk2 <- lgamma(0.5) + lgamma(tau/2) - lgamma((tau+1)/2) + 0.5*log(tau)
    lk2 <- ifelse(tau < 1000000, lk2, 0.5*log(2*pi))
       lk <- lk2 - lk1
        k <- exp(lk)       
    dldtb <- -0.5*log(1+dsq/tau)+(w2*dsq)/(2*tau)
     dldt <- ifelse(y < mu, 0 , dldtb)
     dldt <- dldt+0.5*(digamma((tau+1)/2)-digamma(tau/2)-(1/tau))*(k/(1+k))
    d2ldt2 <- -dldt*dldt
    d2ldt2 <- ifelse(d2ldt2 < -1e-15, d2ldt2,-1e-15)                                    
    d2ldt2
                                     } ,
              d2ldmdd = function(y,mu,sigma,nu,tau)  
                     {
        s <- sigma
     dsq  <- ((y-mu)/s)^2
    w1 <- ifelse(nu < 1000000, (nu+1)/(nu+dsq),1)
    w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq),1)
     dldm <- ifelse(y < mu, (w1*(y-mu))/(s^2) , (w2*(y-mu))/(s^2))
     dldd <- ifelse(y < mu, (w1*dsq-1)/(sigma) , (w2*dsq-1)/(sigma) )
  d2ldmdd <- -(dldm*dldd) 
  d2ldmdd                 
                      },
              d2ldmdv = function(y,mu,sigma,nu,tau)
                      {
        s <- sigma
     dsq  <- ((y-mu)/s)^2
    w1 <- ifelse(nu < 1000000, (nu+1)/(nu+dsq),1)
    w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq),1)
     dldm <- ifelse(y < mu, (w1*(y-mu))/(s^2) , (w2*(y-mu))/(s^2))
      lk1 <- lgamma(0.5) + lgamma(nu/2) - lgamma((nu+1)/2) + 0.5*log(nu)
    lk1 <- ifelse(nu < 1000000, lk1, 0.5*log(2*pi))
      lk2 <- lgamma(0.5) + lgamma(tau/2) - lgamma((tau+1)/2) + 0.5*log(tau)
    lk2 <- ifelse(tau < 1000000, lk2, 0.5*log(2*pi))
       lk <- lk2 - lk1
        k <- exp(lk)       
    dldva <- -0.5*log(1+dsq/nu)+(w1*dsq)/(2*nu)
     dldv <- ifelse(y < mu, dldva , 0)
     dldv <- dldv+0.5*(digamma((nu+1)/2)-digamma(nu/2)-(1/nu))*(k/(1+k))
  d2ldmdv <- -(dldm*dldv)
  d2ldmdv 
                      } ,
              d2ldmdt = function(y,mu,sigma,nu,tau)
                     {  
        s <- sigma
     dsq  <- ((y-mu)/s)^2
    w1 <- ifelse(nu < 1000000, (nu+1)/(nu+dsq),1)
    w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq),1)
     dldm <- ifelse(y < mu, (w1*(y-mu))/(s^2) , (w2*(y-mu))/(s^2))
      lk1 <- lgamma(0.5) + lgamma(nu/2) - lgamma((nu+1)/2) + 0.5*log(nu)
    lk1 <- ifelse(nu < 1000000, lk1, 0.5*log(2*pi))
      lk2 <- lgamma(0.5) + lgamma(tau/2) - lgamma((tau+1)/2) + 0.5*log(tau)
    lk2 <- ifelse(tau < 1000000, lk2, 0.5*log(2*pi))
       lk <- lk2 - lk1
        k <- exp(lk)       
    dldtb <- -0.5*log(1+dsq/tau)+(w2*dsq)/(2*tau)
     dldt <- ifelse(y < mu, 0 , dldtb)
     dldt <- dldt+0.5*(digamma((tau+1)/2)-digamma(tau/2)-(1/tau))*(k/(1+k))
  d2ldmdt <- -(dldm*dldt)
  d2ldmdt
                      },
              d2ldddv = function(y,mu,sigma,nu,tau)  
                      {
        s <- sigma
     dsq  <- ((y-mu)/s)^2
    w1 <- ifelse(nu < 1000000, (nu+1)/(nu+dsq),1)
    w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq),1)
     dldd <- ifelse(y < mu, (w1*dsq-1)/(sigma) , (w2*dsq-1)/(sigma) )
      lk1 <- lgamma(0.5) + lgamma(nu/2) - lgamma((nu+1)/2) + 0.5*log(nu)
    lk1 <- ifelse(nu < 1000000, lk1, 0.5*log(2*pi))
      lk2 <- lgamma(0.5) + lgamma(tau/2) - lgamma((tau+1)/2) + 0.5*log(tau)
    lk2 <- ifelse(tau < 1000000, lk2, 0.5*log(2*pi))
       lk <- lk2 - lk1
        k <- exp(lk)       
    dldva <- -0.5*log(1+dsq/nu)+(w1*dsq)/(2*nu)
     dldv <- ifelse(y < mu, dldva , 0)
     dldv <- dldv+0.5*(digamma((nu+1)/2)-digamma(nu/2)-(1/nu))*(k/(1+k))
   d2ldmdt<- -(dldd*dldv)
   d2ldmdt
                      },
              d2ldddt = function(y,mu,sigma,nu,tau)  
                      {
        s <- sigma
     dsq  <- ((y-mu)/s)^2
    w1 <- ifelse(nu < 1000000, (nu+1)/(nu+dsq),1)
    w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq),1)
     dldd <- ifelse(y < mu, (w1*dsq-1)/(sigma) , (w2*dsq-1)/(sigma) )
      lk1 <- lgamma(0.5) + lgamma(nu/2) - lgamma((nu+1)/2) + 0.5*log(nu)
    lk1 <- ifelse(nu < 1000000, lk1, 0.5*log(2*pi))
      lk2 <- lgamma(0.5) + lgamma(tau/2) - lgamma((tau+1)/2) + 0.5*log(tau)
    lk2 <- ifelse(tau < 1000000, lk2, 0.5*log(2*pi))
       lk <- lk2 - lk1
        k <- exp(lk)       
    dldtb <- -0.5*log(1+dsq/tau)+(w2*dsq)/(2*tau)
     dldt <- ifelse(y < mu, 0 , dldtb)
     dldt <- dldt+0.5*(digamma((tau+1)/2)-digamma(tau/2)-(1/tau))*(k/(1+k))
  d2ldddt <- -(dldd*dldt) 
  d2ldddt                  
                      },
              d2ldvdt = function(y,mu,sigma,nu,tau)  
                      {
        s <- sigma
     dsq  <- ((y-mu)/s)^2
    w1 <- ifelse(nu < 1000000, (nu+1)/(nu+dsq),1)
    w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq),1)
      lk1 <- lgamma(0.5) + lgamma(nu/2) - lgamma((nu+1)/2) + 0.5*log(nu)
    lk1 <- ifelse(nu < 1000000, lk1, 0.5*log(2*pi))
      lk2 <- lgamma(0.5) + lgamma(tau/2) - lgamma((tau+1)/2) + 0.5*log(tau)
    lk2 <- ifelse(tau < 1000000, lk2, 0.5*log(2*pi))
       lk <- lk2 - lk1
        k <- exp(lk)       
    dldva <- -0.5*log(1+dsq/nu)+(w1*dsq)/(2*nu)
     dldv <- ifelse(y < mu, dldva , 0)
     dldv <- dldv+0.5*(digamma((nu+1)/2)-digamma(nu/2)-(1/nu))*(k/(1+k))
    dldtb <- -0.5*log(1+dsq/tau)+(w2*dsq)/(2*tau)
     dldt <- ifelse(y < mu, 0 , dldtb)
     dldt <- dldt+0.5*(digamma((tau+1)/2)-digamma(tau/2)-(1/tau))*(k/(1+k))
  d2ldvdt <- -(dldv*dldt)
  d2ldvdt    
                      }, 
          G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
                                  -2*dST4(y,mu,sigma,nu,tau,log=TRUE), 
                 rqres = expression(rqres(pfun="pST4", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)) ,
            mu.initial = expression(mu <- (y+mean(y))/2), 
         sigma.initial = expression(sigma<- rep(sd(y)/4, length(y))),
            nu.initial = expression(nu <- rep(5, length(y))), 
           tau.initial = expression(tau <-rep(5, length(y))), 
              mu.valid = function(mu) TRUE, 
           sigma.valid = function(sigma)  all(sigma > 0),
              nu.valid = function(nu) all(nu > 0), 
             tau.valid = function(tau) all(tau > 0), 
               y.valid = function(y)  TRUE
          ),
            class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------  
dST4 <- function(x, mu=0, sigma=1, nu=1, tau=10, log=FALSE)
 {
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
          if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))  
      lk1 <- lgamma(0.5) + lgamma(nu/2) - lgamma((nu+1)/2) + 0.5*log(nu)
    lk1 <- ifelse(nu < 1000000, lk1, 0.5*log(2*pi))
      lk2 <- lgamma(0.5) + lgamma(tau/2) - lgamma((tau+1)/2) + 0.5*log(tau)
    lk2 <- ifelse(tau < 1000000, lk2, 0.5*log(2*pi))
          lk <- lk2 - lk1
          k <- exp(lk)       
           loglikb <- dNO(((x-mu)/sigma), mu=0, sigma=1, log=TRUE)
         loglik1 <- dt((x-mu)/sigma, df=nu, log=TRUE)
     if (length(nu)>1) loglik1 <- ifelse(nu<1000000, 
                      loglik1, 
                      loglikb) 
    else loglik1 <- if (tau<1000000) loglik1
                else  loglikb

          loglik2 <- dt((x-mu)/sigma, df=tau, log=TRUE)
     if (length(tau)>1) loglik2 <- ifelse(tau<1000000, 
                      loglik2, 
                      loglikb) 
    else loglik2 <- if (tau<1000000) loglik2
                else  loglikb

          loglik <- ifelse(x < mu, loglik1, loglik2 + lk)
          loglik <- loglik + log(2/(1+k)) - log(sigma)

       if(log==FALSE) ft  <- exp(loglik) else ft <- loglik 
       ft
  }    
#-----------------------------------------------------------------  
pST4 <- function(q, mu=0, sigma=1, nu=1, tau=10, lower.tail = TRUE, log.p = FALSE)
 {  
         if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
         if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
         if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
      lk1 <- lgamma(0.5) + lgamma(nu/2) - lgamma((nu+1)/2) + 0.5*log(nu)
    lk1 <- ifelse(nu < 1000000, lk1, 0.5*log(2*pi))
      lk2 <- lgamma(0.5) + lgamma(tau/2) - lgamma((tau+1)/2) + 0.5*log(tau)
    lk2 <- ifelse(tau < 1000000, lk2, 0.5*log(2*pi))
         lk <- lk2 - lk1
         k <- exp(lk)       
         cdf1 <- 2*pt((q-mu)/sigma, df=nu)
         cdf2 <- 1 + 2*k*(pt( (q-mu)/sigma, df=tau) - 0.5)
         cdf <- ifelse(q < mu, cdf1, cdf2)
         cdf <- cdf/(1+k)
         if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf 
         if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf) 
         cdf    
 }
#-----------------------------------------------------------------  
qST4 <- function(p, mu=0, sigma=1, nu=1, tau=10, lower.tail = TRUE, log.p = FALSE)
 { 
    if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
    if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
    if (lower.tail==TRUE) p <- p else p <- 1-p
      lk1 <- lgamma(0.5) + lgamma(nu/2) - lgamma((nu+1)/2) + 0.5*log(nu)
    lk1 <- ifelse(nu < 1000000, lk1, 0.5*log(2*pi))
      lk2 <- lgamma(0.5) + lgamma(tau/2) - lgamma((tau+1)/2) + 0.5*log(tau)
    lk2 <- ifelse(tau < 1000000, lk2, 0.5*log(2*pi))
    lk <- lk2 - lk1
    k <- exp(lk)       
    suppressWarnings(q1 <- mu + sigma*qt(p*(1+k)/2, df=nu))
    suppressWarnings(q2 <- mu + sigma*qt((p*(1+k)-1)/(2*k) + 0.5, df=tau))
    q <- ifelse(p < (1/(1+k)), q1, q2)
    q
 }
#-----------------------------------------------------------------  
rST4 <- function(n, mu=0, sigma=1, nu=1, tau=10)
  {
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
    if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))  
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))    
    n <- ceiling(n)
    p <- runif(n)
    r <- qST4(p,mu=mu,sigma=sigma,nu=nu,tau=tau)
    r
  }
