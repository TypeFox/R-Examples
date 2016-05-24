#   28_11_2007
SEP3 <- function (mu.link="identity", sigma.link="log", nu.link ="log", tau.link="log")
{
    mstats <- checklink("mu.link", "skew exponential power type 2", substitute(mu.link), 
                         c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "skew exponential power type 2", substitute(sigma.link), 
                         c("inverse", "log", "identity", "own"))
    vstats <- checklink("nu.link", "skew exponential power type 2",substitute(nu.link), 
                         c("inverse", "log", "identity", "own"))
    tstats <- checklink("tau.link", "skew exponential power type 2 ",substitute(tau.link), 
                         c("inverse", "log", "identity", "own")) 
    structure(
          list(family = c("SEP3",  "skew exponential power type 3"),
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
    dldma <- sign(z)*(nu*tau/(2*sigma))*((nu*abs(z))^(tau-1))
    dldmb <- sign(z)*(tau/(2*sigma*nu))*((abs(z)/nu)^(tau-1))
     dldm <- ifelse(y < mu, dldma , dldmb)
     dldm
                                    },
               d2ldm2 = function(y,mu,sigma,nu,tau){
        z <- (y-mu)/sigma
    dldma <- sign(z)*(nu*tau/(2*sigma))*((nu*abs(z))^(tau-1))
    dldmb <- sign(z)*(tau/(2*sigma*nu))*((abs(z)/nu)^(tau-1))
     dldm <- ifelse(y < mu, dldma , dldmb)
    d2ldm2 <- -dldm*dldm
   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)  
   d2ldm2
                                   },
                 dldd = function(y,mu,sigma,nu,tau) {
        z <- (y-mu)/sigma
    dldda <- (nu*abs(z))^tau
    dlddb <- (abs(z)/nu)^tau
     dldd <- ifelse(y < mu, dldda , dlddb )
     dldd <- dldd*(tau/(2*sigma))-(1/sigma)
     dldd
                                     } ,
               d2ldd2 = function(y,mu,sigma,nu,tau) {
        z <- (y-mu)/sigma
    dldda <- (nu*abs(z))^tau
    dlddb <- (abs(z)/nu)^tau
     dldd <- ifelse(y < mu, dldda , dlddb )
     dldd <- dldd*(tau/(2*sigma))-(1/sigma)
    d2ldd2 <- -dldd*dldd
   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
   d2ldd2
                                     },
                 dldv = function(y,mu,sigma,nu,tau) {
        z <- (y-mu)/sigma  
    dldva <- sign(z)*((nu*abs(z))^tau)
    dldvb <- sign(z)*((abs(z)/nu)^tau)
     dldv <- ifelse(y < mu, dldva , dldvb )
     dldv <- dldv*(tau/(2*nu)) + (1/nu) - (2*nu)/(1+(nu^2))
     dldv
                                     } ,
               d2ldv2 = function(y,mu,sigma,nu,tau) { 
        z <- (y-mu)/sigma  
    dldva <- sign(z)*((nu*abs(z))^tau)
    dldvb <- sign(z)*((abs(z)/nu)^tau)
     dldv <- ifelse(y < mu, dldva , dldvb )
     dldv <- dldv*(tau/(2*nu)) + (1/nu) - (2*nu)/(1+(nu^2))
    d2ldv2 <- -dldv*dldv
   d2ldv2 <- ifelse(d2ldv2 < -1e-4, d2ldv2,-1e-4)  
   d2ldv2
                                     },
                 dldt = function(y,mu,sigma,nu,tau) { 
        z <- (y-mu)/sigma  
    dldta <- -(log(nu*abs(z)))*((nu*abs(z))^tau) 
    dldtb <- -(log(abs(z)/nu))*((abs(z)/nu)^tau) 
     dldt <- ifelse(y < mu, dldta , dldtb)
     dldt <- dldt/2 + ((digamma(1+(1/tau)))/(tau^2)) + log(2)/(tau^2)
     dldt
                                     } ,
               d2ldt2 = function(y,mu,sigma,nu,tau) {
        z <- (y-mu)/sigma  
    dldta <- -(log(nu*abs(z)))*((nu*abs(z))^tau) 
    dldtb <- -(log(abs(z)/nu))*((abs(z)/nu)^tau) 
     dldt <- ifelse(y < mu, dldta , dldtb)
     dldt <- dldt/2 + ((digamma(1+(1/tau)))/(tau^2)) + log(2)/(tau^2)
    d2ldt2 <- -dldt*dldt
    d2ldt2 <- ifelse(d2ldt2 < -1e-4, d2ldt2,-1e-4)                                    
    d2ldt2
                                     } ,


              d2ldmdd = function(y,mu,sigma,nu,tau)
                      {  
        z <- (y-mu)/sigma
    dldma <- sign(z)*(nu*tau/(2*sigma))*((nu*abs(z))^(tau-1))
    dldmb <- sign(z)*(tau/(2*sigma*nu))*((abs(z)/nu)^(tau-1))
     dldm <- ifelse(y < mu, dldma , dldmb)
    dldda <- (nu*abs(z))^tau
    dlddb <- (abs(z)/nu)^tau
     dldd <- ifelse(y < mu, dldda , dlddb )
     dldd <- dldd*(tau/(2*sigma))-(1/sigma)
  d2ldmdd <- -(dldm*dldd)
  d2ldmdd                  
                      },
              d2ldmdv = function(y,mu,sigma,nu,tau)
                      {
        z <- (y-mu)/sigma
    dldma <- sign(z)*(nu*tau/(2*sigma))*((nu*abs(z))^(tau-1))
    dldmb <- sign(z)*(tau/(2*sigma*nu))*((abs(z)/nu)^(tau-1))
     dldm <- ifelse(y < mu, dldma , dldmb)
    dldva <- sign(z)*((nu*abs(z))^tau)
    dldvb <- sign(z)*((abs(z)/nu)^tau)
     dldv <- ifelse(y < mu, dldva , dldvb )
     dldv <- dldv*(tau/(2*nu)) + (1/nu) - (2*nu)/(1+(nu^2))
  d2ldmdv <- -(dldm*dldv)
  d2ldmdv
                      } ,
              d2ldmdt = function(y,mu,sigma,nu,tau)
                     {  
        z <- (y-mu)/sigma
    dldma <- sign(z)*(nu*tau/(2*sigma))*((nu*abs(z))^(tau-1))
    dldmb <- sign(z)*(tau/(2*sigma*nu))*((abs(z)/nu)^(tau-1))
     dldm <- ifelse(y < mu, dldma , dldmb)
    dldta <- -(log(nu*abs(z)))*((nu*abs(z))^tau) 
    dldtb <- -(log(abs(z)/nu))*((abs(z)/nu)^tau) 
     dldt <- ifelse(y < mu, dldta , dldtb)
     dldt <- dldt/2 + ((digamma(1+(1/tau)))/(tau^2)) + log(2)/(tau^2)
  d2ldmdt <- -(dldm*dldt)
  d2ldmdt
                      },
              d2ldddv = function(y,mu,sigma,nu,tau)  
                      {
        z <- (y-mu)/sigma
    dldda <- (nu*abs(z))^tau
    dlddb <- (abs(z)/nu)^tau
     dldd <- ifelse(y < mu, dldda , dlddb )
     dldd <- dldd*(tau/(2*sigma))-(1/sigma)
    dldva <- sign(z)*((nu*abs(z))^tau)
    dldvb <- sign(z)*((abs(z)/nu)^tau)
     dldv <- ifelse(y < mu, dldva , dldvb )
     dldv <- dldv*(tau/(2*nu)) + (1/nu) - (2*nu)/(1+(nu^2))
  d2ldddv <- -(dldd*dldv)
  d2ldddv
                      },
              d2ldddt = function(y,mu,sigma,nu,tau)  
                      {
        z <- (y-mu)/sigma
    dldda <- (nu*abs(z))^tau
    dlddb <- (abs(z)/nu)^tau
     dldd <- ifelse(y < mu, dldda , dlddb )
     dldd <- dldd*(tau/(2*sigma))-(1/sigma)
    dldta <- -(log(nu*abs(z)))*((nu*abs(z))^tau) 
    dldtb <- -(log(abs(z)/nu))*((abs(z)/nu)^tau) 
     dldt <- ifelse(y < mu, dldta , dldtb)
     dldt <- dldt/2 + ((digamma(1+(1/tau)))/(tau^2)) + log(2)/(tau^2)
  d2ldddt <- -(dldd*dldt)
  d2ldddt                   
                      },
              d2ldvdt = function(y,mu,sigma,nu,tau)  
                      {
        z <- (y-mu)/sigma  
    dldva <- sign(z)*((nu*abs(z))^tau)
    dldvb <- sign(z)*((abs(z)/nu)^tau)
     dldv <- ifelse(y < mu, dldva , dldvb )
     dldv <- dldv*(tau/(2*nu)) + (1/nu) - (2*nu)/(1+(nu^2))
    dldta <- -(log(nu*abs(z)))*((nu*abs(z))^tau) 
    dldtb <- -(log(abs(z)/nu))*((abs(z)/nu)^tau) 
     dldt <- ifelse(y < mu, dldta , dldtb)
     dldt <- dldt/2 + ((digamma(1+(1/tau)))/(tau^2)) + log(2)/(tau^2)
  d2ldvdt <- -(dldv*dldt)
  d2ldvdt    
                      }, 
          G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
                                  -2*dSEP3(y,mu,sigma,nu,tau,log=TRUE), 
                 rqres = expression(rqres(pfun="pSEP3", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)) ,
            mu.initial = expression(mu <- (y+mean(y))/2), 
         sigma.initial = expression(sigma<- rep(sd(y), length(y))),
            nu.initial = expression(nu <- rep(1, length(y))), 
           tau.initial = expression(tau <-rep(2, length(y))), 
              mu.valid = function(mu) TRUE, 
           sigma.valid = function(sigma)  all(sigma > 0),
              nu.valid = function(nu) all(nu > 0), 
             tau.valid = function(tau) all(tau > 0), 
               y.valid = function(y)  TRUE
          ),
            class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------  
dSEP3 <- function(x, mu=0, sigma=1, nu=2, tau=2, log=FALSE)
 {
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
          if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))    
          z <- (x-mu)/sigma
           suppressWarnings(loglik1 <- -0.5*((nu*abs(z))^tau))
           suppressWarnings(loglik2 <- -0.5*((abs(z)/nu)^tau))
#    if (length(mu)>1) loglik <- ifelse(x < mu, loglik1, loglik2)
#    else loglik <- if (x < mu)   loglik1  else  loglik2 
          loglik <- ifelse(x < mu, loglik1, loglik2)
          loglik <- loglik-log(sigma)+log(nu)-log(1+(nu^2))-(1/tau)*log(2)-lgamma(1+(1/tau))
          fy <- if(log==FALSE) exp(loglik) else loglik 
       fy
  }    
#-----------------------------------------------------------------  
pSEP3 <- function(q, mu=0, sigma=1, nu=2, tau=2, lower.tail = TRUE, log.p = FALSE)
 {  
         if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
         if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
         if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
        k <- nu^2       
       z1 <- nu*(q-mu)/(sigma*(2^(1/tau)))
       z2 <- (q-mu)/(sigma*nu*(2^(1/tau)))
       s1 <- (abs(z1)^tau)
       s2 <- (abs(z2)^tau)
     cdf1 <- 1-pgamma(s1,shape=1/tau,scale=1)
     cdf2 <- 1+k*pgamma(s2,shape=1/tau,scale=1)
#    if (length(mu)>1) cdf <- ifelse(q < mu, cdf1, cdf2)
#    else cdf <- if (q < mu)   cdf1  else  cdf2
      cdf <- ifelse(q < mu, cdf1, cdf2)
      cdf <- cdf/(1+k)
      
    if (length(tau)>1) cdf <- ifelse(tau>10000, 
                                   (q-(mu-(sigma/nu)))/(sigma*((1/nu)+nu)),
                                    cdf)
    else cdf <- if (tau>10000)  (q-(mu-(sigma/nu)))/(sigma*((1/nu)+nu))  else  cdf
      
      if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf  
      if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf) 
      cdf    
 }
#-----------------------------------------------------------------  
qSEP3 <- function(p, mu=0, sigma=1, nu=2, tau=2, lower.tail = TRUE, log.p = FALSE)
 { 
    if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
    if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
    if (lower.tail==TRUE) p <- p else p <- 1-p
    k <- nu^2       
    suppressWarnings(q1 <- mu - 
    (sigma*(2^(1/tau))/nu)*((qgamma( 1-p*(1+k), shape=1/tau, scale=1))^(1/tau)))
    suppressWarnings(q2 <- mu + 
    (sigma*nu*(2^(1/tau)))*((qgamma( (-1/k)*(1-p*(1+k)), shape=1/tau, scale=1))^(1/tau)))
    q <- ifelse(p < (1/(1+k)), q1, q2)
    q
 }
#-----------------------------------------------------------------  
rSEP3 <- function(n, mu=0, sigma=1, nu=2, tau=2)
  {
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
    if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))  
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))    
    n <- ceiling(n)
    p <- runif(n)
    r <- qSEP3(p,mu=mu,sigma=sigma,nu=nu,tau=tau)
    r
  }
