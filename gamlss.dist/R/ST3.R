# updated 27_11_2007
ST3 <- function (mu.link="identity", sigma.link="log", nu.link ="log", tau.link="log")
{
    mstats <- checklink("mu.link", "skew t type 2", substitute(mu.link), 
                         c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "skew t type 2", substitute(sigma.link), 
                         c("inverse", "log", "identity", "own"))
    vstats <- checklink("nu.link", "skew t type 2",substitute(nu.link), 
                         c("inverse", "log", "identity", "own"))
    tstats <- checklink("tau.link", "skew t type 2 ",substitute(tau.link), 
                         c("inverse", "log", "identity", "own")) 
    structure(
          list(family = c("ST3",  "skew t type 3"),
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
       s1 <- sigma/nu
       s2 <- sigma*nu
     dsq1 <- ((y-mu)/s1)^2
     dsq2 <- ((y-mu)/s2)^2
       w1 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq1),1)
       w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq2),1)
     dldm <- ifelse(y < mu, (w1*(y-mu))/(s1^2) , (w2*(y-mu))/(s2^2))
     dldm
                                    },
               d2ldm2 = function(y,mu,sigma,nu,tau){
       s1 <- sigma/nu
       s2 <- sigma*nu
     dsq1 <- ((y-mu)^2)/(s1^2)
     dsq2 <-((y-mu)^2)/(s2^2) 
       w1 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq1),1)
       w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq2),1)
     dldm <- ifelse(y < mu, (w1*(y-mu))/(s1^2) , (w2*(y-mu))/(s2^2) )
   d2ldm2 <- -dldm*dldm
   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)  # NOTE this was added afrer testing
   d2ldm2
                                   },
                 dldd = function(y,mu,sigma,nu,tau) {
       s1 <- sigma/nu
       s2 <- sigma*nu
     dsq1 <- ((y-mu)^2)/(s1^2)
     dsq2 <- ((y-mu)^2)/(s2^2)
       w1 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq1),1)
       w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq2),1)
     dldd <- ifelse(y < mu, (w1*dsq1-1)/(sigma) , (w2*dsq2-1)/(sigma) )
     dldd
                                     } ,
               d2ldd2 = function(y,mu,sigma,nu,tau) {
       s1 <- sigma/nu
       s2 <- sigma*nu
     dsq1 <- ((y-mu)^2)/(s1^2)
     dsq2 <- ((y-mu)^2)/(s2^2)
       w1 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq1),1)
       w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq2),1)
     dldd <- ifelse(y < mu, (w1*dsq1-1)/(sigma) , (w2*dsq2-1)/(sigma) )
   d2ldd2 <- -dldd*dldd
   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
   d2ldd2
                                     },
                 dldv = function(y,mu,sigma,nu,tau) {
       s1 <- sigma/nu
       s2 <- sigma*nu
     dsq1 <- ((y-mu)^2)/(s1^2)
     dsq2 <- ((y-mu)^2)/(s2^2)
       w1 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq1),1)
       w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq2),1)
     dldv <- ifelse(y < mu, -(w1*dsq1-1)/(nu) , (w2*dsq2+1)/(nu) )
     dldv <- dldv - 2*nu/(1+nu^2)
     dldv
                                     } ,
               d2ldv2 = function(y,mu,sigma,nu,tau) { 
       s1 <- sigma/nu
       s2 <- sigma*nu
     dsq1 <- ((y-mu)^2)/(s1^2)
     dsq2 <- ((y-mu)^2)/(s2^2)
    w1 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq1),1)
    w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq2),1)
     dldv <- ifelse(y < mu, -(w1*dsq1-1)/(nu) , (w2*dsq2+1)/(nu) )
     dldv <- dldv - 2*nu/(1+nu^2)
    d2ldv2 <- -dldv*dldv
   d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)  
   d2ldv2
    
                                     },
                 dldt = function(y,mu,sigma,nu,tau) { 
       s1 <- sigma/nu
       s2 <- sigma*nu
     dsq1 <- ((y-mu)^2)/(s1^2)
     dsq2 <- ((y-mu)^2)/(s2^2)
    w1 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq1),1)
    w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq2),1)
    dldta <- -0.5*log(1+dsq1/tau)+(w1*dsq1-1)/(2*tau)
    dldtb <- -0.5*log(1+dsq2/tau)+(w2*dsq2-1)/(2*tau)
     dldt <- ifelse(y < mu, dldta , dldtb)
     dldt <- dldt+0.5*digamma((tau+1)/2)-0.5*digamma(tau/2)
     dldt
                                     } ,
               d2ldt2 = function(y,mu,sigma,nu,tau) {
       s1 <- sigma/nu
       s2 <- sigma*nu
     dsq1 <- ((y-mu)^2)/(s1^2)
     dsq2 <- ((y-mu)^2)/(s2^2)
    w1 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq1),1)
    w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq2),1)
    dldta <- -0.5*log(1+dsq1/tau)+(w1*dsq1-1)/(2*tau)
    dldtb <- -0.5*log(1+dsq2/tau)+(w2*dsq2-1)/(2*tau)
     dldt <- ifelse(y < mu, dldta , dldtb)
     dldt <- dldt+0.5*digamma((tau+1)/2)-0.5*digamma(tau/2)
    d2ldt2 <- -dldt*dldt
    d2ldt2 <- ifelse(d2ldt2 < -1e-15, d2ldt2,-1e-15)                                    
    d2ldt2
                                     } ,
              d2ldmdd = function(y,mu,sigma,nu,tau)  
                     {
       s1 <- sigma/nu
       s2 <- sigma*nu
     dsq1 <- ((y-mu)/s1)^2
     dsq2 <- ((y-mu)/s2)^2
    w1 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq1),1)
    w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq2),1)
     dldm <- ifelse(y < mu, (w1*(y-mu))/(s1^2) , (w2*(y-mu))/(s2^2))
     dldd <- ifelse(y < mu, (w1*dsq1-1)/(sigma) , (w2*dsq2-1)/(sigma) )                
  d2ldmdd <- -(dldm*dldd)
  d2ldmdd                  
                      },
              d2ldmdv = function(y,mu,sigma,nu,tau)
                      {
       s1 <- sigma/nu
       s2 <- sigma*nu
     dsq1 <- ((y-mu)/s1)^2
     dsq2 <- ((y-mu)/s2)^2
    w1 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq1),1)
    w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq2),1)
     dldm <- ifelse(y < mu, (w1*(y-mu))/(s1^2) , (w2*(y-mu))/(s2^2))
     dldv <- ifelse(y < mu, -(w1*dsq1-1)/(nu) , (w2*dsq2+1)/(nu) )
     dldv <- dldv - 2*nu/(1+nu^2)                       
  d2ldmdv <- -(dldm*dldv)
  d2ldmdv
                      } ,
              d2ldmdt = function(y,mu,sigma,nu,tau)
                     {  
       s1 <- sigma/nu
       s2 <- sigma*nu
     dsq1 <- ((y-mu)/s1)^2
     dsq2 <- ((y-mu)/s2)^2
    w1 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq1),1)
    w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq2),1)
     dldm <- ifelse(y < mu, (w1*(y-mu))/(s1^2) , (w2*(y-mu))/(s2^2))                      
    dldta <- -0.5*log(1+dsq1/tau)+(w1*dsq1-1)/(2*tau)
    dldtb <- -0.5*log(1+dsq2/tau)+(w2*dsq2-1)/(2*tau)
     dldt <- ifelse(y < mu, dldta , dldtb)
     dldt <- dldt+0.5*digamma((tau+1)/2)-0.5*digamma(tau/2)                   
  d2ldmdt <- -(dldm*dldt)
   d2ldmdt
                      },
              d2ldddv = function(y,mu,sigma,nu,tau)  
                      {
       s1 <- sigma/nu
       s2 <- sigma*nu
     dsq1 <- ((y-mu)^2)/(s1^2)
     dsq2 <- ((y-mu)^2)/(s2^2)
    w1 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq1),1)
    w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq2),1)
     dldd <- ifelse(y < mu, (w1*dsq1-1)/(sigma) , (w2*dsq2-1)/(sigma) )
     dldv <- ifelse(y < mu, -(w1*dsq1-1)/(nu) , (w2*dsq2+1)/(nu) )
     dldv <- dldv - 2*nu/(1+nu^2)         
  d2ldddv <- -(dldd*dldv)
  d2ldddv
                      },
              d2ldddt = function(y,mu,sigma,nu,tau)  
                      {
       s1 <- sigma/nu
       s2 <- sigma*nu
     dsq1 <- ((y-mu)^2)/(s1^2)
     dsq2 <- ((y-mu)^2)/(s2^2)
    w1 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq1),1)
    w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq2),1)
     dldd <- ifelse(y < mu, (w1*dsq1-1)/(sigma) , (w2*dsq2-1)/(sigma) )      
    dldta <- -0.5*log(1+dsq1/tau)+(w1*dsq1-1)/(2*tau)
    dldtb <- -0.5*log(1+dsq2/tau)+(w2*dsq2-1)/(2*tau)
     dldt <- ifelse(y < mu, dldta , dldtb)
     dldt <- dldt+0.5*digamma((tau+1)/2)-0.5*digamma(tau/2)            
  d2ldddt <- -(dldd*dldt) 
  d2ldddt                  
                      },
              d2ldvdt = function(y,mu,sigma,nu,tau)  
                      {
        s1 <- sigma/nu
       s2 <- sigma*nu
     dsq1 <- ((y-mu)^2)/(s1^2)
     dsq2 <- ((y-mu)^2)/(s2^2)
    w1 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq1),1)
    w2 <- ifelse(tau < 1000000, (tau+1)/(tau+dsq2),1)
     dldv <- ifelse(y < mu, -(w1*dsq1-1)/(nu) , (w2*dsq2+1)/(nu) )
     dldv <- dldv - 2*nu/(1+nu^2)
    dldta <- -0.5*log(1+dsq1/tau)+(w1*dsq1-1)/(2*tau)
    dldtb <- -0.5*log(1+dsq2/tau)+(w2*dsq2-1)/(2*tau)
     dldt <- ifelse(y < mu, dldta , dldtb)
     dldt <- dldt+0.5*digamma((tau+1)/2)-0.5*digamma(tau/2)          
  d2ldvdt <- -(dldv*dldt)
  d2ldvdt
      
                      }, 
          G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
                                  -2*dST3(y,mu,sigma,nu,tau,log=TRUE), 
                 rqres = expression(rqres(pfun="pST3", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)) ,
            mu.initial = expression(mu <- (y+mean(y))/2), 
         sigma.initial = expression(sigma<- rep(sd(y), length(y))),
            nu.initial = expression(nu <- rep(1, length(y))), 
           tau.initial = expression(tau <-rep(10, length(y))), 
              mu.valid = function(mu) TRUE, 
           sigma.valid = function(sigma)  all(sigma > 0),
              nu.valid = function(nu) all(nu > 0), 
             tau.valid = function(tau) all(tau > 0), 
               y.valid = function(y)  TRUE
          ),
            class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------  
dST3 <- function(x, mu=0, sigma=1, nu=1, tau=10, log=FALSE)
 {
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
          if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))  
          loglik1a <- dt((nu*(x-mu)/sigma), df=tau, log=TRUE)
          loglik2a <- dt((x-mu)/(sigma*nu), df=tau, log=TRUE)
          loglika <- ifelse(x < mu, loglik1a, loglik2a)
          loglika <- loglika+log(2*nu/(1+nu^2)) - log(sigma)
           loglik1b <- dNO((nu*(x-mu)/sigma), mu=0, sigma=1, log=TRUE)
          loglik2b <- dNO((x-mu)/(sigma*nu), mu=0, sigma=1, log=TRUE)
          loglikb <- ifelse(x < mu, loglik1b, loglik2b)
          loglikb <- loglikb+log(2*nu/(1+nu^2)) - log(sigma)
     if (length(tau)>1) loglik <- ifelse(tau<1000000, 
                      loglika, 
                      loglikb) 
    else loglik <- if (tau<1000000) loglika
                else  loglikb
       if(log==FALSE) ft  <- exp(loglik) else ft <- loglik 
       ft
  }    
#-----------------------------------------------------------------  
pST3 <- function(q, mu=0, sigma=1, nu=1, tau=10, lower.tail = TRUE, log.p = FALSE)
 {  
         if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
         if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
         if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
         cdf1 <- 2*pt(nu*(q-mu)/sigma, df=tau)
         cdf2 <- 1 + 2*nu*nu*(pt((q-mu)/(sigma*nu), df=tau) - 0.5)
         cdf <- ifelse(q < mu, cdf1, cdf2)
         cdf <- cdf/(1+nu^2)
         if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf 
         if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf) 
         cdf    
 }
#-----------------------------------------------------------------  
qST3 <- function(p, mu=0, sigma=1, nu=1, tau=10, lower.tail = TRUE, log.p = FALSE)
 { 
    if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
    if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
    if (lower.tail==TRUE) p <- p else p <- 1-p
    suppressWarnings(q1 <- mu+(sigma/nu)*qt(p*(1+nu^2)/2, df=tau))
    suppressWarnings(q2 <- mu+(sigma*nu)*qt((p*(1+nu^2)-1)/(2*nu^2) + 0.5, df=tau))
    q <- ifelse(p < (1/(1+nu^2)), q1, q2)
    q
 }
#-----------------------------------------------------------------  
rST3 <- function(n, mu=0, sigma=1, nu=1, tau=10)
  {
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
    if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))  
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))    
    n <- ceiling(n)
    p <- runif(n)
    r <- qST3(p,mu=mu,sigma=sigma,nu=nu,tau=tau)
    r
  }
