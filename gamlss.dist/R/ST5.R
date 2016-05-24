# amended 28_11_2007
# the first derivatives squares have been used here
#----------------------------------------------------------------------------------------
ST5<- function (mu.link="identity", sigma.link="log", nu.link ="identity", tau.link="log")
{
    mstats <- checklink(   "mu.link", "Skew t, type 3", substitute(mu.link), 
                           c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Skew t, type 3", substitute(sigma.link), 
                           c("inverse", "log", "identity", "own"))
    vstats <- checklink(   "nu.link", "Skew t, type 3", substitute(nu.link),    
                           c("1/nu^2", "log", "identity", "own"))
    tstats <- checklink(  "tau.link", "Skew t, type 3", substitute(tau.link),   
                           c("inverse", "log", "identity", "own")) 
    structure(
          list(family = c("ST5", "Skew t, type 5"),
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
       u <- (((2*sigma^2)+tau*((y-mu)^2))^(-0.5))*sqrt(tau)*(y-mu)
       r <- sqrt((nu^2)+2*tau)
      su <- sqrt(1-u^2)
    dldm <- (1-nu/r+tau/2)*su*(1+u)-(1+nu/r+tau/2)*su*(1-u)
    dldm <- dldm/(sigma*sqrt(2*tau))
                    },
   d2ldm2 = function(y,mu,sigma,nu,tau){
       u <- (((2*sigma^2)+tau*((y-mu)^2))^(-0.5))*sqrt(tau)*(y-mu)
       r <- sqrt((nu^2)+2*tau)
      su <- sqrt(1-u^2)
    dldm <- (1-nu/r+tau/2)*su*(1+u)-(1+nu/r+tau/2)*su*(1-u)
    dldm <- dldm/(sigma*sqrt(2*tau))
  d2ldm2 <- -dldm*dldm
   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)  
   d2ldm2
                      },
    dldd = function(y,mu,sigma,nu,tau) {  
       u <- (((2*sigma^2)+tau*((y-mu)^2))^(-0.5))*sqrt(tau)*(y-mu)
       r <- sqrt((nu^2)+2*tau)
    dldd <- (1-nu/r+tau/2)*u*(1+u)-(1+nu/r+tau/2)*u*(1-u)
    dldd <- -(1/sigma) + dldd/(sigma*tau) 
      } ,
   d2ldd2 = function(y,mu,sigma,nu,tau){
       u <- (((2*sigma^2)+tau*((y-mu)^2))^(-0.5))*sqrt(tau)*(y-mu)
       r <- sqrt((nu^2)+2*tau)
    dldd <- (1-nu/r+tau/2)*u*(1+u)-(1+nu/r+tau/2)*u*(1-u)
    dldd <- -(1/sigma) + dldd/(sigma*tau) 
  d2ldd2 <- -dldd*dldd
  d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
  d2ldd2
                     },   
     dldv = function(y,mu,sigma,nu,tau){ 
        u <- (((2*sigma^2)+tau*((y-mu)^2))^(-0.5))*sqrt(tau)*(y-mu)
        r <- sqrt((nu^2)+2*tau)
     dldv <- -digamma((1+nu/r)/tau) + digamma((1-nu/r)/tau)+log(1+u)-log(1-u)
     dldv <- 2*dldv/(r^3)
                       } ,
    d2ldv2 = function(y,mu,sigma,nu,tau){ 
        u <- (((2*sigma^2)+tau*((y-mu)^2))^(-0.5))*sqrt(tau)*(y-mu)
        r <- sqrt((nu^2)+2*tau)
     dldv <- -digamma((1+nu/r)/tau) + digamma((1-nu/r)/tau) +log(1+u)-log(1-u)
     dldv <- 2*dldv/(r^3) 
   d2ldv2 <- -dldv*dldv             
   d2ldv2 <- ifelse(d2ldv2 < -1e-4, d2ldv2,-1e-4)  
   d2ldv2
                       },
      dldt = function(y,mu,sigma,nu,tau){
        u <- (((2*sigma^2)+tau*((y-mu)^2))^(-0.5))*sqrt(tau)*(y-mu)
        r <- sqrt((nu^2)+2*tau)
     dldt <- 4*digamma(2/tau)-tau-4*log(2)
     dldt <- dldt + (1-nu/r+tau/2)*u*(1+u)-(1+nu/r+tau/2)*u*(1-u)
     dldt <- dldt + 2*(1+nu*(r^2+tau)/r^3)*(log(1+u)-digamma((1+nu/r)/tau))
     dldt <- dldt + 2*(1-nu*(r^2+tau)/r^3)*(log(1-u)-digamma((1-nu/r)/tau))
     dldt <- -dldt/(2*tau*tau)
                       } ,
    d2ldt2 = function(y,mu,sigma,nu,tau){ 
        u <- (((2*sigma^2)+tau*((y-mu)^2))^(-0.5))*sqrt(tau)*(y-mu)
        r <- sqrt((nu^2)+2*tau)
     dldt <- 4*digamma(2/tau)-tau-4*log(2)
     dldt <- dldt + (1-nu/r+tau/2)*u*(1+u)-(1+nu/r+tau/2)*u*(1-u)
     dldt <- dldt + 2*(1+nu*(r^2+tau)/r^3)*(log(1+u)-digamma((1+nu/r)/tau))
     dldt <- dldt + 2*(1-nu*(r^2+tau)/r^3)*(log(1-u)-digamma((1-nu/r)/tau))
     dldt <- -dldt/(2*tau*tau)
   d2ldt2 <-  -dldt*dldt   
   d2ldt2 <- ifelse(d2ldt2 < -1e-4, d2ldt2,-1e-4)  
   d2ldt2
                       } ,
  d2ldmdd = function(y,mu,sigma,nu,tau) {
       u <- (((2*sigma^2)+tau*((y-mu)^2))^(-0.5))*sqrt(tau)*(y-mu)
       r <- sqrt((nu^2)+2*tau)
      su <- sqrt(1-u^2)
    dldm <- (1-nu/r+tau/2)*su*(1+u)-(1+nu/r+tau/2)*su*(1-u)
    dldm <- dldm/(sigma*sqrt(2*tau))
    dldd <- (1-nu/r+tau/2)*u*(1+u)-(1+nu/r+tau/2)*u*(1-u)
    dldd <- -(1/sigma) + dldd/(sigma*tau)   
 d2ldmdd <- -(dldm*dldd)
 d2ldmdd
                       },
  d2ldmdv = function(y,mu,sigma,nu,tau) {
       u <- (((2*sigma^2)+tau*((y-mu)^2))^(-0.5))*sqrt(tau)*(y-mu)
       r <- sqrt((nu^2)+2*tau)
      su <- sqrt(1-u^2)
    dldm <- (1-nu/r+tau/2)*su*(1+u)-(1+nu/r+tau/2)*su*(1-u)
    dldm <- dldm/(sigma*sqrt(2*tau))
    dldv <- -digamma((1+nu/r)/tau) + digamma((1-nu/r)/tau)+log(1+u)-log(1-u)
    dldv <- 2*dldv/(r^3)   
 d2ldmdv <- -(dldm*dldv)
 d2ldmdv
                       },
  d2ldmdt = function(y,mu,sigma,nu,tau) {
       u <- (((2*sigma^2)+tau*((y-mu)^2))^(-0.5))*sqrt(tau)*(y-mu)
       r <- sqrt((nu^2)+2*tau)
      su <- sqrt(1-u^2)
    dldm <- (1-nu/r+tau/2)*su*(1+u)-(1+nu/r+tau/2)*su*(1-u)
    dldm <- dldm/(sigma*sqrt(2*tau))
    dldt <- 4*digamma(2/tau)-tau-4*log(2)
    dldt <- dldt + (1-nu/r+tau/2)*u*(1+u)-(1+nu/r+tau/2)*u*(1-u)
    dldt <- dldt + 2*(1+nu*(r^2+tau)/r^3)*(log(1+u)-digamma((1+nu/r)/tau))
    dldt <- dldt + 2*(1-nu*(r^2+tau)/r^3)*(log(1-u)-digamma((1-nu/r)/tau))
    dldt <- -dldt/(2*tau*tau)  
 d2ldmdt <- -(dldm*dldt)
 d2ldmdt
                       },
  d2ldddv = function(y,mu,sigma,nu,tau) {
       u <- (((2*sigma^2)+tau*((y-mu)^2))^(-0.5))*sqrt(tau)*(y-mu)
       r <- sqrt((nu^2)+2*tau)
    dldd <- (1-nu/r+tau/2)*u*(1+u)-(1+nu/r+tau/2)*u*(1-u)
    dldd <- -(1/sigma) + dldd/(sigma*tau)
    dldv <- -digamma((1+nu/r)/tau) + digamma((1-nu/r)/tau)+log(1+u)-log(1-u)
    dldv <- 2*dldv/(r^3) 
 d2ldddv <- -(dldd*dldv)
 d2ldddv
                       },
  d2ldddt = function(y,mu,sigma,nu,tau) {
       u <- (((2*sigma^2)+tau*((y-mu)^2))^(-0.5))*sqrt(tau)*(y-mu)
       r <- sqrt((nu^2)+2*tau)
    dldd <- (1-nu/r+tau/2)*u*(1+u)-(1+nu/r+tau/2)*u*(1-u)
    dldd <- -(1/sigma) + dldd/(sigma*tau)
    dldt <- 4*digamma(2/tau)-tau-4*log(2)
    dldt <- dldt + (1-nu/r+tau/2)*u*(1+u)-(1+nu/r+tau/2)*u*(1-u)
    dldt <- dldt + 2*(1+nu*(r^2+tau)/r^3)*(log(1+u)-digamma((1+nu/r)/tau))
    dldt <- dldt + 2*(1-nu*(r^2+tau)/r^3)*(log(1-u)-digamma((1-nu/r)/tau))
    dldt <- -dldt/(2*tau*tau) 
 d2ldddt <- -(dldd*dldt)
  d2ldddt  
                       },
  d2ldvdt = function(y,mu,sigma,nu,tau) {
       u <- (((2*sigma^2)+tau*((y-mu)^2))^(-0.5))*sqrt(tau)*(y-mu)
       r <- sqrt((nu^2)+2*tau)
    dldv <- -digamma((1+nu/r)/tau) + digamma((1-nu/r)/tau)+log(1+u)-log(1-u)
    dldv <- 2*dldv/(r^3)
    dldt <- 4*digamma(2/tau)-tau-4*log(2)
    dldt <- dldt + (1-nu/r+tau/2)*u*(1+u)-(1+nu/r+tau/2)*u*(1-u)
    dldt <- dldt + 2*(1+nu*(r^2+tau)/r^3)*(log(1+u)-digamma((1+nu/r)/tau))
    dldt <- dldt + 2*(1-nu*(r^2+tau)/r^3)*(log(1-u)-digamma((1-nu/r)/tau))
    dldt <- -dldt/(2*tau*tau)
 d2ldvdt <- -(dldv*dldt)
 d2ldvdt  
                       },
 G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
                       { 
         -2*dST5(y,mu,sigma,nu,tau,log=TRUE)
                        } ,                     
         rqres = expression(
            rqres(pfun="pST5", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)   
                           ) ,
    mu.initial = expression(mu <- (y+mean(y))/2) ,     
 sigma.initial = expression(sigma <- rep(sd(y)/4, length(y))),
    nu.initial = expression(nu <- rep(0.03, length(y))), 
   tau.initial = expression(tau <-rep(3, length(y))), 
      mu.valid = function(mu) TRUE, 
   sigma.valid = function(sigma)  all(sigma > 0),
      nu.valid = function(nu) TRUE , 
     tau.valid = function(tau) all(tau > 0), 
       y.valid = function(y)  TRUE
          ),
            class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------
dST5<- function(x, mu = 0, sigma = 1, nu = 0, tau = 1, log = FALSE)
 {
       if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
       if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
      z <- (x-mu)/sigma
     ve <- 2/tau
    lam <- 2*nu/(tau*sqrt(2*tau+nu*nu))
      a <- (ve+lam)/2
      b <- (ve-lam)/2
 loglik <- (a+0.5)*log(1+z/(sqrt(a+b+z*z)))+(b+0.5)*log(1-z/(sqrt(a+b+z*z)))-(a+b-1)*log(2)-0.5*log(a+b)-lbeta(a, b)-log(sigma)
       if(log==FALSE) ft  <- exp(loglik) else ft <- loglik 
       ft
  }    
#----------------------------------------------------------------------------------------
pST5<- function(q, mu = 0, sigma = 1, nu = 0, tau = 1, lower.tail = TRUE, log.p = FALSE)
 {  
      if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
      if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))           
      z <- (q-mu)/sigma
     ve <- 2/tau
#    if (length(tau)>1) lam <- ifelse(tau<1000000,2*nu/(tau*sqrt(2*tau+nu*nu)),2/tau) 
#     else lam <- if (tau<1000000) 2*nu/(tau*sqrt(2*tau+nu*nu)) else 2/tau
    lam <- 2*nu/(tau*sqrt(2*tau+nu*nu))
      a <- (ve+lam)/2
      b <- (ve-lam)/2
  alpha <- (1+z/(sqrt(a+b+z*z)))/2
# we need the incomplete beta function ratio, i.e. cdf of Beta(a,b) at alpha
      p <- pbeta(alpha,a,b)
      if(lower.tail==TRUE) p  <- p else  p <- 1-p 
      if(log.p==FALSE) p  <- p else  p <- log(p) 
      p
 }
#----------------------------------------------------------------------------------------
qST5<-  function(p, mu=0, sigma=1, nu=0, tau=1, lower.tail = TRUE, log.p = FALSE)
 {   
    if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p <= 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
    if (lower.tail==TRUE) p <- p else p <- 1-p
# we need the inverse cdf of Beta(a,b) corresponding to probability p (i.e. our alpha)
    ve <- 2/tau
#    if (length(tau)>1) lam <- ifelse(tau<1000000,2*nu/(tau*sqrt(2*tau+nu*nu)),2/tau) 
#     else lam <- if (tau<1000000) 2*nu/(tau*sqrt(2*tau+nu*nu)) else 2/tau
    lam <- 2*nu/(tau*sqrt(2*tau+nu*nu))
    a <- (ve+lam)/2
    b <- (ve-lam)/2
    Balpha <- qbeta(p,a,b) 
#    lzalpha <- 0.5*log(a+b) + log(2*Balpha-1) -log(2) - 0.5*log(Balpha) - 0.5*log(1-Balpha)
#    zalpha <- exp(lzalpha)
    zalpha <- (sqrt(a+b))*(2*Balpha-1)/(2*sqrt(Balpha*(1-Balpha)))
    q <- mu + sigma*zalpha  
    q
 }
#----------------------------------------------------------------------------------------
rST5<- function(n, mu=0, sigma=1, nu=0, tau=1)
  {
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    n <- ceiling(n)
    p <- runif(n)
    r <- qST5(p,mu=mu,sigma=sigma,nu=nu,tau=tau)
    r
  }
#----------------------------------------------------------------------------------------
