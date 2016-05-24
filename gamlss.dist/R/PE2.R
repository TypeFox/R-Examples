# amended 27_11_2007
PE2 <- function (mu.link="identity", sigma.link="log", nu.link ="log") 
{

    mstats <- checklink(   "mu.link", "Power Exponential type 2", substitute(mu.link),    
                                                         c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Power Exponential type 2", substitute(sigma.link), 
                                                         c("inverse", "log", "identity", "own"))
    vstats <- checklink(   "nu.link", "Power Exponential type 2", substitute(nu.link),    
                                                         c("logshiftto1", "log", "identity", "own"))
        
    structure(
          list(family = c("PE2", "Power Exponential type 2"),
           parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE), 
                nopar = 3,
                 type = "Continuous", 
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
                 dldm = function(y,mu,sigma,nu) {
                               z <- (y-mu)/sigma
                           dldm <- (sign(z)*nu*(abs(z)^(nu-1)))/sigma
                           dldm
                                    },
               d2ldm2 = function(y,mu,sigma,nu) {
                               z <- (y-mu)/sigma
                           dldm <- (sign(z)*nu*(abs(z)^(nu-1)))/sigma
                         d2ldm2 <- -(nu*nu*gamma(2-(1/nu))*gamma(1/nu))/(sigma^2) 
                         d2ldm2 <- if (any(nu<1.05)) -dldm*dldm else d2ldm2
                         d2ldm2
                                   },
                 dldd = function(y,mu,sigma,nu) {
                              z <- (y-mu)/sigma
                           dldd <- (nu*(abs(z)^nu)-1)/(sigma)
                           dldd   
                                    },
               d2ldd2 = function(sigma,nu) { 
                                  d2ldd2 <- -nu/(sigma^2)
                                  d2ldd2                                   
                                   },
                 dldv = function(y,mu,sigma,nu) {
                              z <- (y-mu)/sigma  
                           dldv <- (1/nu)-((log(abs(z)))*((abs(z))^nu)) 
                           dldv <- dldv+(digamma(1/nu))/(nu^2)
                           dldv 
                                    }, 
               d2ldv2 = function(y,mu,sigma,nu) { 
                               p <- (1+nu)/nu                      
                          part1 <- p*trigamma(p)+(digamma(p))^2
                          part2 <- 2*(digamma(p))
                         d2ldv2 <- part1+part2
                         d2ldv2 <- -d2ldv2/nu^3    
                         d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)                                    
                         d2ldv2 
                                   },
              d2ldmdd = function(y) rep(0,length(y)),
              d2ldmdv = function(y) rep(0,length(y)),
              d2ldddv = function(y,mu,sigma,nu) (1+nu+digamma(1/nu))/(sigma*nu), 
          G.dev.incr  = function(y,mu,sigma,nu,...)  -2*dPE2(y,mu,sigma,nu,log=TRUE),
                 rqres = expression(
                    rqres(pfun="pPE2", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu) 
                                    ),
            mu.initial = expression( mu <- (y+mean(y))/2 ),
         sigma.initial = expression( sigma <- (abs(y-mean(y))+sd(y))/2 ),
            nu.initial = expression( nu <- rep(1.8, length(y))), 
              mu.valid = function(mu) TRUE , 
           sigma.valid = function(sigma)  all(sigma > 0),
              nu.valid = function(nu) all(nu > 0), 
               y.valid = function(y)  TRUE
          ),
            class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------
dPE2<-function(x, mu=0, sigma=1, nu=2, log=FALSE)
 {
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))  
      z <- (x-mu)/sigma
   log.lik <- -log(sigma)+log(nu)-log(2)-((abs(z)^nu))-lgamma(1/nu)
     if(log==FALSE) fy  <- exp(log.lik) else fy <- log.lik
      fy 
  }
#----------------------------------------------------------------------------------------  
pPE2<- function(q, mu=0, sigma=1, nu=2, lower.tail = TRUE, log.p = FALSE)
  {
   if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
   if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))  
       z <- (q-mu)/sigma
       s <- (abs(z)^nu)
     cdf <- 0.5*(1+pgamma(s,shape=1/nu,scale=1)*sign(z))

    if (length(nu)>1) cdf <- ifelse(nu>10000, 
                                   (q-(mu-sigma))/(2*sigma),
                                    cdf)
    else cdf <- if (nu>10000)    (q-(mu-sigma))/(2*sigma) 
                   else cdf

    if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf 
    if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf) 
     cdf
   }
#----------------------------------------------------------------------------------------  
 qPE2<- function(p, mu=0, sigma=1, nu=2, lower.tail = TRUE, log.p = FALSE)
  {
   if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
   if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))  
   if (log.p==TRUE) p <- exp(p) else p <- p
   if (lower.tail==TRUE) p <- p else p <- 1-p
   if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))     
        suppressWarnings(s <- qgamma((2*p-1)*sign(p-0.5),shape=(1/nu),scale=1))
       z <- sign(p-0.5)*((s)^(1/nu))
      ya <- mu + sigma*z  
      ya
   }                               
#----------------------------------------------------------------------------------------
rPE2 <- function(n, mu=0, sigma=1, nu=2)
  {
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))    
    n <- ceiling(n)
    p <- runif(n)
    r <- qPE2(p,mu=mu,sigma=sigma,nu=nu)
    r
  }
#----------------------------------------------------------------------------------------
