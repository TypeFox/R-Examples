# amended 27_11_07
PE <- function (mu.link="identity", sigma.link="log", nu.link ="log") 
{

    mstats <- checklink(   "mu.link", "Power Exponential", substitute(mu.link),    
                                                         c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "Power Exponential", substitute(sigma.link), 
                                                         c("inverse", "log", "identity", "own"))
    vstats <- checklink(   "nu.link", "Power Exponential", substitute(nu.link),    
                                                         c("logshiftto1", "log", "identity", "own"))
        
    structure(
          list(family = c("PE", "Power Exponential"),
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
                          log.c <- 0.5*(-(2/nu)*log(2)+lgamma(1/nu)-lgamma(3/nu))
                              c <- exp(log.c)
                              z <- (y-mu)/sigma
                           dldm <- (sign(z)*nu)/(2*sigma*abs(z)) 
                           dldm <- dldm*((abs(z/c))^nu) 
                                    },
               d2ldm2 = function(y,mu,sigma,nu) {
                          log.c <- 0.5*(-(2/nu)*log(2)+lgamma(1/nu)-lgamma(3/nu))
                              c <- exp(log.c)
                              z <- (y-mu)/sigma
                           dldm <- (sign(z)*nu)/(2*sigma*abs(z)) 
                           dldm <- dldm*((abs(z/c))^nu) 
                         d2ldm2 <- -(nu*nu*gamma(2-(1/nu))*gamma(3/nu))/((sigma*gamma(1/nu))^2) 
                         d2ldm2 <- if (any(nu<1.05)) -dldm*dldm else d2ldm2
                                   },
                 dldd = function(y,mu,sigma,nu) {
                          log.c <- 0.5*(-(2/nu)*log(2)+lgamma(1/nu)-lgamma(3/nu))
                              c <- exp(log.c)
                              z <- (y-mu)/sigma
                           dldd <- ((nu/2)*((abs(z/c))^nu)-1)/sigma
                           dldd   
                                    },
               d2ldd2 = function(sigma,nu) { 
                                  d2ldd2 <- -nu/(sigma^2)
                                  d2ldd2                                   
                                   },
                 dldv = function(y,mu,sigma,nu) {
                          log.c <- 0.5*(-(2/nu)*log(2)+lgamma(1/nu)-lgamma(3/nu))
                              c <- exp(log.c)
                              z <- (y-mu)/sigma  
                       dlogc.dv <- (1/(2*nu^2))*(2*log(2)-digamma(1/nu)+3*digamma(3/nu))
                           dldv <- (1/nu)-0.5*((log(abs(z/c)))*((abs(z/c))^nu)) 
                           dldv <- dldv+log(2)/(nu^2)+digamma(1/nu)/(nu^2)
                           dldv <- dldv+(-1+(nu/2)*((abs(z/c))^nu))*dlogc.dv
                                    }, 
               d2ldv2 = function(y,mu,sigma,nu) { 
                       dlogc.dv <- (1/(2*nu^2))*(2*log(2)-digamma(1/nu)+3*digamma(3/nu))
                              p <- (1+nu)/nu                      
                          part1 <- p*trigamma(p)+2*(digamma(p))^2
                          part2 <- digamma(p)*(log(2)+3-3*digamma(3/nu)-nu)
                          part3 <- -3*(digamma(3/nu))*(1+log(2))    
                          part4 <- -(nu+log(2))*log(2)
                          part5 <- -nu+(nu^4)*(dlogc.dv)^2
                         d2ldv2 <- part1+part2+part3+part4+part5
                         d2ldv2 <- -d2ldv2/nu^3    
                         d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)                                    
                         d2ldv2 
                                   },
              d2ldmdd = function(y) rep(0,length(y)),
              d2ldmdv = function(y) rep(0,length(y)),
              d2ldddv = function(y,mu,sigma,nu) (1/(2*sigma))*((3/nu)*(digamma(1/nu)
                                    -digamma(3/nu))+2+(2/nu)), #
          G.dev.incr  = function(y,mu,sigma,nu,...)  -2*dPE(y,mu,sigma,nu,log=TRUE),
                 rqres = expression(
                    rqres(pfun="pPE", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu) 
                                    ),
            mu.initial = expression( mu <- (y+mean(y))/2 ),#y+0.00001
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
dPE<-function(x, mu=0, sigma=1, nu=2, log=FALSE)
 {
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))  
  log.c <- 0.5*(-(2/nu)*log(2)+lgamma(1/nu)-lgamma(3/nu))
      c <- exp(log.c)
      z <- (x-mu)/sigma
   log.lik <- -log(sigma)+log(nu)-log.c-(0.5*(abs(z/c)^nu))-(1+(1/nu))*log(2)-lgamma(1/nu)
     if(log==FALSE) fy  <- exp(log.lik) else fy <- log.lik
      fy 
  }
#----------------------------------------------------------------------------------------  
pPE<- function(q, mu=0, sigma=1, nu=2, lower.tail = TRUE, log.p = FALSE)
  {
   if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
   if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))  
   log.c <- 0.5*(-(2/nu)*log(2)+lgamma(1/nu)-lgamma(3/nu))
       c <- exp(log.c)
       z <- (q-mu)/sigma
       s <- 0.5*((abs(z/c))^nu)
     cdf <- 0.5*(1+pgamma(s,shape=1/nu,scale=1)*sign(z))

    if (length(nu)>1) cdf <- ifelse(nu>10000, 
                                   (q-(mu-sqrt(3)*sigma))/(sqrt(12)*sigma),
                                    cdf)
    else cdf <- if (nu>10000)   (q-(mu-sqrt(3)*sigma))/(sqrt(12)*sigma) 
                   else cdf
    if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf 
    if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf) 
     cdf
   }
#----------------------------------------------------------------------------------------  
 qPE<- function(p, mu=0, sigma=1, nu=2, lower.tail = TRUE, log.p = FALSE)
  {
   if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
   if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))  
   if (log.p==TRUE) p <- exp(p) else p <- p
   if (lower.tail==TRUE) p <- p else p <- 1-p
   if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))     
   log.c <- 0.5*(-(2/nu)*log(2)+lgamma(1/nu)-lgamma(3/nu))
       c <- exp(log.c)
        suppressWarnings(s <- qgamma((2*p-1)*sign(p-0.5),shape=(1/nu),scale=1))
       z <- sign(p-0.5)*((2*s)^(1/nu))*c
      ya <- mu + sigma*z  
      ya
   }                               
#----------------------------------------------------------------------------------------
rPE <- function(n, mu=0, sigma=1, nu=2)
  {
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))    
    n <- ceiling(n)
    p <- runif(n)
    r <- qPE(p,mu=mu,sigma=sigma,nu=nu)
    r
  }
#----------------------------------------------------------------------------------------
