#   amended 27_11_2007  This is working very well
GIG <- function (mu.link="log", sigma.link="log", nu.link ="identity") 
{
    mstats <- checklink("mu.link", "GIG", substitute(mu.link), 
                         c("1/mu^2", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "GIG", substitute(sigma.link), #
                         c("inverse", "log", "identity", "own"))
    vstats <- checklink("nu.link", "GIG",substitute(nu.link), 
                         c("inverse", "log", "identity", "own"))  
    
    structure(
          list(family = c("GIG", "Generalised Inverse Gaussian"),
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
                 dldm = function(y,mu,sigma,nu) 
                        {
       c <- exp(log(besselK(1/(sigma^2),nu+1))-log(besselK(1/(sigma^2),nu)))  
    dldm <- -(nu/mu)+((c*y)/(2*(sigma^2)*mu^2))-1/(2*(sigma^2)*c*y)
    dldm
                        },
        
           d2ldm2 = function(y,mu,sigma,nu) 
                        {
         c <- exp(log(besselK(1/(sigma^2),nu+1))-log(besselK(1/(sigma^2),nu)))  
      d2ldm2 <- (nu-(c/(sigma^2)))/(mu^2)
      d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15) 
      d2ldm2
                        },

                 dldd = function(y,mu,sigma,nu) 
                        {
       c <- exp(log(besselK(1/(sigma^2),nu+1))-log(besselK(1/(sigma^2),nu)))  
    dcdd <- (c*(sigma^2)*(2*nu+1)+1-c*c)/((sigma^2)*(sigma^2))
    dldd <- (1/(sigma^2))*(nu-(c/(sigma^2))+(1/(2*(sigma^2)))*((c*y/mu)+(mu/(c*y)))+dcdd*((sigma^2)*nu/c-(1/2)*((y/mu)-(mu/(c^2*y)))))    
    dldd <- dldd*(2*sigma)
    dldd
                         },
               d2ldd2 = function(y,mu,sigma,nu) 
                         {
#     this uses the squared first derivative
        c <- exp(log(besselK(1/(sigma^2),nu+1))-log(besselK(1/(sigma^2),nu)))  
     dcdd <- (c*(sigma^2)*(2*nu+1)+1-c*c)/((sigma^2)*(sigma^2))
    dldd <- (1/(sigma^2))*(nu-(c/(sigma^2))+(1/(2*(sigma^2)))*((c*y/mu)+(mu/(c*y)))+dcdd*((sigma^2)*nu/c-(1/2)*((y/mu)-(mu/(c^2*y)))))    
    dldd <- dldd*(2*sigma)
    d2ldd2 <- -dldd*dldd
   d2ldd2 <- ifelse(d2ldd2 < -1e-6, d2ldd2,-1e-6)  
   d2ldd2
                          },
                 dldv = function(y,mu,sigma,nu) 
                          {
       nd <- numeric.deriv(dGIG(y, mu, sigma, nu, log=TRUE), "nu", delta=0.01)
     dldv <- as.vector(attr(nd, "gradient"))                
     dldv
                          },
               d2ldv2 = function(y,mu,sigma,nu) 
                          {
        nd <- numeric.deriv(dGIG(y, mu, sigma, nu, log=TRUE), "nu", delta=0.01)
      dldv <- as.vector(attr(nd, "gradient"))                
    d2ldv2 <- -dldv*dldv
   d2ldv2 <- ifelse(d2ldv2 < -1e-6, d2ldv2,-1e-6)  
   d2ldv2
                          },
              d2ldmdd = function(y,mu,sigma,nu) 
                          {
        c <- exp(log(besselK(1/(sigma^2),nu+1))-log(besselK(1/(sigma^2),nu)))  
     dldm <- -nu/mu+c*y/(2*(sigma^2)*mu^2)-1/(2*(sigma^2)*c*y)
     dcdd <- (c*(sigma^2)*(2*nu+1)+1-c*c)/((sigma^2)*(sigma^2))
    dldd <- (1/(sigma^2))*(nu-(c/(sigma^2))+(1/(2*(sigma^2)))*((c*y/mu)+(mu/(c*y)))+dcdd*((sigma^2)*nu/c-(1/2)*((y/mu)-(mu/(c^2*y)))))    
    dldd <- dldd*(2*sigma)
  d2ldmdd <- -dldm*dldd
  d2ldmdd
                        },
              d2ldmdv = function(y,mu,sigma,nu) 
                        { 
       c <- exp(log(besselK(1/(sigma^2),nu+1))-log(besselK(1/(sigma^2),nu)))   
    dldm <- -nu/mu+c*y/(2*(sigma^2)*mu^2)-1/(2*(sigma^2)*c*y)
      nd <- numeric.deriv(dGIG(y, mu, sigma, nu, log=TRUE), "nu", delta=0.01)
    dldv <- as.vector(attr(nd, "gradient"))  
 d2ldmdv <- -dldm*dldv      
 d2ldmdv        
                        },
              d2ldddv = function(y,mu,sigma,nu) 
                        {
       c <- exp(log(besselK(1/(sigma^2),nu+1))-log(besselK(1/(sigma^2),nu)))   
    dcdd <- (c*(sigma^2)*(2*nu+1)+1-c*c)/((sigma^2)*(sigma^2))
    dldd <- (1/(sigma^2))*(nu-(c/(sigma^2))+(1/(2*(sigma^2)))*((c*y/mu)+(mu/(c*y)))+dcdd*((sigma^2)*nu/c-(1/2)*((y/mu)-(mu/(c^2*y)))))    
    dldd <- dldd*(2*sigma)
      nd <- numeric.deriv(dGIG(y, mu, sigma, nu, log=TRUE), "nu", delta=0.01)
    dldv <- as.vector(attr(nd, "gradient"))
 d2ldddv <- -dldv*dldd  
 d2ldddv
                        },
          G.dev.incr  = function(y,mu,sigma,nu,...) 
                     -2*dGIG(y,mu=mu,sigma=sigma,nu=nu,log=TRUE), 
             rqres = expression(rqres(pfun="pGIG", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)),
        mu.initial = expression( mu <- (y+mean(y))/2), 
     sigma.initial = expression( sigma <- sd(y)/mean(y)), 
        nu.initial = expression( nu <- rep(-0.5, length(y))),  
          mu.valid = function(mu) TRUE , 
       sigma.valid = function(sigma)  all(sigma > 0),
          nu.valid = function(nu) TRUE , 
           y.valid = function(y) all(y>0)
          ),
            class = c("gamlss.family","family"))
}
#--------------------------------------------------------------
dGIG <- function(x, mu=1, sigma=1, nu=1,  log = FALSE)
 {
          if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(x < 0))  stop(paste("x must be positive", "\n", "")) 
               c <- exp(log(besselK(1/(sigma^2),nu+1))-log(besselK(1/(sigma^2),nu)))  
          loglik <- nu*log(c)-nu*log(mu)+(nu-1)*log(x)-log(2)-log(besselK(1/(sigma^2),nu))-1/(2*(sigma^2))*(c*x/mu+mu/(c*x))
          if(log==FALSE) ft  <- exp(loglik) else ft <- loglik 
          ft
  }    
#--------------------------------------------------------------  
pGIG <- function(q, mu=1, sigma=1, nu=1,  lower.tail = TRUE, log.p = FALSE)
 {  
          if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(q < 0))  stop(paste("q must be positive", "\n", ""))  
         lq <- length(q)       
      sigma <- rep(sigma, length = lq)
         mu <- rep(mu, length = lq)   
         nu <- rep(nu, length = lq) 
        cdf <-rep(0, lq)
       for (i in 1:lq)
          {
        cdf[i] <- integrate(function(x) 
                 dGIG(x, mu = 1, sigma = sigma[i], nu = nu[i]), 0.001, q[i]/mu[i] )$value #md br 7-10-11
          }    
    if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf 
    if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf) 
    cdf
 }
#--------------------------------------------------------------
qGIG <- function(p,  mu=1, sigma=1, nu=1, lower.tail = TRUE, log.p = FALSE) 
 {
    #---functions--------------------------------------------   
       h1 <- function(q)
       { 
     pGIG(q , mu = mu[i], sigma = sigma[i], nu = nu[i])-p[i]   
       }
       h <- function(q)
       { 
     pGIG(q , mu = mu[i], sigma = sigma[i], nu = nu[i])   
       }
     #-------------------------------------------------------
    if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))      
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (lower.tail==TRUE) p <- p else p <- 1-p
    if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))     
         lp <-  max(length(p),length(mu),length(sigma),length(nu))
          p <- rep(p, length = lp)                                                                     
      sigma <- rep(sigma, length = lp)
         mu <- rep(mu, length = lp)
         nu <- rep(nu, length = lp)
          q <- rep(0,lp)      
         for (i in  seq(along=p)) 
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
           interval <-  interval <- c(.Machine$double.xmin, mu[i])
           }
        q[i] <- uniroot(h1, interval)$root
         }
    q
   }
#--------------------------------------------------------------
rGIG <- function(n, mu=1, sigma=1, nu=1, ...)
  {
    if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))    
    n <- ceiling(n)
    p <- runif(n)
    r <- qGIG(p,mu=mu,sigma=sigma,nu=nu, ...)
    r
  }
#--------------------------------------------------------------
