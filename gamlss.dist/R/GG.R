# last change 29 Jul 2011 DS
GG <- function (mu.link="log", sigma.link="log", nu.link ="identity") 
{
    mstats <- checklink("mu.link", "GG", substitute(mu.link), 
                         c("1/mu^2", "log", "identity"))
    dstats <- checklink("sigma.link", "GG", substitute(sigma.link), #
                         c("inverse", "log", "identity"))
    vstats <- checklink("nu.link", "GG",substitute(nu.link), 
                         c("1/nu^2", "log", "identity"))  
    structure(
          list(family = c("GG", "generalised Gamma Lopatatsidis-Green"),
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
                                   z <- (y/mu)^nu
                               theta <- 1/(sigma^2*abs(nu)^2)
                                dldm <- ifelse(abs(nu)>1e-06, (z-1)*theta*nu/mu, (1/(mu*(sigma^2))*(log(y)-log(mu))))
                                dldm
                                   },  
               d2ldm2 = function(mu,sigma,nu) 
                                   {
                             d2ldm2 <- ifelse(abs(nu)>1e-06, -1/((mu^2)*(sigma^2)), -(1/(mu^2*sigma^2)))
                             d2ldm2
                                   },
                 dldd = function(y,mu,sigma,nu) 
                                   {
                                  z <- (y/mu)^nu
                              theta <- 1/(sigma^2*abs(nu)^2)
                               dldd <- ifelse(abs(nu) > 1e-06, -2*theta*(log(theta)+1+log(z)-z-digamma(theta))/sigma,
                                             -(1/sigma)+(1/sigma^3)*(log(y)-log(mu))^2) #DS 29-7-11
                               dldd
                                   },
               d2ldd2 = function(y,mu,sigma,nu) 
                                   {
                              theta <- 1/(sigma^2*abs(nu)^2)
                             d2ldd2 <- ifelse(abs(nu) > 1e-06, 4*(theta/(sigma^2))*(1-theta*trigamma(theta)),-2/sigma^2) # DS 29-7-11
                             d2ldd2
                                   },
                 dldv = function(y,mu,sigma,nu) 
                                   {
                                  z <- (y/mu)^nu
                              theta <- 1/(sigma^2*abs(nu)^2)
                               dldv <- (1/nu)*(1+2*theta*(digamma(theta)+z-log(theta)-1- ((z+1)/2)*log(z)))                 
                               dldv
                                   },
               d2ldv2 = function(y,mu,sigma,nu) 
                                   {
                                 # z <- (y/mu)^nu
                              theta <- 1/(sigma^2*abs(nu)^2)
                             d2ldv2 <- -(theta/nu^2)*(trigamma(theta)*(1+4*theta)-(4+3/theta)-log(theta)*(2/theta-log(theta))+digamma(theta)*(digamma(theta)+(2/theta)-2*log(theta)))
                             d2ldv2
                                   },
              d2ldmdd = function(y) rep(0,length(y)),
              d2ldmdv = function(y,mu,sigma,nu) 
                                   { 
                              theta <- 1/(sigma^2*abs(nu)^2)
                                ddd <- (theta/mu)*(digamma(theta)+(1/theta)-log(theta))
                                ddd
                                   },
              d2ldddv = function(y,mu,sigma,nu) 
                                   {
                              theta <- 1/(sigma^2*abs(nu)^2)
                            d2ldddv <- -2*sign(nu)*theta^(3/2)*(2*theta*trigamma(theta)-(1/theta)-2)
                            d2ldddv
                                   },
          G.dev.incr  = function(y,mu,sigma,nu,...) -2*dGG(y,mu=mu,sigma=sigma,nu=nu,log=TRUE), 
             rqres = expression(rqres(pfun="pGG", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)),
        mu.initial = expression( mu <- (y+mean(y))/2), 
     sigma.initial = expression( sigma <- rep(1,length(y))), 
        nu.initial = expression( nu <- rep(1, length(y))), 
          mu.valid = function(mu) all(mu > 0), 
       sigma.valid = function(sigma)  all(sigma > 0),
          nu.valid = function(nu) TRUE , 
           y.valid = function(y) all(y>0)
          ),
            class = c("gamlss.family","family"))
}
#--------------------------------------------------------------
dGG <- function(x, mu=1, sigma=0.5, nu=1,  log = FALSE)
 {
          if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(x < 0))  stop(paste("x must be positive", "\n", "")) 
               z <- (x/mu)^nu 
           theta <- 1/(sigma^2*nu^2) 
         # loglik <- theta*log(theta)+theta*log(z)-theta*z-lgamma(theta)+log(abs(nu))-log(x)
           if(length(nu)>1)  loglik <- ifelse(abs(nu)>1e-06,dGA(z, mu=1, sigma=sigma*abs(nu), log=TRUE)+ log(abs(nu)*z/x),
                                    -log(x)-.5*log(2*pi)-log(sigma)-(1/(2*sigma^2))*(log(x)-log(mu))^2)
          else  if (abs(nu)>1e-06) loglik <-  dGA(z, mu=1, sigma=sigma*abs(nu), log=TRUE)+ log(abs(nu)*z/x) 
                 else loglik <- -log(x)-.5*log(2*pi)-log(sigma)-(1/(2*sigma^2))*(log(x)-log(mu))^2
          if(log==FALSE) ft  <- exp(loglik) else ft <- loglik 
          ft
  }    
#--------------------------------------------------------------  
pGG <- function(q, mu=1, sigma=0.5, nu=1,  lower.tail = TRUE, log.p = FALSE)
 {  
          if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(q < 0))  stop(paste("q must be positive", "\n", ""))  
            z <- (q/mu)^nu
         if(length(nu)>1)  cdf <- ifelse(abs(nu)>1e-06,
                           pGA(z, mu=1, sigma=sigma*abs(nu), lower.tail = lower.tail, log.p = log.p),
                           pNO(log(z), mu=log(mu), sigma=sigma))
          else  if (abs(nu)>1e-06) cdf <- pGA(z, mu=1, sigma=sigma*abs(nu), lower.tail = lower.tail, log.p = log.p)
                 else cdf <- pNO(log(q), mu=log(mu), sigma=sigma)
       cdf
 }

#--------------------------------------------------------------
qGG <- function(p, mu=1, sigma=0.5, nu=1,  lower.tail = TRUE, log.p = FALSE )
 { 
    if (any(mu < 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))  
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (lower.tail==TRUE) p <- p else p <- 1-p
    if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))       

    
    if(length(nu)>1) 
      {
       p <- ifelse(nu>0, p, 1-p)
       z <- ifelse(abs(nu)>1e-06, qGA(p, mu=1, sigma=sigma*abs(nu)),qNO(p, mu=log(mu), sigma=sigma,))
       y <- ifelse(abs(nu)>1e-06, mu*z^(1/nu), exp(z))
      }
    else  if (abs(nu)>1e-06) 
              {
               p <- if(nu>0)  p else 1-p
              z <- qGA(p, mu=1, sigma=sigma*abs(nu))
              y <- mu*z^(1/nu)
              }
          else 
              {
             #  p <- if(nu>0)  p else 1-p
              z <- qNO(p, mu=log(mu), sigma=sigma)
              y <- exp(z)
              }
    y
 }
#--------------------------------------------------------------
rGG <- function(n, mu=1, sigma=0.5, nu=1)
  {
    if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))     
    n <- ceiling(n)
    p <- runif(n)
    r <- qGG(p,mu=mu,sigma=sigma,nu=nu)
    r
  }
#--------------------------------------------------------------
