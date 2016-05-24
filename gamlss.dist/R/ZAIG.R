# The inflated inverse Gaussian distribution
# created by Bob Rigby and Mikis Stasinopoulos, suggested by Gillian Heller 
# ---------------------------------------------------------------------------------------
ZAIG <- function (mu.link ="log", sigma.link="log", nu.link ="logit")
{
    mstats <- checklink("mu.link", "ZAIG", substitute(mu.link), c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "ZAIG", substitute(sigma.link), c("inverse", "log", "identity", "own"))
    vstats <- checklink("nu.link", "ZAIG", substitute(nu.link),  c("logit", "probit", "cloglog", "cauchit", "log", "own"))
    structure(
          list(family = c("ZAIG", "Zero adjusted IG"),
           parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE), 
                nopar = 3, 
                 type = "Mixed",
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
                 dldm = function(y,mu,sigma) ifelse(y==0,0,(y-mu)/((sigma^2)*(mu^3))),
               d2ldm2 = function(y,mu,sigma) ifelse(y==0,0,-1/((mu^3)*(sigma^2))),
                 dldd = function(y,mu,sigma) ifelse(y==0,0,(-1/sigma) +((y-mu)^2)/(y*(sigma^3)*(mu^2))),
               d2ldd2 = function(y,sigma) ifelse(y==0,0,-2/(sigma^2)),
                 dldv = function(y,nu) ifelse(y==0,1/nu,-1/(1-nu)),
               d2ldv2 = function(nu)  -1/(nu*(1-nu)) ,
              d2ldmdd = function(y)  rep(0,length(y)),
              d2ldmdv = function(y)  rep(0,length(y)),
              d2ldddv = function(y)   rep(0,length(y)),
          G.dev.incr  = function(y,mu,sigma,nu,...) 
                                 -2*dZAIG(y,mu,sigma,nu,log=TRUE),                           
                rqres = expression(rqres(pfun="pZAIG", type="Mixed",  mass.p=0,  
                             prob.mp=nu, y=y, mu=mu, sigma=sigma, nu=nu)),
           mu.initial =  expression(mu <- (y+mean(y))/2), 
         sigma.initial =  expression(sigma <- rep(.5,length(y))),
            nu.initial =  expression(   nu <- rep(0.5, length(y))), 
              mu.valid = function(mu) TRUE , 
           sigma.valid = function(sigma)  all(sigma > 0),
              nu.valid = function(nu) all(nu > 0) && all(nu < 1), 
               y.valid = function(y)  all(y>=0)
          ),
            class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------
dZAIG<-function(x, mu=1, sigma=1, nu=.1, log=FALSE)
 {        if (any(mu < 0))  stop(paste("mu must be positive", "\n", "")) 
          if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(x < 0))  stop(paste("x must be positive", "\n", ""))  
 log.lik <- ifelse(x==0, log(nu), log(1-nu)+(-0.5*log(2*pi)-log(sigma)-(3/2)*log(x)-((x-mu)^2)/(2*sigma^2*(mu^2)*x)))
     if(log==FALSE) fy  <- exp(log.lik) else fy <- log.lik
      fy 
  }
#----------------------------------------------------------------------------------------
pZAIG <- function(q, mu=1, sigma=1, nu=0.1, lower.tail = TRUE, log.p = FALSE)
  {       if (any(mu < 0))  stop(paste("mu must be positive", "\n", "")) 
          if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(q < 0))  stop(paste("y must be positive", "\n", ""))  
    cdf1 <- pnorm(((q/mu)-1)/(sigma*sqrt(q))) 
    cdf2 <- exp(2/(mu*sigma^2))*pnorm((-((q/mu)+1))/(sigma*sqrt(q)))
     cdf <- cdf1+ cdf2
     ## the problem with this approximation is that it is not working with 
     ## small sigmas and produce NA's. So here it is a solution
    if (any(is.na(cdf)))
      {
       index <- seq(along=q)[is.na(cdf)]
       for (i in index)
          {
        cdf[i] <- integrate(function(x) 
                 dIG(x, mu = mu[i], sigma = sigma[i], log=FALSE), 0.001, q[i] )$value
          }    
      }
     cdf <- ifelse((q==0), nu, nu+(1-nu)*cdf)    
    if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf 
    if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf) 
    cdf
   }
#---------------------------------------------------------------------------------------- 
qZAIG <- function(p, mu=1, sigma=1,  nu=0.1, lower.tail = TRUE, log.p = FALSE, 
                upper.limit = mu+10*sqrt(sigma^2*mu^3))
  { 
    #---function--------------------------------------------   
       h1 <- function(q)
       { 
     pZAIG(q , mu = mu[i], sigma = sigma[i], nu = nu[i]) - p[i]  
       }
     #-----------------------------------------------------------------
     #-----------------------------------------------------------------
    if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))
    if (any(nu < 0)|any(nu > 1))  stop(paste("nu must be between 0 and 1", "\n", ""))     
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (lower.tail==TRUE) p <- p else p <- 1-p
    if (any(p < 0)|any(p >= 1))  stop(paste("p must be between 0 and 1", "\n", ""))     
    #     lp <- length(p) 
         lp <- max(length(p),length(mu),length(sigma),length(nu)) 
          p <- rep(p, length = lp)
      sigma <- rep(sigma, length = lp)
         mu <- rep(mu, length = lp)
         nu <- rep(nu, length = lp)
      upper <- rep(upper.limit, length = lp )
      lower <- rep(0, length = lp )
          q <- rep(0,lp)    
         for (i in seq(along=p))
          {
           q[i] <- if (nu[i]>=p[i]) 0
                   else  uniroot(h1, c(lower[i], upper[i]))$root                
          if (q[i]>=upper[i]) warning("q is at the upper limit, increase the upper.limit")
         # if (q[i]<=lower[i]) warning("q is at the lower limit, decrease the lower.limit")
          }                                                                               
    q
   }
#-----------------------------------------------------------------------------------------
rZAIG <- function(n, mu=1, sigma=1, nu=0.1, ...)
  { 
    if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(nu < 0)|any(nu > 1))  stop(paste("nu must be between 0 and 1", "\n", ""))   
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
    n <- ceiling(n)
    p <- runif(n)
    r <- qZAIG(p,mu=mu,sigma=sigma, nu = nu, ...)
    r
  }
   
#----------------------------------------------------------------------------------------
plotZAIG <- function( mu =5 , sigma=1, nu = 0.1, from = 0, to=10, n = 101, ...)
 {
  y = seq(from=0.001, to=to, length.out=n )
  pdf<- dZAIG(y, mu = mu ,sigma = sigma, nu = nu) 
  pr0<-c(dZAIG(0, mu=mu ,sigma=sigma,  nu=nu))
  po<-c(0)
  plot(pdf~y, main="Zero adjusted IG", ylim=c(0,max(pdf,pr0)), type="l")
  points(po,pr0,type="h")
  points(po,pr0,type="p", col="blue")
 }
#----------------------------------------------------------------------------------------
meanZAIG <- function(obj)
  {
  if ( obj$family[1]!="ZAIG") stop("the object do not have a ZAIG distribution")
  meanofY<-(1-fitted(obj,"nu"))*fitted(obj,"mu")
  meanofY
  }
#---------------------------------------------------------------------------------------- 
