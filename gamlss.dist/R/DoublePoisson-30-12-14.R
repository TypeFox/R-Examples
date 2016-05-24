#-------------------------------------------------------------------------------
# Bob Rigby  Mikis Stasinopoulos and Marco Enea
# last change Monday, 30  Dec 2014
# the Double Poisson distribution 
#-------------------------------------------------------------------------------
# TO DO
# To put in gamlss.dist
# ------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# numerical derivatives for logC 
DPO <- function (mu.link = "log", sigma.link = "log") 
{
  mstats <- checklink("mu.link", "Double Poisson", substitute(mu.link), 
                      c("inverse", "log", "identity", "sqrt"))
  dstats <- checklink("sigma.link", "Double Poisson", substitute(sigma.link), 
                      c("inverse", "log", "identity", "sqrt"))
  structure(
    list(    family = c("DPO", "Double Poisson"),
         parameters = list(mu = TRUE,sigma = TRUE), 
              nopar = 2, 
               type = "Discrete", 
            mu.link = as.character(substitute(mu.link)),
         sigma.link = as.character(substitute(sigma.link)), 
         mu.linkfun = mstats$linkfun, 
      sigma.linkfun = dstats$linkfun, 
         mu.linkinv = mstats$linkinv, 
      sigma.linkinv = dstats$linkinv,
              mu.dr = mstats$mu.eta, 
           sigma.dr = dstats$mu.eta, 
               dldm = function(y, mu, sigma)
                 {
                  -(1/sigma)+y/(mu*sigma)+
                   as.vector(attr(numeric.deriv(get_C(y, mu, sigma), "mu"),"gradient"))
                 }, 
            d2ldm2 = function(y, mu, sigma) 
                 {
               dldm <- -(1/sigma)+y/(mu*sigma) +
                   as.vector(attr(numeric.deriv(get_C(y, mu, sigma), "mu"),"gradient"))
           d2ldm2 <- -dldm^2
           d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)    
           d2ldm2  
         },
         dldd = function(y,mu,sigma)
         { -(1/(2*sigma)) + (mu/(sigma^2)) -
             (y*log(mu)/(sigma^2)) - (y/(sigma^2)) +
             (y*ifelse(y==0,1,log(y)))/(sigma^2)  +
             as.vector(attr(numeric.deriv(get_C(y, mu, sigma), "sigma"),"gradient"))
         },
         d2ldd2 = function(y,mu,sigma) {# eval.parent(quote(-dldp*dldp))
           dldd <-       -(0.5/sigma) + (mu/(sigma^2)) -
             (y*log(mu)/(sigma^2)) - (y/(sigma^2)) +
             (y*ifelse(y==0,1,log(y)))/(sigma^2)  +
             as.vector(attr(numeric.deriv(get_C(y, mu, sigma), "sigma"),"gradient"))
           d2ldd2 <- -dldd^2
           d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)       
           d2ldd2               
         }, #change this
         d2ldmdd = function(y) rep(0,length(y)),     
         G.dev.incr  = function(y,mu,sigma,...) -2*dDPO(y, mu = mu, sigma = sigma, log = TRUE), 
         rqres = expression(
           rqres(pfun="pDPO", type="Discrete", ymin=0, y=y, mu=mu, sigma=sigma)
         ), 
         mu.initial = expression(mu <-  (y+mean(y))/2),
         sigma.initial = expression(
           sigma <- rep(1,length(y))),
         mu.valid = function(mu) all(mu > 0) , 
         sigma.valid = function(sigma)  all(sigma > 0), 
         y.valid = function(y)  all(y >= 0)
    ),
    class = c("gamlss.family","family"))
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# approximate derivative for logC 
DPO1 <- function (mu.link = "log", sigma.link = "log") 
{
    mstats <- checklink("mu.link", "Double Poisson", substitute(mu.link), 
                        c("inverse", "log", "identity", "sqrt"))
    dstats <- checklink("sigma.link", "Double Poisson", substitute(sigma.link), 
                        c("inverse", "log", "identity", "sqrt"))
    structure(
          list(family = c("DPO1", "Double Poisson type 1"),
           parameters = list(mu = TRUE,sigma = TRUE), 
                nopar = 2, 
                 type = "Discrete", 
              mu.link = as.character(substitute(mu.link)),
           sigma.link = as.character(substitute(sigma.link)), 
           mu.linkfun = mstats$linkfun, 
        sigma.linkfun = dstats$linkfun, 
           mu.linkinv = mstats$linkinv, 
        sigma.linkinv = dstats$linkinv,
                mu.dr = mstats$mu.eta, 
             sigma.dr = dstats$mu.eta, 
                 dldm = function(y, mu, sigma)
                   {
                   -(1/sigma)+y/(mu*sigma)+      # approx:
                    (1+ (sigma-1)/(12*mu) + (sigma-1)*sigma/(12*mu^2))*
                    (24*mu * (12*mu^2 + (sigma-1)*mu - sigma^2 + sigma) - 
                    12*mu^2*(24*mu +sigma-1 ))/(12*mu^2+(sigma-1)*mu-sigma^2+sigma)^2
                    }, 
               d2ldm2 = function(y, mu, sigma) 
                    {
                 dldm <- (1/sigma)+y/(mu*sigma)+      # approx:
                   (1+ (sigma-1)/(12*mu) + (sigma-1)*sigma/(12*mu^2))*
                   (24*mu * (12*mu^2 + (sigma-1)*mu - sigma^2 + sigma) - 
                      12*mu^2*(24*mu +sigma-1 ))/(12*mu^2+(sigma-1)*mu-sigma^2+sigma)^2 
                   d2ldm2 <- -dldm^2
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)    
                   d2ldm2  
               },
                 dldd = function(y,mu,sigma)
                                 { -(0.5/sigma) + (mu/(sigma^2)) -
                                     (y*log(mu)/(sigma^2)) - (y/(sigma^2)) +
                                     (y*ifelse(y==0,1,log(y)))/(sigma^2)  + # approx:
                                     (1+(sigma-1)/(12*mu)+(sigma-1)*sigma/(12*mu^2))*
                                     (-12*mu^2*(mu-2*sigma+1))/
                                     (12*mu^2+(sigma-1)*mu-sigma^2+sigma)^2
                                     #numeric.deriv(gC(mu,sigma), "sigma", delta=0.01)
                                     
                                    # numeric.deriv(get_C(y, mu, sigma), "sigma", delta=0.01)
                                   },
               d2ldd2 = function(y,mu,sigma) {# eval.parent(quote(-dldp*dldp))
                     dldd <-       -(0.5/sigma) + (mu/(sigma^2)) -
                                    (y*log(mu)/(sigma^2)) - (y/(sigma^2)) +
                       (y*ifelse(y==0,1,log(y)))/(sigma^2)  +# approx:
                                    (1+(sigma-1)/(12*mu)+(sigma-1)*sigma/(12*mu^2))*
                                    (-12*mu^2*(mu-2*sigma+1))/
                                    (12*mu^2+(sigma-1)*mu-sigma^2+sigma)^2
                                    #numeric.deriv(gC(mu,sigma), "sigma", delta=0.01)
                                     #altenatively
                                    #numeric.deriv(get_C(y, mu, sigma), "sigma", delta=0.01)
                  d2ldd2 <- -dldd^2
                  d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)       
                  d2ldd2               
                                   }, #change this
              d2ldmdd = function(y) rep(0,length(y)),     
          G.dev.incr  = function(y,mu,sigma,...) -2*dDPO(y, mu = mu, sigma = sigma, log = TRUE), 
               rqres = expression(
                          rqres(pfun="pDPO", type="Discrete", ymin=0, y=y, mu=mu, sigma=sigma)
                                 ), 
            mu.initial = expression(mu <-  (y+mean(y))/2),
         sigma.initial = expression(
                      sigma <- rep(1,length(y))),
              mu.valid = function(mu) all(mu > 0) , 
           sigma.valid = function(sigma)  all(sigma > 0), 
               y.valid = function(y)  all(y >= 0)
          ),
            class = c("gamlss.family","family"))
}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#gC <-  function(mu, sigma) 1/(1+(sigma-1)/(12*mu)+(sigma-1)*sigma/(12*mu^2))
get_C <- function(x, mu, sigma)
{ 
  maxV <- max(max(x)*3, 500)
  y <- 0:maxV
  ly <- max(length(x),length(mu),length(sigma)) 
  theC <- .C("dDPOgetC5_C",as.double(mu),as.double(sigma),as.integer(ly),as.integer(maxV+1),ans=double(ly))$ans
  log(theC)          
}
#-------------------------------------------------------------------------------
# the d function
#-------------------------------------------------------------------------------
dDPO <- function(x, mu = 1, sigma = 1, log = FALSE)
{ 
  if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
  if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
  if (any(x < 0) )  stop(paste("x must be >=0", "\n", ""))  
  ly <- max(length(x),length(mu),length(sigma)) 
  x <- rep(x, length = ly)      
  sigma <- rep(sigma, length = ly)
  mu <- rep(mu, length = ly) 
  maxV <- max(max(x)*3,500)
  y <- 0:maxV
  #theC <- .C("dDPOgetC5_C",as.double(mu),as.double(sigma),as.integer(ly),as.integer(maxV+1),ans=double(ly))$ans
  logofx <- ifelse(x==0,1,log(x))
  lh <- -0.5*log(sigma)-(mu/sigma)-lgamma(x+1)+x*logofx-x+
    (x*log(mu))/sigma+x/sigma-(x*logofx)/sigma+get_C(x, mu, sigma)#log(theC)        
  if(log==FALSE) fy <- exp(lh) else fy <- lh 
  fy
}

#-------------------------------------------------------------------------------
#  The p function  
#-------------------------------------------------------------------------------
pDPO<-function(q, mu = 1, sigma = 1, lower.tail = TRUE, log.p = FALSE)
{ 
  if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
  if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
  if (any(q < 0) )  stop(paste("q must be >=0", "\n", ""))  
  ly <- max(length(q),length(mu),length(sigma)) 
  q <- rep(q, length = ly)      
  sigma <- rep(sigma, length = ly)
  mu <- rep(mu, length = ly) 
  maxV <- max(max(q)*3,500)
  y <- 0:maxV
  den <- unlist(lapply(1:ly,function(x)
    .C("dDPOgetC5_C",as.double(mu[x]),as.double(sigma[x]),as.integer(1),as.integer(q[x]+1),ans=double(1))$ans))
  num <- .C("dDPOgetC5_C",as.double(mu),as.double(sigma),as.integer(ly),as.integer(maxV+1),ans=double(ly))$ans;
  cdf <- num/den
  cdf <- if(lower.tail==TRUE) cdf else 1-cdf
  cdf <- if(log.p==FALSE) cdf else log(cdf)                                                                    
  cdf
}
#-------------------------------------------------------------------------------
# the q function
#-------------------------------------------------------------------------------
qDPO <- function(p, mu=1, sigma=1,  lower.tail = TRUE, log.p = FALSE,  
                 max.value = 10000)
{      
  if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
  if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
  if (any(p < 0) | any(p > 1.0001))  stop(paste("p must be between 0 and 1", "\n", "")) 
  if (log.p==TRUE) p <- exp(p) else p <- p
  if (lower.tail==TRUE) p <- p else p <- 1-p    
  ly <- length(p)                                                       
  QQQ <- rep(0,ly)                         
  nsigma <- rep(sigma, length = ly)
  nmu <- rep(mu, length = ly)                
  for (i in seq(along=p))                                                          
  {
    cumpro <- 0                                                                         
    if (p[i]+0.000000001 >= 1) QQQ[i] <- Inf
    else  
    {  
      for (j in seq(from = 0, to = max.value))
      {
        cumpro <-  pDPO(j, mu = nmu[i], sigma = nsigma[i], log.p = FALSE) 
        # else  cumpro+dSICHEL(j, mu = nmu[i], sigma = nsigma[i], nu = nnu[i], log = FALSE)# the above is faster 
        QQQ[i] <- j 
        if  (p[i] <= cumpro ) break 
      } 
    }
  }          
  QQQ   
} 

#-------------------------------------------------------------------------------
# the r function 
#-------------------------------------------------------------------------------
rDPO <- function(n, mu=1, sigma=1, max.value = 10000)
{ 
  if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
  if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", ""))  
  if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
  n <- ceiling(n)
  p <- runif(n)
  r <- qDPO(p, mu=mu, sigma=sigma, max.value = max.value )
  r
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------



#dyn.load("/Volumes/Data/Users/stasinom/Dropbox/Cpp/dDPO/dDPOgetC5_C.so")
#is.loaded("dDPOgetC5_C")




