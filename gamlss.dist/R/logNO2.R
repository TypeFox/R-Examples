# LOGNO2 with Re-parameretised log(mu) Mikis 1-7-12
#---------------------------------------------------------
# probably OK needs further testing 
#---------------------------------------------------------
#source('/Volumes/Data/Users/stasinom/Documents/gamlss/projects/DISTRIBUTIONS/testingDistributions/testContDist.R')
# testContDist("LOGNO2", y.range=c(0,Inf), mu.range = c(0,Inf),sigma.range=c(0, Inf), mu.val = c(1,.5,10), sigma.val=c(1, 2, 5))
LOGNO2 <-function (mu.link ="log", sigma.link="log") 
{
    mstats <- checklink("mu.link", "LOGNO2", substitute(mu.link), c("inverse", "log", "identity"))
    dstats <- checklink("sigma.link", "LOGNO2", substitute(sigma.link), c("inverse", "log", "identity"))        
    structure(
          list(family = c("LOGNO2", "Log Normal 2"),
           parameters = list(mu=TRUE,sigma=TRUE),
                nopar = 2, 
                 type = "Continuous",
              mu.link = as.character(substitute(mu.link)), 
           sigma.link = as.character(substitute(sigma.link)), 
           mu.linkfun = mstats$linkfun, 
        sigma.linkfun = dstats$linkfun, 
           mu.linkinv = mstats$linkinv, 
        sigma.linkinv = dstats$linkinv,
                mu.dr = mstats$mu.eta, 
             sigma.dr = dstats$mu.eta, 
                 dldm = function(y, mu, sigma){
                 	      dldm <- (1/sigma^2)*(log(y)-log(mu))/mu #(log(y)-log(mu))/((sigma^2)*mu) 
                          dldm
                           },
               d2ldm2 = function(mu, sigma)-(1/(sigma^2)*(mu^2)),
                 dldd = function(y,mu,sigma) {
                 	           dldd <- (1/(sigma^3))*((log(y)-log(mu))^2-sigma^2)
                              dldd}, 
               d2ldd2 = function(sigma) -(2/(sigma^2)),
              d2ldmdd = function(y) rep(0,length(y)),
          G.dev.incr  = function(y,mu,sigma,...)  -2*dLOGNO2(y,mu,sigma,log=TRUE),                         
                rqres = expression(rqres(pfun="pLOGNO2", type="Continuous", y=y, mu=mu, sigma=sigma)),
           mu.initial = expression({   mu <- exp((log(y)+mean(log(y)))/2 ) }),
        sigma.initial = expression({sigma <- rep(sd(log(y)),length(y))}), 
             mu.valid = function(mu) all (mu > 0) , 
          sigma.valid = function(sigma) all(sigma > 0), 
              y.valid = function(y)  all(y>0)
          ),
            class = c("gamlss.family","family")
          )
}

#---------------------------------------------
#---------------------------------------------
dLOGNO2<-function(x, mu=1, sigma=1, log=FALSE)
 { 
         if ( any(mu <= 0))    stop(paste("mu must be positive", "\n", ""))
		 if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
           z <- log(x)
    loglik <- dnorm(z, mean=log(mu), sd=sigma, log=TRUE)-z
    if(log==FALSE) fy  <- exp(loglik) else fy <- loglik 
    fy
  }
#---------------------------------------------
#---------------------------------------------
pLOGNO2 <- function(q, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE)
  {     
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    qz <- log(q)
    cdf <- pnorm(qz, mean=log(mu), sd=sigma, lower.tail = lower.tail, log.p = log.p)
    cdf
   }
#---------------------------------------------
#---------------------------------------------
qLOGNO2 <- function(p, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE)
  { if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", "")) 
    qz <- qnorm(p, mean=log(mu), sd=sigma, lower.tail = lower.tail )
    qy=exp(qz)
    qy
   }
#---------------------------------------------
#---------------------------------------------
rLOGNO2 <- function(n, mu=1, sigma=1)
  { 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    rz <- rnorm(n, mean=log(mu), sd=sigma)
    ry=exp(rz)
    ry
  }
#---------------------------------------------
#---------------------------------------------
