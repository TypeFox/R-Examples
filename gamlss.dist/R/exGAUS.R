# the ex-Gaussian distribution
# Mikis Stasinopoulos and Bob Rigby (suggested by Jonathan Williams) 
# 28_11_07
exGAUS <- function (mu.link="identity", sigma.link="log", nu.link ="log") 
{
    mstats <- checklink("mu.link", "ex-Gaussian", substitute(mu.link), 
                         c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "ex-Gaussian", substitute(sigma.link), #
                         c("inverse", "log", "identity", "own"))
    vstats <- checklink("nu.link", "ex-Gaussian",substitute(nu.link), 
                        c("logshifted", "log", "identity", "own"))
    
    structure(
          list(family = c("exGAUS", "ex-Gaussian"),
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
        z <- y-mu-((sigma^2)/nu)
     pphi <- (dnorm(z/sigma)/(pnorm(z/sigma)))
     dldm <- 1/nu-(1/sigma)*pphi
     dldm
                                    },
               d2ldm2 = function(y,mu,sigma,nu) {
        z <- y-mu-((sigma^2)/nu)
     pphi <- (dnorm(z/sigma)/(pnorm(z/sigma)))
     dldm <- 1/nu-(1/sigma)*pphi
   d2ldm2 <- -dldm*dldm
      d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15) 
      d2ldm2
                                    },
                 dldd = function(y,mu,sigma,nu) {
      z <- y-mu-((sigma^2)/nu) 
   pphi <- (dnorm(z/sigma)/(pnorm(z/sigma)))
   dldd <- (sigma/(nu^2))-((z/sigma^2)+(2/nu))*pphi
   dldd
                                    },
               d2ldd2 = function(y,mu,sigma,nu) {
      z <- y-mu-((sigma^2)/nu) 
   pphi <- (dnorm(z/sigma)/(pnorm(z/sigma)))
   dldd <- (sigma/(nu^2))-((z/sigma^2)+(2/nu))*pphi
 d2ldd2 <- -dldd*dldd
   d2ldd2 <- ifelse(d2ldd2 < -1e-10, d2ldd2,-1e-10)  
   d2ldd2
                                    },
                 dldv = function(y,mu,sigma,nu) {
       z <- y-mu-((sigma^2)/nu)
    pphi <- (dnorm(z/sigma)/(pnorm(z/sigma))) 
    dldv <- -(1/nu)+(z/nu^2)+(sigma/nu^2)*pphi
    dldv
                                    },
               d2ldv2 = function(y,mu,sigma,nu) {
       z <- y-mu-((sigma^2)/nu)
    pphi <- (dnorm(z/sigma)/(pnorm(z/sigma))) 
    dldv <- -(1/nu)+(z/nu^2)+(sigma/nu^2)*pphi
  d2ldv2 <- -dldv*dldv
   d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)  
   d2ldv2
                                    },
              d2ldmdd = function(y,mu,sigma,nu) {
         z <- y-mu-((sigma^2)/nu)
      pphi <- (dnorm(z/sigma)/(pnorm(z/sigma)))
      dldm <- 1/nu-(1/sigma)*pphi
      dldd <- (sigma/(nu^2))-((z/sigma^2)+(2/nu))*pphi
   d2ldmdd <- -dldm *dldd          
   d2ldmdd                          },
              d2ldmdv = function(y,mu,sigma,nu) {
          z <- y-mu-((sigma^2)/nu)
       pphi <- (dnorm(z/sigma)/(pnorm(z/sigma))) 
       dldm <- 1/nu-(1/sigma)*pphi
       dldv <- -(1/nu)+(z/nu^2)+(sigma/nu^2)*pphi
    d2ldmdv <- -dldm *dldv  
    d2ldmdv         
                                    },
              d2ldddv = function(y,mu,sigma,nu) {
       z <- y-mu-((sigma^2)/nu)
    pphi <- (dnorm(z/sigma)/(pnorm(z/sigma))) 
   dldd <- (sigma/(nu^2))-((z/sigma^2)+(2/nu))*pphi
    dldv <- -(1/nu)+(z/nu^2)+(sigma/nu^2)*pphi
     d2ldddv <- -dldd *dldv
     d2ldddv         
                                   },
          G.dev.incr  = function(y,mu,sigma,nu,...) {
        -2*dexGAUS(y,mu=mu,sigma=sigma,nu=nu,log=TRUE)
                                                    }, 
             rqres = expression(rqres(pfun="pexGAUS", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)),
        mu.initial = expression( mu <- (y+mean(y))/2), 
     sigma.initial = expression( sigma <- rep(sd(y)/2,length(y))), 
        nu.initial = expression( nu <- rep(max(2*(mean(y)-median(y)),0.1), length(y))), # 2*max((mean(y)-median(y),0.1) # (mean(y)+sd(y))/4
          mu.valid = function(mu) TRUE , 
       sigma.valid = function(sigma)  all(sigma > 0),
          nu.valid = function(nu) all(nu > 0), 
           y.valid = function(y) TRUE
          ),
            class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------
dexGAUS<-function(x, mu=5, sigma=1, nu=1, log=FALSE)
  { 
   if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
   if (any(nu <= 0) )  stop(paste("nu must be greater than 0 ", "\n", "")) 
    ly <- length(x)       
 sigma <- rep(sigma, length = ly)
    mu <- rep(mu, length = ly)   
    nu <- rep(nu, length = ly) 
     z <- x-mu-((sigma^2)/nu)
logfy <- ifelse(nu>0.05*sigma,  
          -log(nu)-(z+(sigma^2/(2*nu)))/nu+log(pnorm(z/sigma)),
           dnorm(x, mean=mu, sd=sigma, log=TRUE)
               )
#logfy <- ifelse(nu > 0.05, logfy, dnorm(x,mean=mu,sd=sigma, log=TRUE))
#logfy <- -log(nu)+((mu-x)/nu)+(sigma^2/(2*nu^2))+log(pnorm(((x-mu)/sigma)-sigma/nu))
  if(log==FALSE) fy <- exp(logfy) else fy <- logfy
  fy
  }
#----------------------------------------------------------------------------------------
pexGAUS<-function(q, mu = 5, sigma = 1,  nu=1,  lower.tail = TRUE, log.p = FALSE)
  { 
   if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
   if (any(nu <= 0) )  stop(paste("nu must be greater than 0 ", "\n", "")) 
    ly <- length(q)       
 sigma <- rep(sigma, length = ly)
    mu <- rep(mu, length = ly)   
    nu <- rep(nu, length = ly) 
 index <- seq(along=q)
     z <- q-mu-((sigma^2)/nu)
   cdf <- ifelse(nu>0.05*sigma,  
  pnorm((q-mu)/sigma)-pnorm(z/sigma)*exp(((mu+(sigma^2/nu))^2-(mu^2)-2*q*((sigma^2)/nu))/(2*sigma^2)),
          pnorm(q, mean=mu, sd=sigma)
               )
    if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf 
    if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf) 
    cdf
 }
#----------------------------------------------------------------------------------------
qexGAUS <- function(p, mu = 5, sigma = 1, nu = 1,  lower.tail = TRUE,  log.p = FALSE) 
  { 
    #---functions--------------------------------------------   
       h1 <- function(q)
       { 
     pexGAUS(q , mu = mu[i], sigma = sigma[i], nu = nu[i]) - p[i]  #???????????????
       }
       h <- function(q)
       { 
     pexGAUS(q , mu = mu[i], sigma = sigma[i], nu = nu[i])  #??????????????????? 
       }
     #-----------------------------------------------------------------
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))      
    if (any(nu <= 0))  stop(paste("nu must be positive", "\n", "")) 
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
           interval <-  c(mu[i]-sigma[i], mu[i])
           j <-2
           while (h(interval[1]) > p[i]) 
              {interval[1]<- mu[i]-j*sigma[i]
              j<-j+1 
              }
           }
        q[i] <- uniroot(h1, interval)$root
        #interval <- c(.Machine$double.xmin, 20)
         }
    q
   }
#----------------------------------------------------------------------------------------
rexGAUS <- function(n, mu=5, sigma=1, nu=1, ...)
  { 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(nu <= 0))  stop(paste("nu must be positive", "\n", "")) 
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
    n <- ceiling(n)
    p <- runif(n)
    r <- qexGAUS(p,mu=mu,sigma=sigma, nu=nu, ...)
    r
  }
#----------------------------------------------------------------------------------------
