# Amended on 19/10/06 the Reverse Generalized Extreme distribution
RGE <- function (mu.link="identity", sigma.link="log", nu.link ="log") 
{
    mstats <- checklink("mu.link", "RGE", substitute(mu.link), 
                         c("1/mu^2", "log", "identity"))
    dstats <- checklink("sigma.link", "RGE", substitute(sigma.link), 
                         c("inverse", "log", "identity"))
    vstats <- checklink("nu.link", "RGE",substitute(nu.link), 
                         c("inverse", "log", "identity"))  
 
     structure(
          list(family = c("RGE", "Reverse-Generalized-Extreme"),
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
                                     toto <- nu*(y-mu)/sigma
                                     dldm <-   (nu-1)/((1+toto) *sigma)
                                     dldm <- dldm + (1/sigma)*((1+toto)^((1/nu)-1))
                                     dldm
                                    },
               d2ldm2 = function(y,mu,sigma,nu) {
                                     toto <- nu*(y-mu)/sigma
                                     dldm <-   (nu-1)/((1+toto) *sigma)
                                     dldm <- dldm + (1/sigma)*((1+toto)^((1/nu)-1))
                                 d2ldm2 <-  -dldm*dldm
                                 d2ldm2
                                    },
               dldd = function(y,mu,sigma,nu) {
                                     toto <- nu*(y-mu)/sigma
                                     dldd <-  (-1-((1/nu)-1)*(toto)/(1+toto)) 
                                     dldd <- dldd +((1+toto)^((1/nu)-1))*toto/nu
                                     dldd <- dldd/sigma
                                     dldd
                                    },
               d2ldd2 = function(y,mu,sigma,nu) {
                                     toto <- nu*(y-mu)/sigma
                                     dldd <-  (-1-((1/nu)-1)*(toto)/(1+toto)) 
                                     dldd <- dldd +((1+toto)^((1/nu)-1))*toto/nu
                                     dldd <- dldd/sigma
                                 d2ldd2 <-  -dldd*dldd
                                 d2ldd2
                                       },
                 dldv = function(y,mu,sigma,nu) {
                                    toto <- nu*(y-mu)/sigma 
                                     dldv <- -(ifelse((1+toto) <= 0,0,log(1+toto)) +(nu-1)*toto/(1+toto))
                                     dldv <- dldv - ((1+toto)^(1/nu))*((toto/(1+toto)) -
                                                                  ifelse((1+toto) <= 0,0,log(1+toto)))
                                     dldv <- dldv/(nu*nu)
                                     dldv
                                    },
               d2ldv2 = function(y,mu,sigma,nu)  {
                                    toto <- nu*(y-mu)/sigma 
                                     dldv <- -(ifelse((1+toto) <= 0,0,log(1+toto)) +(nu-1)*toto/(1+toto))
                                     dldv <- dldv - ((1+toto)^(1/nu))*((toto/(1+toto)) -
                                                                  ifelse((1+toto) <= 0,0,log(1+toto)))
                                     dldv <- dldv/(nu*nu)
                                    d2ldv2 <-  -dldv*dldv
                                    d2ldv2
                                      },

              d2ldmdd = function(y,mu,sigma,nu)  rep(0,length(y)),
              d2ldmdv = function(y,mu,sigma,nu)  rep(0,length(y)),
              d2ldddv = function(y,mu,sigma,nu)  rep(0,length(y)),
          G.dev.incr  = function(y,mu,sigma,nu,...) {
                                   -2*dRGE(y,mu=mu,sigma=sigma,nu=nu,log=TRUE)
                                                    }, 
                rqres = expression(rqres(pfun="pRGE", type="Continuous", y=y, mu=mu, 
                                              sigma=sigma, nu=nu) ),
             mu.initial = expression( mu <- y+0.45*sd(y)), #rep(mean(y)+0.45*sd(y),length(y))), 
          sigma.initial = expression( sigma <- rep(0.78*sd(y),length(y))),
             nu.initial = expression( nu <- rep(0.1,length(y))), 
              mu.valid = function(mu) TRUE , 
           sigma.valid = function(sigma)  all(sigma > 0),
              nu.valid = function(nu) all(nu > 0) , 
               y.valid = function(y) TRUE     
          ),
            class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------
#--------------------------------------------------------------
#----------------------------------------------------------------------------------------
dRGE <- function(x, mu=1, sigma=0.1, nu=1,  log = FALSE)
 {
      #   if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(nu <= 0))  stop(paste("nu must be positive", "\n", "")) 
          if (any(x < mu-(sigma/nu)))  stop(paste("x must greater than mu-sigma/nu", "\n", ""))  
         y.sca <- nu*(x-mu)/sigma 
       loglik1 <- (-log(sigma)+((1/nu)-1)*ifelse((1+y.sca) <= 0,0,log(1+y.sca)))-(1+y.sca)^(1/nu)
       loglik2 <- -log(sigma)+(x-mu)/sigma-exp((x-mu)/sigma)
       if(length(nu)>1)  loglik <- ifelse(abs(nu)<0.001,loglik2,loglik1)
       else   loglik <- if(abs(nu)<0.001) loglik2 else loglik1 
       ft <- if(log==FALSE) exp(loglik) else loglik 
       ft
  }    
#--------------------------------------------------------------  
pRGE <- function(q, mu=1, sigma=0.1, nu=1,  lower.tail = TRUE, log.p = FALSE)
 {  
      #   if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
       if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
       if (any(nu <= 0))  stop(paste("nu must be positive", "\n", "")) 
      #   if (any(q < 0))  stop(paste("q must be positive", "\n", ""))  
          y.sca <- nu*(q-mu)/sigma 
           FYy1 <-  1- exp(-((1+y.sca)^(1/nu)))
           FYy2 <-   1-exp(-exp((q-mu)/(sigma)))
       if(length(nu)>1)  FYy <- ifelse(abs(nu)<0.001,FYy2,FYy1)
       else   FYy <- if(abs(nu)<0.001) FYy2 else FYy1 
       if(lower.tail==TRUE) FYy  <- FYy else  FYy <- 1-FYy 
       if(log.p==FALSE) FYy  <- FYy else  FYy<- log(FYy) 
       FYy     
 }
#--------------------------------------------------------------
qRGE <- function(p, mu=1, sigma=0.1, nu=1,  lower.tail = TRUE, log.p = FALSE )
 { 
      #   if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(nu <= 0))  stop(paste("nu must be positive", "\n", "")) 
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (lower.tail==TRUE) p <- p else p <- 1-p
    if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
           ya1 <-   mu + (sigma/nu)*(((-log(1-p))^nu)-1)
           ya2 <-   mu + sigma*log(-log(1-p))
    if(length(nu)>1)  ya <-  ifelse(abs(nu)<0.001,ya2,ya1)
    else   ya <- if(abs(nu)<0.001) ya2 else ya1 
    ya
 }
#--------------------------------------------------------------
rRGE <- function(n, mu=1, sigma=0.1, nu=1)
  {
      #   if (any(mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(nu <= 0))  stop(paste("nu must be positive", "\n", "")) 
    if (any(n <= 0))  stop(paste("n must be a positive inteRGEr", "\n", ""))    
    
    n <- ceiling(n)
    p <- runif(n)
    r <- qRGE(p,mu=mu,sigma=sigma,nu=nu)
    r
  }
#--------------------------------------------------------------
