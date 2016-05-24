#----------------------------------------------------------------------------------------
# MS + BR last change Thursday, April 13, 2006
TF <- function (mu.link="identity", sigma.link="log", nu.link ="log")
{
    mstats <- checklink("mu.link", "t Family", substitute(mu.link), c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "t Family", substitute(sigma.link), c("inverse", "log", "identity", "own"))
    vstats <- checklink("nu.link", "t Family", substitute(nu.link), c("inverse", "log", "identity", "own"))
    
    structure(
          list(family = c("TF", "t Family"),
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
                                       s2 <- sigma^2
                                     dsq <- ((y-mu)^2)/s2
                                   omega <- (nu+1)/(nu+dsq) 
                                    dldm <- (omega*(y-mu))/s2
                                    dldm
                                    },
               d2ldm2 = function(sigma,nu) -(nu+1)/((nu+3)*(sigma^2)),
                 dldd = function(y,mu,sigma,nu) {
                                      s2 <- sigma^2
                                     dsq <- ((y-mu)^2)/s2
                                   omega <- (nu+1)/(nu+dsq) 
                                    dldd <- (omega*dsq-1)/sigma 
                                    dldd
                                    },
               d2ldd2 = function(sigma,nu) {
                                      s2 <- sigma^2
                                  d2ldd2 <- -(2*nu)/((nu+3)*s2)
                                  d2ldd2 
                                    },
                 dldv = function(y,mu,sigma,nu) {
                                      s2 <- sigma^2
                                     dsq <- ((y-mu)^2)/s2
                                   omega <- (nu+1)/(nu+dsq)
                                    dsq3 <- 1+(dsq/nu)
                                      v2 <- nu/2
                                      v3 <- (nu+1)/2 
                                    dldv <- -log(dsq3)+(omega*dsq-1)/nu +digamma(v3)-digamma(v2) 
                                    dldv <- dldv/2
                                    },
               d2ldv2 = function(y,mu,sigma,nu)  {
                                      v2 <- nu/2
                                      v3 <-(nu+1)/2
                                  d2ldv2 <- trigamma(v3)-trigamma(v2)+(2*(nu+5))/(nu*(nu+1)*(nu+3))
                                  d2ldv2 <- d2ldv2/4
                                  d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)  
                                    },
              d2ldmdd = function(y)  rep(0,length(y)),
              d2ldmdv = function(y)  rep(0,length(y)),
              d2ldddv = function(sigma,nu)  2/(sigma*(nu+3)*(nu+1)),
          G.dev.incr  = function(y,mu,sigma,nu,...) -2*dTF(y,mu,sigma,nu,log=TRUE),                           
                rqres = expression(
                rqres(pfun="pTF", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu) 
                                   ),
           mu.initial =  expression(mu <- (y+mean(y))/2), 
         sigma.initial =  expression(sigma <- rep(sd(y),length(y))),
            nu.initial =  expression(   nu <- rep(10, length(y))), 
              mu.valid = function(mu) TRUE , 
           sigma.valid = function(sigma)  all(sigma > 0),
              nu.valid = function(nu) all(nu > 0), 
               y.valid = function(y)  TRUE
          ),
            class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------
dTF<-function(x, mu=0, sigma=1, nu=10, log=FALSE)
 {  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))  
    if (length(nu)>1) lik <- ifelse(nu>1000000, 
                       dNO(x,mu=mu,sigma=sigma,log=FALSE), 
                       (1/sigma)*dt((x-mu)/sigma, df=nu, log =FALSE)) 
    else lik <- if (nu>1000000) dNO(x,mu=mu,sigma=sigma,log=FALSE)
                else  (1/sigma)*dt((x-mu)/sigma, df=nu, log =FALSE)
    fy <- if(log==FALSE) lik else log(lik)
    fy 
  }
#---------------------------------------------------------------------------------------- 
pTF <- function(q, mu=0, sigma=1, nu=10, lower.tail = TRUE, log.p = FALSE)
  { if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))  
    if (length(nu)>1) cdf <- ifelse(nu>1000000, 
                                   pNO(q, mu=mu, sigma=sigma, lower.tail=lower.tail, log.p=log.p),
                                    pt((q-mu)/sigma, df=nu,   lower.tail=lower.tail, log.p=log.p))
    else cdf <- if (nu>1000000)    pNO(q, mu=mu, sigma=sigma, lower.tail=lower.tail, log.p=log.p)
                   else pt((q-mu)/sigma, df=nu, lower.tail = lower.tail, log.p = log.p)
    cdf
   }
#----------------------------------------------------------------------------------------
qTF <- function(p, mu=0, sigma=1, nu=10, lower.tail = TRUE, log.p = FALSE)
  { if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))  
    if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
    if (length(nu)>1) q <- ifelse(nu>1000000, 
                                 qNO(p, mu=mu, sigma=sigma, lower.tail=lower.tail, log.p=log.p),
                                  mu+sigma*qt(p,df=nu, lower.tail = lower.tail))
    else q <- if (nu>1000000)    qNO(p, mu=mu, sigma=sigma, lower.tail=lower.tail, log.p=log.p)
              else               mu+sigma*qt(p,df=nu, lower.tail = lower.tail)
    q
   }


rTF <- function(n, mu=0, sigma=1, nu=10)
  { if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))  
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
    n <- ceiling(n)
    p <- runif(n)
    r <- qTF(p, mu=mu, sigma=sigma, nu=nu)
    r
  }
