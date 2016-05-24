# ---------------------------------------------------------------------------------------
# zero altered negative binomial type I (with probability y=0 is nu) 01/03/10
# ---------------------------------------------------------------------------------------
ZANBI = function (mu.link = "log", sigma.link = "log", nu.link = "logit") 
{
    mstats <- checklink("mu.link", "ZANBI", substitute(mu.link), 
        c("inverse", "log", "identity"))
    dstats <- checklink("sigma.link", "ZANBI", substitute(sigma.link), 
        c("inverse", "log", "identity"))
    vstats <- checklink("nu.link", "ZANBI", substitute(nu.link), 
        c("logit", "probit", "cloglog", "log", "own"))

    structure(list(family = c("ZANBI", "Zero altered negative binomial type I"),
               parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE),
                    nopar = 3,
                     type = "Discrete",
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
                  dldm = function(y,mu,sigma,nu) {dldm0 <- NBI()$dldm(y,mu,sigma) + dNBI(0,mu,sigma)*NBI()$dldm(0,mu,sigma)/(1-dNBI(0,mu,sigma))
                         dldm <- ifelse(y==0, 0 , dldm0)
                         dldm}, 
               d2ldm2 = function(y,mu,sigma,nu) {dldm0 <- NBI()$dldm(y,mu,sigma) + dNBI(0,mu,sigma)*NBI()$dldm(0,mu,sigma)/(1-dNBI(0,mu,sigma))
                         dldm <- ifelse(y==0, 0 , dldm0)
                        d2ldm2 <- -dldm*dldm
                         d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)    
                        d2ldm2},
                 dldd = function(y,mu,sigma,nu) {dldd0 <- NBI()$dldd(y,mu,sigma) + dNBI(0,mu,sigma)*NBI()$dldd(0,mu,sigma)/(1-dNBI(0,mu,sigma))
                         dldd <- ifelse(y==0, 0 , dldd0)
                         dldd}, 
               d2ldd2 = function(y,mu,sigma,nu) {dldd0 <- NBI()$dldd(y,mu,sigma) + dNBI(0,mu,sigma)*NBI()$dldd(0,mu,sigma)/(1-dNBI(0,mu,sigma))
                         dldd <- ifelse(y==0, 0 , dldd0)
                  d2ldd2 <- -dldd^2  
                  d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
                  d2ldd2}, 
                 dldv = function(y,mu,sigma,nu) {dldv <- ifelse(y==0, 1/nu, -1/(1-nu))
                         dldv}, 
               d2ldv2 = function(y,mu,sigma,nu) {d2ldv2 <- -1/(nu*(1-nu))
              d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)  
              d2ldv2},
        d2ldmdd = function(y) {dldm0 <- NBI()$dldm(y,mu,sigma) + dNBI(0,mu,sigma)*NBI()$dldm(0,mu,sigma)/(1-dNBI(0,mu,sigma))
                         dldm <- ifelse(y==0, 0 , dldm0)
                         dldd0 <- NBI()$dldd(y,mu,sigma) + dNBI(0,mu,sigma)*NBI()$dldd(0,mu,sigma)/(1-dNBI(0,mu,sigma))
                         dldd <- ifelse(y==0, 0 , dldd0)
                         d2ldm2<--dldm*dldd
                         d2ldm2}, 
                 d2ldmdv = function(y) {d2ldmdv=0  
            d2ldmdv},
        d2ldddv = function(y) {d2ldddv=0
            d2ldddv},        
        G.dev.incr  = function(y,mu,sigma,...) -2*dZANBI(y, mu = mu, sigma = sigma, nu=nu, log = TRUE), 
        rqres = expression(
                          rqres(pfun="pZANBI", type="Discrete", ymin=0, y=y, mu=mu, sigma=sigma)
                                 ), 
           mu.initial = expression(mu <- (y + mean(y))/2),            
##           mu.initial = expression(mu <-  y+0.5),
         sigma.initial = expression(
         sigma <- rep( max( ((var(y)-mean(y))/(mean(y)^2)),0.1),length(y))),
           nu.initial = expression(nu <- rep(0.3, length(y))),
              mu.valid = function(mu) all(mu > 0) , 
           sigma.valid = function(sigma)  all(sigma > 0), 
             nu.valid = function(nu) all(nu > 0 & nu < 1),           
               y.valid = function(y)  all(y >= 0)
          ),
            class = c("gamlss.family","family"))
}
#-------------------------------------------------------------------------------------------
dZANBI<-function(x, mu = 1, sigma = 1, nu = 0.3, log = FALSE)
 { 
        if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
        if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
        if (any(nu <= 0)|any(nu >= 1))  stop(paste("nu must be between 0 and 1 ", "\n", ""))
        if (any(x < 0) )  stop(paste("x must be >=0", "\n", ""))  
        if (length(sigma)>1) fy <- ifelse(sigma>0.0001, dnbinom(x, size=1/sigma, mu = mu, log = T), 
                                          dPO(x, mu = mu, log = T) )
        else fy <- if (sigma<0.0001) dPO(x, mu = mu, log = T) 
                   else dnbinom(x, size=1/sigma, mu = mu, log = T)
         if (length(sigma)>1) fy0 <- ifelse(sigma>0.0001, dnbinom(0, size=1/sigma, mu = mu, log = T), 
                                          dPO(0, mu = mu, log = T) )
        else fy0 <- if (sigma<0.0001) dPO(0, mu = mu, log = T) 
                   else dnbinom(0, size=1/sigma, mu = mu, log = T)
                  
          logfy <- rep(0, length(x))
          logfy <- ifelse((x==0), log(nu), log(1-nu) + fy - log(1-exp(fy0)))          

          if(log == FALSE) fy2 <- exp(logfy) else fy2 <- logfy
          fy2
  }
#------------------------------------------------------------------------------------------
pZANBI <- function(q, mu = 1, sigma = 1, nu = 0.3, lower.tail = TRUE, log.p = FALSE)
  {     
        if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
        if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
        if (any(nu <= 0)|any(nu >= 1))  #In this parametrization  nu = alpha
                    stop(paste("nu must be between 0 and 1 ", "\n", ""))
        if (any(q < 0) )  stop(paste("y must be >=0", "\n", ""))
        if (length(sigma)>1) cdf1 <- ifelse(sigma>0.0001, pnbinom(q, size=1/sigma, mu=mu), 
                                          ppois(q, lambda = mu) )
        else cdf1 <- if (sigma<0.0001) ppois(q, lambda = mu)
                   else pnbinom(q, size=1/sigma, mu=mu)
 
        if (length(sigma)>1) cdf0 <- ifelse(sigma>0.0001, pnbinom(0, size=1/sigma, mu=mu), 
                                          ppois(0, lambda = mu) )
        else cdf0 <- if (sigma<0.0001) ppois(0, lambda = mu)
                   else pnbinom(0, size=1/sigma, mu=mu)
#         cdf <- rep(0,length(q))
#         cdf1 <- ppois(q, lambda = mu, lower.tail = TRUE, log.p = FALSE)
#         cdf2 <- ppois(0, lambda = mu, lower.tail = TRUE, log.p = FALSE)
         cdf3 <- nu+((1-nu)*(cdf1-cdf0)/(1-cdf0))
         cdf <- ifelse((q==0),nu,  cdf3)
         if(lower.tail == TRUE) cdf <- cdf else cdf <-1-cdf
         if(log.p==FALSE) cdf <- cdf else cdf <- log(cdf)    
         cdf
   }
#------------------------------------------------------------------------------------------
qZANBI <- function(p, mu = 1, sigma = 1, nu = 0.3, lower.tail = TRUE, log.p = FALSE)
  {      
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", ""))
          if (any(nu <= 0)|any(nu >= 1))  #In this parametrization  nu = alpha
                     stop(paste("nu must be between 0 and 1 ", "\n", ""))
          if (any(p < 0) | any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
          if (log.p == TRUE) p <- exp(p)   else p <- p
          if (lower.tail == TRUE)  p <- p  else p <- 1 - p
          pnew <- (p-nu)/(1-nu)
        if (length(sigma)>1) cdf0 <- ifelse(sigma>0.0001, pnbinom(0, size=1/sigma, mu=mu), 
                                          ppois(0, lambda = mu) )
        else cdf0 <- if (sigma<0.0001) ppois(0, lambda = mu)
                   else pnbinom(0, size=1/sigma, mu=mu)
          pnew2 <- cdf0*(1-pnew) + pnew           
          pnew2 <- ifelse((pnew2 > 0 ),pnew2, 0)
          if (length(sigma)>1) q <- ifelse(sigma>0.0001,  qnbinom(pnew2, size=1/sigma, mu=mu), 
                                          qpois(pnew2, lambda = mu) )
          else q <- if (sigma<0.0001) qpois(pnew2, lambda = mu)
                   else qnbinom(pnew2, size=1/sigma, mu=mu)
           # q2 <- suppressWarnings(ifelse((pnew > 0 ), q, 0))
#          suppressWarnings(q <- ifelse((pnew > 0 ), qpois(pnew, lambda = mu, ), 0))
          q
   }
#------------------------------------------------------------------------------------------
rZANBI <- function(n, mu = 1, sigma = 1, nu = 0.3)
  { 
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
          if (any(nu <= 0)|any(nu >= 1))  #In this parametrization  nu = alpha
          stop(paste("nu must be between 0 and 1 ", "\n", ""))
          if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
          n <- ceiling(n)
          p <- runif(n)
          r <- qZANBI(p, mu=mu, sigma=sigma, nu=nu)
          r
  }
#------------------------------------------------------------------------------------------
