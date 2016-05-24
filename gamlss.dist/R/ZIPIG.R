# 3/3/10  
# ---------------------------------------------------------------------------------------
# zero inflated Poisson inverse Gaussian (with probability y=0 is nu) 20/04/10
# ---------------------------------------------------------------------------------------
ZIPIG = function (mu.link = "log", sigma.link = "log", nu.link = "logit") 
{
    mstats <- checklink("mu.link", "ZIPIG", substitute(mu.link), 
        c("inverse", "log", "identity"))
    dstats <- checklink("sigma.link", "ZIPIG", substitute(sigma.link), 
        c("inverse", "log", "identity"))
    vstats <- checklink("nu.link", "ZIPIG", substitute(nu.link), 
        c("logit", "probit", "cloglog", "log", "own"))

    structure(list(family = c("ZIPIG", "Zero inflated Poisson inverse Gaussian"),
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
                  dldm = function(y,mu,sigma,nu) {dldm0 <- (1-nu)*((nu+(1-nu)*dPIG(0,mu,sigma))^(-1))*dPIG(0,mu,sigma)*PIG()$dldm(0,mu,sigma)
                         dldm <- ifelse(y==0, dldm0, PIG()$dldm(y,mu,sigma))
                         dldm}, 
               d2ldm2 = function(y,mu,sigma,nu) {dldm0 <- (1-nu)*((nu+(1-nu)*dPIG(0,mu,sigma))^(-1))*dPIG(0,mu,sigma)*PIG()$dldm(0,mu,sigma)
                          dldm <- ifelse(y==0, dldm0, PIG()$dldm(y,mu,sigma))
                        d2ldm2 <- -dldm*dldm
                         d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)    
                        d2ldm2},
                 dldd = function(y,mu,sigma,nu) {dldd0 <- (1-nu)*((nu+(1-nu)*dPIG(0,mu,sigma))^(-1))*dPIG(0,mu,sigma)*PIG()$dldd(0,mu,sigma)
                         dldd <- ifelse(y==0, dldd0, PIG()$dldd(y,mu,sigma))
                         dldd}, 
               d2ldd2 = function(y,mu,sigma,nu) {dldd0 <- (1-nu)*((nu+(1-nu)*dPIG(0,mu,sigma))^(-1))*dPIG(0,mu,sigma)*PIG()$dldd(0,mu,sigma)
                         dldd <- ifelse(y==0, dldd0, PIG()$dldd(y,mu,sigma))
                  d2ldd2 <- -dldd^2  
                  d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
                  d2ldd2}, 
                 dldv = function(y,mu,sigma,nu) {dldv0 <- ((nu+(1-nu)*dPIG(0,mu,sigma))^(-1))*(1-dPIG(0,mu,sigma))
                         dldv <- ifelse(y==0, dldv0, -1/(1-nu))
                         dldv}, 
               d2ldv2 = function(y,mu,sigma,nu) {dldv0 <- ((nu+(1-nu)*dPIG(0,mu,sigma))^(-1))*(1-dPIG(0,mu,sigma))
                         dldv <- ifelse(y==0, dldv0, -1/(1-nu))      
              d2ldv2 <- -dldv^2
              d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)  
              d2ldv2},
        d2ldmdd = function(y,mu,sigma,nu) {dldm0 <- (1-nu)*((nu+(1-nu)*dPIG(0,mu,sigma))^(-1))*dPIG(0,mu,sigma)*PIG()$dldm(0,mu,sigma)
                         dldm <- ifelse(y==0, dldm0, PIG()$dldm(y,mu,sigma))
                         dldd0 <- (1-nu)*((nu+(1-nu)*dPIG(0,mu,sigma))^(-1))*dPIG(0,mu,sigma)*PIG()$dldd(0,mu,sigma)
                         dldd <- ifelse(y==0, dldd0, PIG()$dldd(y,mu,sigma))
                         d2ldm2<--dldm*dldd
                         d2ldm2}, 
                 d2ldmdv = function(y,mu,sigma,nu) {  #partial derivate of log-density respect to mu and alpha
        dldm0 <- (1-nu)*((nu+(1-nu)*dPIG(0,mu,sigma))^(-1))*dPIG(0,mu,sigma)*PIG()$dldm(0,mu,sigma)
                         dldm <- ifelse(y==0, dldm0, PIG()$dldm(y,mu,sigma))
               dldv0 <- ((nu+(1-nu)*dPIG(0,mu,sigma))^(-1))*(1-dPIG(0,mu,sigma))
                         dldv <- ifelse(y==0, dldv0, -1/(1-nu))          
            d2ldmdv <- -dldm*dldv
            d2ldmdv
        },
        
        d2ldddv = function(y,mu,sigma,nu) {   #partial derivate of log-density respect to sigma and alpha
dldd0 <- (1-nu)*((nu+(1-nu)*dPIG(0,mu,sigma))^(-1))*dPIG(0,mu,sigma)*PIG()$dldd(0,mu,sigma)
                         dldd <- ifelse(y==0, dldd0, PIG()$dldd(y,mu,sigma))
               dldv0 <- ((nu+(1-nu)*dPIG(0,mu,sigma))^(-1))*(1-dPIG(0,mu,sigma))
                         dldv <- ifelse(y==0, dldv0, -1/(1-nu))            
            d2ldddv <- -dldd*dldv
            d2ldddv
        },        
        G.dev.incr  = function(y,mu,sigma,nu,...) -2*dZIPIG(y, mu = mu, sigma = sigma, nu=nu, log = TRUE), 
        rqres = expression(
                          rqres(pfun="pZIPIG", type="Discrete", ymin=0, y=y, mu=mu, sigma=sigma, nu=nu)
                                 ), 
           mu.initial = expression(mu <- (y + mean(y))/2),            
##           mu.initial = expression(mu <-  y+0.5),
         sigma.initial = expression(
                      sigma <- rep( max( ((var(y)-mean(y))/(mean(y)^2)),0.1),length(y))),
           nu.initial = expression(nu <- rep(0.1, length(y))),#expression(nu <- rep(((sum(y==0)/length(y))+0.01)/2, length(y))),#expression(nu <- rep(0.3, length(y))),
              mu.valid = function(mu) all(mu > 0) , 
           sigma.valid = function(sigma)  all(sigma > 0), 
             nu.valid = function(nu) all(nu > 0 & nu < 1),           
               y.valid = function(y)  all(y >= 0)
          ),
            class = c("gamlss.family","family"))
}
#-------------------------------------------------------------------------------------------
dZIPIG<-function(x, mu = 1, sigma = 1, nu = 0.3, log = FALSE)
 { 
        if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
        if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
        if (any(nu <= 0)|any(nu >= 1))  stop(paste("nu must be between 0 and 1 ", "\n", ""))
        if (any(x < 0) )  stop(paste("x must be >=0", "\n", ""))  
        if (length(sigma)>1) fy <- ifelse(sigma>0.0001, dPIG(x, mu = mu, sigma=sigma, log = T), 
                                          dPO(x, mu = mu, log = T) )
        else fy <- if (sigma<0.0001) dPO(x, mu = mu, log = T) 
                   else dPIG(x, mu = mu, sigma=sigma, log = T)
          logfy <- rep(0, length(x))
          logfy <- ifelse((x==0), log(nu+(1-nu)*exp(fy)), (log(1-nu) + fy ))          
          if(log == FALSE) fy2 <- exp(logfy) else fy2 <- logfy
          fy2
  }
#------------------------------------------------------------------------------------------
pZIPIG <- function(q, mu = 1, sigma = 1, nu = 0.3, lower.tail = TRUE, log.p = FALSE)
  {     
        if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
        if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
        if (any(nu <= 0)|any(nu >= 1))  #In this parametrization  nu = alpha
                    stop(paste("nu must be between 0 and 1 ", "\n", ""))
        if (any(q < 0) )  stop(paste("y must be >=0", "\n", ""))
        if (length(sigma)>1) cdf <- ifelse(sigma>0.0001, pPIG(q, mu = mu, sigma=sigma, lower.tail=T, log.p = F), 
                                          ppois(q, lambda = mu) )
        else cdf <- if (sigma<0.0001) ppois(q, lambda = mu)
                   else pPIG(q, mu = mu, sigma=sigma, lower.tail=T, log.p = F)
         cdf <- nu + (1-nu)*cdf
         if(lower.tail == TRUE) cdf <- cdf else cdf <-1-cdf
         if(log.p==FALSE) cdf <- cdf else cdf <- log(cdf)    
         cdf
   }
#------------------------------------------------------------------------------------------
qZIPIG <- function(p, mu = 1, sigma = 1, nu = 0.3, lower.tail = TRUE, log.p = FALSE)
  {      
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", ""))
          if (any(nu <= 0)|any(nu >= 1))  #In this parametrization  nu = alpha
                     stop(paste("nu must be between 0 and 1 ", "\n", ""))
          if (any(p < 0) | any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
          if (log.p == TRUE) p <- exp(p)   else p <- p
          if (lower.tail == TRUE)  p <- p  else p <- 1 - p
          pnew <- (p-nu)/(1-nu)-1e-10
          pnew <- ifelse((pnew > 0 ),pnew, 0)
          if (length(sigma)>1) q <- ifelse(sigma>0.0001,  qPIG(pnew, mu = mu, sigma=sigma, lower.tail=T, log.p = F), 
                                          qpois(pnew, lambda = mu) )
          else q <- if (sigma<0.0001) qpois(pnew, lambda = mu)
                   else qPIG(pnew, mu = mu, sigma=sigma, lower.tail=T, log.p = F)
           # q2 <- suppressWarnings(ifelse((pnew > 0 ), q, 0))
#          suppressWarnings(q <- ifelse((pnew > 0 ), qpois(pnew, lambda = mu, ), 0))
          q
   }
#------------------------------------------------------------------------------------------
rZIPIG <- function(n, mu = 1, sigma = 1, nu = 0.3)
  { 
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
          if (any(nu <= 0)|any(nu >= 1))  #In this parametrization  nu = alpha
          stop(paste("nu must be between 0 and 1 ", "\n", ""))
          if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
          n <- ceiling(n)
          p <- runif(n)
          r <- qZIPIG(p, mu=mu, sigma=sigma, nu=nu)
          r
  }
#------------------------------------------------------------------------------------------
