############################################
##### STANDARDIZED SKEW T DISTRIBUTION #####
############################################
####   working  09-07-2012
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
# for testing
#library(gamlss)
#source('/Volumes/Data/Users/stasinom/Documents/gamlss/projects/DISTRIBUTIONS/testingDistributions/testContDist.R')
#testContDist("SST", y.range=c(-Inf,Inf), mu.range = c(-Inf,Inf), sigma.range=c(0, Inf),  nu.range = c(0, Inf), tau.range=c(2, Inf), mu.val = c(0), sigma.val=c(1), nu=c(0.5, 1, 2), tau=c(2.5, 4, 10))
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
SST<-function (mu.link = "identity", sigma.link = "log", nu.link = "log", 
               tau.link = "logshiftto2") 
{
    mstats <- checklink("mu.link", "SST", substitute(mu.link), 
        c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "SST", substitute(sigma.link), 
        c("inverse", "log", "identity", "own"))
    vstats <- checklink("nu.link", "SST", substitute(nu.link), 
        c("inverse", "log", "identity", "own"))
    tstats <- checklink("tau.link", "SST ", substitute(tau.link), 
                        c("logshiftto2", "log", "identity", "own"))
    structure(
           list(family = c("SST", "SST"), 
            parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE, tau = TRUE), 
                 nopar = 4, 
                  type = "Continuous", 
               mu.link = as.character(substitute(mu.link)), 
            sigma.link = as.character(substitute(sigma.link)), 
               nu.link = as.character(substitute(nu.link)), 
              tau.link = as.character(substitute(tau.link)), 
            mu.linkfun = mstats$linkfun, 
         sigma.linkfun = dstats$linkfun, 
            nu.linkfun = vstats$linkfun, 
           tau.linkfun = tstats$linkfun, 
            mu.linkinv = mstats$linkinv, 
         sigma.linkinv = dstats$linkinv, 
            nu.linkinv = vstats$linkinv, 
           tau.linkinv = tstats$linkinv, 
                 mu.dr = mstats$mu.eta, 
              sigma.dr = dstats$mu.eta, 
                 nu.dr = vstats$mu.eta, 
                tau.dr = tstats$mu.eta, 
                  dldm = function(y, mu, sigma, nu, tau) {
                                 m1 <- (2*tau^(1/2)*(nu^2-1))/((tau-1)*beta(1/2, tau/2)*nu)
                                 m2 <- (tau*(nu^3+(1/nu^3)))/((tau-2)*(nu+(1/nu)))
                                 s1 <- sqrt(m2-m1^2)
                                mu1 <- mu- ((sigma*m1)/s1)
                             sigma1 <- sigma/s1
                               dldm <- ST3()$dldm(y, mu1, sigma1, nu, tau)
                               dldm
                                }, 
                  d2ldm2 = function(y, mu, sigma, nu, tau) {
                                 m1 <- (2*tau^(1/2)*(nu^2-1))/((tau-1)*beta(1/2, tau/2)*nu)
                                 m2 <- (tau*(nu^3+(1/nu^3)))/((tau-2)*(nu+(1/nu)))
                                 s1 <- sqrt(m2-m1^2)
                                mu1 <- mu- ((sigma*m1)/s1)
                             sigma1 <- sigma/s1        
                               dldm <- ST3()$dldm(y, mu1, sigma1, nu, tau)
                             d2ldm2 <- -dldm*dldm
                             d2ldm2
                               }, 
                     dldd = function(y, mu, sigma, nu, tau) {
                                m1 <- (2*tau^(1/2)*(nu^2-1))/((tau-1)*beta(1/2, tau/2)*nu)
                                m2 <- (tau*(nu^3+(1/nu^3)))/((tau-2)*(nu+(1/nu)))
                                s1 <- sqrt(m2-m1^2)
                               mu1 <- mu- ((sigma*m1)/s1)
                            sigma1 <- sigma/s1                  
                              dldd <- -(m1/s1)*ST3()$dldm(y, mu1, sigma1, nu, tau) + 
                                        (1/s1)*ST3()$dldd(y, mu1, sigma1, nu, tau)
                              dldd
                              }, 
                   d2ldd2 = function(y, mu, sigma, nu, tau) {
                                m1 <- (2*tau^(1/2)*(nu^2-1))/((tau-1)*beta(1/2, tau/2)*nu)
                                m2 <- (tau*(nu^3+(1/nu^3)))/((tau-2)*(nu+(1/nu)))
                                s1 <- sqrt(m2-m1^2)
                               mu1 <- mu- ((sigma*m1)/s1)
                            sigma1 <- sigma/s1                  
                              dldd <- -(m1/s1)*ST3()$dldm(y, mu1, sigma1, nu, tau) + 
                                       (1/s1)*ST3()$dldd(y, mu1, sigma1, nu, tau)
                            d2ldd2 <- -dldd*dldd
                            d2ldd2
                             }, 
                      dldv = function(y, mu, sigma, nu, tau) {
                                m1 <- ((2*tau^(1/2)*(nu^2-1))/((tau-1)*beta(1/2, tau/2)*nu))
                                m2 <- ((tau*(nu^3+(1/nu^3)))/((tau-2)*(nu+(1/nu))))
                                s1 <- sqrt(m2-m1^2)
                               mu1 <- mu- ((sigma*m1)/s1)
                            sigma1 <- sigma/s1                  
                           dl1dmu1 <- ST3()$dldm(y, mu1, sigma1, nu, tau)
                            dl1dd1 <- ST3()$dldd(y, mu1, sigma1, nu, tau)
                             dl1dv <- ST3()$dldv(y, mu1, sigma1, nu, tau)
                           dmu1dm1 <- -sigma/s1
                           dmu1ds1 <- (sigma*m1)/(s1^2)
                            dd1ds1 <- -sigma/(s1^2)
                             dm1dv <- ((2*tau^(1/2))/((tau-1)*beta(1/2, tau/2)))*((nu^2+1)/(nu^2))
                             dm2dv <- m2*((6*nu^5/(nu^6+1))-(2/nu)-(2*nu/(nu^2+1)))
                             ds1dv <- (dm2dv - 2*m1*dm1dv)/(2*s1)
                              dldv <- dl1dmu1*dmu1dm1*dm1dv + dl1dmu1*dmu1ds1*ds1dv + 
                                      dl1dd1*dd1ds1*ds1dv + dl1dv
                              dldv
                             }, 
                      d2ldv2 = function(y, mu, sigma, nu, tau) {
                                m1 <- (2*tau^(1/2)*(nu^2-1))/((tau-1)*beta(1/2, tau/2)*nu)
                                m2 <- (tau*(nu^3+(1/nu^3)))/((tau-2)*(nu+(1/nu)))
                                s1 <- sqrt(m2-m1^2)
                               mu1 <- mu- ((sigma*m1)/s1)
                            sigma1 <- sigma/s1                  
                           dl1dmu1 <- ST3()$dldm(y, mu1, sigma1, nu, tau)
                            dl1dd1 <- ST3()$dldd(y, mu1, sigma1, nu, tau)
                             dl1dv <- ST3()$dldv(y, mu1, sigma1, nu, tau)
                           dmu1dm1 <- -sigma/s1
                           dmu1ds1 <- (sigma*m1)/(s1^2)
                            dd1ds1 <- -sigma/(s1^2) 
                             dm1dv <- ((2*tau^(1/2))/((tau-1)*beta(1/2, tau/2)))*((nu^2+1)/(nu^2))
                             dm2dv <- m2*((6*nu^5/(nu^6+1))-(2/nu)-(2*nu/(nu^2+1)))
                             ds1dv <- (dm2dv - 2*m1*dm1dv)/(2*s1)
                              dldv <- dl1dmu1*dmu1dm1*dm1dv + dl1dmu1*dmu1ds1*ds1dv + 
                                      dl1dd1*dd1ds1*ds1dv + dl1dv
                            d2ldv2 <- -dldv*dldv
                            d2ldv2
                             }, 
                      dldt = function(y, mu, sigma, nu, tau) {
                               m1 <- (2*tau^(1/2)*(nu^2-1))/((tau-1)*beta(1/2, tau/2)*nu)
                               m2 <- (tau*(nu^3+(1/nu^3)))/((tau-2)*(nu+(1/nu)))
                               s1 <- sqrt(m2-m1^2)
                              mu1 <- mu- ((sigma*m1)/s1)
                           sigma1 <- sigma/s1                 
                          dl1dmu1 <- ST3()$dldm(y, mu1, sigma1, nu, tau)
                           dl1dd1 <- ST3()$dldd(y, mu1, sigma1, nu, tau)
                            dl1dt <- ST3()$dldt(y, mu1, sigma1, nu, tau)
                          dmu1dm1 <- -sigma/s1
                          dmu1ds1 <- (sigma*m1)/(s1^2)
                           dd1ds1 <- -sigma/(s1^2)  
                            dm1dt <- m1*((1/(2*tau))-(1/(tau-1))- 0.5*(digamma(tau/2)) + 0.5*(digamma((tau+1)/2)))
                            dm2dt <- -m2*(2/(tau*(tau-2)))                 
                            ds1dt <- (dm2dt - 2*m1*dm1dt)/(2*s1)         
                             dldt <- dl1dmu1*dmu1dm1*dm1dt + dl1dmu1*dmu1ds1*ds1dt + dl1dd1*dd1ds1*ds1dt + dl1dt
                             dldt
                             }, 
                       d2ldt2 = function(y, mu, sigma, nu, tau) {
                               m1 <- (2*tau^(1/2)*(nu^2-1))/((tau-1)*beta(1/2, tau/2)*nu)
                               m2 <- (tau*(nu^3+(1/nu^3)))/((tau-2)*(nu+(1/nu)))
                               s1 <- sqrt(m2-m1^2)
                              mu1 <- mu- ((sigma*m1)/s1)
                           sigma1 <- sigma/s1                 
                          dl1dmu1 <- ST3()$dldm(y, mu1, sigma1, nu, tau)
                           dl1dd1 <- ST3()$dldd(y, mu1, sigma1, nu, tau)
                            dl1dt <- ST3()$dldt(y, mu1, sigma1, nu, tau)
                          dmu1dm1 <- -sigma/s1
                          dmu1ds1 <- (sigma*m1)/(s1^2)
                           dd1ds1 <- -sigma/(s1^2)
                            dm1dt <- m1*((1/(2*tau))-(1/(tau-1))- 0.5*(digamma(tau/2)) + 0.5*(digamma((tau+1)/2)))
                            dm2dt <- -m2*(2/(tau*(tau-2)))                 
                            ds1dt <- (dm2dt - 2*m1*dm1dt)/(2*s1)
                             dldt <- dl1dmu1*dmu1dm1*dm1dt + dl1dmu1*dmu1ds1*ds1dt + dl1dd1*dd1ds1*ds1dt + dl1dt
                           dl2dt2 <- -dldt*dldt
                           dl2dt2
                            }, 
                     d2ldmdd = function(y, mu, sigma, nu, tau) {
                              m1 <- (2*tau^(1/2)*(nu^2-1))/((tau-1)*beta(1/2, tau/2)*nu)
                              m2 <- (tau*(nu^3+(1/nu^3)))/((tau-2)*(nu+(1/nu)))
                              s1 <- sqrt(m2-m1^2)
                             mu1 <- mu- ((sigma*m1)/s1)
                          sigma1 <- sigma/s1                  
                            dldm <- ST3()$dldm(y, mu1, sigma1, nu, tau)
                            dldd <- -(m1/s1)*ST3()$dldm(y, mu1, sigma1, nu, tau) + (1/s1)*ST3()$dldd(y, mu1, sigma1, nu, tau)
                         d2ldmdd <- -dldm*dldd
                        d2ldmdd
                            }, 
                     d2ldmdv = function(y, mu, sigma, nu, tau) {
          m1 <- (2*tau^(1/2)*(nu^2-1))/((tau-1)*beta(1/2, tau/2)*nu)
          m2 <- (tau*(nu^3+(1/nu^3)))/((tau-2)*(nu+(1/nu)))
          s1 <- sqrt(m2-m1^2)
         mu1 <- mu- ((sigma*m1)/s1)
      sigma1 <- sigma/s1                  
          dl1dmu1<- ST3()$dldm(y, mu1, sigma1, nu, tau)
          dl1dd1 <- ST3()$dldd(y, mu1, sigma1, nu, tau)
          dl1dv <- ST3()$dldv(y, mu1, sigma1, nu, tau)
          dmu1dm1 <- -sigma/s1
          dmu1ds1 <- (sigma*m1)/(s1^2)
          dd1ds1 <- -sigma/(s1^2)
          dm1dv <- ((2*tau^(1/2))/((tau-1)*beta(1/2, tau/2)))*((nu^2+1)/(nu^2))
          dm2dv <- m2*((6*nu^5/(nu^6+1))-(2/nu)-(2*nu/(nu^2+1)))
          ds1dv <- (dm2dv - 2*m1*dm1dv)/(2*s1)
          dldv <- dl1dmu1*dmu1dm1*dm1dv + dl1dmu1*dmu1ds1*ds1dv + 
                  dl1dd1*dd1ds1*ds1dv + dl1dv
          dldm <- ST3()$dldm(y, mu1, sigma1, nu, tau)
          d2ldmdv <- -dldm*dldv
          d2ldmdv
        }, 
        d2ldmdt = function(y, mu, sigma, nu, tau) {
          m1 <- (2*tau^(1/2)*(nu^2-1))/((tau-1)*beta(1/2, tau/2)*nu)
          m2 <- (tau*(nu^3+(1/nu^3)))/((tau-2)*(nu+(1/nu)))
          s1 <- sqrt(m2-m1^2)
         mu1 <- mu- ((sigma*m1)/s1)
      sigma1 <- sigma/s1                  
          dl1dmu1 <- ST3()$dldm(y, mu1, sigma1, nu, tau)
          dl1dd1 <- ST3()$dldd(y, mu1, sigma1, nu, tau)
          dl1dt <- ST3()$dldt(y, mu1, sigma1, nu, tau)
          dmu1dm1 <- -sigma/s1
          dmu1ds1 <- (sigma*m1)/(s1^2)
          dd1ds1 <- -sigma/(s1^2)
           dm1dt <- m1*((1/(2*tau))-(1/(tau-1))- 0.5*(digamma(tau/2)) + 0.5*(digamma((tau+1)/2)))
          dm2dt <- -m2*(2/(tau*(tau-2)))                 
          ds1dt <- (dm2dt - 2*m1*dm1dt)/(2*s1)
          dldt <- dl1dmu1*dmu1dm1*dm1dt + dl1dmu1*dmu1ds1*ds1dt +
                  dl1dd1*dd1ds1*ds1dt + dl1dt
          dldm <- ST3()$dldm(y, mu1, sigma1, nu, tau)
          d2ldmdt <- -dldm*dldt
          d2ldmdt
        }, 
        d2ldddv = function(y, mu, sigma, nu, tau) {
          m1 <- (2*tau^(1/2)*(nu^2-1))/((tau-1)*beta(1/2, tau/2)*nu)
          m2 <- (tau*(nu^3+(1/nu^3)))/((tau-2)*(nu+(1/nu)))
          s1 <- sqrt(m2-m1^2)
         mu1 <- mu- ((sigma*m1)/s1)
      sigma1 <- sigma/s1                  
          dldd <- -(m1/s1)*ST3()$dldm(y, mu1, sigma1, nu, tau) + 
                    (1/s1)*ST3()$dldd(y, mu1, sigma1, nu, tau)
          dl1dmu1<- ST3()$dldm(y, mu1, sigma1, nu, tau)
          dl1dd1 <- ST3()$dldd(y, mu1, sigma1, nu, tau)
          dl1dv <- ST3()$dldv(y, mu1, sigma1, nu, tau)
          dmu1dm1 <- -sigma/s1
          dmu1ds1 <- (sigma*m1)/(s1^2)
          dd1ds1 <- -sigma/(s1^2)
          dm1dv <- ((2*tau^(1/2))/((tau-1)*beta(1/2, tau/2)))*((nu^2+1)/(nu^2))
          dm2dv <- m2*((6*nu^5/(nu^6+1))-(2/nu)-(2*nu/(nu^2+1)))
          ds1dv <- (dm2dv - 2*m1*dm1dv)/(2*s1)
          dldv <- dl1dmu1*dmu1dm1*dm1dv + dl1dmu1*dmu1ds1*ds1dv + 
                  dl1dd1*dd1ds1*ds1dv + dl1dv
          d2ldddv <- -dldd*dldv
          d2ldddv
        }, 
        d2ldddt = function(y, mu, sigma, nu, tau) {
          m1 <- (2*tau^(1/2)*(nu^2-1))/((tau-1)*beta(1/2, tau/2)*nu)
          m2 <- (tau*(nu^3+(1/nu^3)))/((tau-2)*(nu+(1/nu)))
          s1 <- sqrt(m2-m1^2)
         mu1 <- mu- ((sigma*m1)/s1)
      sigma1 <- sigma/s1                  
          dldd <- -(m1/s1)*ST3()$dldm(y, mu1, sigma1, nu, tau) + 
                    (1/s1)*ST3()$dldd(y, mu1, sigma1, nu, tau)
          dl1dmu1 <- ST3()$dldm(y, mu1, sigma1, nu, tau)
          dl1dd1 <- ST3()$dldd(y, mu1, sigma1, nu, tau)
          dl1dt <- ST3()$dldt(y, mu1, sigma1, nu, tau)
          dmu1dm1 <- -sigma/s1
          dmu1ds1 <- (sigma*m1)/(s1^2)
          dd1ds1 <- -sigma/(s1^2)
          dm1dt <- m1*((1/(2*tau))-(1/(tau-1))- 0.5*(digamma(tau/2)) + 0.5*(digamma((tau+1)/2)))
          dm2dt <- -m2*(2/(tau*(tau-2)))                 
          ds1dt <- (dm2dt - 2*m1*dm1dt)/(2*s1)
          dldt <- dl1dmu1*dmu1dm1*dm1dt + dl1dmu1*dmu1ds1*ds1dt +
                  dl1dd1*dd1ds1*ds1dt + dl1dt
          d2ldddt <- -dldd*dldt
          d2ldddt
        }, 
        d2ldvdt = function(y, mu, sigma, nu, tau) {
          m1 <- (2*tau^(1/2)*(nu^2-1))/((tau-1)*beta(1/2, tau/2)*nu)
          m2 <- (tau*(nu^3+(1/nu^3)))/((tau-2)*(nu+(1/nu)))
          s1 <- sqrt(m2-m1^2)
         mu1 <- mu- ((sigma*m1)/s1)
      sigma1 <- sigma/s1                  
          dl1dmu1<- ST3()$dldm(y, mu1, sigma1, nu, tau)
          dl1dd1 <- ST3()$dldd(y, mu1, sigma1, nu, tau)
          dl1dv <- ST3()$dldv(y, mu1, sigma1, nu, tau)
          dl1dt <- ST3()$dldt(y, mu1, sigma1, nu, tau)
          dmu1dm1 <- -sigma/s1
          dmu1ds1 <- (sigma*m1)/(s1^2)
          dd1ds1 <- -sigma/(s1^2)
          dm1dv <- ((2*tau^(1/2))/((tau-1)*beta(1/2, tau/2)))*((nu^2+1)/(nu^2))
          dm2dv <- m2*((6*nu^5/(nu^6+1))-(2/nu)-(2*nu/(nu^2+1)))
          ds1dv <- (dm2dv - 2*m1*dm1dv)/(2*s1)

          dm1dt <- m1*((1/(2*tau))-(1/(tau-1))- 0.5*(digamma(tau/2)) + 0.5*(digamma((tau+1)/2)))
          dm2dt <- -m2*(2/(tau*(tau-2)))                 
          ds1dt <- (dm2dt - 2*m1*dm1dt)/(2*s1)

          dldv <- dl1dmu1*dmu1dm1*dm1dv + dl1dmu1*dmu1ds1*ds1dv + 
                  dl1dd1*dd1ds1*ds1dv + dl1dv
          dldt <- dl1dmu1*dmu1dm1*dm1dt + dl1dmu1*dmu1ds1*ds1dt +
                  dl1dd1*dd1ds1*ds1dt + dl1dt
          d2ldvdt <- -dldv*dldt
          d2ldvdt
        }, 
        G.dev.incr = function(y, mu, sigma, nu, tau, ...) -2 * dSST(y, mu, sigma, nu, tau, log = TRUE), 
             rqres = expression(rqres(pfun = "pSST", 
              type = "Continuous", y = y, mu = mu, sigma = sigma, nu = nu, tau = tau)), 
        mu.initial = expression(mu <- (y + mean(y))/2), 
     sigma.initial = expression(sigma <- rep(sd(y), length(y))), 
        nu.initial = expression(nu <- rep(1, length(y))), 
       tau.initial = expression(tau <- rep(4, length(y))), 
          mu.valid = function(mu) TRUE, 
       sigma.valid = function(sigma) all(sigma > 0), 
          nu.valid = function(nu) all(nu > 0), 
         tau.valid = function(tau) all(tau > 0), 
           y.valid = function(y) TRUE), 
             class = c("gamlss.family", "family"))
}
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#Probability density function
dSST <- function(x, mu=0, sigma=1, nu=0.8, tau=7, log = FALSE){
#  if (any(tau <= 2)) 
#      stop(paste("tau must be greater than 2", "\n", ""))
  
  m1 <- (2*tau^(1/2)*(nu^2-1))/((tau-1)*beta(1/2, tau/2)*nu)
  m2 <- (tau*(nu^3+(1/nu^3)))/((tau-2)*(nu+(1/nu)))
 # cat(tau[1], m1[1],m2[1], "\n")
  s1 <- sqrt(m2-m1^2)
  mu1 <- mu- ((sigma*m1)/s1)
  sigma1 <- sigma/s1
  fy <- dST3(x, mu1, sigma1, nu, tau, log = log)
  fy 
}
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#Cumulative density function
pSST <- function(q, mu=0, sigma=1, nu=0.8, tau=7, lower.tail = TRUE, 
                 log.p = FALSE){
#  if (any(tau <= 2)) 
#      stop(paste("tau must be greater than 2", "\n", ""))
  m1 <- (2*tau^(1/2)*(nu^2-1))/((tau-1)*beta(1/2, tau/2)*nu)
  m2 <- (tau*(nu^3+(1/nu^3)))/((tau-2)*(nu+(1/nu)))
  s1 <- sqrt(m2-m1^2)
  mu1 <- mu- ((sigma*m1)/s1)
  sigma1 <- sigma/s1
  fy <- pST3(q, mu1, sigma1, nu, tau, lower.tail = lower.tail, log.p = log.p)
  fy 
}
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#Quantile function
qSST <- function(p, mu=0, sigma=1, nu=0.8, tau=7, lower.tail = TRUE, 
                 log.p = FALSE){
  if (any(tau <= 2)) 
      stop(paste("tau must be greater than 2", "\n", ""))
  m1 <- (2*tau^(1/2)*(nu^2-1))/((tau-1)*beta(1/2, tau/2)*nu)
  m2 <- (tau*(nu^3+(1/nu^3)))/((tau-2)*(nu+(1/nu)))
  s1 <- sqrt(m2-m1^2)
  mu1 <- mu- ((sigma*m1)/s1)
  sigma1 <- sigma/s1
  fy <- qST3(p, mu1, sigma1, nu, tau, lower.tail = lower.tail, 
             log.p = log.p)
  fy 
}
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#Random generating function
rSST <- function(n, mu=0, sigma=1, nu=0.8, tau=7){
  if (any(tau <= 2)) 
      stop(paste("tau must be greater than 2", "\n", ""))
  if (any(n <= 0)) 
        stop(paste("n must be a positive integer", "\n", ""))  
  n <- ceiling(n)
  p <- runif(n)
  r <- qSST(p, mu = mu, sigma = sigma, nu=nu, tau=tau)
  r
}
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
