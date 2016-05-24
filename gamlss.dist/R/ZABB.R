# the zero 
ZABB <- function (mu.link ="logit", sigma.link = "log", nu.link = "logit")
{
    mstats <- checklink("mu.link", "ZABB", substitute(mu.link),
                        c("logit", "probit", "cloglog", "cauchit", "log", "own"))
    dstats <- checklink("sigma.link", "Beta Binomial", substitute(sigma.link), 
                        c("inverse", "log", "identity", "sqrt", "own"))   
    vstats <- checklink("nu.link", "ZABB", substitute(nu.link),
                        c("logit", "probit", "cloglog", "cauchit", "log", "own"))
    structure(
            list(family = c("ZABB", "Zero Adjusted Beta Binomial"),
             parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE),
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
#========================First Derivatives========================
     dldm = function(y,mu,sigma,nu,bd) {dldm <- ifelse(y==0, 0, BB()$dldm(y,mu,sigma,bd) + (dBB(0,mu,sigma,bd)*BB()$dldm(0,mu,sigma,bd))/(1-dBB(0,mu,sigma,bd)))
                         dldm}, 
     dldd = function(y,mu,sigma,nu,bd) {dldd <- ifelse(y==0, 0, BB()$dldd(y,mu,sigma,bd) + (dBB(0,mu,sigma,bd)*BB()$dldd(0,mu,sigma,bd))/(1-dBB(0,mu,sigma,bd)))
                         dldd}, 
     dldv = function(y,mu,sigma,nu,bd) {dldv <- ifelse(y==0, 1/nu, -1/(1-nu))
                         dldv},
#========================Second Derivatives 1========================
   d2ldm2 = function(y,mu,sigma,nu,bd) {dldm <- ifelse(y==0, 0, BB()$dldm(y,mu,sigma,bd) + (dBB(0,mu,sigma,bd)*BB()$dldm(0,mu,sigma,bd))/(1-dBB(0,mu,sigma,bd)))
                                      d2ldm2 <- ifelse(y==0, 0, -(dldm)^2)
                                      d2ldm2},
   d2ldd2 = function(y,mu,sigma,nu,bd) {dldd <- ifelse(y==0, 0, BB()$dldd(y,mu,sigma,bd) + (dBB(0,mu,sigma,bd)*BB()$dldd(0,mu,sigma,bd))/(1-dBB(0,mu,sigma,bd)))
                                      d2ldd2 <- ifelse(y==0, 0, -(dldd)^2)
                                      d2ldd2},
   d2ldv2 = function(y,mu,sigma,nu,bd) {d2ldv2 <- -1/(nu*(1-nu))
                                      d2ldv2},
#========================Second Derivatives 2========================
    d2ldmdd = function(y,mu,sigma,nu,bd) {dldm <- ifelse(y==0, 0, BB()$dldm(y,mu,sigma,bd) + (dBB(0,mu,sigma,bd)*BB()$dldm(0,mu,sigma,bd))/(1-dBB(0,mu,sigma,bd)))
                        dldd <- ifelse(y==0, 0, BB()$dldd(y,mu,sigma,bd) + (dBB(0,mu,sigma,bd)*BB()$dldd(0,mu,sigma,bd))/(1-dBB(0,mu,sigma,bd)))
                     d2ldmdd <- -(dldm*dldd)
                     d2ldmdd},
    d2ldmdv = function(y,mu,sigma,nu,bd) {dldm <- ifelse(y==0, 0, BB()$dldm(y,mu,sigma,bd) + (dBB(0,mu,sigma,bd)*BB()$dldm(0,mu,sigma,bd))/(1-dBB(0,mu,sigma,bd)))
                        dldv <- ifelse(y==0, 1/nu, -1/(1-nu))
                     d2ldmdv <- -(dldm*dldv)
                     d2ldmdv},
    d2ldddv = function(y,mu,sigma,nu,bd) {dldd <- ifelse(y==0, 0, BB()$dldd(y,mu,sigma,bd) + (dBB(0,mu,sigma,bd)*BB()$dldd(0,mu,sigma,bd))/(1-dBB(0,mu,sigma,bd)))
                        dldv <- ifelse(y==0, 1/nu, -1/(1-nu))
                     d2ldddv <- -(dldd*dldv)
                     d2ldddv},
    G.dev.incr = function(y,mu,sigma,nu,bd,...) -2*dZABB(y,mu,sigma,nu,bd,log=TRUE),  # the error was here nu and sigma wrong way round!
         rqres = expression(rqres(pfun="pZABB", type="Discrete", ymin=0, y=y, mu=mu, sigma=sigma, nu=nu, bd=bd)),
#            mu.initial = expression (   mu <- (y+0.5)/(bd+1)),      
#         sigma.initial = expression (sigma <- rep(1,length(y))),
    mu.initial = expression(mu <- rep(0.5,length(y))),
 sigma.initial = expression(sigma <- rep(0.2,length(y))),
    nu.initial = expression(nu <-rep(0.3, length(y))),
      mu.valid = function(mu) all(mu > 0 & mu < 1),
      sigma.valid = function(sigma)  all(sigma > 0), 
      nu.valid = function(nu) all(nu > 0 & nu < 1), 
       y.valid = function(y)  all(y >= 0)),
        class = c("gamlss.family","family"))
}
#------------------------------------------------------------------------------------------
dZABB <- function(x, mu = 0.5, sigma = 0.1, nu = 0.1, bd = 1, log = FALSE)
 { 
          if (any(mu <= 0) |  any(mu >= 1) )  stop(paste("mu must be between 0 and 1", "\n", "")) 
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
          if (any(nu <= 0) |  any(nu >= 1) )  stop(paste("nu must be between 0 and 1", "\n", ""))
          if (any(x < 0) )  stop(paste("x must be 0 or greater than 0", "\n", ""))   
          logfy <- rep(0, length(x))
          logfy <- ifelse((x==0), log(nu), log(1-nu)+dBB(x,mu,sigma,bd,log=TRUE)-log(1-dBB(0,mu,sigma,bd)))          
          if(log == FALSE) fy <- exp(logfy) else fy <- logfy
          fy

  }
#------------------------------------------------------------------------------------------
pZABB <- function(q, mu = 0.5, sigma = 0.1, nu = 0.1, bd = 1, lower.tail = TRUE, log.p = FALSE)
  {     
         if (any(mu <= 0) |  any(mu >= 1) )  stop(paste("mu must be between 0 and 1", "\n", "")) 
         if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
         if (any(nu <= 0) |  any(nu >= 1) )  stop(paste("nu must be between 0 and 1", "\n", "")) 
         if (any(q < 0) )  stop(paste("y must be 0 or greater than 0", "\n", ""))  
         cdf <- rep(0,length(q))
         cdf1 <- pBB(q, mu, sigma, bd, lower.tail = TRUE, log.p = FALSE)
         cdf2 <- pBB(0, mu, sigma, bd, lower.tail = TRUE, log.p = FALSE)
         cdf3 <- nu+((1-nu)*(cdf1-cdf2)/(1-cdf2))
         cdf <- ifelse((q==0), nu,  cdf3)
         cdf <- ifelse(cdf>1L, 1L , cdf)
         if(lower.tail == TRUE) cdf <- cdf else cdf <-1-cdf
         if(log.p==FALSE) cdf <- cdf else cdf <- log(cdf)    
         cdf
   }
#------------------------------------------------------------------------------------------
qZABB <- function(p, mu = 0.5, sigma = 0.1, nu = 0.1, bd = 1, lower.tail = TRUE, log.p = FALSE)
  {      
         if (any(mu <= 0) |  any(mu >= 1) )  stop(paste("mu must be between 0 and 1", "\n", ""))   
         if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
         if (any(nu <= 0) |  any(nu >= 1) )  stop(paste("nu must be between 0 and 1", "\n", ""))        
         if (any(p < 0) | any(p > 1))  stop(paste("p must be between 0 and 1", "\n", "")) 
         if (log.p == TRUE) p <- exp(p)   else p <- p
         if (lower.tail == TRUE)  p <- p  else p <- 1 - p
          pnew  <- (p-nu)/(1-nu)- (1e-7)
          pnew <- ifelse((pnew<0),0,pnew)  #  corrected 16-04-10
          pnew2 <- pBB(0, mu, sigma, bd, lower.tail = TRUE, log.p = FALSE)*(1-pnew) + pnew 
          suppressWarnings(q <- ifelse((pnew > 0 ), qBB(pnew2, mu, sigma, bd, ), 0))
          q
   }
#------------------------------------------------------------------------------------------
rZABB <- function(n, mu = 0.5, sigma = 0.1, nu = 0.1, bd = 1)
  { 
    if (any(mu <= 0) |  any(mu >= 1) )  stop(paste("mu must be between 0 and 1", "\n", ""))    
    if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
    if (any(nu <= 0) |  any(nu >= 1) )  stop(paste("nu must be between 0 and 1", "\n", ""))    
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
          n <- ceiling(n)
          p <- runif(n)
          r <- qZABB(p, mu = mu, sigma = sigma, nu = nu, bd = bd)
          r
  }
#-----------------------------------------------------------------------------------------
