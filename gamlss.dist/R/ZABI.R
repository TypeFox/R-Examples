# Wednesday, May 5, 2010
# ZABI
ZABI <- function (mu.link = "logit", sigma.link = "logit")
{
    mstats <- checklink("mu.link", "ZABI", substitute(mu.link),    
                           c("logit", "probit", "cloglog", "cauchit", "log", "own"))   
    dstats <- checklink("sigma.link", "ZABI", substitute(sigma.link), 
                           c("logit", "probit", "cloglog", "cauchit", "log", "own"))   
    structure(
          list(family = c("ZABI", "Zero Adjusted Binomial"),
           parameters = list(mu=TRUE, sigma=TRUE), 
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
     dldm = function(y,mu,sigma,bd) {dldm <- ifelse(y==0, 0, BI()$dldm(y,mu,bd) + (dBI(0,bd,mu)*BI()$dldm(0,mu,bd))/(1-dBI(0,bd,mu)))
                         dldm}, 
    d2ldm2 = function(y,mu,sigma,bd) {dldm <- ifelse(y==0, 0, BI()$dldm(y,mu,bd) + (dBI(0,bd,mu)*BI()$dldm(0,mu,bd))/(1-dBI(0,bd,mu)))
                          d2ldm2 <- ifelse(y==0, 0, -(dldm)^2)
                         d2ldm2},
     dldd = function(y,mu,sigma,bd) {dldd <- ifelse(y==0, 1/sigma, -1/(1-sigma))
                         dldd}, 
    d2ldd2 = function(y,mu,sigma,bd) {d2ldd2 <- -1/(sigma*(1-sigma))
                         d2ldd2},
  d2ldmdd = function(y,mu,sigma,bd) {d2ldmdd <- 0
                      d2ldmdd},
 G.dev.incr  = function(y,mu,sigma,bd,...) -2*dZABI(y,bd,mu,sigma,log=TRUE),                      
         rqres = expression(rqres(pfun="pZABI", type="Discrete", ymin=0, y=y, bd=bd, mu=mu, sigma=sigma)) ,
    mu.initial = expression(mu <- rep(0.5,length(y))),         #(((y + 0.5)/(bd +1))+0.5)/2}), 
   sigma.initial = expression(sigma <-rep(0.3, length(y))), 
      mu.valid = function(mu) all(mu > 0 & mu < 1),
   sigma.valid = function(sigma)  all(sigma > 0 & sigma < 1), 
       y.valid = function(y)  all(y >= 0)
          ),
            class = c("gamlss.family","family"))
}
#------------------------------------------------------------------------------------------
dZABI<-function(x, bd = 1, mu = 0.5, sigma = 0.1, log = FALSE)
 { 
 #browser()
          if (any(mu <= 0) |  any(mu >= 1) )  stop(paste("mu must be between 0 and 1", "\n", "")) 
          if (any(sigma <= 0) | any(sigma >= 1) )  stop(paste("sigma must be between 0 and 1", "\n", "")) 
          if (any(x < 0) )  stop(paste("x must be 0 or greater than 0", "\n", ""))   
          logfy <- rep(0, length(x))
          logfy <- ifelse((x==0), log(sigma), log(1-sigma)+dBI(x,bd,mu,log=TRUE)-log(1-dBI(0,bd,mu)))          
          if(log == FALSE) fy <- exp(logfy) else fy <- logfy
          fy

  }
#------------------------------------------------------------------------------------------
pZABI <- function(q,bd = 1, mu = 0.5, sigma = 0.1, lower.tail = TRUE, log.p = FALSE)
  {     
         if (any(mu <= 0) |  any(mu >= 1) )  stop(paste("mu must be between 0 and 1", "\n", ""))  
         if (any(sigma <= 0) | any(sigma >= 1) )  stop(paste("sigma must be between 0 and 1", "\n", "")) 
         if (any(q < 0) )  stop(paste("y must be 0 or greater than 0", "\n", ""))  
         cdf <- rep(0,length(q))
         cdf1 <- pbinom(q, size = bd, prob = mu, lower.tail = TRUE, log.p = FALSE)
         cdf2 <- pbinom(0, size = bd, prob = mu, lower.tail = TRUE, log.p = FALSE)
         cdf3 <- sigma+((1-sigma)*(cdf1-cdf2)/(1-cdf2))
         cdf <- ifelse((q==0),sigma,  cdf3)
         if(lower.tail == TRUE) cdf <- cdf else cdf <-1-cdf
         if(log.p==FALSE) cdf <- cdf else cdf <- log(cdf)    
         cdf
   }
#-----------------------------------------------------------------------------------------
qZABI <- function(p, bd = 1, mu = 0.5, sigma = 0.1, lower.tail = TRUE, log.p = FALSE)
  {      
         if (any(mu <= 0) |  any(mu >= 1) )  stop(paste("mu must be between 0 and 1", "\n", ""))           
         if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0", "\n", "")) 
         if (any(p < 0) | any(p > 1))  stop(paste("p must be between 0 and 1", "\n", "")) 
         if (log.p == TRUE) p <- exp(p)   else p <- p
         if (lower.tail == TRUE)  p <- p  else p <- 1 - p
          pnew  <- (p-sigma)/(1-sigma)-1e-10
          pnew2 <- pbinom(0, size = bd, prob = mu, lower.tail = TRUE, log.p = FALSE)*(1-pnew) + pnew 
          suppressWarnings(q <- ifelse((pnew > 0 ), qbinom(pnew2, size = bd, prob = mu, ), 0))
          q
   }
#-----------------------------------------------------------------------------------------
rZABI <- function(n, bd = 1, mu = 0.5, sigma=0.1)
  { 
    if (any(mu <= 0) |  any(mu >= 1) )  stop(paste("mu must be between 0 and 1", "\n", ""))          
    if (any(sigma <= 0) )  stop(paste("sigma must greated than 0", "\n", "")) 
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
          n <- ceiling(n)
          p <- runif(n)
          r <- qZABI(p, mu = mu, sigma = sigma, bd=bd)
          r
  }
#-----------------------------------------------------------------------------------------
