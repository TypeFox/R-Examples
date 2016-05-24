LNO <- function (mu.link="identity", sigma.link="log") 
{
    mstats <- checklink("mu.link", "Log Normal", substitute(mu.link), c("1/mu^2", "log", "identity"))
    dstats <- checklink("sigma.link", "Log Normal", substitute(sigma.link), c("inverse", "log", "identity"))
    
    structure(
          list(family = c("LNO", "Box-Cox"),
           parameters = list(mu=TRUE, sigma=TRUE, nu=FALSE),
                nopar = 3, 
                 type = "Continuous",
              mu.link = as.character(substitute(mu.link)), 
           sigma.link = as.character(substitute(sigma.link)), 
           mu.linkfun = mstats$linkfun, 
        sigma.linkfun = dstats$linkfun, 
           mu.linkinv = mstats$linkinv, 
        sigma.linkinv = dstats$linkinv,
                mu.dr = mstats$mu.eta, 
             sigma.dr = dstats$mu.eta, 
                  dldm = function(y,mu,sigma,nu) {
                                        cc <- ifelse((nu != 0),((y^nu-1)/nu),log(y))
                                        dldm <- (cc-mu)/sigma^2
                                        dldm
                                    },
               d2ldm2 = function(sigma) -1/sigma^2,
                 dldd = function(y,mu,sigma,nu)  {
                                         cc <- ifelse((nu != 0), ((y^nu-1)/nu),log(y))
                                        dldd <- (1/(sigma^3))*((cc-mu)^2-sigma^2)
                                        dldd
                                    },
               d2ldd2 = function(sigma) -2/sigma^2,
              d2ldmdd = function(y)  rep(0,length(y)),
          G.dev.incr  = function(y,mu,sigma,nu,...) { -2*dLNO(y,mu,sigma,nu,log=TRUE)
                                    }, 
                rqres = expression( rqres(pfun="pLNO", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)),
           mu.initial = expression({if (is.null(nu.start))
                                        nus<-rep(0,length(y)) 
                                    else
                                        {if(length(nu.start)>1) {nus <- nu.start}  else {nus <- rep(nu.start,length(y))}}
                                     cc <- ifelse((nus != 0), ((y^nus-1)/nus),log(y)) 
                                     mu <- (cc+mean(cc))/2  }),
         sigma.initial = expression({if (is.null(nu.start))
                                        {nus<-rep(0,length(y))}
                                     else
                                        {if(length(nu.start)>1) {nus <- nu.start}  else {nus <- rep(nu.start,length(y))}}
                                     cc <- ifelse((nus != 0), ((y^nus-1)/nus),log(y))
                                     sigma <- rep(sd(cc),length(y)) }), 
            nu.initial = expression({  nu <- rep(0,length(y)) }), 
              mu.valid = function(mu) all(mu > 0), 
           sigma.valid = function(sigma)  all(sigma > 0),
              nu.valid = function(nu) TRUE,  
               y.valid = function(y)  all(y > 0)
          ),
            class = c("gamlss.family","family"))
}
#-------------------------------------------
dLNO <- function(x, mu=1, sigma=0.1, nu=0,  log = FALSE)
 {
          if (any(nu!=0 & mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(x < 0))  stop(paste("x must be positive", "\n", ""))  
          if(length(nu)>1)  z <- ifelse(nu != 0,((x^nu)-1)/nu,log(x))
          else   if (nu != 0) z <- ((x^nu)-1)/nu else z <- log(x)
      loglik <- -0.5*log(2*pi)-log(sigma)-0.5*((z-mu)^2)/(sigma^2)+(nu-1)*log(x)
       if(log==FALSE) ft  <- exp(loglik) else ft <- loglik 
       ft
  }    
#--------------------------------------------------------------  
pLNO <- function(q, mu=1, sigma=0.1, nu=0,  lower.tail = TRUE, log.p = FALSE)
 {  
          if (any(nu!=0 & mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
          if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(q < 0))  stop(paste("q must be positive", "\n", ""))  
          if(length(nu)>1)  z <- ifelse(nu != 0,((q^nu)-1)/nu,log(q))
          else   if (nu != 0) z <- ((q^nu)-1)/nu else z <- log(q)
          Fy <- pnorm(z, mean=mu, sd=sigma)
         if(lower.tail==TRUE) Fy  <- Fy else  Fy <- 1-Fy 
         if(log.p==FALSE) Fy  <- Fy else  Fy<- log(Fy) 
         Fy     
 }
#--------------------------------------------------------------
qLNO <- function(p, mu=1, sigma=0.1, nu=0,  lower.tail = TRUE, log.p = FALSE )
 { 
    if (any(nu!=0 & mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (log.p==TRUE) p <- exp(p) else p <- p
    if (lower.tail==TRUE) p <- p else p <- 1-p
    if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))       
     z <- mu+sigma*qnorm(p) 
    if(length(nu)>1) 
         {
    ya <- ifelse((nu!=0),(nu*z+1)^(1/nu),exp(z)) 
         }
    else {
    ya <- if (nu!=0)  (nu*z+1)^(1/nu)
         else        exp(z)                    
         }                 
    ya
 }
#-----------------------------------------------------------------  
rLNO <- function(n, mu=1, sigma=0.1, nu=0)
  {
    if (any(nu!=0 & mu <= 0))  stop(paste("mu must be positive", "\n", "")) 
    if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
    if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))    
    n <- ceiling(n)
    p <- runif(n)
    r <- qLNO(p,mu=mu,sigma=sigma,nu=nu)
    r
  }
#----------------------------------------------------------------------------------------
