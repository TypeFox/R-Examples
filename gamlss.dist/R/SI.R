as.double# as a function of mean=ksi and sigma=tau=1/omega
# Modified by Mikis Monday, March 21, 2005 
# Modified by Bob 1/12/2004
# Wednesday, October 27, 2004 at 14:34
# Modified by Kalliope Akantziliotou 
# Uses the Dean and Lawless iterative algorithm and for the derivative of nu
# uses numerical derivative 
#----------------------------------------------------------------------------------------
SI <-function (mu.link ="log", sigma.link="log", nu.link="identity") 
{
     mstats <- checklink("mu.link", "Sichel", substitute(mu.link), 
                         c("1/mu^2", "log", "identity"))
    dstats <- checklink("sigma.link", "Sichel", substitute(sigma.link), #
                         c("inverse", "log", "identity"))
    vstats <- checklink("nu.link", "Sichel",substitute(nu.link), 
                         c("1/nu^2", "log", "identity"))    
    structure(
          list(family = c("SI", "Sichel"),
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
                 dldm = function(y,mu,sigma,nu) {
                         ty <- tofyS(y=y, mu=mu, sigma=sigma, nu=nu, what=1)
                       dldm <- (y-ty)/mu
                       dldm}, 
               d2ldm2 = function(y,mu,sigma,nu) {
                        ty <- tofyS(y=y, mu=mu, sigma=sigma,nu=nu, what=1)
                      dldm <- (y-ty)/mu
                    d2ldm2 <- - dldm * dldm
                    d2ldm2
                                    },
                 dldd = function(y,mu,sigma,nu) {
                        ty <- tofyS(y=y, mu=mu, sigma=sigma,nu=nu, what=1)
                     lbes2 <- log(besselK((1/sigma),nu+1))-log(besselK((1/sigma),nu))
                      dldd <- ((ty*(1+sigma*mu)/mu) - (sigma*y)-exp(lbes2))/(sigma^2)
                                 dldd},
               d2ldd2 = function(y,mu,sigma,nu){
                       ty  <- tofyS(y=y, mu=mu, sigma=sigma,nu=nu, what=1)
                     lbes2 <- log(besselK((1/sigma),nu+1))-log(besselK((1/sigma),nu))
                      dldd <- ((ty*(1+sigma*mu)/mu) - (sigma*y)-exp(lbes2))/(sigma^2)
                    d2ldd2 <- -dldd*dldd
                    d2ldd2
                                    },
              #d2ldmdd = function() 0,
              d2ldmdd = function(y,mu,sigma,nu) {
                         ty <- tofyS(y=y, mu=mu, sigma=sigma, nu=nu, what=1)
                       dldm <- (y-ty)/mu
                      lbes2 <- log(besselK((1/sigma),nu+1))-log(besselK((1/sigma),nu))
                       dldd <- ((ty*(1+sigma*mu)/mu) - (sigma*y)-exp(lbes2))/(sigma^2)
                    d2ldmdd <- -dldm *dldd
                    d2ldmdd
                                    }, 
              d2ldmdv = function(y,mu,sigma,nu) {
                   ty <- tofyS(y=y, mu=mu, sigma=sigma, nu=nu, what=1)
                 dldm <- (y-ty)/mu
                 #calculates the dldv                 
                   nd <- numeric.deriv(dSI(y, mu, sigma, nu, log=TRUE), "nu", delta=0.0001)
                 dldv <- as.vector(attr(nd, "gradient"))
                #calculates the d2ldmdv
              d2ldmdv <- -dldm *dldv
              d2ldmdv
                                    }, 
              d2ldddv = function(y,mu,sigma,nu) {
                  ty <- tofyS(y=y, mu=mu, sigma=sigma,nu=nu, what=1)
               lbes2 <- log(besselK((1/sigma),nu+1))-log(besselK((1/sigma),nu))
                dldd <- ((ty*(1+sigma*mu)/mu) - (sigma*y)-exp(lbes2))/(sigma^2)
                #calculates the dldv
                  nd <- numeric.deriv(dSI(y, mu, sigma, nu, log=TRUE), "nu", delta=0.0001)
                dldv <- as.vector(attr(nd, "gradient"))
                #calculates the d2ldddv 
             d2ldddv <- -dldd *dldv
             d2ldddv
                                 },               
                 dldv = function(y,mu,sigma,nu) {                           
                   nd <- numeric.deriv(dSI(y, mu, sigma, nu, log=TRUE), "nu", delta=0.0001)
                 dldv <- as.vector(attr(nd, "gradient"))
                 dldv
                                    },
               d2ldv2 = function(y,mu,sigma,nu){
                #  delta  <- 0.01
                #    dldv <- (dSI(y=y, mu=mu, sigma=sigma, nu=nu+delta, log = TRUE)- 
                #            dSI(y=y, mu=mu, sigma=sigma, nu=nu, log= TRUE))/ delta   
                   nd <- numeric.deriv(dSI(y, mu, sigma, nu, log=TRUE), "nu", delta=0.0001)
                 dldv <- as.vector(attr(nd, "gradient"))
               d2ldv2 <- -dldv*dldv
               d2ldv2
                                  } ,
          G.dev.incr  = function(y,mu,sigma,nu, pw=1,..) -2*dSI(y, mu, sigma, nu, log=TRUE),
                rqres = expression(
                  rqres(pfun="pSI", type="Discrete", ymin=0, y=y, mu=mu, sigma=sigma, nu=nu)
                                   ), #
            mu.initial = expression(mu<- y+0.5),
         sigma.initial = expression(
                      sigma <- rep( max( ((var(y)-mean(y))/(mean(y)^2)),0.1),length(y))),
            nu.initial = expression({  nu <- rep(-0.5,length(y)) }), 
              mu.valid = function(mu) all(mu > 0) , 
           sigma.valid = function(sigma)  all(sigma > 0), 
              nu.valid = function(nu) TRUE,  
               y.valid = function(y)  all(y >= 0)
          ),
            class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------
# the old tofyS replaced by new tofySS
#tofyS <- function (y, mu, sigma, nu, bsum=FALSE, ...)
#{
#    ty <- rep(0,length(y))
#    sumlty <- rep(0,length(y)) 
#    for (z1 in length(y):1)
#    {
#    #browser()
#         lyp1 <- Dum <- swi <- dum1 <- alpha <- lbes <- 0
#         lyp1 <- y[z1]+1 
#        tynew <- rep(0,lyp1)
#        alpha <- sqrt(1+2*sigma[z1]*mu[z1])/sigma[z1]
#        lbes <- log(besselK(alpha,nu[z1]+1))-log(besselK((alpha),nu[z1]))
#        tynew[1] <- dum1 <- mu[z1]*((1+2*sigma[z1]*mu[z1])^(-0.5))*exp(lbes) 
#        dum <- ifelse(lyp1==1, 1,2)
#        j <- dum:lyp1 
#        for (i in j)
#        {
#            if (i !=1)
#            tynew[i] <- (sigma[z1]*(2*(i-1+nu[z1])/mu[z1])+(1/tynew[i-1]))*(mu[z1]/(sigma[z1]*alpha))**2
#            tynew[1] <- dum1 
#           ty[z1] <- tynew[lyp1] 
#        }
#        if (bsum ) 
#            sumlty[z1] <- ifelse(lyp1==1, 0 , sum(log(tynew))-log(tynew[i]))
#    }
#    result <- cbind(ty, sumlty)
#}
#----------------------------------------------------------------------------------------
tofyS <- function (y, mu, sigma, nu, what=1)
   {
         ly <- length(y)       
      sigma <- rep(sigma, length = ly)
         mu <- rep(mu, length = ly)   
         nu <- rep(nu, length = ly) 
      alpha <- sqrt(1+2*sigma*mu)/sigma
       lbes <-  log(besselK(alpha,nu+1))-log(besselK((alpha),nu))
     sumlty <- as.double(.C("tofySI1", as.double(y), as.double(mu),as.double(sigma), as.double(nu),  
                            as.double(lbes),ans=double(ly), as.integer(ly), as.integer(max(y)+1), PACKAGE="gamlss.dist")$ans)
   sumlty
   }
#----------------------------------------------------------------------------------------
dSI<-function(x, mu=0.5, sigma=0.02, nu=-0.5, log=FALSE)
  { 
   if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
   if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
   if (any(x < 0) )  stop(paste("x must be >=0", "\n", ""))  
    ly <- max(length(x),length(mu),length(sigma),length(nu)) 
     x <- rep(x, length = ly)        
 sigma <- rep(sigma, length = ly)
    mu <- rep(mu, length = ly)   
    nu <- rep(nu, length = ly) 
 alpha <- sqrt(1+2*sigma*mu)/sigma
  lbes <-  log(besselK(alpha,nu+1))-log(besselK((alpha),nu))
sumlty <- as.double(.C("tofySI2", as.double(x), as.double(mu), as.double(sigma), as.double(nu), 
                       as.double(lbes), ans=double(ly),  as.integer(ly),as.integer(max(x)+1), PACKAGE="gamlss.dist")$ans)
logfy <- -lgamma(x+1)-nu*log(sigma*alpha)+sumlty+log(besselK(alpha,nu))-log(besselK((1/sigma),nu))
  if(log==FALSE) fy <- exp(logfy) else fy <- logfy
  fy
  }
#----------------------------------------------------------------------------------------  
pSI <- function(q, mu=0.5, sigma=0.02, nu=-0.5, lower.tail = TRUE, log.p = FALSE)
  {     
  #--------------------------------------------------
  tocdfS <- function (y, mu, sigma, nu, bsum=TRUE, ...)
  {
      ly <- length(y)       
   sigma <- rep(sigma, length = ly)
      mu <- rep(mu, length = ly)   
      nu <- rep(nu, length = ly) 
      ty <- rep(0,length(y))
     cdf <- rep(0,length(y))
   alpha <- sqrt(1+2*sigma*mu)/sigma
   lbes <-  log(besselK(alpha,nu+1))-log(besselK((alpha),nu)) 
    for (i in 1:length(y))
    {
         lyp1 <- y[i]+1 
        tynew <- lpnew <- rep(0,lyp1)
        tynew[1] <- mu[i]*((1+2*sigma[i]*mu[i])^(-0.5))*exp(lbes[i])
        lpnew[1] <- -nu[i]*log(sigma[i]*alpha[i])+log(besselK(alpha[i],nu[i]))-
                     log(besselK(1/sigma[i],nu[i])) 
        dum <- ifelse(lyp1==1, 1,2)
        for (j in dum:lyp1)
        {
            if (j !=1)
             {
            tynew[j] <- (sigma[i]*(2*(j-1+nu[i])/mu[i])+(1/tynew[j-1]))*
                         (mu[i]/(sigma[i]*alpha[i]))**2
            lpnew[j] <- lpnew[j-1] + log(tynew[j-1]) - log(j-1)
             }           
            ty[i] <- tynew[lyp1] 
        }
        if (bsum ) 
            cdf[i] <- sum(exp(lpnew))
    }
    cdf
  }
  #-----------------------------------------------
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
          if (any(q < 0) )  stop(paste("y must be >=0", "\n", ""))  
          ly <- max(length(q),length(mu),length(sigma),length(nu)) 
           q <- rep(q,length = ly)         
       sigma <- rep(sigma, length = ly)
          mu <- rep(mu, length = ly)   
          nu <- rep(nu, length = ly)   
         cdf <- tocdfS(y=q, mu=mu, sigma=sigma, nu=nu, bsum=TRUE)
          if(lower.tail==TRUE) cdf <- cdf else cdf=1-cdf
          if(log.p==FALSE) cdf <- cdf else cdf <- log(cdf)                                                                    
          cdf
   }
#----------------------------------------------------------------------------------------
qSI <- function(p, mu=0.5, sigma=0.02, nu=-0.5,  lower.tail = TRUE, log.p = FALSE,  
                 max.value = 10000)
  {      
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
          if (any(p < 0) | any(p > 1.0001))  stop(paste("p must be between 0 and 1", "\n", "")) 
          if (log.p==TRUE) p <- exp(p) else p <- p
          if (lower.tail==TRUE) p <- p else p <- 1-p    
           ly <- max(length(p),length(mu),length(sigma),length(nu)) 
            p <- rep(p, length = ly)                                                           
          QQQ <- rep(0, length = ly)                         
       nsigma <- rep(sigma, length = ly)
          nmu <- rep(mu, length = ly)                
          nnu <- rep(nu, length = ly) 
        for (i in seq(along=p))                                                          
      {
       cumpro <- 0                                                                         
     if (p[i]+0.000000001 >= 1) QQQ[i] <- Inf
     else  
        {  
            for (j in seq(from = 0, to = max.value))
            {
            cumpro <-   pSI(j, mu = nmu[i], sigma = nsigma[i], nu = nnu[i], log.p = FALSE)  
            QQQ[i] <- j 
       if  (p[i] <= cumpro ) break 
            } 
        }
      }          
          QQQ   
   }
#----------------------------------------------------------------------------------------
rSI <- function(n, mu=0.5, sigma=0.02, nu=-0.5)
  { 
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
          if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
          n <- ceiling(n)
          p <- runif(n)
          r <- qSI(p, mu=mu, sigma=sigma, nu=nu)
          r
  }
#----------------------------------------------------------------------------------------
