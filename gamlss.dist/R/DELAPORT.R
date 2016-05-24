 DEL <- function (mu.link ="log", sigma.link="log", nu.link="logit") 
{
     mstats <- checklink("mu.link", "DEL", substitute(mu.link), 
                         c("1/mu^2", "log", "identity"))
    dstats <- checklink("sigma.link", "DEL", substitute(sigma.link), #
                         c("inverse", "log", "identity"))
    vstats <- checklink("nu.link", "DEL", substitute(nu.link),
                         c("logit", "probit", "cloglog", "cauchit", "log"))  
    structure(
          list(family = c("DEL", "Delaporte"),
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
                 dldm = function(y,mu,sigma,nu) 
                       {
                     logty <-log(y+1)+dDEL(y+1, mu=mu, sigma=sigma, nu=nu, log=TRUE)-
                                      dDEL(y, mu=mu, sigma=sigma, nu=nu, log=TRUE)
                       ty  <- exp(logty)
                      dldm <- (y-ty)/mu
                      dldm
                       }, 
               d2ldm2 = function(y,mu,sigma,nu) 
                       {
                     logty <-log(y+1)+dDEL(y+1, mu=mu, sigma=sigma, nu=nu, log=TRUE)-
                                      dDEL(y, mu=mu, sigma=sigma, nu=nu, log=TRUE)
                       ty  <- exp(logty)
                      dldm <- (y-ty)/mu
                    d2ldm2 <- - dldm * dldm
                    d2ldm2
                        },
                 dldd = function(y,mu,sigma,nu) 
                        {
                     #    S <- calcS(y = y, mu = mu, sigma = sigma, nu = nu)
                     #    U <- calcU(y = y, mu = mu, sigma = sigma, nu = nu)
                     # dldd <- digamma(1/sigma)/(sigma*sigma)-(mu*(1-nu))/(sigma*(1+mu*sigma*(1-nu)))
                     # dldd <- dldd+log(1+mu*sigma*(1-nu))/(sigma^2)+U/S
                        nd <- numeric.deriv(dDEL(y, mu, sigma, nu, log=TRUE), "sigma", delta=0.01)
                      dldd <- as.vector(attr(nd, "gradient"))
                      dldd
                        },
               d2ldd2 = function(y,mu,sigma,nu)
                        {
                      #S <- calcS(y = y, mu = mu, sigma = sigma, nu = nu)
                      #   U <- calcU(y = y, mu = mu, sigma = sigma, nu = nu)
                      #dldd <- digamma(1/sigma)/(sigma*sigma)-(mu*(1-nu))/(sigma*(1+mu*sigma*(1-nu)))
                      #dldd <- dldd+log(1+mu*sigma*(1-nu))/(sigma^2)+U/S
                     nd <- numeric.deriv(dDEL(y, mu, sigma, nu, log=TRUE), "sigma", delta=0.01)
                   dldd <- as.vector(attr(nd, "gradient"))
                 d2ldd2 <- -dldd*dldd
                 d2ldd2
                        },
              #d2ldmdd = function() 0,
              d2ldmdd = function(y,mu,sigma,nu) 
                       {
                       
                       logty <-log(y+1)+dDEL(y+1, mu=mu, sigma=sigma, nu=nu, log=TRUE)-
                                      dDEL(y, mu=mu, sigma=sigma, nu=nu, log=TRUE)
                       ty  <- exp(logty)
                      dldm <- (y-ty)/mu
                        nd <- numeric.deriv(dDEL(y, mu, sigma, nu, log=TRUE), "sigma", delta=0.01)
                      dldd <- as.vector(attr(nd, "gradient"))
                     #     c <- exp(log(besselK((1/sigma),nu+1))-log(besselK((1/sigma),nu)))
                     #  dcdd <- (c*sigma*(2*nu+1)+1-c*c)/(sigma*sigma)   
                     #  dldd <- (((ty*(c+sigma*mu)/mu) - (sigma*y)-c)/(sigma^2))+dcdd*(ty-y)/c
                    d2ldmdd <- -dldm *dldd
                    d2ldmdd
                       }, 
              d2ldmdv = function(y,mu,sigma,nu) 
                       {
                  logty <-log(y+1)+dDEL(y+1, mu=mu, sigma=sigma, nu=nu, log=TRUE)-
                                      dDEL(y, mu=mu, sigma=sigma, nu=nu, log=TRUE)
                       ty  <- exp(logty)
                      dldm <- (y-ty)/mu
                 #calculates the dldv                 
                   nd <- numeric.deriv(dDEL(y, mu, sigma, nu, log=TRUE), "nu", delta=0.01)
                 dldv <- as.vector(attr(nd, "gradient"))
                #calculates the d2ldmdv
              d2ldmdv <- -dldm *dldv
              d2ldmdv
                       }, 
              d2ldddv = function(y,mu,sigma,nu) 
                       {
                  # ty <- tofyDEL(y=y, mu=mu, sigma=sigma,nu=nu, what=1)
                  #  c <- exp(log(besselK((1/sigma),nu+1))-log(besselK((1/sigma),nu)))
                 #dcdd <- (c*sigma*(2*nu+1)+1-c*c)/(sigma*sigma)   
                 #dldd <- (((ty*(c+sigma*mu)/mu) - (sigma*y)-c)/(sigma^2))+dcdd*(ty-y)/c
                #calculates the dldv
                   nd <- numeric.deriv(dDEL(y, mu, sigma, nu, log=TRUE), "sigma", delta=0.01)
                 dldd <- as.vector(attr(nd, "gradient"))
                   nd <- numeric.deriv(dDEL(y, mu, sigma, nu, log=TRUE), "nu", delta=0.01)
                 dldv <- as.vector(attr(nd, "gradient"))
                #calculates the d2ldddv 
              d2ldddv <- -dldd *dldv
              d2ldddv
                       },               
                 dldv = function(y,mu,sigma,nu) 
                       {                           
                   nd <- numeric.deriv(dDEL(y, mu, sigma, nu, log=TRUE), "nu", delta=0.01)
                 dldv <- as.vector(attr(nd, "gradient"))
                 dldv
                       },
               d2ldv2 = function(y,mu,sigma,nu)
                       {
                   nd <- numeric.deriv(dDEL(y, mu, sigma, nu, log=TRUE), "nu", delta=0.01)
                 dldv <- as.vector(attr(nd, "gradient"))
               d2ldv2 <- -dldv*dldv
               d2ldv2
                       } ,
          G.dev.incr  = function(y,mu,sigma,nu, pw=1,..) -2*dDEL(y, mu, sigma, nu, log=TRUE),
                rqres = expression(
                 rqres(pfun="pDEL", type="Discrete", ymin=0, y=y, mu=mu, sigma=sigma, nu=nu)
                                   ),  #
            mu.initial = expression(mu<- (y+mean(y)/2)),
         sigma.initial = expression(
                      sigma <- rep( max( ((var(y)-mean(y))/(mean(y)^2)),0.1),length(y))),
            nu.initial = expression({  nu <- rep(0.5,length(y)) }), 
              mu.valid = function(mu) all(mu > 0) , 
           sigma.valid = function(sigma)  all(sigma > 0), 
              nu.valid = function(nu) all(nu > 0) && all(nu < 1),    
               y.valid = function(y)  all(y >= 0)
          ),
            class = c("gamlss.family","family"))
}
#----------------------------------------------------------------------------------------
dDEL<-function(x, mu=1, sigma=1, nu=.5, log=FALSE)
  {
   if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
   if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
   if (any(nu <= 0) | any(nu >= 1))  stop(paste("nu must be between 0 and 1", "\n", "")) 
   if (any(x < 0) )  stop(paste("x must be >=0", "\n", ""))  
 ly <- max(length(x),length(mu),length(sigma),length(nu)) 
     x <- rep(x, length = ly)      
 sigma <- rep(sigma, length = ly)
    mu <- rep(mu, length = ly)   
    nu <- rep(nu, length = ly) 
 logpy0 <- -mu*nu-(1/sigma)*(log(1+mu*sigma*(1-nu)))
   S <- tofyDEL2(x, mu, sigma, nu)
  # S <- tofyDELPORT(x, mu, sigma, nu)[,2]
 logfy <-  logpy0-lgamma(x+1)+S
 if(log==FALSE) fy <- exp(logfy) else fy <- logfy
 
 if (length(sigma)>1) fy <- ifelse(sigma>0.0001, fy, 
                                          dPO(x, mu = mu, log = log) )
        else fy <- if (sigma<0.0001) dPO(x, mu = mu, log = log) 
                   else fy 
  fy
  }
#----------------------------------------------------------------------------------------
# tofyDEL <- function (y, mu, sigma, nu, what=1)
#    {
#   ly <- max(length(y),length(mu),length(sigma),length(nu)) 
#   y <- rep(y, length = ly)    
#   sigma <- rep(sigma, length = ly)
#   mu <- rep(mu, length = ly)   
#   nu <- rep(nu, length = ly) 
#   sumlty <- as.double(.C(ifelse(what==1,"tofydel1","tofydel2"), 
#               as.double(y), as.double(mu), 
#               as.double(sigma), as.double(nu), 
#               ans=double(ly), as.integer(length(y)),
#               as.integer(max(y)+1), PACKAGE="gamlss.dist")$ans)
#   
#   sumlty
#    }
#-------------------------------------------------------------------------------
 tofyDEL1 <- function (y, mu, sigma, nu)
 {
   ly <- max(length(y),length(mu),length(sigma),length(nu)) 
   y <- rep(y, length = ly)    
   sigma <- rep(sigma, length = ly)
   mu <- rep(mu, length = ly)   
   nu <- rep(nu, length = ly) 
   sumlty <- as.double(.C("tofydel1", 
                          as.double(y), as.double(mu), 
                          as.double(sigma), as.double(nu), 
                          ans=double(ly), as.integer(length(y)),
                          as.integer(max(y)+1), PACKAGE="gamlss.dist")$ans)
   
   sumlty
 }
#-------------------------------------------------------------------------------
 tofyDEL2 <- function (y, mu, sigma, nu)
 {
   ly <- max(length(y),length(mu),length(sigma),length(nu)) 
   y <- rep(y, length = ly)    
   sigma <- rep(sigma, length = ly)
   mu <- rep(mu, length = ly)   
   nu <- rep(nu, length = ly) 
   sumlty <- as.double(.C("tofydel2", 
                          as.double(y), as.double(mu), 
                          as.double(sigma), as.double(nu), 
                          ans=double(ly), as.integer(length(y)),
                          as.integer(max(y)+1), PACKAGE="gamlss.dist")$ans)
   
   sumlty
 }
 #-------------------------------------------------------------------------------
# tofyDELPORT <- function (y, mu, sigma, nu, what=1) 
#   {
#         ly <- max(length(y),length(mu),length(sigma),length(nu)) 
#          y <- rep(y, length = ly)    
#      sigma <- rep(sigma, length = ly)
#         mu <- rep(mu, length = ly)   
#         nu <- rep(nu, length = ly) 
#         ty <- sumlty <- rep(0,length(y)) 
#    for (i in 1:length(y))
#    {
#         iy <- dum <- 0
#         iy <- y[i]+1 
#      tofyn <- rep(0,iy)
#       tofyn[1] <- mu[i]*nu[i]+mu[i]*(1-nu[i])/(1+mu[i]*sigma[i]*(1-nu[i])) 
#        if (iy<2)
#        {
#          ty[i] <- tofyn[iy] 
#      sumlty[i] <- 0
#        }
#        else
#        {
#        for (j in 2:iy)
#         {
#           dum = (1+(1/(mu[i]*sigma[i]*(1-nu[i]))))
#             tofyn[j] <- ((j-1)+mu[i]*nu[i]+(1/(sigma[i]*(1-nu[i])))-(mu[i]*nu[i]*(j-1))/ tofyn[j-1])/dum      
#         }
#              
#            ty[i] <- tofyn[iy] 
#            sumlty[i] <- sum(log(tofyn[1:iy-1]))
#       }     
#    }
#   result <- cbind(ty, sumlty)
#   result        
#   }   
#----------------------------------------------------------------------------------------
pDEL <- function(q, mu = 1, sigma = 1, nu = .5, lower.tail = TRUE, log.p = FALSE)
  {     
   if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
   if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
   if (any(nu <= 0) | any(nu >= 1))  stop(paste("nu must be between 0 and 1", "\n", "")) 
   if (any(q < 0) )  stop(paste("q must be >=0", "\n", ""))    
        ly <- length(q)                                                       
       FFF <- rep(0,ly)                         
    nsigma <- rep(sigma, length = ly)
       nmu <- rep(mu, length = ly) 
       nnu <- rep(nu, length = ly)                                                        
         j <- seq(along=q) 
   for (i in j)                                                          
      {                                                                 
        y.y <- q[i]                                                   
         nn <- nnu[i]                                                  
         mm <- nmu[i]
       nsig <- nsigma[i]                                                     
     allval <- seq(0,y.y)
     pdfall <- dDEL(allval, mu = mm, sigma = nsig, nu = nn, log = FALSE)
     FFF[i] <- sum(pdfall)                                             
      }  
      cdf <- FFF
      cdf <- if(lower.tail==TRUE) cdf else 1-cdf
      cdf <- if(log.p==FALSE) cdf else log(cdf)                                                                    
      cdf
  }
#----------------------------------------------------------------------------------------
qDEL <- function(p, mu=1, sigma=1, nu=0.5,  lower.tail = TRUE, log.p = FALSE,  
                 max.value = 10000)
  {      
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
          if (any(nu <= 0) | any(nu >= 1))  stop(paste("nu must be between 0 and 1", "\n", "")) 
          if (any(p < 0) | any(p > 1.0001))  stop(paste("p must be between 0 and 1", "\n", "")) 
          if (log.p==TRUE) p <- exp(p) else p <- p
          if (lower.tail==TRUE) p <- p else p <- 1-p    
           ly <- length(p)                                                       
          QQQ <- rep(0,ly)                         
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
            cumpro <-  pDEL(j, mu = nmu[i], sigma = nsigma[i], nu = nnu[i], log.p = FALSE) 
                       # else  cumpro+dSICHEL(j, mu = nmu[i], sigma = nsigma[i], nu = nnu[i], log = FALSE)# the above is faster 
           QQQ[i] <- j 
       if  (p[i] <= cumpro ) break 
            } 
        }
      }          
          QQQ   
   } 
#----------------------------------------------------------------------------------------
rDEL <- function(n, mu=1, sigma=1, nu=0.5, max.value = 10000)
  { 
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
          if (any(nu <= 0) | any(nu >= 1))  stop(paste("nu must be between 0 and 1", "\n", "")) 
          if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
          n <- ceiling(n)
          p <- runif(n)
          r <- qDEL(p, mu=mu, sigma=sigma, nu=nu, max.value = max.value )
          r
  }
  
  
#dyn.load("c:/gamlss/fortran/tofydel.dll")
#gfortran -O3 -c tofydel.f -o tofydel.0
#gcc -std-gnu99 -shared -s -o tofydel.dll tofydel.o -lgfortran 
