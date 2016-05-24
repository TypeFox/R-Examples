## the Poison Inverse Gaussian distribution
## KA RAR and MS
## will not work with CG()
##-----------------------------------------------------------------------------------------
PIG <- function (mu.link = "log", sigma.link = "log") 
{
    mstats <- checklink("mu.link", "Beta Binomial", substitute(mu.link),
                         c("inverse", "log", "identity", "sqrt"))   
    dstats <- checklink("sigma.link", "Beta Binomial", substitute(sigma.link), 
                        c("inverse", "log", "identity", "sqrt"))   
    structure(
          list(family = c("PIG","Poisson.Inverse.Gaussian"),
           parameters = list(mu = TRUE, sigma = TRUE),
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
                 dldm = function(y,mu,sigma) {
                   ty <- as.double(.C("tofyPIG1", as.double(y), as.double(mu), as.double(sigma),
                            ans=double(length(y)), as.integer(length(y)), as.integer(max(y)+1),
                            PACKAGE="gamlss.dist")$ans)
              dldm <- (y-ty)/mu
              dldm
                                    }, 
               d2ldm2 = function(y,mu,sigma) { 
                                    #d2ldm2 <- eval.parent(quote(-dldp*dldp))
                 ty <- as.double(.C("tofyPIG1", as.double(y), as.double(mu), as.double(sigma),
                            ans=double(length(y)), as.integer(length(y)), as.integer(max(y)+1),
                            PACKAGE="gamlss.dist")$ans)
                       dldm <- (y-ty)/mu
                     d2ldm2 <- -dldm*dldm
                     d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15) 
                     d2ldm2
                                    },
                 dldd = function(y,mu,sigma) {
                   ty <- as.double(.C("tofyPIG1", as.double(y), as.double(mu), as.double(sigma),
                            ans=double(length(y)), as.integer(length(y)), as.integer(max(y)+1),
                            PACKAGE="gamlss.dist")$ans)
             dldd <- ((ty*(1+sigma*mu)/mu)-(1+sigma*y))/(sigma^2)
              dldd
                                    },
               d2ldd2 = function(y,mu,sigma){
                                    #d2ldd2 <- eval.parent(quote(-dldp*dldp))
                 ty <- as.double(.C("tofyPIG1", as.double(y), as.double(mu), as.double(sigma),
                              ans=double(length(y)), as.integer(length(y)), as.integer(max(y)+1),
                              PACKAGE="gamlss.dist")$ans)
                     dldd <- ((ty*(1+sigma*mu)/mu)-(1+sigma*y))/(sigma^2)
                   d2ldd2 <- -dldd*dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
                   d2ldd2
                                   },
              d2ldmdd = function(y,mu,sigma) {
                ty <- as.double(.C("tofyPIG1", as.double(y), as.double(mu), as.double(sigma),
                                ans=double(length(y)), as.integer(length(y)), as.integer(max(y)+1),
                                PACKAGE="gamlss.dist")$ans)
                      dldm <- (y-ty)/mu
                      dldd <- ((ty*(1+sigma*mu)/mu)-(1+sigma*y))/(sigma^2)
                   d2ldmdd <- -dldm*dldd
                    d2ldmdd
                                    },
              
           G.dev.incr  = function(y,mu,sigma,pw=1,..) -2*dPIG(y, mu, sigma, log=TRUE),
                rqres = expression(
                  rqres(pfun="pPIG", type="Discrete", ymin=0, y=y, mu=mu, sigma=sigma)
                                   ), 
           mu.initial = expression(mu <- (y+mean(y))/2),
        sigma.initial = expression( sigma <- rep( max( ((var(y)-mean(y))/(mean(y)^2)),0.1),
                                        length(y)) ),
             mu.valid = function(mu) all(mu > 0) , 
          sigma.valid = function(sigma)  all(sigma > 0), 
              y.valid = function(y)  all(y >= 0)
          ),
                class = c("gamlss.family","family"))
}
##-----------------------------------------------------------------------------------------
# 12-3-06 it seems that this function do work with y=0
#t.of.y <- function (y, mu, sigma, bsum = FALSE, ...)
#{
#  tofy <- sumlty <- rep(0,length(y))
#  lofyp1 <- 0
# for (i in seq(along=y))
#    {
#         lofyp1 <- y[i]+1 
#        tofynew <- rep(0,lofyp1)
#       tofynew[1] <- mu[i]*((1+2*sigma[i]*mu[i])^(-0.5))       
#        if (lofyp1==1)
#         {
#          tofy[i] <- tofynew[lofyp1]
#          if (bsum)  sumlty[i] <-  0           
#         }
#        else
#         { 
#         for (j in 2:lofyp1)
#           {
#       tofynew[j] <- ((sigma[i]*(2*(j-1)-1)/mu[i])+(1/tofynew[j-1]))*(tofynew[1])^2
#           }
#          tofy[i] <- tofynew[lofyp1]
#        if (bsum)             sumlty[i] <- sum(log(tofynew))-log(tofynew[j])
#         }
#    }
# result <- cbind(tofy, sumlty)
# result
#}
#-----------------------------------------------------------------------------------------
dPIG<-function(x, mu = 0.5, sigma = 0.02 , log = FALSE)
 { 
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
          if (any(x < 0) )  stop(paste("x must be >=0", "\n", ""))  
          ly <- length(x)                                                                    
      nsigma <- rep(sigma, length = ly)
         nmu <- rep(mu, length = ly)
      sumlty <- as.double(.C("tofyPIG2", as.double(x), as.double(nmu), as.double(nsigma),
                             ans=double(ly), as.integer(length(x)), 
                             as.integer(max(x)+1), PACKAGE="gamlss.dist")$ans)
       logfy <- -lgamma(x+1)+(1-sqrt(1+2*sigma*mu))/sigma +sumlty
          if(log==FALSE) fy <- exp(logfy) else fy <- logfy
          fy
  }
##-----------------------------------------------------------------------------------------  
pPIG <- function(q, mu=0.5, sigma=0.02, lower.tail = TRUE, log.p = FALSE)
  {     
  ## function to calculate the cdf
  ## function to calculate the cdf
  if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
  if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
  if (any(q < 0) )  stop(paste("y must be >=0", "\n", ""))  
  lq <- length(q)                                                                    
  nsigma <- rep(sigma, length = lq)
  nmu <- rep(mu, length = lq)      
  cdf <-as.double(.C("tocdf", as.double(q), as.double(nmu), as.double(nsigma),
           ans=double(lq), as.integer(lq), PACKAGE="gamlss.dist")$ans)
  if(lower.tail==TRUE) cdf <- cdf else cdf=1-cdf
  if(log.p==FALSE) cdf <- cdf else cdf <- log(cdf)                                                                    
  cdf
   }
##-----------------------------------------------------------------------------------------
qPIG <- function(p, mu = 0.5, sigma = 0.02,  lower.tail = TRUE, log.p = FALSE, 
                 max.value = 10000)
  {      
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
          if (any(p < 0) | any(p > 1.0001))  stop(paste("p must be between 0 and 1", "\n", "")) 
          if (log.p==TRUE) p <- exp(p) else p <- p
          if (lower.tail==TRUE) p <- p else p <- 1-p    
           ly <- length(p)                                                       
          QQQ <- rep(0,ly)                         
       nsigma <- rep(sigma, length = ly)
          nmu <- rep(mu, length = ly)       
        for (i in seq(along=p))                                                          
      { 
       cumpro <- 0                                                                        
     if (p[i]+0.000000001 >= 1) QQQ[i] <- Inf
     else  
        {  
            for (j in seq(from = 0, to = max.value))
            {
           cumpro <-  pPIG(j, mu = nmu[i], sigma = nsigma[i] , log.p = FALSE)  
           QQQ[i] <- j 
       if  (p[i] <= cumpro ) break 
            } 
        }
      }
          QQQ   
   }
##-----------------------------------------------------------------------------------------
rPIG <- function(n, mu = 0.5, sigma = 0.02)
  { 
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
          if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
          n <- ceiling(n)
          p <- runif(n)
          r <- qPIG(p, mu=mu, sigma=sigma)
          r
  }
##-----------------------------------------------------------------------------------------
