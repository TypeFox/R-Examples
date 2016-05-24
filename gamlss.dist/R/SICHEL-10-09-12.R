#as a function of mean=ksi and sigma=tau=1/omega
# Modified by Mikis Monday, March 21, 2005 
# Modified by Bob 1/12/2004
# Wednesday, October 27, 2004 at 14:34
# Modified by Kalliope Akantziliotou 
#Uses the Dean and Lawless iterative algorithm and for the derivative of nu
#uses numerical derivative 
#----------------------------------------------------------------------------------------
SICHEL <-function (mu.link ="log", sigma.link="log", nu.link="identity") 
{
     mstats <- checklink("mu.link", "Sichel", substitute(mu.link), 
                         c("1/mu^2", "log", "identity"))
    dstats <- checklink("sigma.link", "Sichel", substitute(sigma.link), 
                         c("inverse", "log", "identity"))
    vstats <- checklink("nu.link", "Sichel",substitute(nu.link), 
                         c("1/nu^2", "log", "identity"))    
    structure(
          list(family = c("SICHEL", "Sichel"),
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
                      sigma <- ifelse(sigma>10000&nu>0, 10000, sigma )
                         ty <- tofySICHEL(y=y, mu=mu, sigma=sigma, nu=nu)
                       dldm <- (y-ty)/mu
                       dldm}, 
               d2ldm2 = function(y,mu,sigma,nu) 
                           {
                     sigma <- ifelse(sigma>10000&nu>0, 10000, sigma )
                        ty <- tofySICHEL(y=y, mu=mu, sigma=sigma,nu=nu)
                      dldm <- (y-ty)/mu
                    d2ldm2 <- - dldm * dldm
                    d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)  
                    d2ldm2
                           },
                 dldd = function(y,mu,sigma,nu) 
                           {
                             sigma <- ifelse(sigma>10000&nu>0, 10000, sigma )
                        ty <- tofySICHEL(y=y, mu=mu, sigma=sigma,nu=nu)
                         c <- exp(log(besselK((1/sigma),nu+1))-log(besselK((1/sigma),nu)))
                      dcdd <- (c*sigma*(2*nu+1)+1-c*c)/(sigma*sigma)   
                      dldd <- (((ty*(c+sigma*mu)/mu) - (sigma*y)-c)/(sigma^2))+dcdd*(ty-y)/c
                      dldd 
                            },
               d2ldd2 = function(y,mu,sigma,nu)
                            {
                              sigma <- ifelse(sigma>10000&nu>0, 10000, sigma )
                        ty <- tofySICHEL(y=y, mu=mu, sigma=sigma,nu=nu)
                         c <- exp(log(besselK((1/sigma),nu+1))-log(besselK((1/sigma),nu)))
                      dcdd <- (c*sigma*(2*nu+1)+1-c*c)/(sigma*sigma)   
                      dldd <- (((ty*(c+sigma*mu)/mu) - (sigma*y)-c)/(sigma^2))+dcdd*(ty-y)/c
                    d2ldd2 <- -dldd*dldd
                    d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
                    d2ldd2
                             },
              #d2ldmdd = function() 0,
              d2ldmdd = function(y,mu,sigma,nu) 
                             {
                               sigma <- ifelse(sigma>10000&nu>0, 10000, sigma )
                         ty <- tofySICHEL(y=y, mu=mu, sigma=sigma, nu=nu)
                       dldm <- (y-ty)/mu
                          c <- exp(log(besselK((1/sigma),nu+1))-log(besselK((1/sigma),nu)))
                       dcdd <- (c*sigma*(2*nu+1)+1-c*c)/(sigma*sigma)   
                       dldd <- (((ty*(c+sigma*mu)/mu) - (sigma*y)-c)/(sigma^2))+dcdd*(ty-y)/c
                    d2ldmdd <- -dldm *dldd
                    d2ldmdd
                              }, 
              d2ldmdv = function(y,mu,sigma,nu) {
                sigma <- ifelse(sigma>10000&nu>0, 10000, sigma )
                   ty <- tofySICHEL(y=y, mu=mu, sigma=sigma, nu=nu)
                 dldm <- (y-ty)/mu
                 #calculates the dldv                 
                   nd <- numeric.deriv(dSICHEL(y, mu, sigma, nu, log=TRUE), "nu", delta=0.001)
                 dldv <- as.vector(attr(nd, "gradient"))
                #calculates the d2ldmdv
              d2ldmdv <- -dldm *dldv
              d2ldmdv
                                    }, 
              d2ldddv = function(y,mu,sigma,nu) {
                   ty <- tofySICHEL(y=y, mu=mu, sigma=sigma,nu=nu)
                    c <- exp(log(besselK((1/sigma),nu+1))-log(besselK((1/sigma),nu)))
                 dcdd <- (c*sigma*(2*nu+1)+1-c*c)/(sigma*sigma)   
                 dldd <- (((ty*(c+sigma*mu)/mu) - (sigma*y)-c)/(sigma^2))+dcdd*(ty-y)/c
                #calculates the dldv
                   nd <- numeric.deriv(dSICHEL(y, mu, sigma, nu, log=TRUE), "nu", delta=0.001)
                 dldv <- as.vector(attr(nd, "gradient"))
                #calculates the d2ldddv 
              d2ldddv <- -dldd *dldv
              d2ldddv
                                 },               
                 dldv = function(y,mu,sigma,nu) {                           
                   nd <- numeric.deriv(dSICHEL(y, mu, sigma, nu, log=TRUE), "nu", delta=0.001)
                 dldv <- as.vector(attr(nd, "gradient"))
                 dldv
                                    },
               d2ldv2 = function(y,mu,sigma,nu){
                #  delta  <- 0.01
                #    dldv <- (dSI(y=y, mu=mu, sigma=sigma, nu=nu+delta, log = TRUE)- 
                #            dSI(y=y, mu=mu, sigma=sigma, nu=nu, log= TRUE))/ delta   
                   nd <- numeric.deriv(dSICHEL(y, mu, sigma, nu, log=TRUE), "nu", delta=0.001)
                 dldv <- as.vector(attr(nd, "gradient"))
               d2ldv2 <- -dldv*dldv
               d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)  
               d2ldv2
                                  } ,
          G.dev.incr  = function(y,mu,sigma,nu, pw=1,..) -2*dSICHEL(y, mu, sigma, nu, log=TRUE),
                rqres = expression(    
                        rqres(pfun="pSICHEL", type="Discrete", ymin=0, y=y, mu=mu, sigma=sigma, nu=nu)
                                  ), 
            mu.initial = expression(mu<- (y+mean(y))/2 ),
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
# this is only for checking with tofySICHEL
#tofySNN <- function (y, mu, sigma, nu, bsum=FALSE, ...)
#{
#    ty <- rep(0,length(y))
#    sumlty <- rep(0,length(y)) 
#    for (z1 in length(y):1)
#    {
#    #browser()
#         lyp1 <- Dum <- swi <- dum1 <- alpha <- lbes <- c <- 0
#         lyp1 <- y[z1]+1 
#        tynew <- rep(0,lyp1)
#            c <- exp(log(besselK(1/sigma[z1],nu[z1]+1))-log(besselK((1/sigma[z1]),nu[z1])))
#        alpha <- sqrt(1+2*sigma[z1]*(mu[z1]/c))/sigma[z1]
#        #cat("c, alphs",c,alpha,"\n")
#         lbes <- log(besselK(alpha,nu[z1]+1))-log(besselK((alpha),nu[z1]))
#     tynew[1] <- dum1 <- (mu[z1]/c)*((1+2*sigma[z1]*(mu[z1]/c))^(-0.5))*exp(lbes) 
#          dum <- ifelse(lyp1==1, 1,2)
#            j <- dum:lyp1 
#        for (i in j)
#        {
#            if (i !=1)
#           tynew[i] <- (c*sigma[z1]*(2*(i-1+nu[z1])/mu[z1])+(1/tynew[i-1]))*(mu[z1]/(sigma[z1]*alpha*c))**2
#            tynew[1] <- dum1 
#           ty[z1] <- tynew[lyp1] 
#        }
#        if (bsum ) 
#            sumlty[z1] <- ifelse(lyp1==1, 0 , sum(log(tynew))-log(tynew[i]))
#    }
#    result <- cbind(ty, sumlty)
#}
#----------------------------------------------------------------------------------------
tofySICHEL <- function (y, mu, sigma, nu)
   {
  ly <- max(length(y),length(mu),length(sigma),length(nu)) 
  y <- rep(y, length = ly)    
  sigma <- rep(sigma, length = ly)
  mu <- rep(mu, length = ly)   
  nu <- rep(nu, length = ly) 
  cvec <- exp(log(besselK((1/sigma),nu+1))-log(besselK((1/sigma),nu)))
  #cvec is  c=R_nu(1/sigma) where R_lambda(t)=K_{lambda+1}(t)/K_{lambda}(t) page 10
  alpha <- sqrt(1+2*sigma*(mu/cvec))/sigma # alpha^2=sigma^-2+(2*mu)/(c*sigma)
  lbes <-  log(besselK(alpha,nu+1))-log(besselK((alpha),nu)) #R_nu(alpha) page 30
  sumlty <- as.double(.C("tofySICHEL1", as.double(y), as.double(mu), 
               as.double(sigma), as.double(nu), as.double(lbes),
               as.double(cvec), ans=double(length(y)),as.integer(length(y)),
               as.integer(max(y)+1), PACKAGE="gamlss.dist")$ans)
   sumlty
   }
#----------------------------------------------------------------------------------------
# this function is using tofySNN
#dSInew<-function(y, mu=1, sigma=1, nu=-0.5, log=FALSE)
#  { 
#   if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
#   if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
#   if (any(y < 0) )  stop(paste("y must be >=0", "\n", ""))  
#    ly <- length(y)       
# sigma <- rep(sigma, length = ly)
#    mu <- rep(mu, length = ly)   
#    nu <- rep(nu, length = ly) 
#  cvec <- exp(log(besselK((1/sigma),nu+1))-log(besselK((1/sigma),nu)))
# alpha <- sqrt(1+2*sigma*(mu/cvec))/sigma
## lbes <-  log(besselK(alpha,nu+1))-log(besselK((alpha),nu))
#sumlty <- as.double(.C("tofys_", as.double(y), as.double(mu), 
#                   as.double(sigma), as.double(nu), as.double(lbes),
#                   as.integer(length(y)), as.integer(max(y)+1))[[2]])
#sumlty<-tofySNN(y=y, mu=mu, sigma=sigma, nu=nu, bsum=TRUE)[,2]
#browser()
#logfy <- -lgamma(y+1)-nu*log(sigma*alpha)+sumlty+log(besselK(alpha,nu))-log(besselK((1/sigma),nu))
#  if(log==FALSE) fy <- exp(logfy) else fy <- logfy
#  fy
#  }
#------------------------------------------------------------------------------------------
dSICHEL<-function(x, mu=1, sigma=1, nu=-0.5, log=FALSE)
  { 
   if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
   if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
   if (any(x < 0) )  stop(paste("x must be >=0", "\n", ""))  
    ly <- max(length(x),length(mu),length(sigma),length(nu)) 
     x <- rep(x, length = ly)      
 sigma <- rep(sigma, length = ly)
    mu <- rep(mu, length = ly)   
    nu <- rep(nu, length = ly) 
  cvec <- exp(log(besselK((1/sigma),nu+1))-log(besselK((1/sigma),nu)))
 alpha <- sqrt(1+2*sigma*(mu/cvec))/sigma
  lbes <-  log(besselK(alpha,nu+1))-log(besselK((alpha),nu))
 #cat("mu, sigma and nu", mu[1], sigma[1], nu[1], "\n")
sumlty <- as.double(.C("tofySICHEL2", as.double(x), as.double(mu), 
                       as.double(sigma), as.double(nu), as.double(lbes),
                       as.double(cvec), ans=double(length(x)),as.integer(length(x)),
                       as.integer(max(x)+1), PACKAGE="gamlss.dist")$ans)
logfy <- -lgamma(x+1)-nu*log(sigma*alpha)+sumlty+log(besselK(alpha,nu))-log(besselK((1/sigma),nu))
  
  if(log==FALSE) fy <- exp(logfy) else fy <- logfy
  if (length(sigma)>1) fy <- ifelse((sigma>10000)&(nu>0), dNBI(x, mu = mu, sigma= 1/nu, log = log) ,fy)
        else fy <- if ((sigma>10000)&(nu>0)) dNBI(x, mu = mu, sigma= 1/nu, log = log)  
                   else  fy
  fy
  }
#----------------------------------------------------------------------------------------     
#--------------------------------------------------
#  tocdfS <- function (y, mu, sigma, nu, bsum=TRUE, ...)
#  {
#      ly <- length(y)       
#   sigma <- rep(sigma, length = ly)
#      mu <- rep(mu, length = ly)   
#      nu <- rep(nu, length = ly) 
#      ty <- rep(0,length(y))
#     cdf <- rep(0,length(y))
#    cvec <- exp(log(besselK((1/sigma),nu+1))-log(besselK((1/sigma),nu)))
#   alpha <- sqrt(1+2*sigma*mu/cvec)/sigma
#    lbes <-  log(besselK(alpha,nu+1))-log(besselK((alpha),nu)) 
#    for (i in 1:length(y))
#    {
#         lyp1 <- y[i]+1 
#        tynew <- lpnew <- rep(0,lyp1)
#        tynew[1] <- (mu[i]/cvec[i])*((1+2*sigma[i]*(mu[i]/cvec[i]))^(-0.5))*exp(lbes[i])
#        lpnew[1] <- -nu[i]*log(sigma[i]*alpha[i])+log(besselK(alpha[i],nu[i]))-
#                     log(besselK(1/sigma[i],nu[i])) 
#        dum <- ifelse(lyp1==1, 1,2)
#        for (j in dum:lyp1)
#        {
#            if (j !=1)
#             {
#            tynew[j] <- (cvec[i]*sigma[i]*(2*(j-1+nu[i])/mu[i])+(1/tynew[j-1]))*
#                         (mu[i]/(sigma[i]*alpha[i]*cvec[i]))**2
#            lpnew[j] <- lpnew[j-1] + log(tynew[j-1]) - log(j-1)
#             }           
#            ty[i] <- tynew[lyp1] 
#        }
#        if (bsum ) 
#            cdf[i] <- sum(exp(lpnew))
#    }
#    cdf
#  }
#-----------------------------------------------
pSICHEL <- function(q, mu=1, sigma=1, nu=-0.5, lower.tail = TRUE, log.p = FALSE)
{  
  if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
  if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
  if (any(q < 0) )  stop(paste("q must be >=0", "\n", ""))  
  ly <- max(length(q),length(mu),length(sigma),length(nu)) 
  q <- rep(q,length = ly)       
  sigma <- rep(sigma, length = ly)
  mu <- rep(mu, length = ly)   
  nu <- rep(nu, length = ly)   
  cdf <- as.double(.C("cdfSICHEL", as.double(q), as.double(mu),as.double(sigma), as.double(nu), 
            ans=double(ly), as.integer(ly), PACKAGE="gamlss.dist")$ans)
  # as.integer(max(q)+1))#
  if(lower.tail==TRUE) cdf <- cdf else cdf=1-cdf
  if(log.p==FALSE) cdf <- cdf else cdf <- log(cdf)                                                                    
  cdf
}
#----------------------------------------------------------------------------------------
qSICHEL <- function(p, mu=1, sigma=1, nu=-0.5,  lower.tail = TRUE, log.p = FALSE,  
                 max.value = 10000)
  {      
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
          if (any(p < 0) | any(p > 1.0001))  stop(paste("p must be between 0 and 1", "\n", "")) 
          if (log.p==TRUE) p <- exp(p) else p <- p
          if (lower.tail==TRUE) p <- p else p <- 1-p  
           ly <- max(length(p),length(mu),length(sigma),length(nu)) 
            p <- rep(p, length = ly)                                                         
          QQQ <- rep(0,length = ly)                         
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
            cumpro <-  pSICHEL(j, mu = nmu[i], sigma = nsigma[i], nu = nnu[i], log.p = FALSE) 
                       # else  cumpro+dSICHEL(j, mu = nmu[i], sigma = nsigma[i], nu = nnu[i], log = FALSE)# the above is faster 
           QQQ[i] <- j 
       if  (p[i] <= cumpro ) break 
            } 
        }
      }          
#        for (i in seq(along=p))                                                          
#      {                                                                 
#     cumpro <- pSI(0, mu = nmu[i], sigma = nsigma[i],  nu = nnu, log.p = FALSE)
#      for (j in seq(from = 1,to = max.value))
#        {
#        if (p[i]+0.000000001 >= 1) QQQ[i] <- Inf
#        else if (p[i]-cumpro > 0 | identical(p[i],cumpro) )#((cumpro+.Machine$double.eps) < p.p) 
#           { cumpro <-  if (fast == FALSE) pSI(j, mu = nmu[i], sigma = nsigma[i], nu = nnu[i], log.p = FALSE) 
#                        else  cumpro+dSI(j, mu = nmu[i], sigma = nsigma[i], nu = nnu[i], log.p = FALSE)# the above is faster 
#             QQQ[i] <- j-1 #if (! identical(all.equal(p[i],1),TRUE)) j-1 else Inf
#           }
#        else break
#        }                                  
#       }   
          QQQ   
   }
#----------------------------------------------------------------------------------------
rSICHEL <- function(n, mu=1, sigma=1, nu=-0.5, max.value = 10000)
  { 
          if (any(mu <= 0) )  stop(paste("mu must be greater than 0 ", "\n", "")) 
          if (any(sigma <= 0) )  stop(paste("sigma must be greater than 0 ", "\n", "")) 
          if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
          n <- ceiling(n)
          p <- runif(n)
          r <- qSICHEL(p, mu=mu, sigma=sigma, nu=nu, max.value = max.value )
          r
  }
#----------------------------------------------------------------------------------------
VSICHEL<- function(obj)
 {
 #if (!is.gamlss(obj)) stop("The object is not a gamlss object")
 if (obj$family[1]!="SICHEL") stop("The distribution is not a Sichel distribution")
     mu <- fitted(obj,"mu")
  sigma <- fitted(obj,"sigma")
     nu <- fitted(obj,"nu")
     cc <- exp(log(besselK(1/sigma,nu+1))-log(besselK((1/sigma),nu)))
 varY <-  mu + mu^2 * (((2*sigma*(nu+1))/cc) + (1/cc^2) -1)
 varY
 }
#Etofy <- function(mu=1, sigma=1, nu=1, maxval=1000)
#  {
#       ly <- maxval+1 
#       yy <- 0:maxval      
#       sigma <- rep(sigma, length = ly)
#          mu <- rep(mu, length = ly)   
#          nu <- rep(nu, length = ly)
#          ty <- tofySICHEL(yy, mu=mu, sigma=sigma, nu=nu)
#          Py <- dSICHEL(yy, mu=mu, sigma=sigma, nu=nu)
#          Ev <-sum(ty*ty*Py)
#          Ev  
#  }
#----------------------------------------------------------------------------------------
#find.corr <- function(mu, sigma,  nu )
# {
# #browser()
#      cc <- exp(log(besselK((1/sigma),nu+1))-log(besselK((1/sigma),nu)))
#   Alpha <- (cc/mu+2*sigma*(nu+1)+(1/cc)-cc)/(sigma*sigma) 
#    Ety2 <- Etofy( mu = mu, sigma = sigma, nu = nu)  
#    dcds <- (cc*sigma*(2*nu+1)+1-cc*cc)/(sigma*sigma)
#     Emu <- -(1/mu)+((((2*sigma*(nu+1))/cc)+(1/(cc*cc))))-Ety2/(mu*mu)
#  Esigma <- -(Alpha*Alpha)*mu+ (Alpha*cc)/(sigma*sigma) +
#            (Alpha*Alpha)*(mu*mu)*((2*sigma*(nu+1))/(cc)+1/(cc*cc))-
#            (Alpha*Alpha)*Ety2
#Emusigma <-  (1/sigma)+(dcds/cc)-Alpha*mu*((2*sigma*(nu+1))/(cc)+1/(cc*cc))+
#             (Alpha*Ety2)/mu 
#    MMM  <- matrix(NA,2,2)
#MMM[1,1] <- Emu
#MMM[2,2] <- Esigma
#MMM[1,2] <- MMM[2,1] <- Emusigma
#MMM <- -MMM
##print(MMM)
##browser()
#cov<- solve(MMM)
##print(cov)
#corr<- cov2cor(cov)
##MIN<- solve(MMM)
##Ecor <-MIN[1,2]/ sqrt(MIN[1,1])*sqrt(MIN[2,2])     
#corr[1,2]
# }
##----------------------------------------------------------------------------------------
## plotting the expected correlation between mu and sigma  
#corr.plot <- function(mu=10, sigma=c(0.5,1,1.5), inter=c(-3,3))
#  {
#  #browser()
#  nu.val <- seq(inter[1],inter[2], length=20)
#  len <- length(nu.val)
#  lens <-length(sigma)
#  cor <- matrix(NA, nrow=len, ncol=lens)
#  #browser()
#  co<-1
#  plot(cor[,1]~nu.val,  ty="n", ylim=c(-1,1), ylab="asymtotic correlation.", xlab=expression(nu))
#   for (j in 1:lens)
#    {
#     for (i in seq(along=nu.val)) 
#      {
#       cor[i,j] <- find.corr(mu=mu, sigma=sigma[j], nu= nu.val[i])
#      }
#    lines(cor[,j]~nu.val, col=co)
#    #browser()
#    text(x=nu.val[5], y=cor[5,j]+0.05, label=expression(sigma))
#    text(x=nu.val[5]+0.4, y=cor[5,j]+0.05, label=paste("=",sigma[j]))
#    co <-co+1
#    }
#  
#  lines(rep(0,len)~nu.val)
#  lines(nu.val~rep(0,len))
#}
##---------------------------------------------------------------------------------------
## plotting the c as a function of sigma ans nu
#c.plot <- function(sigma=c(0.001,2), nu=c(-3,3))
# {
#eee <- expand.grid(sigma=seq(sigma[1],sigma[2],length=20),nu=seq(nu[1],nu[2],length=20)) 
#ccc <- exp(log(besselK((1/eee$sigma),eee$nu+1))-log(besselK((1/eee$sigma),eee$nu)))
#eee <- data.frame(eee, c=ccc)
#sss <- seq(sigma[1],sigma[2], length=20)
#nnn <- seq(nu[1],nu[2], length=20)
#op <- par(mfrow=c(1,1))
#contour(sss,nnn,matrix(ccc,nrow=20),nlevels=50, ylab=expression(nu), xlab=expression(sigma), xlim=c(0,1)) ##bar(x) == sum(frac(x[i], n)
##image(Fln,An,matrix(newrent$pred,nrow=length(Fln)),col=topo.colors(40))
#library(lattice)
##wireframe(c~sigma*nu, eee, aspect=c(1,0.5), drape=TRUE, colorkey=list(space="right", height=0.6))
## par(op)
# }
##---------------------------------------------------------------------------------------
##dyn.load("c:/gamlss/fortran/tofySICHEL.dll")##
