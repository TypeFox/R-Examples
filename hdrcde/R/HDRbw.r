########## R function: hdrbw ##########

# Obtain an estimate for the 'HDRlevel' highest
# density region optimal bandwidth.

# Last changed: 06 AUG 2009

hdrbw <- function(x,HDRlevel,gridsize=801,nMChdr=1000000,
                  graphProgress=FALSE)
{

    if(HDRlevel > 1 & HDRlevel < 100)
        HDRlevel <- HDRlevel/100
    if(is.list(x))
        x <- x[[1]]

   # Obtain pilot bandwidths:

   hHat0 <- densDPI(x)
   hHat1 <- densDPI(x,drv=1)
   hHat2 <- densDPI(x,drv=2)

   # Obtain pilot estimate of the density function
   # and use this to estimate the "fTau" value:

   n <- length(x)
   Ivec <- sample(1:n,nMChdr,replace=TRUE)
   xHat <- x[Ivec] + hHat0*rnorm(nMChdr)
   fhat <- KernSmooth::bkde(x,bandwidth=hHat0,gridsize=gridsize,
                range.x=c(min(xHat),max(xHat)))
   fhatx <- approx(fhat$x,fhat$y,xHat,rule=2)$y
   fTauHat <- quantile(fhatx,HDRlevel)
   names(fTauHat) <- NULL
   
   if (graphProgress)
   {
      par(mfrow=c(2,1))
      plot(fhat,type="l",bty="l",main=expression(paste("Pilot estimate of ",f[HDRlevel])),
           xlab="x",ylab="density estimate")
      lines(range(fhat$x),rep(fTauHat,2),col="red")
   }  

   # Now use density estimate to approximate the "xcuts" vector:

   discrim <- diff(sign(c(fhat$y-fTauHat)))
   discrim <- c(0,discrim)
   xcutsInds <- (1:length(fhat$x))[discrim!=0]
   xcutsHat <- (fhat$x[xcutsInds] + fhat$x[xcutsInds-1])/2

   if (graphProgress)
      points(xcutsHat,rep(fTauHat,length(xcutsHat)),col="green4") 
           
   # Compute asymptotic approximation:

   fdxcutsHat <- rep(NA,length(xcutsHat))
   fddxcutsHat <- rep(NA,length(xcutsHat))
   for (j in 1:length(xcutsHat))
   {
      fdxcutsHat[j] <- drvKDE(x,xcutsHat[j],hHat1,1)
      fddxcutsHat[j] <- drvKDE(x,xcutsHat[j],hHat2,2)
   }   
   Bhats <- hdrBvals(fTauHat,fdxcutsHat,fddxcutsHat)
   B1hat <- Bhats$B1 ; B2hat <- Bhats$B2 ; B3hat <- Bhats$B3

   # Search for minimiser of estimated asymptotic risk:

   riskFuns <- function(s,B1,B2,B3,drv=0)
   {
      farg <- B2*s^5
      if (drv==0)
         return(sum((1/s)*B1*dnorm(farg)+B3*(s^4)*(2*pnorm(farg)-1)))
      if (drv==1)
         return(sum(((-1/s^2)*B1+10*B2*B3*(s^8)-5*B1*(B2^2)*(s^8))*dnorm(farg)
                 +4*B3*(s^3)*(2*pnorm(farg)-1)))
      if (drv==2)
         return(sum((2*B1*(1/s^3)+5*B2*(24*B3-7*B1*B2)*(s^7)
                -25*(B2^3)*(2*B3-B1*B2)*(s^17))*dnorm(farg)
                 +12*B3*(s^2)*(2*pnorm(farg)-1)))
   }

   hHat <- hHat0
   sMid <- n^(1/10)*sqrt(hHat)
   sLow <- sMid/5 ; sUpp <- 5*sMid
   rootCaptured <- FALSE
   while(!rootCaptured)
   { 
      fdLow <- riskFuns(sLow,B1hat,B2hat,B3hat,1)
      fdUpp <- riskFuns(sUpp,B1hat,B2hat,B3hat,1)
      if (fdLow>0) sLow <- sLow/2
      if (fdUpp<0) sUpp <- 2*sUpp
      if ((fdLow<0)&(fdUpp>0)) rootCaptured <- TRUE
   }
   numBisect <- 100
   for (ib in 1:numBisect)
   {
      sMid <- (sLow + sUpp)/2
      fdLow <- riskFuns(sLow,B1hat,B2hat,B3hat,1)
      fdUpp <- riskFuns(sUpp,B1hat,B2hat,B3hat,1)
      fdMid <- riskFuns(sMid,B1hat,B2hat,B3hat,1)
      if (fdMid<0) sLow <- sMid
      if (fdMid>0) sUpp <- sMid
   }
   s <- sMid
   maxit <- 100 ; tolerance <- 0.0000001
   finished <- FALSE
   itnum <- 0 

   # Do some Newtown-Raphson updates:
   
   while (!finished)
   {
      itnum <- itnum + 1

      farg <- B2hat*s^5
      fds <- sum(((-1/s^2)*B1hat+10*B2hat*B3hat*(s^8)-5*B1hat*(B2hat^2)*(s^8))*dnorm(farg)
              +4*B3hat*(s^3)*(2*pnorm(farg)-1))

      fdds <- sum((2*B1hat*(1/s^3)+5*B2hat*(24*B3hat-7*B1hat*B2hat)*(s^7)
                -25*(B2hat^3)*(2*B3hat-B1hat*B2hat)*(s^17))*dnorm(farg)
                 +12*B3hat*(s^2)*(2*pnorm(farg)-1))

      sNew <- s - fds/fdds
      if (sNew<0) stop("Newton-Raphson failure (out of range).")
      relErr <- abs((sNew-s)/s)
      if (relErr<tolerance) finished <- TRUE
      if (itnum>maxit)
         stop("Newton-Raphson failure (maximum number of iterations exceeded).")
      s <- sNew 
   }

   hHat <- (s^2)*n^(-1/5)

   if (graphProgress)
   {
      hgrid <- exp(seq(log(hHat/5),log(5*hHat),length=101))
      AriskGrid <- rep(0,length(hgrid))
      for (j in 1:length(B1hat))
      {
         AriskGrid <- AriskGrid + B1hat[j]*dnorm(B2hat[j]*sqrt(n*hgrid^5))/sqrt(n*hgrid)
         AriskGrid <- AriskGrid + B3hat[j]*(hgrid^2)*(2*pnorm(B2hat[j]*sqrt(n*hgrid^5))-1)
      }
      plot(log(hgrid),AriskGrid,col="blue",type="l",bty="l",xlab="h",ylab="estimated risk")
      lines(rep(log(hHat),2),c(0,100),col="DarkOrange",lwd=2)
   }
     
   return(hHat)
}

############ End of hdrbw ###########

########## R function: hdrBvals ##########

# For computing the B_{1j}, B_{2,j} and B_{3,j}
# values that arise in the asymptotic approximation
# to the HDR (highest density region) estimation
# risk.

# Last changed: 06 AUG 2009

hdrBvals <- function(fTau,fdxcuts,fddxcuts)
{

   RK <- 1/(2*sqrt(pi))
   rVal <- length(fdxcuts)/2
   EveInds <- seq(2,(2*rVal),by=2)
   OddInds <- seq(1,(2*rVal-1),by=2)
  
   fdHarSum <- 1/sum(1/abs(fdxcuts))

   D1sum1 <- sum(fddxcuts/abs(fdxcuts))
   D1sum2 <- sum(fdxcuts[EveInds]-fdxcuts[OddInds])/fTau
   D1 <- 0.5*fdHarSum*(D1sum1 + D1sum2)

   D2 <- RK*fTau*(fdHarSum^2)*sum(1/(fdxcuts^2))

   D3 <- RK*fTau*fdHarSum/abs(fdxcuts)

   B1 <- 2*fTau*sqrt(RK*fTau-2*D3+D2)/abs(fdxcuts)
   B2 <- fTau*abs(0.5*fddxcuts-D1)/sqrt(RK*fTau-2*D3+D2)
   B3 <- fTau*abs(0.5*fddxcuts-D1)/abs(fdxcuts)

   return(list(B1=B1,B2=B2,B3=B3))
}

############ End of hdrBvals ############

########## R function: drvKDE ##########

# For obtaining an exact kernel density derivative
# estimate at a point.

# Last changed: 16 SEP 2008

drvKDE <- function(x,x0,bandwidth,drv)
{
   h <- bandwidth
   n <- length(x)

   if (drv==(-1))
      return(sum(pnorm((x0-x)/h))/n)
   
   if (drv==0)
      return(sum(dnorm((x0-x)/h))/(n*h))

   if (drv==1)
      return(sum(-((x0-x)/h)*dnorm((x0-x)/h))/(n*h^2))

   if (drv==2)
      return(sum((((x0-x)/h)^2-1)*dnorm((x0-x)/h))/(n*h^3))
}

############ End of drvKDE ############

########## R function: densDPI ##########

# For obtaining direct plug-in bandwidths
# for densities and low-order derivatives.

# Last changed: 20 SEP 2008

densDPI <- function(x,drv=0,gridsize=401,range.x=range(x),
                    truncate=TRUE)
{
   # Set preliminary values:
    
   n <- length(x)
   a <- range.x[1]
   b <- range.x[2]
   gpoints <- seq(a,b,length = gridsize)
   gcounts <- linbin(x,gpoints,truncate)
   scalest <- min(sd(x),(quantile(x,3/4)-quantile(x,1/4))/1.349)

   # Obtain estimated bandwidth via two-stage plug-in strategy
   # (depending on derivative being estimated):

   if (drv==(-1))
   {
      alpha <- (2*(sqrt(2)*scalest)^7/(9*n))^(1/7)
      psi4hat <- KernSmooth::bkfe(gcounts,4,alpha,
                                   range.x=c(a,b),binned=TRUE)
      alpha <- (sqrt(2/pi)/(psi4hat*n))^(1/5)
      psi2hat <- KernSmooth::bkfe(gcounts,2,alpha,
                                   range.x=c(a,b),binned=TRUE)
      hHat <- (-1/(sqrt(pi)*psi2hat*n))^(1/3)
   }
   
   if (drv==0)
   {
      alpha <- (2*(sqrt(2)*scalest)^9/(7*n))^(1/9)
      psi6hat <- KernSmooth::bkfe(gcounts,6,alpha,
                                   range.x=c(a,b),binned=TRUE)
      alpha <- (-3*sqrt(2/pi)/(psi6hat*n))^(1/7)
      psi4hat <- KernSmooth::bkfe(gcounts,4,alpha,
                                   range.x=c(a,b),binned=TRUE)
      hHat <- (1/(2*sqrt(pi)*psi4hat*n))^(1/5)
   }

   if (drv==1)
   {
      alpha <- (2*(sqrt(2)*scalest)^11/(9*n))^(1/11)
      psi8hat <- KernSmooth::bkfe(gcounts,8,alpha,
                                   range.x=c(a,b),binned=TRUE)
      alpha <- (15*sqrt(2/pi)/(psi8hat*n))^(1/9)
      psi6hat <- KernSmooth::bkfe(gcounts,6,alpha,
                                   range.x=c(a,b),binned=TRUE)
      hHat <- (-3/(4*sqrt(pi)*psi6hat*n))^(1/7)
   }

   if (drv==2)
   {
      alpha <- (2*(sqrt(2)*scalest)^13/(11*n))^(1/13)
      psi10hat <- KernSmooth::bkfe(gcounts,10,alpha,
                                    range.x=c(a,b),binned=TRUE)
      alpha <- (-105*sqrt(2/pi)/(psi10hat*n))^(1/11)
      psi8hat <- KernSmooth::bkfe(gcounts,8,alpha,
                                   range.x=c(a,b),binned=TRUE)
      hHat <- (15/(8*sqrt(pi)*psi8hat*n))^(1/9)
   }
   
   return(hHat)
}

############ End of densDPI ############
