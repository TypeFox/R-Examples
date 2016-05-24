### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### noKernel: remove Gaussian kernel smoother from data (SRM: May 2, 2013; May 20, 2013).
###
###########################################################################

noKernel <- function (dat,smooth=0.1,sort=F,output=1,genplot=T,verbose=T)
{

if(verbose) cat("\n----- REMOVING GAUSSIAN KERNEL SMOOTHER FROM STRATIGRAPHIC SERIES -----\n")
if(smooth <= 0) {stop("**** ERROR, FUNCTION TERMINATED: smooth must be > 0")}

# ensure we have a data frame
dat <- data.frame(dat)
npts <- length(dat[,1])
if(verbose) cat(" * Number of data points input=", npts,"\n")
# save a version for later plotting
dat1 <- data.frame(dat)

# ensure data is sorted to increasing xID order (for ksmooth)
#  missing depths (NA) are removed during sort
if (sort) 
  {
    if(verbose) cat(" * Sorting into increasing order\n")
    dat <- dat[order(dat[1],na.last=NA,decreasing=F),]
  }

# if local smoothing for detrending involved
if (smooth != 0)
 {
  if(verbose) cat(" * Smoothing/detrending with Gaussian kernel\n")
# factor of 2 necessary, because by default kernels are scaled so that their quartiles 
#  (viewed as prob densities) are at +/- 0.25*bandwidth
     smoothScaled= abs( dat[1,1]-dat[npts,1] ) * smooth * 2 
# if x.points = dat[,1], will evalute at original sample locations only
     smooth2 <- ksmooth(dat[,1],dat[,2],kernel=c("normal"),bandwidth=smoothScaled,x.points=dat[,1])
# remove smoother from data series
     resid<-data.frame(cbind(smooth2$x,(dat[,2]-smooth2$y)))
# save smoother
     smoother<-data.frame( cbind(smooth2$x,smooth2$y) )
# plot smoother
     if(genplot) 
       {
         par(mfrow=c(2,2))
         plot(dat,cex=0.5, main="Stratigraphic Data Series",xlab="Location",ylab="Value")
         lines(smoother,col="red")
         plot(resid, cex=0.5,main="Gaussian Kernel Detrended Series",xlab="Location",ylab="Detrended Value"); lines(resid)
### plot the denisty and the histogram together
         hist(resid[,2],freq=F, xlab="Value", main="Distribution of Residual Values"); lines(density(resid[,2], bw="nrd"),col="red"); grid()
### boxplot
         boxplot(resid[,2], ylab="Value", main="Boxplot of Residual Values")
        } 
 }
    
if (output==1) return(resid)
if (output==2) return(smoother)

### END function noKernel
}