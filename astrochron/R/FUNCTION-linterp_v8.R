### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### linterp: linearly interpolate data - (SRM: January 24, 2012; Feb. 10 2012; 
###                                            July 31, 2012; Dec. 14, 2012;
###                                            May 20, 2013; June 5-7, 2013;
###                                            June 19-21, 2013; July 29, 2013;
###                                            August 9, 2013; October 16-17, 2015)
###########################################################################

linterp <- function (dat,dt=NULL,start=NULL,genplot=T,verbose=T)
{

  if(verbose) cat("\n----- APPLYING PIECEWISE-LINEAR INTERPOLATION TO STRATIGRAPHIC SERIES -----\n")
  dat <- data.frame(dat)
  dt2 <- dat[2,1]-dat[1,1]
      
  if(dt2<0)
     { 
       if (verbose) cat("\n * Sorting data into increasing height/depth/time, removing empty entries\n")
       dat <- dat[order(dat[1], na.last = NA, decreasing = F), ]
       npts <- length(dat[,1])
     }
     
  npts <- length(dat[,1]) 
  if(verbose) cat("\n * Number of samples=", npts,"\n")

### first define points for interpolation
  if(is.null(start)) start <- dat[1,1]
  end <- dat[length(dat[,1]),1]

# interpolated to median dx by default
if(is.null(dt))
  {
     if(verbose) cat(" * Determining median sampling interval for series\n")
 ### now evaluate sampling statistics
     x1<-dat[1:(npts-1),1]
     x2<-dat[2:(npts),1]
     dx=x2-x1    
     dt=median(dx)
     if(verbose) cat(" * Will interpolate to median sampling interval of",dt,"\n")
   } 

# error checking 
   if(dt==0) 
    {
      if (verbose) cat("\n**** ERROR: dt is ZERO!\n")
      stop("**** TERMINATING NOW!")
    }  

   if(dt<0) 
    {
      if (verbose) cat("\n**** WARNING: dt must be positive. Will take absolute value.\n")
      dt=abs(dt)
    }  

 
### then interpolate
  xout <- seq(start,end,by=dt)
### redefine npts
  npts <- length(xout)
  interp <- approx(dat[,1],dat[,2],xout,method="linear",n=npts)

### assign interpolated data to data.frame d
  d <- as.data.frame(interp)
  nnpts=length(d[,1])
  if(verbose) cat(" * New number of samples=",nnpts,"\n")
  colnames(d)[1] <- colnames(dat[1])
  colnames(d)[2] <- colnames(dat[2])

### plot data
  if(genplot)
   {
    par(mfrow=c(2,2))
    plot(d,cex=0.5, xlab="Location",ylab=colnames(dat[2]),main="Raw (black) and Interpolated (red) Data",col="red"); lines(dat)
### plot the denisty and the histogram together
    hist(interp$y,freq=F,xlab="Interpolated Value",main="Distribution of Interpolated Values"); lines(density(interp$y, bw="nrd"),col="red");grid()
### boxplot
    boxplot(interp$y,ylab="Interpolated Value",main="Boxplot of Interpolated Values")
### Normal probabilty plot (Normal Q-Q Plot)
    qqnorm(interp$y); qqline(interp$y, col="red"); grid()
   } 
    
  return(d)

### END function linterp
}
