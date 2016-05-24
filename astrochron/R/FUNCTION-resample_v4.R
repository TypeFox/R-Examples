### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### resample: resample a time series using a new time axis
###           with variable sample intervals. values are
###           piecewise linearly interpolated from original data. 
###               (SRM: April 4, 2013; May 20, 2013; June 18-19, 2013;
###                     February 19, 2014; April 14, 2014)
###
###########################################################################

resample <- function (dat,xout,genplot=T,verbose=T)
{

if(verbose) cat("\n----- RESAMPLING STRATIGRAPHIC SERIES -----\n")
   
  dat <- data.frame(dat) 
  datpts = length(dat[,1])
  xout<-data.frame(xout)
  npts <- length(xout[,1])
  
  if(verbose) cat(" * Number of data points in series prior to resampling:",datpts,"\n")
  if(verbose) cat(" * Number of resampling sites:",npts,"\n")

### sort to ensure increasing depth/height/time
  if(verbose) cat(" * Sorting data, removing empty entries\n")
  dat <- dat[order(dat[1],na.last=NA,decreasing=F),]
  xout <- xout[order(xout[1],na.last=NA,decreasing=F),]
# note: xout is no longer a data.frame, but dat is.
 
### cull xout if falls outside of the range of dat[,1]
  xout=xout[xout>=dat[1,1]]
  xout=xout[xout<=dat[datpts,1]]
  
  if(length(xout) < npts)
   {
     if(verbose) cat("\n**** WARNING: Some of the specified sample locations fall outside of the data set and will be ignored.\n")
   }
  
  interp <- approx(dat[,1],dat[,2],xout,method="linear")
### assign interpolated data to data.frame d
  d <- as.data.frame(interp)
### define new npts
  npts <- length(d[,1])

  if(verbose) cat(" * Number of data points following resampling=",npts,"\n")

### plot data
  if(genplot)
   {
    par(mfrow=c(2,2))
    plot(dat,type="l",xlab="Resampled Location",ylab="Resampled Value",main="Raw (black), Resampled (red)"); points(interp,cex=0.5,col="red")
### plot the denisty and the histogram together
    hist(interp$y,freq=F,xlab="Resampled Value",main="Distribution of Resampled Values"); lines(density(interp$y, bw="nrd"),col="red")
### boxplot
    boxplot(interp$y,ylab="Resampled Value",main="Boxplot of Resampled Values")
### Normal probabilty plot (Normal Q-Q Plot)
    qqnorm(interp$y); qqline(interp$y, col="red")
   } 
    
  return(d)

### END function resample
}
