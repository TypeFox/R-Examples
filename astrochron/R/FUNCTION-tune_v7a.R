### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### Read tune - (SRM: April 26, 2012; October 29, 2012; May 20, 2013; 
###                   June 5, 2013; June 13, 2013; July 1, 2013; Nov. 22, 2013;
###                   January 14, 2014; June 27, 2014)
###
### Tune spatial series to time, using time control points
###########################################################################

tune <- function (dat,controlPts,extrapolate=F,genplot=T,verbose=T)
{

if(verbose) cat("\n----- TUNING STRATIGRAPHIC SERIES -----\n")

# make sure the input series are data frames
   dat=data.frame(dat)
   controlPts=data.frame(controlPts)

   npts <- length(dat[,1]) 
   if(verbose) cat(" * Number of data points=", npts,"\n")

   ictrl <- length(controlPts[,1])
   if(verbose) cat(" * Number of time control points=", ictrl,"\n")
   
### sort to ensure increasing depth/height/time
   if(verbose) cat(" * Sorting datasets into ensure increasing order, removing empty entries\n")
   dat <- dat[order(dat[1],na.last=NA,decreasing=F),]
   controlPts <- controlPts[order(controlPts[1],na.last=NA,decreasing=F),]
   
### check for duplicate depths/heights in dat
   dx1=dat[2:npts,1]-dat[1:(npts-1),1]
   if(min(dx1) == 0)
     {
       cat("\n**** WARNING: duplicate depth/height datum found in dat\n")
     }  

### check for duplicate depths/heights in controlPts
   dx2=controlPts[2:ictrl,1]-controlPts[1:(ictrl-1),1]
   if(min(dx2) == 0)
     {
       cat("\n**** ERROR: duplicate depth/height datum found in controlPts\n")
       stop("**** TERMINATING NOW!")
     }  


tuneit <- function (npts,x,ictrl,ctrl,t)
 {
    F_dat = .Fortran( 'tune_r',PACKAGE='astrochron',
    npts=as.integer(npts),x=as.double(x),ictrl=as.integer(ictrl),
    ctrl=as.double(ctrl),t=as.double(t),
    
    tuned=double(npts) 
    )

# return the results
    return(F_dat)
 }
 
   tuneout <- tuneit(npts,dat[,1],ictrl,controlPts[,1],controlPts[,2])

   out <- data.frame (cbind (tuneout$tuned,dat[,2]) )

# now remove points from both start and end if needed
   dtDir=controlPts[1,2]-controlPts[ictrl,2]
   if(!extrapolate && dtDir <0) out = subset(out, (out[1] >= controlPts[1,2]) & (out[1] <= controlPts[ictrl,2]) ) 
   if(!extrapolate && dtDir >0) out = subset(out, (out[1] <= controlPts[1,2]) & (out[1] >= controlPts[ictrl,2]) ) 

   ipts=length(out[,1])

### now evaluate sampling statistics
     t1<-out[1:(ipts-1),1]
     t2<-out[2:ipts,1]
     dt=t2-t1
     dtMin=min(dt) 
     dtMax=max(dt)
     dtMean=mean(dt)     
     dtMedian=median(dt)

     if(verbose)
      {
       cat("\n * Mean sampling interval=", dtMean,"\n")
       cat(" * Median sampling interval=",dtMedian,"\n")
       cat(" * Maximum sampling interval=",dtMax,"\n")
       cat(" * Minimum sampling interval=", dtMin,"\n") 
      }
     
     if(genplot)
      {
### plot data series. Note, cex is the factor by which to increase or decrease default symbol size
       par(mfrow=c(2,1))
       plot(dat, xlab="Location",ylab="Value",main="Data Series",bty="n",lwd=2,cex.axis=1.3,cex.lab=1.3,cex.main=1.4,cex=0.5)
       lines(dat)
       plot(out, xlab="Tuned",ylab="Value",main="Tuned Data Series",bty="n",lwd=2,cex.axis=1.3,cex.lab=1.3,cex.main=1.4,cex=0.5)
       lines(out)
     }
     
     return(out)

### END function tune
}
