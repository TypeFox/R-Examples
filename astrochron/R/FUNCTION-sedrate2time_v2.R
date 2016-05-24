### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### sedrate2time: integrate sedimentation rate curve to obtain 
###               time-space map - (SRM: June 12, 2013; February 20, 2014)
###
###########################################################################

sedrate2time <- function (sedrates,timedir=1,genplot=T,verbose=T)
{

  if (verbose) cat("\n----- INTEGRATING SEDIMENTATION RATE CURVE-----\n")

  sedrates<-data.frame(sedrates)

# sort to ensure increasing depth/height
  if (verbose) 
    {
       cat("\n * Sorting sedrates into increasing depth/height order.\n")
       cat("     Will remove empty entries.\n")
    }   
  sedrates <- sedrates[order(sedrates[1],na.last=NA,decreasing=F),]
  ipts <- length(sedrates[,1])
  if (verbose) cat(" * Number of sedimentation rates=", ipts,"\n")

# check for duplicate values at a given depth
   dx1=sedrates[2:ipts,1]-sedrates[1:(ipts-1),1]
   if(min(dx1) == 0)
     {
       cat("\n**** ERROR: duplicate depth/height datum found\n")
       stop("**** TERMINATING NOW!")
     }  

# interpolate to even sampling interval
  dxmin = min(sedrates[2:ipts,1] - sedrates[1:(ipts-1),1])
  sedrates <- linterp(sedrates,dt=dxmin,genplot=F,verbose=F)
# new number of points
  npts <- length(sedrates[,1])
  
# integrate using mid-point algorithm
# convert depth/height to cm for integration purposes
  sedrates[1] = sedrates[1]*100
# convert sedrate to ka/cm for integration purposes
  sedrates[2] = 1/sedrates[2]
# calculate midpoints
  dx = sedrates[2,1]-sedrates[1,1]
  midptx = ( sedrates[2:npts,1] + sedrates[1:(npts-1),1] ) / 2
  slope = ( sedrates[2:npts,2] - sedrates[1:(npts-1),2] ) /dx
  yint = sedrates[2:npts,2]-(slope*sedrates[2:npts,1])
  midpty = (slope*midptx) + yint
# now add up midpts to get integration estimate
  hsum = cumsum(midpty*dx)
  hsum = append(0,hsum)
 
  if(timedir == 2) hsum = hsum[npts]-hsum
  
  out = data.frame(cbind(sedrates[,1]/100,hsum))
  colnames(out)[1] <- 'meters'
  colnames(out)[2] <- 'ka'

if(genplot)
  {
### plots
    par(mfrow=c(1,1))
    plot(out,cex=0.8,xlab="meters",ylab="Time (ka)",main="Time-Space Map"); lines(out)
  }
  
  return(out)

### END function sedrate2time
}
