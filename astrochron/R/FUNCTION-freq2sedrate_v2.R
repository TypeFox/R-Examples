### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### freq2sedrate: convert local spatial frequency to sedimentation rate - 
###                 (SRM: June 12, 2013; February 20, 2014)
###
###########################################################################

freq2sedrate <- function (freqs,period=NULL,ydir=1,genplot=T,verbose=T)
{

  if (verbose) cat("\n----- CONVERTING RECORD OF SPATIAL FREQUENCY TO SEDIMENTATION RATE CURVE-----\n")
# error checking
  if(is.null(period))
   { 
     cat("\n**** ERROR: you must assign a temporal period\n")
     stop("**** TERMINATING NOW!")
   }

  freqs<-data.frame(freqs)
  
# sort to ensure increasing depth/height
  if (verbose) 
    {
       cat("\n * Sorting into increasing depth/height order.\n")
       cat("     Will remove empty entries.\n")
    }   
  freqs <- freqs[order(freqs[1],na.last=NA,decreasing=F),]
  ipts <- length(freqs[,1])
  if (verbose) cat(" * Number of control points =", ipts,"\n")

# check for duplicate values at a given depth
   dx1=freqs[2:ipts,1]-freqs[1:(ipts-1),1]
   if(min(dx1) == 0)
     {
       cat("\n**** ERROR: duplicate depth/height datum found\n")
       stop("**** TERMINATING NOW!")
     }  

  sedrate=1/( freqs[,2]*period*0.01 )
  out = data.frame(cbind(freqs[,1],sedrate))
  colnames(out)[1] <- 'Depth/Height'
  colnames(out)[2] <- 'Sedrate'

if(genplot)
  {
### plots
    par(mfrow=c(1,2))
    if (ydir == 1) ylimset=c( min(freqs[,1]),max(freqs[,1]) )
    if (ydir == -1) ylimset=c( max(freqs[,1]),min(freqs[,1]))
    plot(freqs[,2],freqs[,1],cex=0.8,ylim=ylimset,xlab="Spatial Frequency (cycles/m)",ylab="Location (m)",main="Local Frequency"); lines(freqs[,2],freqs[,1])
    plot(out[,2],out[,1],cex=0.8,ylim=ylimset,xlab="Sedimentation Rate (cm/ka)",ylab="Location (m)",main="Sedimentation Rate"); lines(out[,2],out[,1])
  }
  
  return(out)

### END function freq2sedrate
}
