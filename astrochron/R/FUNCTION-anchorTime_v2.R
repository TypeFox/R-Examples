### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### anchorTime function - (SRM: November 23, 2012; June 10, 2013)
###
### This script will anchor a floating time series to a radioisotopic age (or otherwise)
###  it centers a time series on a given (floating) time datum and 
###  assigns the anchored age.
###########################################################################

anchorTime <- function (dat,time,age,timeDir=1,flipOut=F,verbose=T,genplot=T)
{
   dat <- data.frame(dat)
   npts <- length(dat[,1]) 
   if(verbose)   { cat("\n----- ANCHOR FLOATING ASTROCHRONOLOGY-----\n") }   
   if(verbose)   { cat(" * Number of data points=", npts,"\n") }

### sort to ensure increasing depth/height/time
   if(verbose)   { cat (" * Sorting data and removing empty entries\n") }
   dat2 <- dat[order(dat[1],na.last=NA,decreasing=flipOut),]
   if (timeDir==1) dat2[,1] = time-dat2[,1] + age
   if (timeDir==2) dat2[,1] = dat2[,1]-time + age

   if(genplot)
    {
      par(mfrow=c(2,1))
      plot(dat, cex=.5,xlab="Location",ylab="Value",main="Floating Series",bty="n")
      lines(dat)
      plot(dat2, cex=.5,xlab="Location",ylab="Value",main="Anchored Series",bty="n")
      lines(dat2)
    }
   
     return(data.frame(dat2))

### END function anchorTime
}
