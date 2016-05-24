### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### flip function - (SRM: June 10, 2013)
###
### flip stratigraphic order, and assign new start datum
###########################################################################

flip <- function (dat,begin=0,genplot=T,verbose=T)
{

   dat <- data.frame(dat)
   npts <- length(dat[,1]) 
# determine the current direction of the data series   
   dt <- dat[2,1]-dat[1,1]
   if(verbose)   { cat("\n----- FLIP STRATIGRAPHIC SERIES-----\n") }   
   if(verbose)   { cat(" * Number of data points=", npts,"\n") }
   if(dt > 0) flipOut = TRUE
   if(dt < 0) flipOut = FALSE
# flip/sort
   if(verbose)   { cat (" * Flipping/sorting data, removing empty entries\n") }
   dat2 <- dat[order(dat[1],na.last=NA,decreasing=flipOut),]
# dx is the sampling interval of the flipped series
   dx = dat2[2:npts,1]-dat2[1:npts-1,1]
   dat2[1,1] = begin
   if(dt > 0) dat2[2:npts,1] = begin - cumsum(dx) 
   if(dt < 0) dat2[2:npts,1] = begin + cumsum(dx)    

# reassign R indices
   dat2 <- data.frame( cbind(dat2[,1],dat2[,2]) )

### plot data series. Note, cex is the factor by which to increase or decrease default symbol size
   if(genplot)
    {
     par(mfrow=c(2,1))
     plot(dat, cex=.5,xlab="Location",ylab="Value",main="Original Data Series",bty="n")
     lines(dat)
     plot(dat2, cex=.5,xlab="Location",ylab="Value",main="Flipped Data Series",bty="n")
     lines(dat2)
    }
    
     return(data.frame(dat2))

### END function flip
}
