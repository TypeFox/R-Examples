### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### sedRamp : apply sedrate model to covert time to 
###  stratigraphy - (SRM: August 2012; June 11, 2013; June 21, 2013
###                       June 30, 2013)
###
###########################################################################

sedRamp<- function (dat,srstart=0.01,srend=0.05,genplot=T,verbose=T)
{
  
   if(verbose) cat("\n----- APPLY RAMPING SEDIMENTATION RATE MODEL TO CONVERT TIME TO STRATIGRAPHY -----\n")
   dat <-  data.frame(dat)
   npts <- length(dat[,1])
   if(verbose) cat(" * Number of points in stratigraphic series=", npts,"\n") 
   dtt <- dat[2:npts,1]-dat[1:(npts-1),1]
# error checking   
   if((min(dtt)-max(dtt)) != 0) 
     {
       cat("\n**** ERROR: sampling interval is not uniform.\n")
       stop("**** TERMINATING NOW!")
     }
     
   dt=dtt[1]  
   if(verbose) cat(" * Sample interval=", dt,"\n") 
   
### sedrate step increase
    srstep=(srend-srstart)/(npts-1)
### compute sedimentation rates
    sedrate <- srstart + srstep*(  (1:npts)-1 )
### convert time to space
    depth <- cumsum(dt*sedrate)
    depth <- depth - min(depth)
### assign to data frame
    out <- as.data.frame(cbind(depth,dat[,2]))

### plot noise model
    if(genplot)
      {
        par(mfrow=c(2,1))
        plot(dat, cex=.5,xlab="Time",ylab="Value",main="Time Series",bty="n")
        lines(dat)
        plot(out, cex=.5,xlab="Space",ylab="Value",main="Modeled Stratigraphic Series",bty="n")
        lines(out)
       }
        
    return(out)
    
### END function sedRamp
}

