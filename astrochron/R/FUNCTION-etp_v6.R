### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### etp: generate combined eccentricity-tilt-precession curve 
###          (SRM: March 23, 2012; May 3, 2012; October 11, 2012; April 23, 2013; April 29, 2013
###                May 20, 2013; July 10, 2013; August 15, 2013; April 7, 2014; January 21, 2015)
###
###########################################################################


etp <- function (tmin=0,tmax=1000,dt=1,eWt=1,oWt=1,pWt=1,esinw=T,solution=NULL,standardize=T,genplot=T,verbose=T)
{

  if(verbose) cat("\n----- GENERATING ECCENTRICITY-TILT-PRECESSION MODEL -----\n")
  if(is.null(solution))
   {
     if(tmax<=249000)
      {
        if(verbose) cat(" * Downloading Laskar et al. (2004) astronomical solution\n")
        tempLaskar <- tempfile()
        download.file("http://www.geology.wisc.edu/~smeyers/astrochron/la04.txt.bz2",tempLaskar)
        if(verbose) cat(" * Decompressing Laskar et al. (2004) astronomical solution\n")
        la <- read.table(bzfile(tempLaskar),header=T)
        unlink(tempLaskar)
        la=data.frame(cbind((1:249001)-1,la))
      }
     if(tmax>249000) 
      {
         cat("\n**** ERROR: Laskar et al. (2004) solution is not available > 249000 ka.\n")
         stop("    TERMINATING NOW!")
      }
   }  

  if(!is.null(solution)) la=solution
    
# replace preccession angle with sin(precession angle)
  la[2] <- sin(la[2])
# replace precession angle with esin(precession angle)
  if(esinw) la[2] <- la[4]*la[2]
# isolate portion of record that is needed
  la <- subset(la, (la[1] >= (tmin-1)) & (la[1] <= (tmax+1)) )   
# calculate mean values for eccentricity, obliquity, precession
  laMean <- colMeans(la)
# calculate std deviations for eccentricity, obliquity, precession
  laStdev <- sapply(la,sd)
# standardize eccentricity, tilt, precession
  if(standardize)
    {
     prec <- ( la[2]-laMean[2] ) / laStdev[2]
     obl <- ( la[3]-laMean[3] ) / laStdev[3]
     ecc <- ( la[4]-laMean[4] ) / laStdev[4]
    }
  if(!standardize)
    {
     prec <- la[2]
     obl <- la[3]
     ecc <- la[4]
    }

# now weight each term as desired
  prec <- prec * pWt
  obl <- obl * oWt
  ecc <- ecc * eWt
  
# combine eccentricity, tilt and precession
  combined <- ecc + obl + prec
  out <- data.frame( cbind(la[1],combined) )

### then interpolate as needed
  xout <- seq(tmin,tmax,by=dt)
### redefine npts
  npts <- length(xout)
  interp <- approx(out[,1],out[,2],xout,method="linear",n=npts)
  interp <- data.frame(interp)

if (genplot)
  {
### plots
   par(mfrow=c(2,2))
   plot(la[,1],la[,4],cex=0.5,xlab="Time (ka BP)",ylab="Eccentricity",main="Eccentricity",type="l")
   plot(la[,1],la[,3],cex=0.5,xlab="Time (ka BP)",ylab="Obliquity (radians)",main="Obliquity",type="l")
   if(esinw) plot(la[,1],la[,2],cex=0.5,xlab="Time (ka BP)",ylab="Eccentricity*sin(angle)",main="Precession",type="l")
   if(!esinw) plot(la[,1],la[,2],cex=0.5,xlab="Time (ka BP)",ylab="sin(angle)",main="Precession",type="l")
   plot(interp,cex=0.25,xlab="Time (ka BP)",ylab="Value",main="ETP");lines(interp,col="red")
  }
  
  return(interp)

### END function etp
}
