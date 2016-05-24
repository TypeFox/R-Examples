### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### demean: remove mean from data - (SRM: January 26, 2012; Oct. 11, 2012;
###                                        Dec 12, 2012; May 20, 2013)
###########################################################################

demean <- function (dat,genplot=T,verbose=T)
{

  if (verbose) cat("\n----- REMOVING MEAN VALUE FROM STRATIGRAPHIC SERIES -----\n")
  dat<-data.frame(dat)
  ipts <- length(dat[,1]) 
  if (verbose) cat(" * Number of data points=", ipts,"\n")
  dave <- colMeans(dat[2])
  dat[2] <- dat[2] - dave
  if (verbose) cat(" * Mean value removed=",dave, "\n")

if(genplot)
  {
### plots
    par(mfrow=c(2,2))
    plot(dat,cex=0.5,xlab="Location",ylab="Value",main="Data Series"); lines(dat)
### plot the denisty and the histogram together
    hist(dat[,2],freq=F); lines(density(dat[,2], bw="nrd"),col="red"); grid()
### boxplot
    boxplot(dat[,2])
### Normal probabilty plot (Normal Q-Q Plot)
    qqnorm(dat[,2]); qqline(dat[,2], col="red")
  }
  
  return(dat)

### END function demean
}
