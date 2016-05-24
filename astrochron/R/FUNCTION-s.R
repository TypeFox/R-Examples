### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### s: standardize variable in stratigraphic data series
###              - subtract mean, and divide by standard devation
###              (SRM: August 5, 2013)
###                                        
###########################################################################

s <- function (dat,genplot=F,verbose=T)
{

  if (verbose) cat("\n----- STANDARDIZING VARIABLE IN STRATIGRAPHIC SERIES -----\n")
  dat<-data.frame(dat)
  ipts <- length(dat[,1]) 
  if (verbose) cat(" * Number of data points=", ipts,"\n")
  dave <- colMeans(dat[2])
  dat[2] <- dat[2] - dave
  dstd <- sd(dat[,2])
  dat[2] <- dat[2]/dstd
  if (verbose) cat(" * Mean value removed=",dave, "\n")
  if (verbose) cat(" * Divided by standard deviation=",dstd, "\n")

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

### END function s
}
