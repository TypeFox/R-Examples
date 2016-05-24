### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### clipIt: clip or mute stratigraphic series below a theshold value to transfer 
###       power to modulator (SRM: December 25, 2013; January 5, 2014)
###
###########################################################################

clipIt <- function (dat,thresh=NULL,clipval=NULL,clipdiv=NULL,genplot=T,verbose=T)
{

  if (verbose) cat("\n----- CLIP STRATIGRAPHIC SERIES -----\n")
  dat<-data.frame(dat)
# make copy for clipping
  dat2 <- dat
  ipts <- length(dat[,1]) 
  if (verbose) cat(" * Number of data points=", ipts,"\n")
  dave <- colMeans(dat[2])
  dmax <- max(dat[2])
  dmin <- min(dat[2])
  if(is.null(thresh)) 
   {
     thresh <- dave
     if(verbose) cat(" * Will use mean value as clipping threshold\n")
    } 
  if(is.null(clipval)) clipval <- thresh
  if(clipval > dmax && verbose) cat(" * WARNING: threshold is above maximum value\n")
  if(clipval < dmin && verbose) cat(" * WARNING: threshold is below minimum value\n")
  
  if(is.null(clipdiv)) 
    {
      dat2[dat2[2]<thresh,2] <- clipval
     } 
  if(!is.null(clipdiv))
    {
      dat2[dat2[2]<thresh,2] <- dat2[dat2[2]<thresh,2]/clipdiv
    }

if(genplot)
  {
### plots
    par(mfrow=c(2,2))
    plot(dat,cex=0.5,xlab="Location",ylab="Value",main="Original (black) and Clipped (red) Data", type="l"); points(dat2, col="red",cex=0.5)
### plot the denisty and the histogram together
    hist(dat2[,2],freq=F); lines(density(dat2[,2], bw="nrd"),col="red"); grid()
### boxplot
    boxplot(dat2[,2])
### Normal probabilty plot (Normal Q-Q Plot)
    qqnorm(dat2[,2]); qqline(dat2[,2], col="red")
  }
  
  return(dat2)

### END function clipIt
}
