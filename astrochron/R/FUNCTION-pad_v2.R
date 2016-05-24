### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### pad: pad time series with zeros - (SRM: July 31, 2013; Sept. 10, 2015)
###########################################################################

pad <- function (dat,zeros=ipts,genplot=T,verbose=T)
{

  if (verbose) cat("\n----- PADDING STRATIGRAPHIC SERIES WITH ZEROS -----\n")
  dat<-data.frame(dat)
  ipts <- length(dat[,1]) 
  dt = dat[2,1]-dat[1,1]

# error checking 
  dtest <- dat[2:ipts,1]-dat[1:(ipts-1),1] 
  epsm=1e-9
  if( (max(dtest)-min(dtest)) > epsm ) 
     {
       cat("\n**** ERROR: sampling interval is not uniform.\n")
       stop("**** TERMINATING NOW!")
     }
  
  if (verbose) cat(" * Number of data points=", ipts,"\n")  
  if(verbose) cat(" * Sample interval=", dt,"\n")
  padded <- as.numeric(dat[,2])
  padded <- append(padded, rep(0,zeros) )
  addtime <- dat[ipts,1] + dt*(1:zeros)
  paddedtime <- append(dat[,1],addtime)

  out <- data.frame(cbind(paddedtime,padded))
  colnames(out)[1] <- colnames(dat[1])
  colnames(out)[2] <- colnames(dat[2])

if(genplot)
  {
### plots
    par(mfrow=c(2,2))
    plot(out,cex=0.5,xlab="Location",ylab="Value",main="Data Series"); lines(out)
### plot the denisty and the histogram together
    hist(out[,2],freq=F); lines(density(out[,2], bw="nrd"),col="red"); grid()
### boxplot
    boxplot(out[,2])
### Normal probabilty plot (Normal Q-Q Plot)
    qqnorm(out[,2]); qqline(out[,2], col="red")
  }
  
  return(out)

### END function pad
}
