### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### replEps: replace values <=0 with smallest postive value (eps)
###            (SRM: Dec. 14, 2012; April 30, 2013; May 20, 2013; June 21, 2013)
###
###########################################################################

replEps <- function (dat, ID=T, genplot=T, verbose=T)
{

  if(verbose) cat("\n----- REPLACING VALUES <=0 WITH SMALLEST POSTIVE MACHINE VALUE -----\n")

  dat <- data.frame(dat)
  ipts <- length(dat[,1]) 
  if(verbose) cat(" * Number of data points=", ipts,"\n")
  cols=length(dat)
  if(verbose) cat(" * Number of columns=", cols,"\n")
 
  numneg = sum(dat[,2:cols] <= 0)
  if(verbose) cat(" * Number of negative or zero values=", numneg,"\n")
  
  repl.eps <- function(x) {x[x <= 0] <- .Machine$double.eps; x}
  dat[,2:cols] <- repl.eps(dat[,2:cols])
  
  if(cols == 2 && genplot)
    {
### plots
     par(mfrow=c(2,2))
     plot(dat,cex=0.5,xlab="Location",ylab="Value",main="Data Series"); lines(dat)
### plot the denisty and the histogram together
     hist(dat[,2],freq=F,xlab="Value",main="Distribution of isoloated values"); lines(density(dat[,2], bw="nrd"),col="red")
### boxplot
     boxplot(dat[,2],ylab="Value",main="Boxplot for isolated values")
### Normal probabilty plot (Normal Q-Q Plot)
     qqnorm(dat[,2]); qqline(dat[,2], col="red")
    }
   
  return(dat)

### END function replEps
}
