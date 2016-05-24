### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### function prewhiteAR - (SRM: October 3, 2012; Oct. 11, 2012; Oct 15, 2012
###                             April 24, 2013; May 20, 2013; June 13, 2013;
###                             September 10, 2015)              
###
### prewhiten with AR filter, as selected via MLE or other methods
###########################################################################

prewhiteAR <- function (dat,order=0,method="mle",aic=T,genplot=T,verbose=T) 
{

if(verbose) cat("\n----- PREWHITENING STRATIGRAPHIC SERIES WITH AUTOREGRESSIVE FILTER -----\n")

if(order<=0 && aic==F) 
  {
    if(verbose) 
      {
        cat("\n**** WARNING: Selecting an order <=0 initiates automatic evaluation out\n")
        cat("              to the maximum possible order permitted by the algorithm,\n")
        cat("              but you have entered aic=F. Now changing to aic=T.\n\n")
      }
    aic=T
  }  


### dat contains a time/depth series (time, value)
   dat<-data.frame(dat)
   npts <- length(dat[,1]) 
   dt <- dat[2,1]-dat[1,1]

# error checking 
   dtest <- dat[2:npts,1]-dat[1:(npts-1),1] 
   epsm=1e-9
   if( (max(dtest)-min(dtest)) > epsm ) 
     {
       cat("\n**** ERROR: sampling interval is not uniform.\n")
       stop("**** TERMINATING NOW!")
     }

   if(verbose) 
     {
       cat(" * Number of data points=", npts,"\n")
       cat(" * Sampling interval (space or time):",dt,"\n")
     }  

# convert to time series for ar routine
    x <- as.ts(dat[,2]) 

# use MLE
if(order<=0) { x.ar = ar(x,aic=aic,method=method) }
if(order>0) { x.ar = ar(x,aic=aic,order.max=order,method=method) }

    if(verbose) 
      {
        cat(" * AR Order=",x.ar$order,"\n")
        cat(" * AR Coefficients:\n")
        print(x.ar$ar)
      }
    
# the prewhitened signal
    out <- data.frame(cbind(dat[(x.ar$order+1):npts,1],x.ar$resid[(x.ar$order+1):npts]))

### plot data
  if(genplot)
   {    
### plot data series. Note, cex is the factor by which to increase or decrease default symbol size
     par(mfrow=c(3,2))
     plot(dat,cex=.5,main="Original data series",xlab="Position",ylab="Value")
     lines(dat)
     plot(out,cex=.5,main="Whitened data series",xlab="Position",ylab="Whitened value")
     lines(out)
     if(aic) 
     {
        aicOrder=0L:(length(x.ar$aic)-1)
        plot(aicOrder,x.ar$aic,log="y",ylab="AIC",main="Model selection")
        abline(v=x.ar$order,col="red")
     }
### plot the denisty and the histogram together
     hist(out[,2],freq=F,main="Distribution of whitened data",xlab="Whitened value"); lines(density(out[,2], bw="nrd"),col="red")
### boxplot
     boxplot(out[,2],main="Boxplot of whitened data")
### Normal probabilty plot (Normal Q-Q Plot)
     qqnorm(out[,2]); qqline(out[,2], col="red")
   }
   
   return(out)

### END function prewhiteAR
}
