### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### armaGen : Generate autoregressive-moving average series - (SRM: October 2-3, 2012; 
###                           April 24, 2013; May 20, 2013; June 21, 2013)
###########################################################################

armaGen <- function (npts=1024,dt=1,m=0,std=1,rhos=c(0.9),thetas=c(0),genplot=T,verbose=T)
{

if(verbose) cat("\n----- GENERATING AUTOREGRESSIVE-MOVING AVERAGE SERIES -----\n")
### note that this creates a ts object.
    ar.sim <- arima.sim(list(ar=rhos,ma=thetas),n=npts+1)
### generate time axis
    ta <- 1:npts
### change time axis from unit spacing of 1 to desired value
    ta <- (ta*dt) - dt
 
    noise=data.frame(cbind(ta,as.numeric(ar.sim[2:(npts+1)])))

### set mean and standard deviation
    noise[2] = noise[2] - colMeans(noise[2])
    noise[2] = std*noise[2]/sapply(noise[2],sd)
    noise[2] = m + noise[2]

### plot noise model
    if(genplot)
      {
        par(mfrow=c(2,2))
        plot(noise, cex=.5,xlab="Location",ylab="Modeled Value", main="AR-MA Data Series"); lines(noise)
### plot the denisty and the histogram together
        hist(noise[,2],freq=F,xlab="Modeled Value",main="Histogram of Modeled Values"); lines(density(noise[,2], bw="nrd"),col="red");grid()
### boxplot
        boxplot(noise[,2],ylab="Modeled Value", main="Boxplot of Modeled Values")
### Normal probabilty plot (Normal Q-Q Plot)
        qqnorm(noise[,2]); qqline(noise[,2], col="red");grid()
       }
 
    return(noise)
 
 ### END function armaGen
}

