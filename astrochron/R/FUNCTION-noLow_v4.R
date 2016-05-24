### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### noLow: remove lowess smoother from data  - (SRM: March 12, 2012; May 31, 2012;
###                                                  April 24, 2013; May 20, 2013;
###                                                  June 21, 2013)
###
###########################################################################

noLow <- function (dat,smooth=.20,output=1,genplot=T,verbose=T)
{

if(verbose) cat("\n----- REMOVING LOWESS SMOOTHER FROM STRATIGRAPHIC SERIES -----\n")

dat<-data.frame(dat)
### fit a LOWESS smoother to data series
smooth.lo <- loess((dat[,2]) ~ dat[,1], degree=1, span= smooth )
y.predict <-predict(smooth.lo, se= TRUE)

  if(verbose) print(summary(smooth.lo))
  
  out <- dat
  out[2] <- y.predict$fit
  resid <- dat
  resid[2] <- dat[2] - y.predict$fit

### plot linearly detrended data
if(genplot)
 {
  par(mfrow=c(2,2))
  plot(dat,cex=0.5, xlab="Location", ylab="Value", main="Data with LOWESS Fit"); lines(dat[,1],y.predict$fit,col="red")
  plot(resid,cex=0.5, xlab="Location", ylab="Value", main="Residuals from LOWESS Fit"); lines(resid)
### plot the denisty and the histogram together
  hist(resid[,2],freq=F, xlab="Value", main="Distribution of Residual Values"); lines(density(resid[,2], bw="nrd"),col="red"); grid()
### boxplot
  boxplot(resid[,2], ylab="Value", main="Boxplot of Residual Values")
### Normal probabilty plot (Normal Q-Q Plot)
#  qqnorm(resid[,2]); qqline(out[,2], col="red")
 }
 
  if(output==1) return(data.frame(resid))   
  if(output==2) return(data.frame(out))

### END function noLow
}
