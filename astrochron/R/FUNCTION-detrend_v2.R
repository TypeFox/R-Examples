### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### detrend: linearly detrend data - (SRM: January 24, 2012; March 29, 2012;
###                                   Oct. 11, 2012; Dec. 12, 2012; May 20, 2013)
###
###########################################################################

detrend <- function (dat,genplot=T,verbose=T)
{

  if(verbose) cat("\n----- SUBTRACTING LINEAR TREND FROM STRATIGRAPHIC SERIES -----\n")
  dat<-data.frame(dat)

### use least-squares fit to remove linear trend
  lm.0 <- lm(dat[,2] ~ dat[,1])
  if(verbose) cat(" * Slope=",lm.0$coeff[2],"\n")
  if(verbose) cat(" * y-intercept=",lm.0$coeff[1],"\n")
  
  y.predict <-predict(lm.0, se= TRUE)
  
  out <- dat
   out[2] <- dat[2] - y.predict$fit

if(genplot)
  {
### plot linearly detrended data
    par(mfrow=c(2,2))
    plot(dat,cex=0.5, xlab="Location", ylab="Value", main="Data with Trend Line"); lines(dat[,1],y.predict$fit,col="red")
    plot(out,cex=0.5, xlab="Location", ylab="Residual Value", main="Residuals from Trend Line"); lines(out)
### plot the denisty and the histogram together
    hist(out[,2],freq=F, xlab="Residual Value", main="Distribution of Residual Values"); lines(density(out[,2], bw="nrd"),col="red");grid()
### boxplot
    boxplot(out[,2], ylab="Residual Value", main="Boxplot of Residual Values")
### Normal probabilty plot (Normal Q-Q Plot)
#  qqnorm(out[,2]); qqline(out[,2], col="red")
   }
   
  return(out)

### END function detrend
}
