### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### divTrend: divide data value by linear trend - (SRM: March 29, 2012; Oct. 11, 2012; 
###                                                May 20, 2013; June 21, 2013)
###
### motivated by Bloomfield (2000), pg 134.
###########################################################################

divTrend <- function (dat,genplot=T,verbose=T)
{

  if(verbose) cat("\n----- DIVIDE STRATIGRAPHIC SERIES BY LINEAR TREND-----\n")
  dat=data.frame(dat)
### use least-squares fit to remove linear trend
  lm.0 <- lm(dat[,2] ~ dat[,1])
  if(verbose) cat(" * Slope=",lm.0$coeff[2],"\n")
  if(verbose) cat(" * y-intercept=",lm.0$coeff[1],"\n")
  
  y.predict <-predict(lm.0, se= TRUE)
  
  out <- dat
   out[2] <- dat[2]/y.predict$fit

if(genplot)
  {
### plot detrended data
    par(mfrow=c(2,2))
    plot(dat,cex=0.5, xlab="Location", ylab="Value", main="Data with Trend Line"); lines(dat[,1],y.predict$fit,col="red")
    plot(out,cex=0.5, xlab="Location", ylab="Value", main="Data/Trend"); lines(out)
### plot the denisty and the histogram together
    hist(out[,2],freq=F, xlab="Value", main="Distribution of Data/Trend Values"); lines(density(out[,2], bw="nrd"),col="red"); grid()
### boxplot
    boxplot(out[,2], ylab="Value", main="Boxplot of Data/Trend Values")
### Normal probabilty plot (Normal Q-Q Plot)
#  qqnorm(out[,2]); qqline(out[,2], col="red")
   }
   
  return(out)

### END function divTrend
}
