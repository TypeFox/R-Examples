### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### gausTaper: apply cosine taper to data - (SRM: January 26, 2012; Oct. 11, 2012;
###                                               May 20-27, 2013; June 5, 2013;
###                                               June 21, 2013; July 29, 2013;
###                                               Aug. 3, 2013; Aug. 7-8, 2013;
###                                               Aug. 13, 2013; Sept. 10, 2015)
###
### Based on Harris (1978), p.70
### alpha is 1/std. dev. (a measure of the width of the Fourier Transform
### of the Dirichlet kernel)
###########################################################################

gausTaper <- function (dat=NULL,alpha=3,rms=T,demean=T,detrend=F,genplot=T,verbose=T)
{

   if (verbose) cat("\n----- APPLYING GAUSSIAN TAPER TO STRATIGRAPHIC SERIES -----\n")
   if(!is.null(dat)) dat <- data.frame(dat)
# generate series to evalute spectral properties of window
   if(is.null(dat)) 
     {
      dat=data.frame( cbind(1:256,rep(1,256)) )
      demean=F
      detrend=F
      genplot=F
     }    

   ipts <- length(dat[,1]) 
   if (verbose) cat(" * Number of data points=", ipts, "\n")

# error checking 
   dtest <- dat[2:ipts,1]-dat[1:(ipts-1),1] 
   epsm=1e-9
   if( (max(dtest)-min(dtest)) > epsm ) 
     {
       cat("\n**** ERROR: sampling interval is not uniform.\n")
       stop("**** TERMINATING NOW!")
     }

if(genplot)
 {
### plots
  par(mfrow=c(2,2))
  plot(dat,cex=0.5,xlab="Location",ylab="Value",main="Data Series"); lines(dat)
 }

###########################################################################
### remove mean and linear trend
###########################################################################
   dave <- colMeans(dat[2])
   if (demean) 
     { 
       dat[2] <- dat[2] - dave
       if(verbose) { cat(" * Mean value removed=",dave,"\n") }
     }
### use least-squares fit to remove linear trend
   if (detrend) 
    {
      lm.0 <- lm(dat[,2] ~ dat[,1])
      dat[2] <- dat[2] - (lm.0$coeff[2]*dat[1] + lm.0$coeff[1])
      if(verbose) {cat(" * Linear trend removed. m=",lm.0$coeff[2],"b=",lm.0$coeff[1],"\n") }
    }

   taper <- double(ipts)
   for (i in 1:ipts)
     {
# Center time series around x=0
         x= as.double(i) - ( ( as.double(ipts)/2 ) + 0.5)
         y= ( alpha * x / ( as.double(ipts)/2 ) )^2
         taper[i] = exp(-0.5*y)
      }

# normalize taper to RMS=1 to preserve power for white process
   if(rms)
     {
       norm = sqrt( sum(taper*taper)/ipts )
       taper = taper/norm
      }

# now multiply taper and data
   dat[2] <- dat[2] * taper
   
# multiply taper by maximum data value for plotting purposes only!
   tt <- cbind(dat[1],(taper/max(taper))*max(dat[2]))

if(genplot)
  {
### more plots
    plot(dat,cex=0.5,xlab="Location",ylab="Tapered Value",main="Tapered Data Series"); lines(dat,lty=2); lines(tt,col="blue")
    hist(dat[,2],freq=F,xlab="Tapered Value",main="Distribution of Tapered Values"); lines(density(dat[,2], bw="nrd"),col="red");grid()
### Normal probabilty plot (Normal Q-Q Plot)
    qqnorm(dat[,2]); qqline(dat[,2], col="red")
   }
   
  return(dat)

### END function gausTaper
}
