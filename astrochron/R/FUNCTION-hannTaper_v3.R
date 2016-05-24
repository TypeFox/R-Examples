### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### hannTaper: apply Hann (a.k.a. Hanning) taper to data 
###                      - (SRM: July 30, 2013; Aug. 3, 2013; Aug. 8, 2013
###                              Aug. 13, 2013; Sept. 10, 2015)
###
###########################################################################

hannTaper <- function (dat=NULL,rms=T,demean=T,detrend=F,genplot=T,verbose=T)
{

   if (verbose) cat("\n----- APPLYING HANN TAPER TO STRATIGRAPHIC SERIES -----\n")
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
   if (verbose) cat(" * Number of data points=", ipts,"\n")

# error checking 
   dtest <- dat[2:ipts,1]-dat[1:(ipts-1),1] 
   epsm=1e-9
   if( (max(dtest)-min(dtest)) > epsm ) 
     {
       cat("\n**** ERROR: sampling interval is not uniform.\n")
       stop("**** TERMINATING NOW!")
     }

### plots
if(genplot)
  {
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

# Generate Hann taper
   taper=.5*(1-cos(2*pi*(1:ipts)/(ipts+1)))

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
   hist(dat[,2],freq=F,xlab="Tapered Value",main="Distribution of Tapered Values"); lines(density(dat[,2], bw="nrd"),col="red")
### Normal probabilty plot (Normal Q-Q Plot)
   qqnorm(dat[,2]); qqline(dat[,2], col="red")
   }

  return(dat)

### END function hannTaper
}
