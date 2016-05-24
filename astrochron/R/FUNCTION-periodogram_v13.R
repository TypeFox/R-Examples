### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### function periodogram - (SRM: January 26, 2012; April 28, 2012; 
###                         October 10, 2012, November 20, 2012; May 20-24, 2013; 
###                         June 5, 2013; June 13, 2013; July 30-31, 2013;
###                         August 3, 2013; August 7-10, 2013; Nov. 26, 2013;
###                         October 18, 2014; October 22, 2014; January 21, 2015;
###                         September 10, 2015)
###
### simple unwindowed periodogram
###########################################################################

periodogram <- function (dat,padfac=2,demean=T,detrend=F,nrm=1,xmin=0,xmax=Nyq,pl=1,output=0,f0=F,genplot=T,verbose=T)
{

if(verbose) cat("\n----- CALCULATING PERIODOGRAM FOR STRATIGRAPHIC SERIES -----\n")

   d <- data.frame(dat)
   npts <- length(d[,1]) 
   dt = d[2,1]-d[1,1]

# error checking 
   if(dt<0)
     { 
       if (verbose) cat("\n * Sorting data into increasing height/depth/time, removing empty entries\n")
       d <- d[order(d[1], na.last = NA, decreasing = F), ]
       dt <- d[2,1]-d[1,1]
       npts <- length(d[,1])
     }
   dtest <- d[2:npts,1]-d[1:(npts-1),1] 
   epsm=1e-9
   if( (max(dtest)-min(dtest)) > epsm ) 
     {
       cat("\n**** ERROR: sampling interval is not uniform.\n")
       stop("**** TERMINATING NOW!")
     }

  if (verbose) 
    {
      cat(" * Number of data points in stratigraphic series:",npts,"\n")
      cat(" * Stratigraphic series length (space or time):",(npts-1)*dt,"\n")
      cat(" * Sampling interval (space or time):",dt,"\n")
    }

###########################################################################
### remove mean and linear trend
###########################################################################
   if (demean) 
     { 
       dave <- colMeans(d[2])
       d[2] <- d[2] - dave
       if(verbose) cat(" * Mean value removed=",dave,"\n")
     }
     
   if (!demean && verbose) cat(" * Mean value NOT subtracted\n")

### use least-squares fit to remove linear trend
    if (detrend) 
      {
      lm.0 <- lm(d[,2] ~ d[,1])
      d[2] <- d[2] - (lm.0$coeff[2]*d[1] + lm.0$coeff[1])
      if(verbose) cat(" * Linear trend removed. m=",lm.0$coeff[2],"b=",lm.0$coeff[1],"\n")
      }

   if (!detrend && verbose) cat(" * Linear trend NOT subtracted\n") 
 
### Calculate Nyquist
   Nyq <- 1/(2*dt)
#### Calculate Rayleigh Frequency
   Ray <- 1/(npts*dt)

### pad with zeros if like (power of 2 not required!)
### also convert from data frame to numeric
  pad <- as.numeric(d[,2])
  if(padfac>1) pad <- append( pad, rep(0,(npts*padfac-npts)) )
  
### add another zero if we don't have an even number of data points, so Nyquist exists.   
   if((npts*padfac)%%2 != 0) pad <- append(pad,0)

if(verbose)
 {
   cat(" * Nyquist frequency:",Nyq,"\n")
   cat(" * Rayleigh frequency:",Ray,"\n")
   cat(" * Padded to",length(pad),"points\n")
 }

### new resulting frequency grid   
   nf = length(pad)
   df <- 1/(nf*dt)
   freq <- double(nf)
 
### set frequency index vector 
    i <- seq(1,nf,by=1)
     
### take fft
   ft <- fft(pad)
### apply normalization   
   if(nrm == 1) ft=ft/npts
### caculate power 
   pwr <- Mod(ft)^2
### caculate amplitude
   amp <- sqrt(pwr)
### calculate phase
   phase <- atan2(Im(ft),Re(ft))
### now make frequency vector (negative frequencies are not assigned correctly here; we will discard)
   freq <- df*(i-1)
### put all results into a common data frame
   fft.out <- data.frame(cbind(freq,amp,pwr,phase))
### extract results from 0 to positive nyquist (recall that the negative frequencies were not 
###   assigned correctly; they are listed in 'freq' as > Nyquist)
   fft.out <- subset(fft.out,(fft.out[,1] <= Nyq))
   if(!f0) fft.out <- subset(fft.out,(fft.out[,1] > 0))
   colnames(fft.out) <- c('Frequency','Amplitude','Power','Phase')

### do the same for the Fourier coefficients if desired
   if (output==2)
     {
      fc.out <- data.frame(cbind(freq,Re(ft),Im(ft)))
      fc.out <- subset(fc.out,(fc.out[,1] <= Nyq))
      if(!f0) fc.out <- subset(fc.out,(fc.out[,1] > 0))
      colnames(fc.out)[1] <- 'Frequency'
      colnames(fc.out)[2] <- 'Real Coeff.'
      colnames(fc.out)[3] <- 'Imag. Coeff.'
     }

   if(genplot)
    {
      par(mfrow=c(2,2))
### plot the results
      plot(d,type="l", col="blue",ylab="Value", xlab="Location", main="Stratigraphic Series",bty="n")
      if (pl == 2) {plot(fft.out[,1],fft.out[,3], type="l",col="red", ylab="Power", xlim=c(xmin,xmax),xlab="Frequency", main="Periodogram Power",bty="n")}
# do not plot f(0) if using log spectrum
      if (pl == 1)  
        { 
         fft.out2 <- subset(fft.out,(fft.out[,1] > 0))
         plot(fft.out2[,1],log(fft.out2[,3]), type="l",col="red", ylab="Log Power", xlim=c(xmin,xmax), xlab="Frequency", main="Log Periodogram Power",bty="n")
        }
     plot(fft.out[,1],fft.out[,2], type="l",col="red", ylab="Amplitude", xlim=c(xmin,xmax), xlab="Frequency", main="Periodogram Amplitude",bty="n")
     plot(fft.out[,1],fft.out[,4], type="l",col="red", ylab="Phase", xlim=c(xmin,xmax), xlab="Frequency", main="Periodogram Phase",bty="n")
   }
   
if (output==1)  return(fft.out)
if (output==2)  return(fc.out)
   
#### END function periodogram
}