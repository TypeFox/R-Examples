### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### function hilbert - (SRM: January 26-February 1, 2012; April 26, 2012; 
###                     May 25, 2012; November 17, 2012; November 23, 2012;
###                     May 20-24, 2013; May 28, 2013; June 5, 2013; 
###                     June 13, 2013; Nov. 26, 2013; January 13, 2014; 
###                     January 27, 2015; January 31, 2015; February 3, 2015;
###                     Sept. 10, 2015)
###
### hilbert transform
###########################################################################

hilbert <- function (dat,padfac=2,demean=T,detrend=F,output=T,addmean=F,genplot=T,verbose=T)
{

   if(verbose) { cat("\n----- PERFORMING HILBERT TRANSFORM ON STRATIGRAPHIC SERIES -----\n") }
   dat <- data.frame(dat)
   npts <- length(dat[,1]) 
   dt = dat[2,1]-dat[1,1]

# error checking 
   if(dt<0)
     { 
       if (verbose) cat(" * Sorting data into increasing height/depth/time, removing empty entries\n")
       dat <- dat[order(dat[1], na.last = NA, decreasing = F), ]
       dt <- dat[2,1]-dat[1,1]
       npts <- length(dat[,1]) 
     }
   dtest <- dat[2:npts,1]-dat[1:(npts-1),1] 
   epsm=1e-9
   if( (max(dtest)-min(dtest)) > epsm ) 
     {
       cat("\n**** ERROR: sampling interval is not uniform.\n")
       stop("**** TERMINATING NOW!")
     }

   if(!demean && addmean)
    {
       cat("\n**** WARNING: will not add mean value to amplitude envelope, since demean=F.\n")
    }

   d <- dat
   if(verbose) { cat(" * Number of data points=", npts,"\n") }
   if(verbose) { cat(" * Sample interval=", dt,"\n") }
   
###########################################################################
### remove mean and linear trend
###########################################################################
   dave <- colMeans(d[2])
   if (demean) 
     { 
       d[2] <- d[2] - dave
       if(verbose) { cat(" * Mean value removed=",dave,"\n") }
     }
### use least-squares fit to remove linear trend
   if (detrend) 
    {
      lm.0 <- lm(d[,2] ~ d[,1])
      d[2] <- d[2] - (lm.0$coeff[2]*d[1] + lm.0$coeff[1])
      if(verbose) {cat(" * Linear trend removed. m=",lm.0$coeff[2],"b=",lm.0$coeff[1],"\n") }
    }
    
### Calculate Nyquist
   Nyq <- 1/(2*dt)
#### Calculate Rayleigh Frequency
   Ray <- 1/(npts*dt)
   
### pad with zeros if desired (power of 2 not required!)
### also convert from data frame to numeric
   pad <- as.numeric(d[,2])
   if(padfac>1) pad <- append( pad, rep(0,(npts*padfac-npts)) )

### add another zero if we don't have an even number of data points, so Nyquist exists.   
   if((npts*padfac)%%2 != 0) pad <- append(pad,0)
    
### new resulting frequency grid   
   nf = length(pad)
   df <- 1/(nf*dt)
   freq <- double(nf)
   
### take fft
   ft <- fft(pad)   
# Locate real components for zero, nyquist, first neg freq., (zero-df)
   izero = 1
   nyqfreq = 0.5*nf + 1
   negfreq = 0.5*nf + 2
   minusdf = nf
   
### set frequency index vector 
    i <- seq(1,nf,by=1)
### assign positive frequencies out to Nyquist
    freq <- df*(i[1:nyqfreq]-1)
## assign negative frequencies
    f2 <- ( (nf/2) - 1:minusdf ) * df * -1
    freq <- append(freq,f2)    
### caculate amplitude
    amp <- Mod(ft[1:nyqfreq])
### put results into a common matrix
    fft.out <- data.frame(cbind(freq[1:nyqfreq],amp))

if(genplot)
 {
### plot some of the results
   par(mfrow=c(3,1))
   plot(fft.out[,1],fft.out[,2], type="l",col="red",ylab="Amplitude", xlab="Frequency", main="Amplitude Spectrum for Data",bty="n")
 }
 
### Hilbert Transform follows the method outlined in Kanasewich (1973).
### The Hilbert Transform (aka quadrature function) will be stored in 
### the vector ht. The Hilbert Transform is a 90 degree phase 
### shift operator.
### For the frequency of 0, multiply the real and imaginary 
### components by 0
      ht <- complex(nf)
      ht[1]=0
### For frequencies > 0, multiply real and imaginary components by
### i.  Note that this results in the real and imaginary components being
### swapped
      ht[2:nyqfreq]=(Im(ft[2:nyqfreq])*-1) + 1i*(Re(ft[2:nyqfreq]))

### For frequencies < 0, multiply real and imaginary components by
### -i.  Note that this results in the real and imaginary components being
### swapped
      ht[negfreq:minusdf]=Im(ft[negfreq:minusdf]) + 1i*(Re(ft[negfreq:minusdf])*-1)
        
### inverse FFT
###  normalize iFFT output by length of record
     ift <- fft(ht,inverse=TRUE)/nf

### Calculate the instantaneous amplitude
### instantaneous amplitude= envelope function = modulus of complex function
### we only need to look at first 'npts', others are from zero padding
      envelope <- double(npts)
### The complex function (aka analytical function) is: (Original function) - i (Hilbert Transform)
### We will calculate amplitude from this complex function:
###        amplitude = (real)**2 + (imag)**2
### Note that the imaginary component of both 'Original function' and 'Hilbert Transform' is zero.
### See Kanasewich (1981) for additional information 
     envelope = sqrt ( d[1:npts,2]^2 + Re(ift[1:npts])^2 )
           
     d2 <- data.frame(cbind(d[,1],envelope))
     d3<- data.frame( cbind(d[,1],(envelope + dave)) )
     d3b<- data.frame( cbind(d[,1],envelope) )
     d4<- data.frame( cbind(d[,1],(dave-envelope)) )
     d4b<- data.frame( cbind(d[,1],(-1*envelope)) )

if(genplot)
  {
     plot(d2, type="l",col="red", ylab="Value", xlab="Location", main="Instantaneous Amplitude",bty="n")
     plot(dat, type="l",col="blue", ylab="Value", xlab="Location", main="Comparison",bty="n")
     if(demean)
      {
        lines(d3, col="red")
        lines(d4, col="red")
      }  
     if(!demean)
      {     
        lines(d3b, col="red")
        lines(d4b, col="red")
      }  
  }   
 if(output == T ) 
       {
         if(!addmean || !demean) return(d2)
         if(addmean && demean) return(d3)
       }

#### END function hilbert
}