### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### function bandpass - (SRM: January 26-27, 2012; March 28, 2012; 
###                      April 28, 2012; May 25, 2012; October 10, 2012; 
###                      November 17, 2012; May 1, 2013; May 20-24, 2013;
###                      June 5, 2013; June 13, 2013; July 10, 2013; 
###                      Aug. 9, 2013; Aug. 13-14, 2013; Aug. 15, 2013
###                      Aug. 16, 2013; Nov. 26, 2013; Sept. 16, 2014;
###                      Feb. 3, 2015; September 10, 2015)
###
### bandpass filter (allows rectangular, Gaussian and tapered cosine windows) 
###########################################################################


bandpass <- function (dat,padfac=2,flow=NULL,fhigh=NULL,win=0,alpha=3,p=0.25,demean=T,detrend=F,addmean=T,output=1,xmin=0,xmax=Nyq,genplot=T,verbose=T)
{

   if(verbose)   { cat("\n----- BANDPASS FILTERING STRATIGRAPHIC SERIES-----\n") } 
   dat <- data.frame(dat)
   npts <- length(dat[,1]) 
   dt = dat[2,1]-dat[1,1]

# error checking 
   if(dt<0)
     { 
       if (verbose) cat("\n * Sorting data into increasing height/depth/time, removing empty entries\n")
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

   d <- dat
   dt = d[2,1]-d[1,1]  
   if(verbose)   { cat(" * Number of data points=", npts,"\n") }
   if(verbose)   { cat(" * Sample interval=", dt,"\n")}
      
### remove mean
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
      if(verbose)  { cat(" * Linear trend removed. m=",lm.0$coeff[2],"b=",lm.0$coeff[1],"\n") }
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

   if(is.null(flow)) flow <- df
   if(is.null(fhigh)) fhigh <- Nyq-df
   if (flow > fhigh)  {fh <- fhigh; fhigh <- flow; flow <- fh}
   if (fhigh> (Nyq-df))  
     {
      fhigh <- (Nyq-df)
      if(verbose) cat(" * ERROR: maximum frequency too high, reset to Nyquist-df\n")
      }
   if (flow < df)  
     {
      flow <- df
      if(verbose) cat(" * ERROR: minimum frequency too low, reset to df\n") 
     }
   
### center of filter
   fcent=(flow+fhigh)/2

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
    f2 <- ( (nf/2) - (1:(minusdf-negfreq+1) ) ) * df * -1
    freq <- append(freq,f2)
### caculate amplitude 
    amp <- Mod(ft[1:nyqfreq])
### put results into a data frame
    fft.out <- data.frame(cbind(freq[1:nyqfreq],amp))

### plot some of the results
  if(genplot)
   {
     par(mfrow=c(2,2))
     plot(d,type="l", col="blue",ylab="Value", xlab="Location", main="Stratigraphic Series",bty="n")
     plot(fft.out[,1],fft.out[,2], type="l",col="red",xlim=c(xmin,xmax),ylim=c(0,max(amp)),ylab="Amplitude", xlab="Frequency", main="Amplitude",bty="n")
     abline(v=flow,col="purple",lty=3)
     abline(v=fhigh,col="purple",lty=3)
   }

### scaling factor for plotting of taper
     epsm=10^-13
     scaleT=max(amp)-min(amp)
     if(scaleT< epsm) scaleT = max(amp)
   
### now bandpass
    bpts = 0
# Apply rectangular window to zero and positive frequencies (first portion of the output)
    for (j in izero:nyqfreq) { if(freq[j] < flow || freq[j] > fhigh) {ft[j] = 0} else bpts=bpts+1 }
# Apply rectangular window to negative frequencies (second portion of the output)
    for (j in negfreq:minusdf) { if(freq[j] < -1*fhigh || freq[j] > -1*flow) {ft[j] = 0} }
    if(verbose)  
      { 
        cat(" * Center of bandpass filter =",fcent,"\n") 
        cat(" *",bpts,"pos/neg frequency pairs will be bandpassed\n") 
      }

### This portion of code is used if we want to apply a non-rectangular window
### For a Gaussian window, based on Harris (1978), p.70. alpha is 1/std. dev. 
### (a measure of the width of the Fourier Transform of the Dirichlet kernel)
### note that this function does not achieve zero, but approaches it for high alpha
### here are some minima that correspond to different alpha:
###  2.5 = 0.04711472; 3 = 0.01228416 ; 3.5 = 0.002508343.
    if(win == 1)
      {
# make taper
        taper <- double(bpts)
        for (j in 1:bpts)
          {
# center time series around x=0
            x= as.double(j) - ( ( as.double(bpts)/2 ) + 0.5)
            y= ( alpha * x / ( as.double(bpts)/2 ) )^2
            taper[j] = exp(-0.5*y)
           }
# apply taper
        k = 1
        ii = 0
        for (j in izero:nyqfreq) 
         { 
           if(freq[j] >= flow && freq[j] <= fhigh) 
             { 
               ft[j] = ft[j]*taper[k]
               k=k+1
               if (ii == 0) 
                 {
                   fhold <- freq[j] 
                   ii = ii +1
                 }
              } 
         }

# now apply Gaussian taper to negative frequencies (note, this taper is symmetric, so the direction
# doesn't matter)
        k = 1
        for (j in negfreq:minusdf) { if(freq[j] >= -1*fhigh && freq[j] <= -1*flow) {ft[j] = ft[j]*taper[k]; k=k+1} }
# add window to spectrum plot
        i <- 1:bpts
        bpfreq <- fhold+df*(i-1)
        if(genplot) 
          { 
            lines(bpfreq,taper*scaleT,col="blue") 
            points(fcent,max(amp),col="blue")
            mtext(round(fcent,digits=5),side=3,line=0,at=fcent,cex=0.5,font=4,col="blue")
          }
       }  

### For a simple cosine-tapered window. Based on Percival and Walden (1993), p. 209
    Pbpts=(bpts*p*.5)
    if(win == 2 && Pbpts >= 2)
      {
# make taper
        taper <- double(bpts)
        taper[1:bpts] <- 1 
        for (j in 1:Pbpts)
          {
           y= (2*pi*as.double(j) ) / ( (p*as.double(bpts) ) +1 )
# First part of taper
           taper[j] = 0.5*(1-cos(y))
# Last part of taper
           taper[bpts+1-j] = taper[j]
          }   
# apply taper
        k = 1
        ii = 0
        for (j in izero:nyqfreq) 
         { 
           if(freq[j] >= flow && freq[j] <= fhigh) 
             { 
               ft[j] = ft[j]*taper[k]
               k=k+1
               if (ii == 0) 
                 {
                   fhold <- freq[j] 
                   ii = ii +1
                 }
              } 
         }
# now apply cosine-taper to negative frequencies (note, this taper is symmetric, so the direction
# doesn't matter)
        k = 1
        for (j in negfreq:minusdf) { if(freq[j] >= -1*fhigh && freq[j] <= -1*flow) {ft[j] = ft[j]*taper[k]; k=k+1} }
# add window to spectrum plot
        i <- 1:bpts
        bpfreq <- fhold+df*(i-1)
        if(genplot) 
          { 
            lines(bpfreq,taper*scaleT,col="blue")
            points(fcent,max(amp),col="blue")
            mtext(round(fcent,digits=5),side=3,line=0,at=fcent,cex=0.5,font=4,col="blue")
          }
       }   

# if fewer than 4 tapered points (2 on each side of window), to not apply taper     
       if(win ==2 && Pbpts < 2) 
        {
           cat("\n**** WARNING: Too few data points to apply cosine-tapered window.\n")
           cat("              Will use rectangular window. Try increasing padfac,\n")
           cat("              and/or increasing p.\n")
        }
        
# inverse FFT
###  convert to real number and normalize iFFT output by length of record, suppress warning
###   about discarding imaginary component (it is zero!).
     ifft <- suppressWarnings(as.double(fft(ft,inverse=TRUE)/nf))
### isolate prepadded portion
     bp <- ifft[1:npts]
     d[2] <- bp
     d3<- data.frame( cbind(d[1],(bp + dave)) )
     colnames(d)[2] <- colnames(dat[2])
     colnames(d3)[2] <- colnames(dat[2])
     
     if(genplot)
       {
         plot(d, type="l",col="red", ylab="Value", xlab="Location", main="Bandpassed Signal",bty="n")
         plot(dat, type="l",col="blue", ylim=c(min(dat[,2],d3[,2]),max(dat[,2],d3[,2])), ylab="Value", xlab="Location", main="Comparison",bty="n")
         lines(d3, col="red")
       }
       
     if(output == 1)
      {
        if(!addmean) {return(d)}
        if(addmean) {return(d3)}
      }
      
     if(output == 2) return(data.frame(cbind(bpfreq,taper)))
     
#### END function bandpass
}