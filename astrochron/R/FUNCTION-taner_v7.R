### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### function taner - (SRM: August 7, 2013; August 9, 2013; August 13, 2013;
###                         Nov. 26, 2013; Sept. 16, 2014; Jan. 8, 2014; 
###                         Feb. 3, 2015; March 3, 2015; July 31, 2015;
###                         Sept. 10, 2015)
###
###########################################################################

taner <- function (dat,padfac=2,flow=NULL,fhigh=NULL,roll=10^3,demean=T,detrend=F,addmean=T,output=1,xmin=0,xmax=Nyq,genplot=T,verbose=T)
{

   if(verbose)   { cat("\n----- TANER BANDPASS FILTERING STRATIGRAPHIC SERIES-----\n") } 
   dat <- data.frame(dat)
   npts <- length(dat[,1]) 
   dt <- dat[2,1]-dat[1,1]

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
### pad with zeros if desired (power of 2 not required)
### also convert from data frame to numeric
  pad <- as.numeric(d[,2])
  if(padfac>1) pad <- append( pad, rep(0,(npts*padfac-npts)) )

### add another zero if we don't have an even number of data points, so Nyquist exists.   
   if((npts*padfac)%%2 != 0) pad <- append(pad,0)

### new resulting frequency grid   
   nf = length(pad)
   df <- 1/(nf*dt)
   freq <- double(nf)

   if (is.null(fhigh)) fhigh = Nyq
# when flow is not specified, automatically perform lowpass filtering
   if (is.null(flow)) 
    {
       if(verbose)  { cat("\n **** NOTE: flow not specified. Will perform lowpass filtering.\n")}    
       flow = -1*fhigh
    }   
       
   if (flow > fhigh)  {fh <- fhigh; fhigh <- flow; flow <- fh}

   if(fhigh>Nyq)
    {
      if(verbose) { cat("\n**** WARNING: fhigh exceeds Nyquist frequency.\n")}
    }

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
### assign negative frequencies
    f2 <- ( (nf/2) - (1:(minusdf-negfreq+1) ) ) * df * -1
    freq <- append(freq,f2)
### caculate amplitude 
    amp <- Mod(ft[1:nyqfreq])
### put results into a common matrix
    fft.out <- data.frame(cbind(freq[1:nyqfreq],amp))

### plot some of the results
  if(genplot)
   {
     par(mfrow=c(2,2))
     plot(d,type="l", col="blue",ylab="Value", xlab="Location", main="Stratigraphic Series",bty="n")
     plot(fft.out[,1],fft.out[,2], type="l",col="red",xlim=c(xmin,xmax),ylab="Amplitude", xlab="Frequency", main="Amplitude",bty="n")
     abline(v=flow,col="purple",lty=3)
     abline(v=fhigh,col="purple",lty=3)
   }
   
### assign parameters for taner filter, then generate filter. this portion modified from 
### L. Hinnov FORTRAN code TANER.FOR
### also see: http://www.rocksolidimages.com/pdf/attrib_revisited.htm#_Toc328470897
### center of filter
    fcent=(flow+fhigh)/2
### convert to angular frequency
    wl = 2*pi*flow
    wc = 2*pi*fcent
    wh = 2*pi*fhigh
    bw = wh - wl
### amp2 is the maximum value of the taper
### set to be similar to Butterworth filter
    amp2 = 1/sqrt(2)
### roll is the roll-off octave
    arg1 = 1 - (roll*log(10)) / (20*log(amp2))
    arg1 = log(arg1)
    arg2 = ((bw+2)/bw)^2
    arg2 = log(arg2)
    twod = 2*arg1/arg2
### generate filter for positive frequencies
    dw = 2*pi/(nf*dt)
    w = ((1:nyqfreq)-1)*dw
    arg = (2*abs(w-wc)/bw)^twod
    darg=-1*arg
    filter_pos =  (amp2 * exp(darg))
### generate filter for negative frequencies
    w = ((negfreq:nf)-nf-1)*dw
    aw = abs(w)
    arg = (2*abs(aw-wc)/bw)^twod
    filter_neg = (amp2 * exp(-arg))

### normalize window maximum to 1 (rather than 1/sqrt(2)
    filter_pos=filter_pos/max(filter_pos)
    filter_neg=filter_neg/max(filter_neg)

    taper=c(filter_pos,filter_neg)   

### apply filter
    ft = ft*taper

    if(genplot) 
        { 
            lines(freq[1:nyqfreq],filter_pos*(max(amp)),col="blue")
            points(fcent,max(amp),col="blue")
            if(flow>=0) points(flow,max(amp),col="blue")
            points(fhigh,max(amp),col="blue")
            mtext(round(fcent,digits=5),side=3,line=0,at=fcent,cex=0.5,font=4,col="blue")
            if(flow>=0) mtext(round(flow,digits=5),side=3,line=0,at=flow,cex=0.5,font=4,col="blue")
            mtext(round(fhigh,digits=5),side=3,line=0,at=fhigh,cex=0.5,font=4,col="blue")
        }

# inverse FFT
###  convert to real number and normalize iFFT output by length of record, supress warning
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
         plot(dat, type="l",col="blue", ylab="Value", xlab="Location", main="Comparison",bty="n")
         lines(d3, col="red")
       }
       
     if(output == 1)
      {
        if(!addmean) {return(d)}
        if(addmean) {return(d3)}
       }
      
     if(output == 2) return(data.frame( cbind(freq[1:nyqfreq],filter_pos) ) )
       
#### END function taner
}