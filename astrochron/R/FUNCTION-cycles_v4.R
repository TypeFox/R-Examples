### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### cycles: make a time series with specified harmonic components and noise 
###    - (SRM: April 23, 2013; May 20, 2013; July 31, 2013; August 1, 2013)
###
###########################################################################

cycles <- function (freqs=NULL,phase=NULL,amp=NULL,start=0,end=499,dt=1,noisevar=0,genplot=T,verbose=T)
{

if(verbose) cat("\n----- GENERATING A TIME SERIES WITH SPECIFIED HARMONIC COMPONENTS AND NOISE -----\n")

# default: create 3-term precession model, from Berger 
if(is.null(freqs))
 {
  if(verbose) cat(" * Generating 3-term Berger precession model by default\n")
  freqs = c(1/23.716,1/22.428,1/18.976)
  phase = c(0.5586799,3.44144,5.440017)
  amp = c(0.0186080,0.0162752,-0.0130066)
 }

if(!is.null(freqs))
 {
   if(is.null(phase)) 
    {
      if(verbose) cat(" * phase not specified, will use 0.\n")
      phase=rep(0,length(freqs)) 
    }
   if(is.null(amp)) 
    {
      if(verbose) cat(" * amplitude not specified, will use 1.\n")
      amp=rep(1,length(freqs)) 
    }
 }
 
if(length(freqs)!=length(amp) || length(amp)!=length(phase))
 {
   if(verbose) cat(" * ERROR: amplitude and phase must be specified for all frequencies\n")
   stop("function halted")
 } 
 
 if(verbose)
   {
     temp=data.frame(cbind(freqs,amp,phase))
     colnames(temp) <- c("Frequency","Amplitude","Phase")
     print(temp)
   }  
 
### generate position axis
    ta <- seq(start,end,by=dt)
    npts=length(ta)
    
### generate cycles
cycle <- function(freqs,phase,amp,n) 
  {
### initialize a matrix (with zeros), which will contain each 'freqs' cycle
###  in a column, with each record of length of 'n' points.
    x <- matrix(0, n, length(freqs))
    for (i in 1:length(freqs)) 
      {
        x[,i] <- amp[i]*sin( (2*pi*freqs[i]) * ta - (phase[i]) )
      }
    return(x)
  }

c <- cycle(freqs,phase,amp,npts)

totsig = rowSums(c) + rnorm(npts,sd=sqrt(noisevar))
dat <- data.frame( cbind(ta,totsig) )
 

if(genplot)
  {  
### plots
  par(mfrow=c(2,2))
  plot(dat,cex=0.5,xlab="Location",ylab="Value",main="Data Series"); lines(dat)
### plot the denisty and the histogram together
  hist(dat[,2],freq=F,xlab="Value",main="Distribution of values"); lines(density(dat[,2], bw="nrd"),col="red")
### boxplot
  boxplot(dat[,2],ylab="Value",main="Boxplot for values")
### Normal probabilty plot (Normal Q-Q Plot)
  qqnorm(dat[,2]); qqline(dat[,2], col="red")
   }
    
    colnames(dat)[1] <- 'time'
    colnames(dat)[2] <- 'value'
    return(dat)

# end function cycles
}