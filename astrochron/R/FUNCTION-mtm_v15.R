### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### MTM function - (SRM: February 28, 2012; March 29, 2012; 
###                      September 18, 2012; Oct 8, 2012; Oct. 12, 2012; Nov. 23, 2012
###                      May 20-21, 2013; May 23, 2013; May 27, 2013; June 5, 2013; 
###                      June 13, 2013; August 5, 2013; August 12, 2013; Nov. 26, 2013;
###                      July 31, 2014; January 31, 2015; February 1-3, 2015; 
###                      February 26, 2015; March 6, 2015; June 30, 2015; Sept. 10, 2015;
###                      December 14, 2015)
###
### uses multitaper library and built in functions from R
###########################################################################

mtm <- function (dat,tbw=3,ntap=NULL,padfac=5,demean=T,detrend=F,siglevel=0.9,ar1=T,output=0,CLpwr=T,xmin=0,xmax=Nyq,pl=1,sigID=F,genplot=T,verbose=T)
{

if(verbose) cat("\n----- PERFORMING Multitaper Spectral Analysis -----\n")

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

if (verbose) 
 {
   cat(" * Number of data points in stratigraphic series:",npts,"\n")
   cat(" * Stratigraphic series length (space or time):",(npts-1)*dt,"\n")
   cat(" * Sampling interval (space or time):",dt,"\n")
 }

numtap=trunc((2*tbw)-1)
if(is.null(ntap)) 
 {
  ntap=numtap
  if (verbose) cat(" * Will use default setting of",ntap,"DPSS tapers\n")
 }

if(ntap > numtap)
 {
  ntap=numtap
  if (verbose) cat("**** WARNING: The number of DPSS tapers specified is too large. ntap reset to default value of",ntap,"\n")
 } 

if(ntap < 2)
 {
  ntap=numtap
  if (verbose) cat("**** WARNING: The number of DPSS tapers specified is too small. ntap reset to default value of",ntap,"\n")
 } 

  
###########################################################################
### MTM Power spectrum using 'multitaper' library
###########################################################################
### this version computes the adaptive multitaper spectrum

# remove mean and linear trend if requested
if (demean) 
  { 
    dave <- colMeans(dat[2])
    dat[2] <- dat[2] - dave
    if(verbose) cat(" * Mean value subtracted=",dave,"\n")
  }

if (!demean && verbose) cat(" * Mean value NOT subtracted\n")

### use least-squares fit to remove linear trend
if (detrend) 
  {
    lm.0 <- lm(dat[,2] ~ dat[,1])
    dat[2] <- dat[2] - (lm.0$coeff[2]*dat[1] + lm.0$coeff[1])
    if(verbose) cat(" * Linear trend subtracted. m=",lm.0$coeff[2],"b=",lm.0$coeff[1],"\n")
  }

    if (!detrend && verbose) cat(" * Linear trend NOT subtracted\n")

### calculate Nyquist freq
Nyq <- 1/(2*dt)
### calculate Rayleigh frequency
Ray <- 1/(dt*npts)
npad=as.integer(npts*padfac)
### add another zero if we don't have an even number of data points, so Nyquist exists.   
if((npad*padfac)%%2 != 0) npad = npad + 1
# padded frequency grid
df = 1/(npad*dt)

if(verbose)
 {
   cat(" * Nyquist frequency:",Nyq,"\n")
   cat(" * Rayleigh frequency:",Ray,"\n")
   cat(" * MTM Power spectrum bandwidth resolution:",tbw/(npts*dt),"\n")
   cat(" * Padded to",npad,"points\n")
 }

# make dat a time series object, here with unit sampling interval
datTS<-as.ts(dat[,2])
spec <- spec.mtm(datTS,nw=tbw,k=ntap,Ftest=T,nFFT=npad,taper=c("dpss"),centre=c("none"),jackknife=F,returnZeroFreq=F,plot=F)

# assign correct frequencies to spec$freq (note: this variable is returned if output = 4)
spec$freq <- spec$freq/dt
# note: no zero frequency present, also remove Nyquist now
nfreq = length(spec$freq) - 1
freq <- spec$freq[1:nfreq]

# normalize power (divided by npts in spec.mtm)
pwrRaw <- spec$spec[1:nfreq]/npts
FtestRaw <- spec$mtm$Ftest[1:nfreq]

# AR(1) noise model spectrum
if(ar1)
 {
### what is the estimated AR1 coefficient?
    lag0 <- dat[1:(npts-1),2]
    lag1 <- dat[2:npts,2]
    rho_raw <- cor(lag0,lag1)
### Calculate Raw red noise spectrum
### "So" is the average power. This can be determined from the white noise variance
###  as So = var/(1-rho^2), where rho is the lag-1 coeff
###  We will determine average power directly from measured spectrum
    So = mean(pwrRaw)
    RawAR = So * (1-(rho_raw^2)) / (  1 - (2*rho_raw*cos(pi*freq/Nyq)) + (rho_raw^2) )
    dofAR = (2*ntap)
    chiRawAR <-  (pwrRaw/RawAR) * dofAR
    chiCLRawAR <- pchisq(chiRawAR, df=dofAR)
### 90, 95 and 99% CL for power plot
    if(CLpwr)
     {
       AR1_90 <- RawAR*qchisq(0.9, df=dofAR)/dofAR
       AR1_95 <- RawAR*qchisq(0.95, df=dofAR)/dofAR
       AR1_99 <- RawAR*qchisq(0.99, df=dofAR)/dofAR
    }
 }

### f-test CL
dof=2*ntap
prob <- pf(FtestRaw,2,dof-2)

### generate plots
if(genplot)
 {
   if(ar1) 
    {
      par(mfrow=c(3,1))
      if(!CLpwr) mtitle=c("MTM Power (black) & AR1 fit (red)")
      if(CLpwr) mtitle=c("MTM Power (black); AR1 fit (red); 90%CL, 95%CL, 99%CL (dotted)")
    }   
   if(!ar1) 
    {
     par(mfrow=c(2,1))
     mtitle=c("MTM Power")
    }
   if(pl == 1)
    {
      plot(freq,log(pwrRaw),type="l", col="black", xlim=c(xmin,xmax), xlab="Frequency",ylab="Log Power",main=mtitle,cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n")
      if(ar1) 
        {
          lines(freq,log(RawAR),xlim=c(xmin,xmax),col="red",lwd=2)
          if(CLpwr) 
            {
              lines(freq,log(AR1_90),xlim=c(xmin,xmax),col="red",lwd=1,lty=3)
              lines(freq,log(AR1_95),xlim=c(xmin,xmax),col="red",lwd=1,lty=3)
              lines(freq,log(AR1_99),xlim=c(xmin,xmax),col="red",lwd=1,lty=3)
             }
         }     
    }
   if(pl == 2)
    {
      plot(freq,pwrRaw,type="l", col="black", xlim=c(xmin,xmax), xlab="Frequency",ylab="Linear Power",main=mtitle,cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n")
      if(ar1) 
        {
          lines(freq,RawAR,xlim=c(xmin,xmax),col="red",lwd=2)
          if(CLpwr) 
            {
              lines(freq,AR1_90,xlim=c(xmin,xmax),col="red",lwd=1,lty=3)
              lines(freq,AR1_95,xlim=c(xmin,xmax),col="red",lwd=1,lty=3)
              lines(freq,AR1_99,xlim=c(xmin,xmax),col="red",lwd=1,lty=3)
             }
        }
     } 
 }
 

if(verbose) 
  {
    cat("\n * Searching for significant Harmonic F-test peaks\n")
    cat("     that achieve",siglevel*100,"% CL:\n") 
  }  
### identify the f-test probability peaks
### initialize variables and arrays
probmax <- 1:nfreq
freqloc <- 1:nfreq
plateau <- 0
ij <- 1 
for ( j in 1:(nfreq-2) )
  {
# Check for plateaus
    if (prob[j] == prob[j+1] || prob[j+1] == prob [j+2] ) 
      {
         plateau=plateau + 1
         if(verbose) cat("**** WARNING: plateau detected at",freq[j+1],"Probability=",prob[j+1],"\n")
      }
    if ( prob[j] < prob[j+1] && prob[j+1] > prob[j+2] )
      {
# save peak
         probmax[ij] <- prob[j+1]
         freqloc[ij] <- j+1
         ij=ij+1
       }  
  }
numpeak <- ij-1

if(verbose) cat("ID  / Frequency / Period / Harmonic_CL\n") 
Frequency<- double(numpeak)
Harmonic_CL<- double(numpeak)
ii = 1
for (j in 1:numpeak)
   {
    if (probmax[j] >= siglevel)
      {
       Frequency[ii] <- freq[freqloc[j]]
       Harmonic_CL[ii] <- probmax[j] 
       if(verbose) cat(ii," ",Frequency[ii], " ", 1/Frequency[ii]," ", Harmonic_CL[ii]*100,"\n")
       ii=ii+1
      }
    }


if (plateau > 0 && verbose) cat("\n * Number of plateaus detected=",plateau,"\n")


# plot confidence level estimates
if(genplot)
 {
   if(sigID && (ii-1) > 0)
    {
### plot "significant" F-test frequencies (on power plot first)
      plfreq=double(ii-1)
      pltext=double(ii-1)
      for (k in 1:(ii-1))
        {
           plfreq[k]=Frequency[k]
           pltext[k]=k
        }
      abline(v=plfreq,col="gray",lty=3)
      mtext(pltext[seq(from=1,to=(ii-1),by=2)], side=3,line=0.25,at=plfreq[seq(from=1,to=(ii-1),by=2)],cex=0.5,font=4)
      if((ii-1) > 1) mtext(pltext[seq(from=2,to=(ii-1),by=2)], side=3,line=-0.25,at=plfreq[seq(from=2,to=(ii-1),by=2)],cex=0.5,font=4)
    }
    
  if(!ar1) 
   {
     plot(freq,prob*100,type="l",col="red",xlim=c(xmin,xmax),ylim=c(80,100),cex.axis=1.1,cex.lab=1.1,xlab="Frequency",ylab="Confidence Level",main="MTM Harmonic F-Test Confidence Level Estimates",bty="n",lwd=2)
     abline(h=c(90,95,99),col="black",lty=3)
     if(sigID && (ii-1) > 0)
      {
### plot "significant" F-test frequencies (on probabilty plot)
         for (k in 1:(ii-1))
          {
           plfreq[k]=Frequency[k]
           pltext[k]=k
          }
        abline(v=plfreq,col="gray",lty=3)
        mtext(pltext[seq(from=1,to=(ii-1),by=2)], side=3,line=0.25,at=plfreq[seq(from=1,to=(ii-1),by=2)],cex=0.5,font=4)
        if((ii-1) > 1) mtext(pltext[seq(from=2,to=(ii-1),by=2)], side=3,line=-0.25,at=plfreq[seq(from=2,to=(ii-1),by=2)],cex=0.5,font=4)
      }
   }  
  if(ar1) 
   { 
     plot(freq,chiCLRawAR*100,type="l",col="red",xlim=c(xmin,xmax),ylim=c(0,100),cex.axis=1.1,cex.lab=1.1,lwd=2,xlab="Frequency",ylab="Confidence Level",main="AR1 Confidence Level Estimates",bty="n")
     abline(h=c(90,95,99),col="black",lty=3)
     if(sigID && (ii-1) > 0)
      {
### plot "significant" F-test frequencies (on probabilty plot)
         for (k in 1:(ii-1))
          {
           plfreq[k]=Frequency[k]
           pltext[k]=k
          }
        abline(v=plfreq,col="gray",lty=3)
        mtext(pltext[seq(from=1,to=(ii-1),by=2)], side=3,line=0.25,at=plfreq[seq(from=1,to=(ii-1),by=2)],cex=0.5,font=4)
        if((ii-1) > 1) mtext(pltext[seq(from=2,to=(ii-1),by=2)], side=3,line=-0.25,at=plfreq[seq(from=2,to=(ii-1),by=2)],cex=0.5,font=4)
      }
     plot(freq,prob*100,type="l",col="red",xlim=c(xmin,xmax),ylim=c(80,100),cex.axis=1.1,cex.lab=1.1,xlab="Frequency",ylab="Confidence Level",main="Harmonic F-Test Confidence Level Estimates",bty="n",lwd=2)
     abline(h=c(90,95,99),col="black",lty=3)
     if(sigID && (ii-1) > 0)
      {
### plot "significant" F-test frequencies (on probabilty plot)
         for (k in 1:(ii-1))
          {
           plfreq[k]=Frequency[k]
           pltext[k]=k
          }
        abline(v=plfreq,col="gray",lty=3)
        mtext(pltext[seq(from=1,to=(ii-1),by=2)], side=3,line=0.25,at=plfreq[seq(from=1,to=(ii-1),by=2)],cex=0.5,font=4)
        if((ii-1) > 1) mtext(pltext[seq(from=2,to=(ii-1),by=2)], side=3,line=-0.25,at=plfreq[seq(from=2,to=(ii-1),by=2)],cex=0.5,font=4)
      }
   } 
 }
  

if (output==1) 
 {
   if(!ar1) 
     { 
       spectrum <- data.frame(cbind(freq,pwrRaw,prob*100))
       colnames(spectrum)[1] <- 'Frequency'
       colnames(spectrum)[2] <- 'Power'
       colnames(spectrum)[3] <- 'Harmonic_CL'
      }
   if(ar1)
     {
       spectrum <- data.frame(cbind(freq,pwrRaw,prob*100,chiCLRawAR*100,RawAR))
       colnames(spectrum)[1] <- 'Frequency'
       colnames(spectrum)[2] <- 'Power'
       colnames(spectrum)[3] <- 'Harmonic_CL'
       colnames(spectrum)[4] <- 'AR1_CL'
       colnames(spectrum)[5] <- 'AR1_fit'
     }  
   return(spectrum)
 }

if (output==2) 
 {
    sigfreq <- data.frame(Frequency[1:(ii-1)])
    colnames(sigfreq) <- 'Frequency'
    return(sigfreq)
 }

if (output==3) 
 {
    sigfreq <- data.frame(Frequency[1:(ii-1)],Harmonic_CL[1:(ii-1)])
    colnames(sigfreq)[1] <- 'Frequency'
    colnames(sigfreq)[2] <- 'Harmonic_CL'
    return(sigfreq)
 }

if (output==4) 
 {
   return(spec)
 }

#### END function mtm
}
