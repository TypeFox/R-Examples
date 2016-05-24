### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### mtmAR function - (SRM: October 14, 2012; May 24, 2013; June 13, 2013;
###                        June 21, 2013; Nov. 26, 2013; July 30-31, 2014;
###                        February 3, 2015; March 6, 2015; September 10, 2015)
###
### compute 'intermediate spectrum' of Thomson (2001): ratio of mtm to AR spectrum
### uses multitaper library and built in functions from R
###########################################################################

### NOTE: AR fit and MTM include f(0) and Nyquist.

mtmAR <- function (dat,tbw=3,ntap=NULL,order=1,method="mle",CItype=1,padfac=5,demean=T,detrend=F,output=1,xmin=0,xmax=Nyq,pl=1,genplot=T,verbose=T)
{

if(verbose) cat("\n----- DETERMINING INTERMEDIATE SPECTRUM AS IN THOMSON ET AL. (2001) -----\n")

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
npad=npts*padfac
### add another zero if we don't have an even number of data points, so Nyquist exists.   
if((npad*padfac)%%2 != 0) npad = npad + 1
# padded frequency grid
df = 1/(npad*dt)

if(verbose)
 {
   cat(" * Nyquist frequency:",Nyq,"\n")
   cat(" * Rayleigh frequency:",Ray,"\n")
   cat(" * MTM Power spectrum bandwidth resolution:",tbw/(npts*dt),"\n")
 }


# make dat a time series object, here with unit sampling interval
datTS<-as.ts(dat[,2])
spec <- spec.mtm(datTS,nw=tbw,k=ntap,Ftest=F,nFFT=npad,taper=c("dpss"),centre=c("none"),jackknife=F,returnZeroFreq=T,plot=F)

### save frequency and power to freq and pwr
freq <- spec$freq/dt
# normalize power (divided by npts in spec.mtm)
pwrRaw <- spec$spec/npts
FtestRaw <- spec$mtm$Ftest

nfreq = length(freq)

ARspec=spec.ar(dat[,2],method=method,order=order,n.freq=nfreq,plot=F)
# normalize power
pwrAR=ARspec$spec/npts

# ensure equivalent normalization
sumMTM=sum(pwrRaw)
sumAR=sum(pwrAR)
if(verbose) cat(" * Total Power: MTM=",sumMTM,"AR=",sumAR,"\n")
MTMnorm=pwrRaw/sumMTM
ARnorm=pwrAR/sumAR

mtmAR1=MTMnorm/ARnorm

dof = (2*ntap)
# generate confidence intervals for plots
if(CItype==1) 
 {
   CL90 <- dof/qchisq(0.2, df=dof)
   CL95 <- dof/qchisq(0.1, df=dof)
   CL995 <- dof/qchisq(0.01, df=dof)
 }

if(CItype==2) 
 {
   CL95 <- dof/qchisq(0.05, df=dof)
   CL05 <- dof/qchisq(0.95, df=dof)
   CL995 <- dof/qchisq(0.005, df=dof)
   CL005 <- dof/qchisq(0.995, df=dof)
 }

# generate confidence levels 
  mtmAR1dof=mtmAR1*dof
  chiCLAR <- pchisq(mtmAR1dof, df=dof)*100
 
 
### generate plots
if(genplot)
 {
  par(mfrow=c(3,1))
  if (pl==1)
   {
    plot(freq,log(pwrRaw),type="l", col="blue", xlim=c(xmin,xmax), xlab="Frequency",ylab="Log Power",main="MTM Power Spectrum",cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n")
    plot(freq,log(pwrAR),type="l", col="red", xlim=c(xmin,xmax), xlab="Frequency",ylab="Log Power",main="AR Power Spectrum",cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n")
    plot(freq,log(mtmAR1),type="l", col="black", xlim=c(xmin,xmax), xlab="Frequency",ylab="MTM/AR Spectrum",main="MTM/AR Power Estimates",cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n")
    abline(h=0,col="red",lty=3)
    abline(h=log(CL95),col="red",lty=3)
    abline(h=log(CL995),col="blue",lty=3)
    mtext("95%", side=4,line=-1,at=(log(CL95)),cex=0.5,font=4,las=1)
    mtext("99.5%", side=4,line=-1,at=(log(CL995)),cex=0.5,font=4,las=1)
    if(CItype==1) 
      {
       abline(h=log(CL90),col="red",lty=3)
       mtext("90%", side=4,line=-1,at=(log(CL90)),cex=0.5,font=4,las=1)
      }
    if(CItype==2) 
     {
      abline(h=log(CL05),col="blue",lty=3)
      abline(h=log(CL005),col="blue",lty=3)
      mtext("5%", side=4,line=-1,at=(log(CL05)),cex=0.5,font=4,las=1)
      mtext("0.005%", side=4,line=-1,at=(log(CL005)),cex=0.5,font=4,las=1)
     }
   }
 if (pl==2)
  {
   plot(freq,pwrRaw,type="l", col="blue", xlim=c(xmin,xmax), xlab="Frequency",ylab="Linear Power",main="MTM Power Spectrum",cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n")
   plot(freq,pwrAR,type="l", col="red", xlim=c(xmin,xmax), xlab="Frequency",ylab="Linear Power",main="AR Power Spectrum",cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n")
   plot(freq,pwrRaw/pwrAR,type="l", col="black", xlim=c(xmin,xmax), xlab="Frequency",ylab="MTM/AR Power Spectrum",main="MTM/AR Power Estimates",cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n")
   abline(h=1,col="red",lty=3)
   abline(h=CL95,col="red",lty=3)
   abline(h=CL995,col="blue",lty=3)
   mtext("95%", side=4,line=-1,at=(CL95),cex=0.5,font=4,las=1)
   mtext("99.5%", side=4,line=-1,at=(CL995),cex=0.5,font=4,las=1)
   if(CItype==1) 
    {
      abline(h=CL90,col="red",lty=3)
      mtext("90%", side=4,line=-1,at=(CL90),cex=0.5,font=4,las=1)
    }
   if(CItype==2) 
    {
      abline(h=CL05,col="blue",lty=3)
      abline(h=CL005,col="blue",lty=3)
      mtext("5%", side=4,line=-1,at=(CL05),cex=0.5,font=4,las=1)
      mtext("0.005%", side=4,line=-1,at=(CL005),cex=0.5,font=4,las=1)
    }
  }
 }

if (output==1) 
{
  spectrum <- data.frame(cbind(freq,mtmAR1,chiCLAR))
  colnames(spectrum)[1] <- 'Frequency'
  colnames(spectrum)[2] <- 'MTM/AR'
  colnames(spectrum)[3] <- 'AR_CL'
  return(spectrum)
}

if (output==2) 
{
  spectrum <- data.frame(cbind(freq,mtmAR1))
  colnames(spectrum)[1] <- 'Frequency'
  colnames(spectrum)[2] <- 'MTM/AR'
  return(spectrum)
}

if (output==3) 
{
  spectrum <- data.frame(cbind(freq,chiCLAR))
  colnames(spectrum)[1] <- 'Frequency'
  colnames(spectrum)[2] <- 'AR_CL'
  return(spectrum)
}

#### END function mtmAR
}
