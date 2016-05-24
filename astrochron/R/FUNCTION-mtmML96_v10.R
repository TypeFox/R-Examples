### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### mtmML96 function - (SRM: December 22, 2013; July 31, 2014; January 31, 2015;
###                          February 1, 2015; February 5, 2015; February 22, 2015;
###                          February 26, 2015; March 6, 2015; September 10, 2015;
###                          December 14-15, 2015)
###
### Perform Mann and Lees (1996) robust red noise mtm analysis, with some modifications.
### Uses multitaper library and built in functions from R.
###########################################################################

mtmML96 <- function (dat,tbw=3,ntap=NULL,padfac=5,demean=T,detrend=F,medsmooth=0.2,opt=1,linLog=2,siglevel=0.9,output=0,CLpwr=T,xmin=0,xmax=Nyq,sigID=F,pl=1,genplot=T,verbose=T)
{

if(verbose) cat("\n----- PERFORMING Mann and Lees (1996) Robust Red Noise Analysis -----\n")

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

# remove mean and linear trend if requested
if (demean) 
  { 
    dave <- colMeans(dat[2])
    dat[2] <- dat[2] - dave
    if(verbose) cat(" * Mean value subtracted=",dave,"\n")
  }

if (!demean && verbose) cat(" * Mean value NOT subtracted\n")

# use least-squares fit to remove linear trend
if (detrend) 
  {
    lm.0 <- lm(dat[,2] ~ dat[,1])
    dat[2] <- dat[2] - (lm.0$coeff[2]*dat[1] + lm.0$coeff[1])
    if(verbose) cat(" * Linear trend subtracted. m=",lm.0$coeff[2],"b=",lm.0$coeff[1],"\n")
  }

if (!detrend && verbose) cat(" * Linear trend NOT subtracted\n")


# what is the estimated AR1 coefficient for the combined signal?
lag0 <- dat[1:(npts-1),2]
lag1 <- dat[2:npts,2]
rho_raw <- cor(lag0,lag1)
if(verbose) cat(" * Estimated Conventional AR1 coefficient =",rho_raw,"\n")

# calculate Nyquist freq
Nyq <- 1/(2*dt)
# calculate rayleigh frequency
Ray <- 1/(dt*npts)
npad=as.integer(npts*padfac)
# add another zero if we don't have an even number of data points, so Nyquist exists.   
if((npad*padfac)%%2 != 0) npad = npad + 1
# padded frequency grid
df = 1/(npad*dt)

if(verbose)
 {
   cat(" * Nyquist frequency:",Nyq,"\n")
   cat(" * Rayleigh frequency:",Ray,"\n")
   cat(" * MTM Power spectrum bandwidth resolution:",tbw/(npts*dt),"\n")
   cat(" * Padded to",npad,"points\n")
   cat(" * Frequency grid spacing:",df,"\n")
 }

# make dat a time series object, here with unit sampling interval
datTS <- as.ts(dat[,2])
spec <- spec.mtm(datTS,nw=tbw,k=ntap,Ftest=T,nFFT=npad,taper=c("dpss"),centre=c("none"),jackknife=F,returnZeroFreq=F,plot=F)

# assign correct frequencies to spec$freq (note: this variable is returned if output = 4)
spec$freq <- spec$freq/dt
# note: no zero frequency present, also remove Nyquist now
nfreq = length(spec$freq) - 1
freq <- spec$freq[1:nfreq]

# normalize power (divided by npts in spec.mtm)
pwrRaw <- spec$spec[1:nfreq]/npts
FtestRaw <- spec$mtm$Ftest[1:nfreq]


###########################################################################
### Robust AR(1) coefficient estimate by analytic fit to median smoothed
###    spectrum
###########################################################################

### CALCULATE MEDIAN SMOOTHED SPECTRUM USING runmed FUNCTION #############
    if(verbose) cat(" * Calculating median smoothed spectrum using Tukey's robust end-point rule and symmetrical medians\n")
    freqSmooth = Nyq * medsmooth
# number of frequencies per each smoothing window
    nptsSmooth = as.integer( round(freqSmooth/df, digits=1) )
# check to see if nptsSmooth is even
    if(nptsSmooth %% 2 == 0) 
      {
         nptsSmooth=nptsSmooth+1
         if(verbose) cat(" * Median smoothing window increased by 1 point\n")
      }
    if(verbose) cat(" * Number of smoothing points=",nptsSmooth,"\n")
    pwrMedian <- runmed(pwrRaw,nptsSmooth,endrule="median")[1:nfreq]
### determine average power (So)
    So <- mean(pwrMedian)
 
  
###########################################################################
### Analytical fit of AR1 to median smoothed spectrum
###
### Use Brent's method, Gauss-Newton, or brute force grid search approach (if 
### other approaches fail).
###########################################################################

# check to ensure logarthim exists
if(min(pwrMedian) < 0) 
  {
      cat("**** ERROR IN ANALYTIC FIT TO MEDIAN SMOOTHER: Spectrum has negative power\n")
      stop("**** TERMINATING NOW!")
  }

if(opt == 1 || opt == 2)
 {
   if(opt==1) 
     {
# define function for Brent's method (option 1)
       rednoise1 <- function (rho)
         {
          rhospec=So * (1-rho^2)/(1 - (2*rho*cos(pi*freq/Nyq)) + rho^2)
          if(linLog == 1) sum((pwrMedian-rhospec)^2)
          if(linLog == 2) sum((log(pwrMedian)-log(rhospec))^2)
         }
       if(verbose) cat(" * Calculating analytic fit of AR1 to median smoothed spectrum using Brent's method\n")
       MLred=optim(par=rho_raw,rednoise1,method="Brent",lower=0,upper=1)    
       rho <- MLred$par
     }  

   if(opt==2) 
     {
# define function for Gauss-Newton algorithm (option 2)
      rednoise2 <- function (freq,rho,So,Nyq)
         {
          So * (1-rho^2)/(1 - (2*rho*cos(pi*freq/Nyq)) + rho^2)
         }
       if(verbose) cat(" * Calculating analytic fit of AR1 to median smoothed spectrum using Newton-Gauss algorithm\n")
       if(verbose) cat(" *   Gauss-Newton Convergence for AR1 fit:\n")
       if(linLog == 1) MLred <- nls( pwrMedian ~ rednoise2(freq,rho,So,Nyq), data=data.frame(cbind(freq,pwrMedian)), start=list(rho=rho_raw), trace=T )
       if(linLog == 2) MLred <- nls( log(pwrMedian) ~ log(rednoise2(freq,rho,So,Nyq)), data=data.frame(cbind(freq,pwrMedian)), start=list(rho=rho_raw), trace=T )
       rho <- summary(MLred)$coefficient[1]
     }  
  }

# Grid search
if(opt == 3)
  {
     if(verbose) cat(" * Calculating analytic fit of AR1 to median smoothed spectrum using brute force grid search\n")
     min_MSE = 10^36
     for (rho1 in seq(0, 0.999, by=0.001)) 
       { 
# calculate mean square error
         rhospec = So * (1-rho1^2)/(1 - (2*rho1*cos(pi*freq/Nyq)) + rho1^2)
         if(linLog == 1) MSE=sum((pwrMedian-rhospec)^2)
         if(linLog == 2) MSE=sum((log(pwrMedian)-log(rhospec))^2)
         MSE=MSE/nfreq
         if (MSE < min_MSE) 
           {
             min_MSE = MSE
             rho = rho1 
           }
       }
  }

if(verbose) cat(" * Estimated Robust AR1 coefficient =",rho,"\n")


#***************************************************************
###  Calculate Theoretical AR(1) spectrum (EQ's 4 and 5 of Mann and Lees, 1996)
#***************************************************************
### pwrMedAR will be the AR1 power spectrum for the composite signal
pwrMedAR = So * (1-(rho^2)) / (  1 - (2*rho*cos(pi*freq/Nyq)) + (rho^2) )   

# calculate significance levels
dof = (2*ntap)
chiMedAR <-  (pwrRaw/pwrMedAR) * dof
chiCLMedAR <- pchisq(chiMedAR, df=dof)

# 90, 95 and 99% CL for power plot
    if(CLpwr)
     {
       MedAR1_90 <- pwrMedAR*qchisq(0.9, df=dof)/dof
       MedAR1_95 <- pwrMedAR*qchisq(0.95, df=dof)/dof
       MedAR1_99 <- pwrMedAR*qchisq(0.99, df=dof)/dof
    }

### f-test CL
dof=2*ntap
prob <- pf(FtestRaw,2,dof-2)

### generate plots
if(genplot)
 {
   par(mfrow=c(3,1))
   if(pl == 1)
    {
      plot(freq,log(pwrRaw),type="l", col="black", xlim=c(xmin,xmax), xlab="Frequency",ylab="Log Power",main="MTM Power (black), Robust AR1 fit (red), smoothed (green)",cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n")
      lines(freq,log(pwrMedAR),xlim=c(xmin,xmax),col="red",lwd=2)
      lines(freq,log(pwrMedian),xlim=c(xmin,xmax),col="seagreen",lwd=2,lty=3)
      if(CLpwr) 
        {
          lines(freq,log(MedAR1_90),xlim=c(xmin,xmax),col="red",lwd=1,lty=3)
          lines(freq,log(MedAR1_95),xlim=c(xmin,xmax),col="red",lwd=1,lty=3)
          lines(freq,log(MedAR1_99),xlim=c(xmin,xmax),col="red",lwd=1,lty=3)
        }
    }
   if(pl == 2)
    {
      plot(freq,pwrRaw,type="l", col="black", xlim=c(xmin,xmax), xlab="Frequency",ylab="Linear Power",main="MTM Power (black), Robust AR1 fit (red), smoothed (green)",cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n")
      lines(freq,pwrMedAR,xlim=c(xmin,xmax),col="red",lwd=2)
      lines(freq,pwrMedian,xlim=c(xmin,xmax),col="seagreen",lwd=2,lty=3)      
      if(CLpwr) 
        {
          lines(freq,MedAR1_90,xlim=c(xmin,xmax),col="red",lwd=1,lty=3)
          lines(freq,MedAR1_95,xlim=c(xmin,xmax),col="red",lwd=1,lty=3)
          lines(freq,MedAR1_99,xlim=c(xmin,xmax),col="red",lwd=1,lty=3)
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

  plot(freq,chiCLMedAR*100,type="l",col="red",xlim=c(xmin,xmax),ylim=c(0,100),cex.axis=1.1,cex.lab=1.1,lwd=2,xlab="Frequency",ylab="Confidence Level",main="Robust AR1 Confidence Level Estimates",bty="n")
  abline(h=c(90,95,99),col="black",lty=3)

if(sigID && (ii-1) > 0)
 {
### plot "significant" F-test frequencies (on red noise CL plot)
      for (k in 1:(ii-1))
        {
           plfreq[k]=Frequency[k]
           pltext[k]=k
        }
      abline(v=plfreq,col="gray",lty=3)
      mtext(pltext[seq(from=1,to=(ii-1),by=2)], side=3,line=0.25,at=plfreq[seq(from=1,to=(ii-1),by=2)],cex=0.5,font=4)
      if((ii-1) > 1) mtext(pltext[seq(from=2,to=(ii-1),by=2)], side=3,line=-0.25,at=plfreq[seq(from=2,to=(ii-1),by=2)],cex=0.5,font=4)
 }

     plot(freq,prob*100,type="l",col="red",xlim=c(xmin,xmax),ylim=c(80,100),cex.axis=1.1,cex.lab=1.1,xlab="Frequency",ylab="Confidence Level",main="MTM Harmonic F-test Confidence Level Estimates",bty="n",lwd=2)
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

if (output==1) 
 {
       spectrum <- data.frame(cbind(freq,pwrRaw,prob*100,chiCLMedAR*100,pwrMedAR,pwrMedian))
       colnames(spectrum)[1] <- 'Frequency'
       colnames(spectrum)[2] <- 'Power'
       colnames(spectrum)[3] <- 'Harmonic_CL'
       colnames(spectrum)[4] <- 'AR1_CL'
       colnames(spectrum)[5] <- 'AR1_fit'
       colnames(spectrum)[6] <- 'Median_Smoothed_Power'
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

#if (output==4) 
# {
#   return(spec)
# }

#### END function mtmML96
}
