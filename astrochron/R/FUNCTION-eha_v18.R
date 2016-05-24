### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### eha function - (SRM: December 7-19, 2012; May 17-26, 2013; June 5-7, 2013;
###                      June 13, 2013; June 26, 2013; July 27, 2013; 
###                      August 5, 2013; Nov. 26, 2013; January 13, 2014;
###                      February 14, 2014; Nov. 10, 2014; Jan. 31, 2015;
###                      February 3, 2015; March 6, 2015; April 8, 2015;
###                      May 24, 2015; Sept. 10, 2015)
###
### uses EHA (FORTRAN) library and built-in functions from R
###########################################################################

eha <- function (dat,tbw=2,pad=defaultPad,fmin=0,fmax=Nyq,step=dt*10,win=dt*100,demean=T,detrend=T,siglevel=0.90,sigID=F,ydir=1,output=0,pl=1,xlab=NULL,ylab=NULL,genplot=2,verbose=T)
{

# uses fields library for plotting, multitaper library to generate dpss tapers

dat <- data.frame(dat)
mpts <- as.integer( length(dat[,1]) )
dt <- dat[2,1]-dat[1,1]

# error checking 
   if(dt<0)
     { 
       if (verbose) cat("\n * Sorting data into increasing height/depth/time, removing empty entries\n")
       dat <- dat[order(dat[1], na.last = NA, decreasing = F), ]
       dt <- dat[2,1]-dat[1,1]
       mpts <- as.integer( length(dat[,1]) )
     }

   y <- dat[,2]

   dtest <- dat[2:mpts,1]-dat[1:(mpts-1),1] 
   epsm=1e-9
   if( (max(dtest)-min(dtest)) > epsm ) 
     {
       cat("\n**** ERROR: sampling interval is not uniform.\n")
       stop("**** TERMINATING NOW!")
     }

   if( tbw > 10 ) 
     {
       cat("\n**** ERROR: maximum time-bandwidth product = 10.\n")
       stop("**** TERMINATING NOW!")
     }

if(verbose)
 {
   cat("\n ----- PERFORMING EVOLUTIVE HARMONIC ANALYSIS -----\n")
   cat(" * Number of data points in stratigraphic series:",mpts,"\n")
   cat(" * Stratigraphic series length (space or time):",(mpts-1)*dt,"\n")
   cat(" * Sampling interval (space or time):",dt,"\n")
  }
  
# number of points per window
winpts=as.integer( floor(win/dt) + 1 )
# number of points to increment window
wininc=as.integer( round(step/dt) )


# error checking 
if(winpts >= mpts) 
  {
     winpts=mpts
     wininc=1
     if(verbose) 
       {
         cat(" * Specified duration longer than permitted for evolutive analysis.\n")
         cat("         Will calculate single MTM spectrum using full data series.\n")
       }
  }   

# number of spectra
nspec= as.integer( floor( (mpts-winpts) / wininc) + 1 )
### calculate Nyquist freq
Nyq <- 1/(2*dt)
### calculate Rayleigh frequency
Ray <- 1/(dt*winpts)

if(verbose)
 {
   cat(" * Number of data points per window:",winpts,"\n")
   cat(" * Moving window size (space or time):",(winpts-1)*dt,"\n")
   cat(" * Window step points:",wininc,"\n")
   cat(" * Window step (space or time):",wininc*dt,"\n")
   cat(" * Number of windows:",nspec,"\n")
   if(demean) cat(" * Mean value for each window will be subtracted\n")
   if(!demean) cat(" * Mean value for each window will NOT be subtracted\n")
   if (detrend) cat(" * Linear trend for each window will be subtracted\n")
   if (!detrend) cat(" * Linear trend for each window will NOT be subtracted\n")
   cat(" * Nyquist frequency:",Nyq,"\n")
   cat(" * Rayleigh frequency:",Ray,"\n")
   cat(" * MTM Power spectrum bandwidth resolution:",tbw/(winpts*dt),"\n")
 }
 
numtap <- trunc((2*tbw)-1)
if (verbose) cat(" * Will use",numtap,"DPSS tapers\n")
# convert to integer for eha function
kmany <- as.integer(numtap)

### number of points to use in paddding is set at function call, 
### either explicitly, or to default value (calculated below).
### pad to 5-10*winpts, if conducting F-test, otherwise, use newpts*2
### using Singleton FFT, npad must not factor into a prime number > 23
### we will ensure compliance by setting default to power of 2
defaultPad= 2*2^ceiling(log2(winpts))
### add another zero if we don't have an even number of data points, so Nyquist exists.   
if((pad)%%2 != 0) pad = pad + 1
newpts=as.integer(pad)
if(verbose) cat(" * Padded to",newpts,"points\n") 

# error checking 
   if( pad > 200000 ) 
     {
       cat("\n**** ERROR: Number of points (post-padding) must be <= 200,000.\n")
       stop("**** TERMINATING NOW!")
     }

   if(pad < winpts)
     {
       cat("\n**** ERROR: Number of points (post-padding) must be >= window points.\n")
       stop("**** TERMINATING NOW!")
     }

# padded frequency grid
df = 1/(newpts*dt)
# determine index of fmin
ifstart = as.integer( round(fmin/df) + 1 )
# determine index of fmax
ifend = as.integer( round(fmax/df) + 1 )
# number of frequencies
nfreq = ifend - ifstart + 1

if(demean) imean=1
if(!demean) imean=2
if(detrend) itrend=1
if(!detrend) itrend=2

# generate the dpss with library(multitaper).  These are not normalized to RMS=1, but will
# be normalized as such in fortran routine.
dpssT=dpss(n=winpts,k=kmany,nw=tbw)
tapers=dpssT$v
evals=dpssT$eigen

ehaF <- function (mpts,winpts,wininc,nspec,nfreq,dt,imean,itrend,ifstart,ifend,y,tbw,kmany,newpts,tapers,evals)
 {
    F_dat = .Fortran('eha_rv6',PACKAGE='astrochron',
                
                mpts=as.integer(mpts),winpts=as.integer(winpts),
                wininc=as.integer(wininc),nspec=as.integer(nspec),nfreq=as.integer(nfreq), 
                dt=as.double(dt),imean=as.integer(imean),itrend=as.integer(itrend),ifstart=as.integer(ifstart),
                ifend=as.integer(ifend),y=as.double(y),tbw=as.double(tbw),kmany=as.integer(kmany),
                newpts=as.integer(newpts),tapers=as.double(tapers),evals=as.double(evals),
                
                f=double(as.integer(nfreq)),amp=matrix(0,nrow=nspec,ncol=nfreq),
                phase=matrix(0,nrow=nspec,ncol=nfreq),
                ftest=matrix(0,nrow=nspec,ncol=nfreq),power=matrix(0,nrow=nspec,ncol=nfreq),
                height=double(as.integer(nspec)),ier=integer(1) 
                )

# return the results
    return(F_dat)
  }

# set up matrices for results
   freq<-double(nfreq)
   height<-double(nspec)
   pwrRaw<-double(nfreq*nspec)
   dim(pwrRaw)<-c(nfreq,nspec)
   FtestRaw<-double(nfreq*nspec)
   dim(FtestRaw)<-c(nfreq,nspec)
   ampRaw<-double(nfreq*nspec)
   dim(ampRaw)<-c(nfreq,nspec)
   prob<-double(nfreq*nspec)
   dim(prob)<-c(nfreq,nspec)

   spec <- ehaF(mpts,winpts,wininc,nspec,nfreq,dt,imean,itrend,ifstart,ifend,y,tbw,kmany,newpts,tapers,evals)

### check for error in FFT. normally err = 0, but is set to 1 if padding is not approporate
   if(spec$ier==1) 
     {
       cat("\n**** ERROR: Number of points (post-padding) must not factor into a prime number greater than 23\n")
       stop("**** TERMINATING NOW!")
     }

### save to arrays (and transpose if needed)
   freq <- spec$f
   pwrRaw <- t(spec$power)
   FtestRaw <- t(spec$ftest)
   ampRaw <- t(spec$amp)
   height <- spec$height + dat[1,1]

### f-test CL
dof=2*numtap
prob <- pf(FtestRaw,2,dof-2)

# plot results
if (is.null(xlab)) xlab = c("Frequency")
if (is.null(ylab)) ylab = c("Location")

if(nspec > 1)
{
 if (genplot == 1)
  {
    par(mfrow=c(2,2))
    if (ydir == -1) 
     {
# in this case, reset ylim range.
# note that useRaster=T is not a viable option, as it will plot the results backwards, even though the
#  y-axis scale has been reversed!  This option will result in a slower plotting time.
       ylimset=c( max(height),min(height) )
       if(pl == 1) image.plot(freq,height,log(pwrRaw),ylim=ylimset,col = tim.colors(100),xlab=xlab,ylab=ylab,main="(a) EPSA: Log Power")
       if(pl == 2) image.plot(freq,height,pwrRaw,ylim=ylimset,col = tim.colors(100),xlab=xlab,ylab=ylab,main="(a) EPSA: Power")
       image.plot(freq,height,ampRaw,ylim=ylimset,col = tim.colors(100),xlab=xlab,ylab=ylab,main="(b) EHA: Amplitude")
       image.plot(freq,height,log10(FtestRaw),ylim=ylimset,col = tim.colors(100),xlab=xlab,ylab=ylab,main="(c) EHA: Log10 F-value")
       image.plot(freq,height,prob,ylim=ylimset,col = tim.colors(100),xlab=xlab,ylab=ylab,main="(d) EHA: Harmonic F-test CL")
     }

    if (ydir == 1) 
     {
# useRaster=T results in a faster plotting time.
       if(pl == 1) image.plot(freq,height,log(pwrRaw),col = tim.colors(100),useRaster=T,xlab=xlab,ylab=ylab,main="(a) EPSA: Log Power",cex.lab=1.1)
       if(pl == 2) image.plot(freq,height,pwrRaw,col = tim.colors(100),useRaster=T,xlab=xlab,ylab=ylab,main="(a) EPSA: Power",cex.lab=1.1)
       image.plot(freq,height,ampRaw,col = tim.colors(100),useRaster=T,xlab=xlab,ylab=ylab,main="EHA: Amplitude")
       image.plot(freq,height,log10(FtestRaw),col = tim.colors(100),useRaster=T,xlab=xlab,ylab=ylab,main="EHA: Log10 F-value")
       image.plot(freq,height,prob,col = tim.colors(100),useRaster=T,xlab=xlab,ylab=ylab,main="EHA: Harmonic F-test CL")
     }
# end genplot = 1
   }

if (genplot == 2)
  {
    if (ydir == -1) 
     {
       dev.new(title=paste("Data series"),height=6.8,width=2)
       plot(dat[,2],dat[,1], ylim = c( max(dat[,1]),min(dat[,1]) ),type="l",xlab="Value",ylab=ylab,main="Data Series",cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n")
       dev.new(title=paste("Time-frequency results"),height=6,width=10)
       par(mfrow=c(1,3))
# mar is defined as: bottom, left, top, right
       par(mar = c(3.1, 4.1, 5.1, 0.7))
# mgp is dfined as: axis title, axis labels and axis line
       par(mgp = c(2.2,1,0))
# in this case, reset ylim range.
# note that useRaster=T is not a viable option, as it will plot the results backwards, even though the
#  y-axis scale has been reversed!  This option will result in a slower plotting time.
       ylimset=c( max(height),min(height) )
       if(pl == 1) image.plot(freq,height,log(pwrRaw),ylim=ylimset,col = tim.colors(100),xlab=xlab,ylab=ylab,main="(a) EPSA: Log Power",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
       if(pl == 2) image.plot(freq,height,pwrRaw,ylim=ylimset,col = tim.colors(100),xlab=xlab,ylab=ylab,main="(a) EPSA: Power",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
       image.plot(freq,height,ampRaw,ylim=ylimset,col = tim.colors(100),xlab=xlab,ylab="",main="(b) EHA: Amplitude",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
       image.plot(freq,height,prob,ylim=ylimset,col = tim.colors(100),xlab=xlab,ylab="",main="(c) EHA: Harmonic F-test CL",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
     }

    if (ydir == 1) 
     {
       dev.new(title=paste("Data series"),height=6.8,width=2)
       plot(dat[,2],dat[,1],type="l",xlab="Value",ylab="Location",main="Data Series",cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n")
       dev.new(title=paste("Time-frequency results"),height=6,width=10)
       par(mfrow=c(1,3))
       par(mar = c(3.1, 4.1, 5.1, 0.7))
       par(mgp = c(2.2,1,0))
# useRaster=T results in a faster plotting time.
       if(pl == 1) image.plot(freq,height,log(pwrRaw),col = tim.colors(100),useRaster=T,xlab=xlab,ylab=ylab,main="(a) EPSA: Log Power",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
       if(pl == 2) image.plot(freq,height,pwrRaw,col = tim.colors(100),useRaster=T,xlab=xlab,ylab=ylab,main="(a) EPSA: Linear Power",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
       image.plot(freq,height,ampRaw,col = tim.colors(100),useRaster=T,xlab=xlab,ylab=ylab,main="(b) EHA: Amplitude",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
       image.plot(freq,height,prob,col = tim.colors(100),useRaster=T,xlab=xlab,ylab=ylab,main="(c) EHA: Harmonic F-test CL",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
     }
# end genplot = 2
   }

if (genplot == 3)
  {
# normalize amplitude spectra to have maxima of unity
    ampNorm<-double(nfreq*nspec)
    dim(ampNorm)<-c(nfreq,nspec)
    ampNorm=t(ampRaw)/(apply(ampRaw,2,max))
    ampNorm=t(ampNorm)

# filter at given significance level
    ampSig<-ampNorm
    ampSig[prob<siglevel] <- NA

    if (ydir == -1) 
     {
       dev.new(title=paste("Data series"),height=6.8,width=2)
       plot(dat[,2],dat[,1], ylim = c( max(dat[,1]),min(dat[,1]) ),type="l",xlab="Value",ylab=ylab,main="Data Series",cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n")
       dev.new(title=paste("Time-frequency results"),height=6,width=10)
       par(mfrow=c(1,3))
       par(mar = c(3.1, 4.1, 5.1, 0.7))
       par(mgp = c(2.2,1,0))
# in this case, reset ylim range.
# note that useRaster=T is not a viable option, as it will plot the results backwards, even though the
#  y-axis scale has been reversed!  This option will result in a slower plotting time.
       ylimset=c( max(height),min(height) )
       if(pl == 1) image.plot(freq,height,log(pwrRaw),ylim=ylimset,col = tim.colors(100),xlab=xlab,ylab=ylab,main="(a) EPSA: Log Power",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
       if(pl == 2) image.plot(freq,height,pwrRaw,ylim=ylimset,col = tim.colors(100),xlab=xlab,ylab=ylab,main="(a) EPSA: Power",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
       image.plot(freq,height,ampNorm,ylim=ylimset,col = tim.colors(100),xlab=xlab,ylab=ylab,main="(b) EHA: Normalized Amplitude",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
       image.plot(freq,height,ampSig,ylim=ylimset,col=tim.colors(100),xlab=xlab,ylab=ylab,main="(c) EHA: Normalized & Filtered Amplitude",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
     }

    if (ydir == 1) 
     {
       dev.new(title=paste("Data series"),height=6.8,width=2)
       plot(dat[,2],dat[,1],type="l",xlab="Value",ylab=ylab,main="Data Series",cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n")
       dev.new(title=paste("Time-frequency results"),height=6,width=10)
       par(mfrow=c(1,3))
       par(mar = c(3.1, 4.1, 5.1, 0.7))
       par(mgp = c(2.2,1,0))
# useRaster=T results in a faster plotting time.
       if(pl == 1) image.plot(freq,height,log(pwrRaw),col = tim.colors(100),useRaster=T,xlab=xlab,ylab=ylab,main="(a) EPSA: Log Power",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
       if(pl == 2) image.plot(freq,height,pwrRaw,col = tim.colors(100),useRaster=T,xlab=xlab,ylab=ylab,main="(a) EPSA: Power",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
       image.plot(freq,height,ampNorm,col = tim.colors(100),useRaster=T,xlab=xlab,ylab=ylab,main="(b) EHA: Normalized Amplitude",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
       image.plot(freq,height,ampSig,col=tim.colors(100),useRaster=T,xlab=xlab,ylab=ylab,main="(c) EHA: Normalized & Filtered Amplitude",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
     }
# end genplot = 3
   }  
   
   if (genplot == 4)
  {
# normalize amplitude spectra to have maxima of unity
    ampNorm<-double(nfreq*nspec)
    dim(ampNorm)<-c(nfreq,nspec)
    ampNorm=t(ampRaw)/(apply(ampRaw,2,max))
    ampNorm=t(ampNorm)

# normalize power spectra to have maxima of unity
    pwrNorm<-double(nfreq*nspec)
    dim(pwrNorm)<-c(nfreq,nspec)
    pwrNorm=t(pwrRaw)/(apply(pwrRaw,2,max))
    pwrNorm=t(pwrNorm)

# filter at given significance level
    ampSig<-ampNorm
    ampSig[prob<siglevel] <- NA

    if (ydir == -1) 
     {
       dev.new(title=paste("Data series"),height=6.8,width=2)
       plot(dat[,2],dat[,1], ylim = c( max(dat[,1]),min(dat[,1]) ),type="l",xlab="Value",ylab=ylab,main="Data Series",cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n")
       dev.new(title=paste("Time-frequency results"),height=6,width=10)
       par(mfrow=c(1,3))
       par(mar = c(3.1, 4.1, 5.1, 0.7))
       par(mgp = c(2.2,1,0))
# in this case, reset ylim range.
# note that useRaster=T is not a viable option, as it will plot the results backwards, even though the
#  y-axis scale has been reversed!  This option will result in a slower plotting time.
       ylimset=c( max(height),min(height) )
       if(pl == 1) image.plot(freq,height,log(pwrNorm),ylim=ylimset,col = tim.colors(100),xlab=xlab,ylab=ylab,main="(a) EPSA: Log Normalized Power (unity)",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
       if(pl == 2) image.plot(freq,height,pwrNorm,ylim=ylimset,col = tim.colors(100),xlab=xlab,ylab=ylab,main="(a) EPSA: Normalized Power (unity)",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
       image.plot(freq,height,ampNorm,ylim=ylimset,col = tim.colors(100),xlab=xlab,ylab=ylab,main="(b) EHA: Normalized Amplitude (unity)",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
       image.plot(freq,height,ampSig,ylim=ylimset,col=tim.colors(100),xlab=xlab,ylab=ylab,main="(c) EHA: Normalized & Filtered Amplitude",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
     }

    if (ydir == 1) 
     {
       dev.new(title=paste("Data series"),height=6.8,width=2)
       plot(dat[,2],dat[,1],type="l",xlab="Value",ylab=ylab,main="Data Series",cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n")
       dev.new(title=paste("Time-frequency results"),height=6,width=10)
       par(mfrow=c(1,3))
       par(mar = c(3.1, 4.1, 5.1, 0.7))
       par(mgp = c(2.2,1,0))
# useRaster=T results in a faster plotting time.
       if(pl == 1) image.plot(freq,height,log(pwrNorm),col = tim.colors(100),useRaster=T,xlab=xlab,ylab=ylab,main="(a) EPSA: Log Normalized Power (unity)",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
       if(pl == 2) image.plot(freq,height,pwrNorm,col = tim.colors(100),useRaster=T,xlab=xlab,ylab=ylab,main="(a) EPSA: Normalized Power (unity)",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
       image.plot(freq,height,ampNorm,col = tim.colors(100),useRaster=T,xlab=xlab,ylab=ylab,main="(b) EHA: Normalized Amplitude (unity)",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
       image.plot(freq,height,ampSig,col=tim.colors(100),useRaster=T,xlab=xlab,ylab=ylab,main="(c) EHA: Normalized & Filtered Amplitude",cex.lab=1.2,horizontal=T,legend.shrink=0.7,legend.mar=2)
     }
# end genplot = 4
   }  
   
# end nspec > 1
}

if(nspec == 1)
{

if(verbose) 
  {
    cat("\n * Searching for significant Harmonic F-test peaks\n")
    cat("     that achive",siglevel*100,"% confidence level:\n") 
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

if(verbose) cat("\nID  / Frequency / Period / Harmonic_CL\n") 
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

if (plateau > 0 && verbose) cat(" * Number of plateaus detected=",plateau,"\n")
 

### generate plots
if(genplot)
  {
    par(mfrow=c(3,1))
# POWER SPECTRUM
    if(pl == 1) plot(freq,log(pwrRaw),type="l", col="black", xlab=xlab,ylab="Log Power",main="MTM Power Estimates",cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n")
    if(pl == 2) plot(freq,pwrRaw,type="l", col="black", xlab=xlab,ylab="Linear Power",main="MTM Power Estimates",cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n")
### plot "significant" frequencies (on power plot first)
    if(sigID)
     {
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
# AMPLITUDE SPECTRUM
   plot(freq,ampRaw,type="l", col="blue", xlab=xlab,ylab="Amplitude",main="MTM Harmonic Amplitude Estimates",cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n")
   if(sigID)
    {
      for (k in 1:(ii-1))
        {
           plfreq[k]=Frequency[k]
           pltext[k]=k
        }
      abline(v=plfreq,col="gray",lty=3)
      mtext(pltext[seq(from=1,to=(ii-1),by=2)], side=3,line=0.25,at=plfreq[seq(from=1,to=(ii-1),by=2)],cex=0.5,font=4)
      if((ii-1) > 1) mtext(pltext[seq(from=2,to=(ii-1),by=2)], side=3,line=-0.25,at=plfreq[seq(from=2,to=(ii-1),by=2)],cex=0.5,font=4)
     }   
# HARMONIC-CL SPECTRUM
   plot(freq,prob*100,type="l",col="red",ylim=c(siglevel*100,100),xlab=xlab,ylab="Confidence Level",main="MTM Harmonic F-test Confidence Level Estimates",cex.axis=1.1,cex.lab=1.1,lwd=2,bty="n")
   abline(h=c(90,95,99),col="black",lty=3)
   if(sigID)
    {          
      for (k in 1:(ii-1))
        {
           plfreq[k]=Frequency[k]
           pltext[k]=k
        }
      abline(v=plfreq,col="gray",lty=3)
      mtext(pltext[seq(from=1,to=(ii-1),by=2)], side=3,line=0.25,at=plfreq[seq(from=1,to=(ii-1),by=2)],cex=0.5,font=4)
      if((ii-1) > 1) mtext(pltext[seq(from=2,to=(ii-1),by=2)], side=3,line=-0.25,at=plfreq[seq(from=2,to=(ii-1),by=2)],cex=0.5,font=4)
     }   
# end genplot
  }
# end nspec = 1
}

# OUTPUT RESULTS
# add column titles to identify each record for output
   colnames(pwrRaw) <- height
   colnames(ampRaw) <- height
   colnames(prob) <- height
   colnames(FtestRaw) <- height
   
if (output==1) 
  {
# add frequency column
   pwrRaw <- data.frame( cbind(freq,pwrRaw[1:nfreq,]) )
   ampRaw <- data.frame( cbind(freq,ampRaw[1:nfreq,]) )
   prob <- data.frame( cbind(freq,prob[1:nfreq,]) )
   FtestRaw <- data.frame( cbind(freq,FtestRaw[1:nfreq,]) )
   return(list(pwrRaw,ampRaw,prob,FtestRaw))
  }

if (output==2) return( data.frame( cbind(freq,pwrRaw[1:nfreq,]) ) )
if (output==3) return( data.frame( cbind(freq,ampRaw[1:nfreq,]) ) )
if (output==4) return( data.frame( cbind(freq,prob[1:nfreq,]) ) )

if (output==5 && nspec == 1) 
 {
    sigfreq <- data.frame(Frequency[1:(ii-1)])
    colnames(sigfreq) <- 'Frequency'
    return(sigfreq)
 }

if (output==6 && nspec == 1) 
 {
    sigfreq <- data.frame(Frequency[1:(ii-1)],Harmonic_CL[1:(ii-1)])
    colnames(sigfreq)[1] <- 'Frequency'
    colnames(sigfreq)[2] <- 'Harmonic_CL'
    return(sigfreq)
 }

#### END function eha
}
