### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### function ar1etp - (SRM: March 21, 2012; May 3, 2012; Aug 2, 2012; 
###                         Oct. 8, 2012; Jan. 23, 2013; April 29, 2013; 
###                         May 20, 2013; June 21, 2013; Sept. 25, 2013; 
###                         Jan. 30, 2015; February 1, 2015)
###
### Run an ensemble of AR1 + signal simulations
###########################################################################

ar1etp <- function(etpdat=NULL,nsim=100,rho=0.9,wtAR=1,sig=90,tbw=2,padfac=5,ftest=F,fmax=0.1,speed=0.5,pl=2,graphfile=0)
{

if(is.null(etpdat)) 
  {
    cat("\n * No astronomical series input. Will use default ETP from:\n")
    cat("    etp(tmin=0,tmax=1000,dt=5,genplot=F)\n")
    etpdat=etp(tmin=0,tmax=1000,dt=5,genplot=F)
  }  

npts <- length(etpdat[,1]) 
dt <- etpdat[2,1]-etpdat[1,1]

### center and standardize etpdat
etpdat[2]=etpdat[2]-colMeans(etpdat[2])
etpdat[2]=etpdat[2]/sapply(etpdat[2],sd)

nran <- npts

cat("\n----- PERFORMING AR1-ETP SIMULATIONS -----\n")
cat("             press <esc> to stop    \n")

#### start simulation loop
for (ii in 1:nsim)
  {

## output iteration
cat(" * SIMULATION:",ii,"/",nsim,"\n")

### turn on PDF driver
if(graphfile==1) 
  {
     myfile <-paste("ar1etp",ii,".pdf")
     pdf(file=myfile,width=6,height=8)
   }  
### turn on JPG driver
if(graphfile==2) 
  {
     myfile <-paste("ar1etp",ii,".jpg")
     jpeg(filename=myfile,width=6,height=8,units="in",res=200)
   }

par(mfrow=c(2,1))


###########################################################################
### 1. Make AR1 noise 
###########################################################################

### Generate normal deviates, mean = 0, sd = 1
### this is "white noise"
white <- rnorm(nran)
### generate AR(1) red noise
ar1 <- 1:(length(white))
### initialize ar1[1]
ar1[1] <- white[1]
for (i in 2:(length(white)))  
### multiply previous value by coeff, then add innovation
    { ar1[i] <- rho*ar1[i-1]+white[i]}
### generate time axis
ta <- 1:(length(white))
### change time axis from unit spacing of 1 to desired value
ta <- (ta*dt) - dt

### center and standardize ar1
ar1=ar1-mean(ar1)
ar1=ar1/sd(ar1)
### set noise level relative to ETP (ETP already standardized to 1)
ar1=ar1*wtAR

noise <- as.data.frame(cbind(ta,ar1))

### what is the estimated AR1 coefficient?
    lag0 <- ar1[1:(npts-1)]
    lag1 <- ar1[2:npts]

rho_raw <- cor(lag0,lag1)

###########################################################################
### 2. sum ETP and AR1 noise
###########################################################################

signal <- as.data.frame(cbind(ta,etpdat[,2] + ar1))

### what is the estimated AR1 coefficient for the combined signal?
    lag0 <- signal[1:(npts-1),2]
    lag1 <- signal[2:npts,2]

rho_raw_sig <- cor(lag0,lag1)

###########################################################################
### 3. MTM Power spectrum using 'multitaper' library
###########################################################################

### calculate Nyquist freq
Nyq <- 1/(2*dt)
### calculate rayleigh frequency
Ray <- 1/(dt*npts)

numtap <- (2*tbw)-1
nf = as.integer(Nyq/Ray)
crit <- 100 * (1.0 - ( (1-(sig*0.01)) / nf) )

npad=npts*padfac
### add another zero if we don't have an even number of data points, so Nyquist exists.   
if((npad*padfac)%%2 != 0) npad = npad + 1
# padded frequency grid
df = 1/(npad*dt)

# make signal a time series object, here with unit sampling interval
signalTS<-as.ts(signal[,2])
spec <- spec.mtm(signalTS,nw=tbw,k=numtap,Ftest=T,nFFT=npad,jackknife=F,returnZeroFreq=F,plot=F)

# note: no zero frequency present, also remove Nyquist now
nfreq = length(spec$freq) - 1
freq <- spec$freq[1:nfreq]/dt
pwrRaw <- spec$spec[1:nfreq]

#***************************************************************
### Calculate Raw red noise spectrum
### "So" is the average power. This can be determined from the white noise variance
###  as So = var/(1-rho^2), where rho is the lag-1 coeff
###  We will determine average power directly from measured spectrum
So = mean(pwrRaw)

RawAR = So * (1-(rho_raw_sig^2)) / (  1 - (2*rho_raw_sig*cos(pi*freq/Nyq)) + (rho_raw_sig^2) )
#***************************************************************

if(pl == 1) 
 {
   plot(freq,log(pwrRaw),type="l",xlab="Frequency (cycles/ka)",ylab="Log(Power)",xlim=c(0,fmax),cex.axis=1.2,cex.lab=1.2,lwd=2)   
   lines(freq,log(RawAR),col="red",lwd=2)
  }
if(pl == 2)
 {
   plot(freq,pwrRaw,type="l",xlab="Frequency (cycles/ka)",ylab="Power",xlim=c(0,fmax),cex.axis=1.2,cex.lab=1.2,lwd=2)   
   lines(freq,RawAR,col="red",lwd=2)
 }
 
abline(v=c(.0024752475,.0080000000,.0105263157,.0185185185,.0243902439,0.04228,0.04474,0.052780),col="blue",lty=3,lwd=2)
mtext(c("e1","e2","e3","o1","o2","p1","p2","p3"),side=3,line=0,at=c(.0024752475,0.0076,.0115263157,.0185185185,.0243902439,.0416666666,.0454545454,.0526315789))

mtext(c(expression(paste(rho,"-noise=")),round(rho_raw,digits=2)),side=1,line=5,col="red",at=c(0.015,0.03),cex=1.2)
mtext(c(expression(paste(rho,"-noise+signal=")),round(rho_raw_sig,digits=2)),side=1,line=5,col="red",at=c(0.065,0.087),cex=1.2)
  
###########################################################################
### 4. Evaluate confidence levels
###########################################################################

dof = (2*numtap)
chiRawAR <-  (pwrRaw/RawAR) * dof
chiCLRawAR <- pchisq(chiRawAR, df=dof)
plot(freq,chiCLRawAR*100,type="l",col="red",xlab="Frequency (cycles/ka)",ylab="Confidence Level",xlim=c(0,fmax),ylim=c(0,100),cex.axis=1.2,cex.lab=1.2,lwd=2)
abline(h=c(sig,crit),col="black",lty=3)
mtext(c(sig),side=4,line=0,at=c(sig))
abline(v=c(.0024752475,.0080000000,.0105263157,.0185185185,.0243902439,0.04228,0.04474,0.052780),col="blue",lty=3,lwd=2)
mtext(c("e1","e2","e3","o1","o2","p1","p2","p3"),side=3,line=0,at=c(.0024752475,0.0076,.0115263157,.0185185185,.0243902439,.0416666666,.0454545454,.0526315789))

### add f-test CL for raw spectrum
fCLRaw <- pf(spec$mtm$Ftest,2,dof-2)
if (ftest) {lines(freq,fCLRaw*100,type="l",col="green")}


if (nsim == 1)
  {
    return(data.frame(signal))
  }
  
if(graphfile==0) Sys.sleep(speed)

### turn off PDF or JPEG driver  
if(graphfile!=0) dev.off()  

  }

### end function ar1etp
}



