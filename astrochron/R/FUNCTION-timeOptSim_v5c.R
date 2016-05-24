### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### function timeOptSim - (SRM: May 28, 2012; Oct. 14, 2014; Oct. 17, 2014;
###                             Oct. 19, 2014; Jan. 13, 2015; March 9, 2015
###                             June 8, 2015; Sept. 30, 2015; 
###                             October 20-21, 2015; November 19, 2015;
###                             December 17, 2015)
###########################################################################

timeOptSim <- function (dat,sedrate=NULL,numsim=1000,rho=NULL,fit=1,flow=NULL,fhigh=NULL,roll=NULL,targetE=NULL,targetP=NULL,output=0,genplot=T,verbose=T)
{

if(verbose) cat("\n----- TimeOpt Monte Carlo Simulation -----\n")

cormethod=1

# prepare data array
   dat = data.frame(dat)      
   npts <- length(dat[,1]) 
   dx <- dat[2,1]-dat[1,1]

# error checking 
   if(dx<0)
     { 
       if (verbose) cat("\n * Sorting data into increasing height/depth/time, removing empty entries\n")
       dat <- dat[order(dat[1], na.last = NA, decreasing = F), ]
       dx <- dat[2,1]-dat[1,1]
       npts <- length(dat[,1])
     }
   dtest <- dat[2:npts,1]-dat[1:(npts-1),1] 
   epsm=1e-10
   if( (max(dtest)-min(dtest)) > epsm ) 
     {
       cat("\n**** ERROR: sampling interval is not uniform.\n")
       stop("**** TERMINATING NOW!")
     }

if (verbose) 
 {
   cat(" * Number of data points in stratigraphic series:",npts,"\n")
   cat(" * Stratigraphic series length (meters):",(npts-1)*dx,"\n")
   cat(" * Sampling interval (meters):",dx,"\n\n")
 }

# standardize data series
   dat[2]=dat[2]-colMeans(dat[2])
   dat[2]=dat[2]/sapply(dat[2],sd)

# convert sedrate from cm/ka to m/ka for processing
   sedrate=sedrate/100

### what is the estimated AR1 coefficient?
   if(is.null(rho))
    {
      lag0 <- dat[1:(npts-1),2]
      lag1 <- dat[2:npts,2]
      rho <- cor(lag0,lag1)
      if(verbose) cat(" * Raw AR1 =",rho,"\n")
    }  

#######################################################################################
# set up default bandpass frequencies and targets
#  first for precession
if(fit == 1)
 {
   if(is.null(flow)) 
    {
      flow = 0.035
      if(verbose) cat(" * Using default flow =",flow,"\n")
    }  
   if(is.null(fhigh)) 
    {
      fhigh = 0.065
      if(verbose) cat(" * Using default fhigh =",fhigh,"\n")
    } 
   if(is.null(roll))
    {
      roll = 10^3
      if(verbose) cat(" * Using default roll =",roll,"\n")
     }

   if(is.null(targetP))
    {
# the four dominant precession peaks, based on spectral analysis of 
#   Laskar et al. (2004), 0-10 Ma
      targetP <- double(4)
      targetP[1] = 23.62069
      targetP[2] = 22.31868
      targetP[3] = 19.06768 
      targetP[4] = 18.91979 
      if(verbose) cat(" * Using default precession target periods =",targetP,"\n")
    }
  }

# next for short eccentricity
if(fit == 2)
 {
   if(is.null(flow))
    {
      flow = 0.007
      if(verbose) cat(" * Using default flow =",flow,"\n")
    }  
   if(is.null(fhigh))
    {
      fhigh = 0.011
      if(verbose) cat(" * Using default fhigh =",fhigh,"\n")
     }
   if(is.null(roll))
    {
      roll = 10^5
      if(verbose) cat(" * Using default roll =",roll,"\n")
     }
 }  
      
if(is.null(targetE))
 {
# the five domintant eccentricity peaks based on spectral analysis of LA10d solution 
#   (Laskar et al., 2011), 0-20 Ma
     targetE <- double(5)
     targetE[1] = 405.6795
     targetE[2] = 130.719
     targetE[3] = 123.839
     targetE[4] = 98.86307
     targetE[5] = 94.87666
     if(verbose) cat(" * Using default eccentricity target periods =",targetE,"\n")     
  }   

if(fit==2 && !is.null(targetP)) 
  {
    if(verbose) cat("\n**** WARNING: targetP is defined but will not be used in fitting!\n")
  }  

# targetTot is for plotting, and fitting if precession modulations assessed
targetTot = c(targetE,targetP)


#######################################################################################
# Definition of FUNCTIONS: genCycles and fitIt.
# function to generate cos (real) and sin (imaginary) terms for each target period, 
#   and convert to spatial cycles, given a particular sed rate in m/ka
genCycles <- function(sedrate1, targetIn, n) 
  {
    x <- matrix(0, n, 2*length(targetIn))
    for (i in 1:length(targetIn)) 
      {
        x[,2*i-1] <- cos( (2*pi)/(targetIn[i]) * (dx/sedrate1) * (1:n))
        x[,2*i] <- sin( (2*pi)/(targetIn[i]) * (dx/sedrate1) * (1:n))
      }
    return(x)
  }
  
# function to perform fitting and calculate r, r-squared.
#   dx passed into function transparently
fitIt <- function(sedrate1,timeSeries,targetIn) 
  {
    xm <- genCycles(sedrate1, targetIn, npts)
    lm.0 <- lm(timeSeries[,2] ~ xm)
    if(cormethod==1) rval = cor(timeSeries[,2],lm.0$fitted,method=c("pearson"))
    if(cormethod==2) rval = cor(timeSeries[,2],lm.0$fitted,method=c("spearman"))
    rsq=rval^2
    return(cbind(sedrate1,rsq,rval))
   } 
  
#######################################################################################

# CALIBRATE DEPTH SERIES (m) TO TIME (ka)
    ts = dat
# create new time vector
# it is the index vector for time
    it <- seq(1,npts,by=1)
    time = (dx/sedrate) * (it-1)
    ts[1] = time

# bandpass precession or short eccentricity band
    bp = taner(ts,padfac=2,flow=flow,fhigh=fhigh,roll=roll,demean=T,detrend=F,addmean=F,genplot=F,verbose=F)

# hilbert transform for instantaneous amplitude
    hil = hilbert(bp,genplot=F,verbose=F)

# execute functions
# for precession modulations
    if(fit == 1) 
      {
         res = fitIt(sedrate,hil,targetE)
         pwrOut = fitIt(sedrate,ts,targetTot)
       }
# for short eccentricity modulations
    if(fit == 2) 
      {
         res = fitIt(sedrate,hil,targetE[1])
         pwrOut = fitIt(sedrate,ts,targetE)
       }

    if(verbose) {if(res[3] < 0) cat("\n**** WARNING: DATA modulation fit correlation <0 \n")} 
    if(verbose) {if(pwrOut[3] < 0) cat("\n**** WARNING: DATA power fit correlation <0 \n")}          
                
    datPwr=pwrOut[,2]
    datCor=res[,2]
    datCorPwr = datPwr*datCor
   
if(verbose)
 {  
  cat("\n * Spectral Power r^2 =", datPwr,"\n")
  cat(" * Envelope r^2 =", datCor,"\n")
  cat(" * (Envelope r^2) x (Spectral Power r^2) =", datCorPwr,"\n")
 }

#######################################################################################
# Monte Carlo Simulation

if(verbose) {  cat("\n * PLEASE WAIT: Performing", numsim,"simulations \n") }

#  create output array, dimension appropriately
#  simres will contain envelope, power, envelope*power
simres <- rep(NA,numsim*3)
dim(simres) <- c(numsim,3)

# begin simulation loop
isim=0
for (isim in 1:numsim) 
  {
# generate AR1 noise
    sim = ar1(npts, dx, mean=0, sdev=1, rho=rho, genplot=F, verbose=F)
# recenter and standardize
    sim[2]=sim[2]-colMeans(sim[2])
    sim[2]=sim[2]/sapply(sim[2],sd)

# CALIBRATE DEPTH SERIES (m) TO TIME (ka)
# create new time vector
# index vector for time
    it <- seq(1,npts,by=1)
    time = (dx/sedrate) * (it-1)
    sim[1] = time

# bandpass precession or short eccentricity band
    bpSim = taner(sim,padfac=2,flow=flow,fhigh=fhigh,roll=roll,demean=T,detrend=F,addmean=F,genplot=F,verbose=F)

# hilbert transform for instantaneous amplitude
    hilSim = hilbert(bpSim,genplot=F,verbose=F)

# execute functions
# for precession modulations
    if(fit == 1) 
      {
         res = fitIt(sedrate,hilSim,targetE)
         pwrOut = fitIt(sedrate,sim,targetTot)
       }
# for short eccentricity modulations
    if(fit == 2) 
      {
         res = fitIt(sedrate,hilSim,targetE[1])
         pwrOut = fitIt(sedrate,sim,targetE)
       }

    if(verbose) {if(res[3] < 0) cat("\n**** WARNING: SIMULATION modulation fit correlation <0 \n")} 
    if(verbose) {if(pwrOut[3] < 0) cat("\n**** WARNING: SIMULATION power fit correlation <0 \n")}          
         
    simres[isim,1] <- res[,2]
    simres[isim,2] <- pwrOut[,2]
    simres[isim,3] <- res[,2]*pwrOut[,2]
# end simulation loop    
 }

# now sort results, determine how many have values > your result
# envelope
    numgt = sum(simres[,1]>datCor)
    pvalCor=numgt/numsim
    if(pvalCor < (10/numsim) && (10/numsim) <=1 ) pvalCor= 10/numsim
    if(pvalCor >= (10/numsim) && (10/numsim) <=1 ) pvalCor=pvalCor   
    if((10/numsim) > 1 ) pvalCor=1  
# spectral power
    numgt = sum(simres[,2]>datPwr)
    pvalPwr=numgt/numsim 
    if(pvalPwr < (10/numsim) && (10/numsim) <=1 ) pvalPwr= 10/numsim
    if(pvalPwr >= (10/numsim) && (10/numsim) <=1 ) pvalPwr=pvalPwr   
    if((10/numsim) > 1 ) pvalPwr=1
# envelope * spectral power
    numgt = sum(simres[,3]>datCorPwr)
    pvalCorPwr=numgt/numsim 
    if(pvalCorPwr < (10/numsim) && (10/numsim) <=1 ) pvalCorPwr= 10/numsim
    if(pvalCorPwr >= (10/numsim) && (10/numsim) <=1 ) pvalCorPwr=pvalCorPwr   
    if((10/numsim) > 1 ) pvalCorPwr=1
     
    if(verbose) cat(" * (Envelope r^2) * (Spectral Power r^2) p-value =",pvalCorPwr, "\n")

     
    if(genplot)
     { 
      par(mfrow=c(1,1))
      plot(density(simres[,3]), col="red",xlim=c(0,1),type="l",main="Simulation (Envelope r^2) * (Spectral Power r^2) Distribution")
      lines(density(simres[,3]))
      abline(v=datCorPwr,col="black",lty=3) 
     }
     
# output = (0) nothing, (1) envelope*spectral power r^2 p-value, (2) output simulation r^2 results
     if(output == 1) return(data.frame(pvalCorPwr))
     if(output == 2) return(data.frame(simres))
     
### END function timeOptSim
}
