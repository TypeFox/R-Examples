### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### function testPrecession - (SRM: July 10, 2013; August 9, 2013; August 14, 2013;
###                         Sept. 16-17, 2014; Jan. 8, 2015; Jan. 20-22, 2015;
###                         Jan. 29, 2015; February 4, 2015; February 23, 2015;
###                         March 11, 2015; September 10, 2015)
###
### Perform astrochonologic testing as in Zeeden et al. (2015)
###########################################################################

testPrecession <- function(dat,nsim=1000,gen=1,rho=NULL,esinw=NULL,output=T,genplot=T,verbose=T)
{

   if(verbose) { cat("\n----- PERFORMING ZEEDEN ET AL. (2015) ASTROCHRONOLOGIC TESTING -----\n") }
   dat <- data.frame(dat)

# sort data to ensure increasing order
   if(verbose) cat(" * Sorting data into ensure increasing order, removing empty entries\n")
   dat <- dat[order(dat[1],na.last=NA,decreasing=F),]
   npts <- length(dat[,1]) 

# remove mean
   dat[2] <- dat[2] - colMeans(dat[2])
        
# error checking 
   maxdat=max(dat[,1])
   if(maxdat > 50000) 
      {
       cat("\n**** ERROR: this approach is only applicable when testing astrochronologies <= 50 Ma!\n")
       stop("**** TERMINATING NOW!")
     }

# error checking 
   mindat=min(dat[,1])
   if(mindat < 0) 
      {
       cat("\n**** ERROR: time (first column of input series) must be positive.\n")
       stop("**** TERMINATING NOW!")
     }


   dtest <- dat[2:npts,1]-dat[1:(npts-1),1] 
   epsm=1e-9
   if( (max(dtest)-min(dtest)) > epsm ) 
     {
       cat("\n**** ERROR: temporal sampling interval is not uniform.\n")
       stop("**** TERMINATING NOW!")
     }

   dt = dat[2,1]-dat[1,1]

   if(dt > 5)
     {
       cat("\n**** ERROR: temporal sampling interval must be <= 5 ka.\n")
       stop("**** TERMINATING NOW!")
     }

   if(verbose) { cat(" * Number of data points=", npts,"\n") }
   if(verbose) { cat(" * Sample interval=", dt,"\n") }
   
   if(!is.null(esinw)) 
    {
# check to ensure that esinw covers the analyzed series    
     min2=min(esinw[,1])
     max2=max(esinw[,1])
     if(min2>mindat || max2<maxdat)
      {
       cat("\n**** ERROR: esinw does not span the entire temporal interval of your input series.\n")
       stop("**** TERMINATING NOW!")
      }     
     xout <- seq(dat[1,1],dat[npts,1],by=dt)
     precT <- approx(esinw[,1],esinw[,2],xout,method="linear",n=npts)
     precT <- data.frame(precT)
    }
     
###########################################################################
### determine eccentricity signal for time interval
###########################################################################
   if(is.null(esinw)) precT=etp(tmin=dat[1,1],tmax=dat[npts,1],dt=dt,eWt=0,oWt=0,pWt=1,esinw=T,solution=NULL,standardize=F,genplot=F)
# bandpass filter to be consistent with tuned data processing
   precT_BP=taner(precT,padfac=2,flow=0.029,fhigh=0.12,roll=10^3,demean=T,detrend=F,genplot=F,verbose=F)
# determine instantaneous amplitude
   precT_BP_Hil=hilbert(precT_BP,padfac=2,demean=T,detrend=F,genplot=F,verbose=F)
# lowpass instantaneous amplitude
   ecc=taner(precT_BP_Hil,padfac=2,fhigh=0.013,flow=-0.013,xmax=.1,roll=10^4,demean=T,detrend=F,genplot=F,verbose=F)
   
###########################################################################
### process tuned data series
###########################################################################    
# bandpass
   prec_BP=taner(dat,padfac=2,flow=0.029,fhigh=0.12,roll=10^3,demean=T,detrend=F,genplot=F,verbose=F)
# determine instantaneous amplitude
   prec_BP_Hil=hilbert(prec_BP,padfac=2,demean=T,detrend=F,genplot=F,verbose=F)
# lowpass instantaneous amplitude
   prec_BP_Hil_Lowpass=taner(prec_BP_Hil,padfac=2,fhigh=0.013,flow=-0.013,xmax=.1,roll=10^4,demean=T,detrend=F,genplot=F,verbose=F)
   
# generate plot comparing data results to eccentricty target
   if(genplot)
    {
     dev.new(height=7,width=8)
     par(mfrow=c(2,1))
     plot(prec_BP,type="l")
     lines(prec_BP_Hil,col="red")
     lines(prec_BP_Hil_Lowpass,col="blue")
     plot(s(prec_BP_Hil_Lowpass),type="l",col="blue",xlab="Time (kiloyears)",ylab="Standardized variables",main="Comparison of theoretical (black) and observed (blue) modulation",lwd=2)
     lines(s(ecc),lwd=2)
    }
   
# calculate spearman rank correlation coefficient, and pearson (for comparison) 
   datcor=cor(prec_BP_Hil_Lowpass[,2],ecc[,2],method=c("spearman"))
   datcor2=cor(prec_BP_Hil_Lowpass[,2],ecc[,2],method=c("pearson"))
   if(verbose) {cat(" * Spearman rank correlation between tuned data envelope and eccentricity=",datcor,"\n") }
   if(verbose) {cat(" * Pearson correlation between tuned data envelope and eccentricity=",datcor2,"\n") }

# check for correlation <= 0   
   if(datcor<=0)
    {
       cat("\n**** ERROR: correlation must be positive to proceed.\n")
       stop("**** TERMINATING NOW!")
    }

if(nsim>0)
{   
###########################################################################
### perform Monte Carlo simulations using phase-randomized or AR1 surrogates
###########################################################################   
   if(gen==1)sur=surrogates(dat,nsim=nsim,preserveMean=T,std=T,genplot=F,verbose=T)
   if(gen==2)
    {
# estimate rho from data
       if(is.null(rho))
        {
          lag0 <- dat[1:(npts-1),2]
          lag1 <- dat[2:npts,2]
          rho = cor(lag0,lag1)
          if(verbose) cat(" * Estimated AR1 coefficient =",rho,"\n")
         } 
       sur=ar1(npts=npts,dt=dt,mean=0,sdev=1,rho=rho,shuffle=F,nsim=nsim,genplot=F,verbose=T)
    }
     
   simcor=double(nsim)
   for (i in 1:nsim)
    {
      surdat=data.frame(cbind(dat[,1],sur[,i]))
      sur_BP=taner(surdat,padfac=2,flow=0.029,fhigh=0.12,roll=10^3,demean=T,detrend=F,genplot=F,verbose=F)
# determine instantaneous amplitude
      sur_BP_Hil=hilbert(sur_BP,padfac=2,demean=T,detrend=F,genplot=F,verbose=F)
# lowpass instantaneous amplitude
      sur_BP_Hil_Lowpass=taner(sur_BP_Hil,padfac=2,fhigh=0.013,flow=-0.013,xmax=.1,roll=10^4,demean=T,detrend=F,genplot=F,verbose=F)

# calculate spearman rank correlation coefficient 
      simcor[i]=cor(sur_BP_Hil_Lowpass[,2],ecc[,2],method=c("spearman"))
      if(verbose) {cat("SIMULATION",i,"Correlation =",simcor[i],"\n") }
     }

     if(genplot)
      {
        dev.new(height=4,width=6)
        par(mfrow=c(1,1))
        plot(density(simcor, bw="nrd"),xlim=c(-1,1),type="l",col="red",main="Monte Carlo Results",lwd=2)
        grid()
        abline(v=datcor,col="blue",lwd=2,lty=3)
        mtext(round(datcor,digits=5),side=3,line=0,at=datcor,cex=1,font=4,col="blue")
      }  
       
# now sort simulation correlation results, determine how many have r > your result
# this follows Ebisuzaki (1997)
    pnumgt = sum(abs(simcor)>abs(datcor))
    cat("\n * Number of simulations with |r<sim>| > |r<dat>| = ", pnumgt,"\n")
    ppvalue=pnumgt/nsim
    if(ppvalue < (10/nsim) && (10/nsim) <=1 ) cat("\n * P-value <", 10/nsim,"\n")
    if(ppvalue >= (10/nsim) && (10/nsim) <=1 ) cat("\n * P-value =", ppvalue,"\n")    
    if((10/nsim) > 1 ) cat("\n * P-value = 1 \n")    

# return simlulation results for plotting
    if(output) return(simcor)

# end nsim>0 section
} 
 
if(nsim==0)
  { 
    out = data.frame(cbind(prec_BP[1],prec_BP[2],prec_BP_Hil[2],prec_BP_Hil_Lowpass[2],ecc[2]))
    colnames(out)[1]<-"Time"
    colnames(out)[2]<-"Precession_BP"
    colnames(out)[3]<-"Precession_BP_Amp"
    colnames(out)[4]<-"Precession_BP_Amp_LP"
    colnames(out)[5]<-"Eccentricity"
    if(output) return(out) 
   }
     
# end function testPrecession     
}     