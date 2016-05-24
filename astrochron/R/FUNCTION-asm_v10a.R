### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### asm function - (SRM: August 23, 2012; Sept 1, 2012; May 20, 2013; 
###                      June 5-7, 2013; July 2, 2013; January 13, 2014
###                      January 15, 2014; June 27, 2014)
###
### conduct ASM analysis using FORTRAN code
###
### NOTE: this code uses asm1.8_R.f
###########################################################################


asm <- function (freq,target,fper=NULL,rayleigh,nyquist,sedmin=1,sedmax=5,numsed=50,linLog=1,iter=100000,output=F,genplot=T)
{

 cat("\n----- PERFORMING AVERAGE SPECTRAL MISFIT ANALYSIS -----\n")

# make sure freq, target and fper are data frames
 freq=data.frame(freq)
 numfreq=nrow(freq)
 target=data.frame(target)
 numtarg=nrow(target)

 specgen=1
 repl=F
 setfreq=2

# ERROR checking  
 if(!is.null(fper)) fper=data.frame(fper)
 if(is.null(fper)) 
   {
     cat("\n**** WARNING: No uncertainty assigned to astronomical target frequencies.\n")
     fper = double(numtarg)
     fper = data.frame(fper)
   }
 if(numtarg != nrow(fper)) 
   {
     cat("\n**** ERROR: target and fper do not have matching entries\n")
     stop("    TERMINATING NOW!")
   }  
 test=sum(abs(freq[,1]-sort(freq[,1])))
 if(test>0)
    {
     cat("\n**** ERROR: candidate astronomical frequencies are not increasing in order\n")
     stop("    TERMINATING NOW!")
   }  
 if(sedmin>sedmax)
    {
     cat("\n**** ERROR: sedmin > sedmax \n")
     stop("    TERMINATING NOW!")
   }  
 if(linLog != 0 && linLog != 1)
    {
     cat("\n**** ERROR: linLog must be set to 0 or 1 \n")
     stop("    TERMINATING NOW!")
   }  
 if(specgen != 1 && specgen != 2)
    {
     cat("\n**** ERROR: specgen must be set to 1 or 2\n")
     stop("    TERMINATING NOW!")
   }  
 if(setfreq != 1 && setfreq != 2 && setfreq != 3)
    {
     cat("\n**** ERROR: setfreq must be set to 1, 2, or 3\n")
     stop("    TERMINATING NOW!")
   }  
 if(iter <1 && iter>100000)
    {
     cat("\n**** ERROR: iter must be set between 1-100,000\n")
     stop("    TERMINATING NOW!")
   }  

  ispecgen=as.integer(specgen)
  isetfreq=as.integer(setfreq)
  
# note: even if specgen = 1, we must set irepl for call to FORTRAN  
  if(repl) irepl=as.integer(2) 
  if(!repl) irepl=as.integer(1)
   
### wrapper for FORTRAN code   
 asm1.8 <- function (freq,numfreq,target,numtarg,fper,rayleigh,nyquist,sedmin,sedmax,numsed,linLog,ispecgen,irepl,isetfreq,iter)
  {
    F_dat = .Fortran('asm18_r', PACKAGE='astrochron',
                freq=as.double(freq),numfreq=as.integer(numfreq),target=as.double(target),
                numtarg=as.integer(numtarg),fper=as.double(fper),rayleigh=as.double(rayleigh),
                nyquist=as.double(nyquist),sedmin=as.double(sedmin),sedmax=as.double(sedmax),
                numsed=as.integer(numsed),linLog=as.integer(linLog),ispecgen=as.integer(ispecgen),
                irepl=as.integer(irepl),isetfreq=as.integer(isetfreq),iter=as.integer(iter),
                
                newNumsed=integer(1),asmSedrate=double(numsed),asm=double(numsed),
                asmTerms=double(numsed),asmCL=double(numsed)
             )

# return the results
            return(F_dat)
  }

 res <- asm1.8(freq[,1],numfreq,target[,1],numtarg,fper[,1],rayleigh,nyquist,sedmin,sedmax,numsed,linLog,ispecgen,irepl,isetfreq,iter)

 newsed=res$newNumsed

# find sedimentation rate with smallest Ho-SL
 fit=which.min(res$asmCL[1:newsed])

# check for multiple minima
 sumMin = sum( res$asmCL[1:newsed] == min(res$asmCL[1:newsed]) )
 if(sumMin > 1) cat("\n**** WARNING: Multiple minima detected, only lowest sedrate plotted.\n") 

### generate plots
 if(genplot)
  {
    dev.new(title=paste("Average Spectral Misfit Results"),height=6,width=8)
    par(mfrow=c(2,2))
    plot(res$asmSedrate[1:newsed],res$asm[1:newsed],type="l", col="blue", xlab="Sedimentation Rate (cm/ka)",ylab="ASM (cycles/ka)",main="(a) Average Spectral Misfit",bty="n",lwd=2,cex.axis=1.1,cex.lab=1.1)
    abline(v=res$asmSedrate[fit],col="red",lwd=2,lty=4)
    plot(res$asmSedrate[1:newsed],res$asmTerms[1:newsed],type="l", col="blue", ylim=c(0,numtarg),xlab="Sedimentation Rate (cm/ka)",ylab="# Astronomical Terms",main="(c) Number of Terms Evaluated",bty="n",lwd=2,cex.axis=1.1,cex.lab=1.1)
    abline(v=res$asmSedrate[fit],col="red",lwd=2,lty=4)
    plot(res$asmSedrate[1:newsed],res$asmCL[1:newsed],log="y",type="l", col="blue",xlab="Sedimentation Rate (cm/ka)",ylab="Ho-SL (%)",main="(b) Null Hypothesis Significance Level",bty="n",lwd=2,cex.axis=1.1,cex.lab=1.1)
    abline(h=(100/newsed),col="black",lty=3)
    abline(v=res$asmSedrate[fit],col="red",lwd=2,lty=4)
    mtext(paste(round(res$asmSedrate[fit],digits=3),"cm/ka"),side=3,line=0,at=res$asmSedrate[fit],cex=0.9,font=4,col="red")

# plot calibrated spectrum versus target
    plot(1,1,xlim=c(0,target[numtarg,1]),ylim=c(0,1), yaxt='n',bty='n',ylab="",xlab="cycles/ka",main="(d) Data (black) vs. Target (red)",cex.axis=1.1,cex.lab=1.1)
# add horizontal line to extend x-axis
    abline(h=-0.04,lwd=2)
    for (i in 1:numtarg)
      {
        abline(v=target[i,1],col="red",lwd=2)
# include uncertainties
        xt0=target[i,1]-(target[i,1]*fper[i,1])
        xt1=target[i,1]+(target[i,1]*fper[i,1])
        segments(xt0,1,xt1,1,col="red",lwd=1.5)
      }  
 
     for (i in 1:numfreq)
      {
        ff=freq[i,1]*res$asmSedrate[fit]/100
        abline(v=ff,col="black",lwd=1,lty=3)
# include uncertainties
        xt0=ff-(0.5*rayleigh*res$asmSedrate[fit]/100)
        xt1=ff+(0.5*rayleigh*res$asmSedrate[fit]/100)
        segments(xt0,1.037,xt1,1.037,col="black",lwd=1)
      }  
  }

 cat("\n * Analysis complete:\n")
 cat("    Optimal Sedimentation Rate (cm/ka) at =",res$asmSedrate[fit],"\n")
 cat("    Ho-SL (%) =",res$asmCL[fit],"\n")
 cat("      or p-value =",res$asmCL[fit]/100,"\n")
 cat("    ASM (cycles/ka) =",res$asm[fit],"\n")
 cat("    Number of Astronomical Terms Fit =",res$asmTerms[fit],"\n")

 if (output) 
  {
    asmresult <- data.frame(cbind(res$asmSedrate[1:newsed],res$asm[1:newsed],res$asmCL[1:newsed],res$asmTerms[1:newsed]))
    return(asmresult)
  }

#### END function asm
}
