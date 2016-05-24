### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### eAsm function - (SRM: January 13, 2014, January 29, 2015; May 18, 2015)
###
### wrapper to conduct Evolutive ASM analysis using FORTRAN code
###
### NOTE: this code uses asm1.8_R.f
###########################################################################

eAsm <- function (spec,siglevel=0.9,target,fper=NULL,rayleigh,nyquist,sedmin=1,sedmax=5,numsed=50,linLog=1,iter=100000,ydir=1,output=4,genplot=F)
{

 cat("\n----- PERFORMING EVOLUTIVE AVERAGE SPECTRAL MISFIT ANALYSIS -----\n")

# ensure we have a data frame
  spec=data.frame(spec)

# assign frequencies from first column of spec
  freq=spec[,1]

  cols=length(spec)
# assign locations for each spectrum (column headers)
  loc=suppressWarnings(as.numeric(substr(colnames(spec[2:cols]),start=2,stop=100)))
# for negative depth/height/time values, "-" has been changed to "."
# this will create NAs. implement modification of fix recommended by Mathieu Martinez
  neg=grepl(".",substr(colnames(spec[2:cols]), start=2,stop=2),fixed=T)
  fixloc=which(neg)
  if(any(neg)) {loc[fixloc]=-1*as.numeric(substr(colnames(spec[(fixloc+1)]),start=3,stop=100))}
# assign specta (amplitude, power, or probability)
  sp=as.matrix( spec[2:cols] )

  numrec=length(loc)
  numfreq=length(freq)
  cat("\n * Number of spectra to analyze =",numrec,"\n")
  cat(" * Number of frequencies per spectrum =",numfreq,"\n")

# isolate portion of the results that you want to analyze (part 1: 'freq')
   ifreq= which(freq <= nyquist)
   freq2=freq[ifreq]

# loop over all windows (note, these must be set to number of sedrates)
   resASM<-double(numsed*numrec)
   dim(resASM)<-c(numsed,numrec)
   resHo<-double(numsed*numrec)
   dim(resHo)<-c(numsed,numrec)

  ii=1
  for (i in 1:numrec)
    {
     cat("\n * PROCESSING WINDOW=",i,"; Location=",loc[i],"\n")
# isolate portion of the results that you want to analyze (part 2: 'sp')
      sp2 = sp[ifreq,i]
# if no peaks acheive specified siglevel, skip
      ij=sum(sp2>=siglevel)
      if(ij==0)
       {
         cat("\n**** WARNING: no peaks achieve specified significance level\n")
         resASM[,i] <- NA
         resHo[,i] <- 100
       }

      if(ij>0)
       {  
# find peaks, store in freqPeak
         freqPeak=peak(cbind(freq2,sp2),level=siglevel,genplot=F,verbose=T)[2]
# calculate ASM
         res <- asm(freq=freqPeak,target=target,fper=fper,rayleigh=rayleigh,nyquist=nyquist,sedmin=sedmin,sedmax=sedmax,numsed=numsed,linLog=linLog,iter=iter,output=TRUE,genplot=genplot)
         if(ii==1)
          {
            resSed <- res[,1]
            resTerms <- res[,4]
            ii=2
          }  
         resASM[,i] <- res[,2]
         resHo[,i] <- res[,3]
       }
    }

# now plot. 
    dev.new(title=paste("Evolutive ASM Results"))
# widths= a vector of values for the widths of columns
# heights= a vector of values for the heights of rows.
# Relative widths are specified with numeric values. Absolute widths (in centimetres) are specified with the lcm() function.
    layout(matrix(c(1,2), 1, 2, byrow = TRUE),widths=c(3.5,1), heights=c(2,2))
    brks <- c(0.00,0.05,0.1,0.2,0.4,0.6,0.8,1,5,10,100)    
    colorScale=append(rainbow(9,start=0,end=0.75),gray(0))
    if(ydir == 1) ylimset=c( min(loc),max(loc) )
    if(ydir == -1) ylimset=c( max(loc),min(loc) )
    image(log(resSed),loc,resHo,ylim=ylimset,breaks=brks,col=colorScale,xlab="Log(Sedimentation Rate, cm/ka)",ylab="Location")
 # add linear sedimentation rates at top
    mtext("Linear Sedimentation Rate, cm/ka",side=3,line=3)
    linSed=seq(from=par("xaxp")[1],to=par("xaxp")[2],by=(par("xaxp")[2]-par("xaxp")[1])/(par("xaxp")[3]))
    axis(side = 3, at = linSed,labels=F)
    for (i in 1:length(linSed))
     {
       mtext(round(exp(1)^linSed[i],digits=2),side=3,line=1,at=linSed[i])
     }
# add scale
    plot(-1,-1,ylim=c(0,1),xlim=c(0,1), yaxt='n',xaxt='n',bty='n',ylab="",xlab="",cex.lab=1.5,main="Ho-SL (%)")
    rect(0,0,1,0.1,col=colorScale[1])
    rect(0,0.1,1,0.2,col=colorScale[2])
    rect(0,0.2,1,0.3,col=colorScale[3])
    rect(0,0.3,1,0.4,col=colorScale[4])
    rect(0,0.4,1,0.5,col=colorScale[5])
    rect(0,0.5,1,0.6,col=colorScale[6])
    rect(0,0.6,1,0.7,col=colorScale[7])
    rect(0,0.7,1,0.8,col=colorScale[8])
    rect(0,0.8,1,0.9,col=colorScale[9])
    rect(0,0.9,1,1,col=colorScale[10])
    text(par("usr")[1]-1, 0, srt=0, adj = 0, labels = "0.0", xpd = TRUE)
    text(par("usr")[1]-1.3, 0.1, srt=0, adj = 0, labels = "0.05", xpd = TRUE)
    text(par("usr")[1]-1, 0.2, srt=0, adj = 0, labels = "0.1", xpd = TRUE)
    text(par("usr")[1]-1, 0.3, srt=0, adj = 0, labels = "0.2", xpd = TRUE)
    text(par("usr")[1]-1, 0.4, srt=0, adj = 0, labels = "0.4", xpd = TRUE)    
    text(par("usr")[1]-1, 0.5, srt=0, adj = 0, labels = "0.6", xpd = TRUE) 
    text(par("usr")[1]-1, 0.6, srt=0, adj = 0, labels = "0.8", xpd = TRUE)     
    text(par("usr")[1]-0.6, 0.7, srt=0, adj = 0, labels = "1", xpd = TRUE)     
    text(par("usr")[1]-0.6, 0.8, srt=0, adj = 0, labels = "5", xpd = TRUE)     
    text(par("usr")[1]-0.9, 0.9, srt=0, adj = 0, labels = "10", xpd = TRUE)     
    text(par("usr")[1]-1.2, 1, srt=0, adj = 0, labels = "100", xpd = TRUE)     

 if (output > 0) 
  {
# add column titles to identify each record for output
    colnames(resASM) <- loc
    colnames(resHo) <- loc
    outHo=data.frame(cbind(resSed,resHo))
    outASM=data.frame(cbind(resSed,resASM))
    outTerms=data.frame(cbind(resSed,resTerms))
    if(output == 1) return(outHo)
    if(output == 2) return(outASM)
    if(output == 3) return(outTerms)
    if(output == 4) return(list(outHo,outASM,outTerms))
  }

#### END function eAsm
}
