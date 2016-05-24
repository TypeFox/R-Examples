CreateTimeAnalysis <-
function(HRVData, size=300, numofbins=NULL, interval=7.8125, verbose=NULL ) {
# ----------------------------------------------------
# Creates a Time analysis associated to the data model
# ----------------------------------------------------
#  	size: size of window (sec.)
#  	interval: width of bins in histogram for TINN and HRV index (msec.)
  	if (!is.null(verbose)) {
	  	cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
		  SetVerbose(HRVData,verbose)
	  } 
	
   	if (HRVData$Verbose) {
      	cat("** Creating time analysis\n")
   	}

   	num=length(HRVData$TimeAnalysis)

   	HRVData$TimeAnalysis[[num+1]]=list()

   	HRVData$TimeAnalysis[[num+1]]$size=size # length size for analysis

    # vecthist contains the bins for the histogram
    minRR=min(HRVData$Beat$RR)
    maxRR=max(HRVData$Beat$RR)
    if (!is.null(numofbins)){
      interval = (maxRR-minRR)/(numofbins-2)
      vecthist = seq(minRR-interval/2,maxRR+interval/2,len=numofbins)
    }else{
      medRR=(min(HRVData$Beat$RR)+max(HRVData$Beat$RR))/2.0
      lowhist=medRR-interval*ceiling((medRR-minRR)/interval)
      longhist=ceiling((maxRR-lowhist)/interval)+1
      vecthist=seq(from=lowhist,by=interval,length.out=longhist)
    }
   	if (HRVData$Verbose) {
      	cat("   Size of window:",size,"seconds \n")
      	cat("   Width of bins in histogram:",interval,"milliseconds \n")
   	}

   	# SDNN
   	HRVData$TimeAnalysis[[num+1]]$SDNN=sd(HRVData$Beat$RR)

   	WindowMin=head(HRVData$Beat$Time,n=1)
   	WindowMax=WindowMin + size
   	WindowIndex=1
   	RRWindowMean=c(0)
   	RRWindowSD=c(0)
   	while (WindowMax < tail(HRVData$Beat$Time,1)) {
      	RRWindow=HRVData$Beat$RR[HRVData$Beat$Time >= WindowMin & HRVData$Beat$Time < WindowMax]
        # check if there is an interval without beats
        if (length(RRWindow) == 0){
          message = paste(sep="", "Interval without beats from ",WindowMin,
                          " to ",WindowMax," seconds!!\n  Returning NA in SDANN and SDNNIDX\n")
          warning(message)
          # introduce the NAs to ensure that the user notices the warning
          RRWindowMean[WindowIndex] = NA
          RRWindowSD[WindowIndex] = NA
          # there is no need to compute more windows
          break;
        }
      	RRWindowMean[WindowIndex]=mean(RRWindow)
      	RRWindowSD[WindowIndex]=sd(RRWindow)
      	WindowMin = WindowMin+size
      	WindowMax = WindowMax+size
      	WindowIndex = WindowIndex+1
   	}

    numberOfWindows = WindowIndex-1
   	if (HRVData$Verbose) {
       cat("   Number of windows:",numberOfWindows,"\n")
   	}

    
    if (numberOfWindows <= 1){
      warning("There is no window or just one window. Cannot compute the standard deviation!! Returning NA in SDANN\n")
    }
   	# SDANN
   	HRVData$TimeAnalysis[[num+1]]$SDANN=sd(RRWindowMean) 

   	# SDNNIDX
   	HRVData$TimeAnalysis[[num+1]]$SDNNIDX=mean(RRWindowSD) 

   	# pNN50
	  NRRs=length(HRVData$Beat$RR)
   	RRDiffs = diff(HRVData$Beat$RR)
   	RRDiffs50=RRDiffs[abs(RRDiffs)>50]
   	HRVData$TimeAnalysis[[num+1]]$pNN50=100.0*length(RRDiffs50)/length(RRDiffs)

    # SDSD
    HRVData$TimeAnalysis[[num+1]]$SDSD = sd(RRDiffs)

   	# rMSSD
   	HRVData$TimeAnalysis[[num+1]]$rMSSD=sqrt(mean(RRDiffs^2))

   	# IRRR
   	RRQuant=quantile(RRDiffs)
   	HRVData$TimeAnalysis[[num+1]]$IRRR=RRQuant[[4]]-RRQuant[[2]]

   	# MADRR
   	HRVData$TimeAnalysis[[num+1]]$MADRR=median(abs(RRDiffs))

   	# TINN and HRV index
   	h = hist(HRVData$Beat$RR, breaks=vecthist, plot=FALSE)
   	area=length(HRVData$Beat$RR)*interval
   	maxhist=max(h$counts)
   	HRVData$TimeAnalysis[[num+1]]$TINN=area/maxhist
   	HRVData$TimeAnalysis[[num+1]]$HRVi=length(HRVData$Beat$RR)/maxhist



   	if (HRVData$Verbose) {
      	cat("   Data has now",num+1,"time analyses\n")
      	cat("      SDNN:",HRVData$TimeAnalysis[[num+1]]$SDNN,"msec. \n")
      	cat("      SDANN:",HRVData$TimeAnalysis[[num+1]]$SDANN,"msec. \n")
     	  cat("      SDNNIDX:",HRVData$TimeAnalysis[[num+1]]$SDNNIDX,"msec. \n")
      	cat("      pNN50:",HRVData$TimeAnalysis[[num+1]]$pNN50,"%\n")
      	cat("      SDSD:",HRVData$TimeAnalysis[[num+1]]$SDSD,"msec.\n")
      	cat("      r-MSSD:",HRVData$TimeAnalysis[[num+1]]$rMSSD,"msec.\n")
      	cat("      IRRR:",HRVData$TimeAnalysis[[num+1]]$IRRR,"msec.\n")
      	cat("      MADRR:",HRVData$TimeAnalysis[[num+1]]$MADRR,"msec.\n")
      	cat("      TINN:",HRVData$TimeAnalysis[[num+1]]$TINN,"msec.\n")
      	cat("      HRV index:",HRVData$TimeAnalysis[[num+1]]$HRVi,"\n")
   	}

   	return(HRVData)
}

