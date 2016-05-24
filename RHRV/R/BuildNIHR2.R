# source('/data/pulse/rhrv/pkg/R/BuildNIHR2.R', chdir = TRUE)
BuildNIDHR <-
function(HRVData, verbose=NULL) {
#------------------------------------------------------ 
# Obtains instantaneous heart rate variation from beats positions
# D for difference
#------------------------------------------------------ 
	if (!is.null(verbose)) {
		cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
		SetVerbose(HRVData,verbose)
	}
	
	if (HRVData$Verbose) {
		cat("** Calculating non-interpolated heart rate differences **\n")
	}

	if (is.null(HRVData$Beat$Time)) {
		cat("   --- ERROR: Beats positions not present... Impossible to calculate Heart Rate!! ---\n")
		return(HRVData)
	}
	
	NBeats=length(HRVData$Beat$Time)
	if (HRVData$Verbose) {
		cat("   Number of beats:",NBeats,"\n");
	}
	

   #using NA, not constant extrapolation as else in RHRV  
   #drr=c(NA,NA,1000.0*	diff(HRVData$Beat$Time, lag=1 , differences=2))
   HRVData$Beat$dRR=c(NA, NA, 
   	1000.0*diff(HRVData$Beat$Time, lag=1, differences=2))

   HRVData$Beat$avRR=(c(NA,HRVData$Beat$RR[-1])+HRVData$Beat$RR)/2

	return(HRVData)
}

