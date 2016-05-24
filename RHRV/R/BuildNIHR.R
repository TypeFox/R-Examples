BuildNIHR <-
function(HRVData, verbose=NULL) {
#------------------------------------------------------ 
# Obtains instantaneous heart rate from beats positions
#------------------------------------------------------ 

	if (!is.null(verbose)) {
		cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
		SetVerbose(HRVData,verbose)
	}
	
	if (HRVData$Verbose) {
		cat("** Calculating non-interpolated heart rate **\n")
	}

	if (is.null(HRVData$Beat$Time)) {
		cat("   --- ERROR: Beats positions not present... Impossible to calculate Heart Rate!! ---\n")
		return(HRVData)
	}
	
	NBeats=length(HRVData$Beat$Time)
	if (HRVData$Verbose) {
		cat("   Number of beats:",NBeats,"\n");
	}
	
	hr=c(0)
hr[2:NBeats]=60.0/diff(HRVData$Beat$Time)
	hr[1]=hr[2] # Not a real data
	HRVData$Beat$niHR = hr

   rr=c(0)
   rr[2:NBeats]=1000.0*diff(HRVData$Beat$Time)
	rr[1]=rr[2] # Not a real data
   HRVData$Beat$RR=rr


	return(HRVData)
}

