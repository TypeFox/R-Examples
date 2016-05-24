CreateFreqAnalysis <-
function(HRVData, verbose=NULL) {
# ---------------------------------------------------------
# Creates a frequency analysis associated to the data model
# ---------------------------------------------------------

	if (!is.null(verbose)) {
		cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
		SetVerbose(HRVData,verbose)
	}
	
	if (HRVData$Verbose) {
		cat("** Creating frequency analysis\n")
	}

   	num=length(HRVData$FreqAnalysis)

   	HRVData$FreqAnalysis[[num+1]]=list()


	if (HRVData$Verbose) {
		cat("   Data has now",num+1,"frequency analysis\n")
	}

   	return (HRVData)
}

