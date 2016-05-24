CreateNonLinearAnalysis <-
function(HRVData, verbose=NULL) {
# ----------------------------------------------------------
# Creates a non linear analysis associated to the data model
# ----------------------------------------------------------

	if (!is.null(verbose)) {
		cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
		SetVerbose(HRVData,verbose)
	}
	
	if (HRVData$Verbose) {
		cat("** Creating non linear analysis\n")
   	}

   	num=length(HRVData$NonLinearAnalysis)

   	HRVData$NonLinearAnalysis[[num+1]]=list()

   	if (HRVData$Verbose) {
      	cat("   Data has now ",num+1," nonlinear analysis\n")
   	}

   	return(HRVData)
}

