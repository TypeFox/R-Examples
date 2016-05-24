ReadFromFile <-
function(name, verbose=FALSE) {
# ---------------------------
# Reads data model from a file
# ---------------------------

  HRVData = CreateHRVData()
	nameext=sprintf("%s.%s",name,HRVData$Ext)
	
	if (verbose) {
		cat("** Reading file:",nameext,"\n")
	}
	
	if (!file.exists(nameext)) {
		stop("  --- File does not exist!! ---\n    --- Quitting now!! ---\n")
	}
	
	HRVData=dget(nameext)
	
	if (HRVData$Verbose) {
			cat("   ",file.info(nameext)$size," bytes read\n",sep="")
	}
  HRVData = SetVerbose(HRVData=HRVData,Verbose=verbose)
  return (HRVData)
}

