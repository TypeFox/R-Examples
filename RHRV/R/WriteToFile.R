WriteToFile <-
function(HRVData, name, overwrite=TRUE, verbose=NULL) {
# ---------------------------
# Writes data model to a file
# ---------------------------
#	overwrite: if true, overwrites previously existing file
	
	if (!is.null(verbose)) {
		cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
		SetVerbose(HRVData,verbose)
	}
	
	nameext=sprintf("%s.%s",name,HRVData$Ext)
	
	if (HRVData$Verbose) {
		cat("** Writing file:",nameext,"\n")
	}
	
	if (file.exists(nameext)) {
		if (HRVData$Verbose) {
			cat("   File ",nameext," already exists\n",sep="")
		}
		if (!overwrite) {
			stop("  --- File exists... No overwriting it!! ---\n    --- Quitting now!! ---\n")
		}
	}
		
	dput(HRVData,file=nameext)
	if (HRVData$Verbose) {
			cat("   ",file.info(nameext)$size," bytes written\n",sep="")
	}
}

