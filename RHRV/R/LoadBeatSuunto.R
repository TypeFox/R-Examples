LoadBeatSuunto <- function(HRVData, RecordName, RecordPath=".", verbose = NULL) {
#-------------------------------
# Loads beats from an ascii file
#-------------------------------
#	RecordName -> record containing RR values
#	RecordPath -> path
#-------------------------------

	dir=getwd()

	if (!is.null(verbose)) {
		cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
		SetVerbose(HRVData, verbose)
	    }
	if (HRVData$Verbose) {
		cat("** Loading beats positions for record:", RecordName,"**\n")
		cat("   Path:",RecordPath,"\n")
	}

	setwd(RecordPath)

	#Date and time information
	aux=scan(RecordName,what=character(0),sep="=",quiet=TRUE)
	date=aux[which(aux=="STARTTIME")+1]
	dateAux = substr(date,1,10)
	dateAux = gsub("\\.","-",dateAux)

	time = substr(date,12,19)
	time = gsub("\\.",":",time)
		
		
	if (HRVData$Verbose) {
		
		cat("   Date: ",dateAux, "\n")
		cat("   Time: ",time, "\n")
	}
	datetimeinfo = paste(dateAux,time,sep = " ")
	HRVData$datetime=datetimeinfo

	aux=scan(RecordName,what=character(0),sep="=",quiet=TRUE)
	HRVData$Beat$RR=as.numeric(aux[-(1:which(aux=="[CUSTOM1]"))])

	HRVData$Beat$Time=cumsum(HRVData$Beat$RR)/1000
	if (HRVData$Verbose) {
        	cat("   Number of beats:", length(HRVData$Beat$Time), "\n")
    	}
	setwd(dir)
	return(HRVData)
}
