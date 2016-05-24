LoadBeatPolar <- function(HRVData, RecordName, RecordPath=".", verbose = NULL) {
#-------------------------------
# Loads beats from an ascii file
#-------------------------------
#	RecordName -> record containing RR values
#	RecordPath -> path
#-------------------------------

	dir=getwd()
  on.exit(setwd(dir))

	if (!is.null(verbose)) {
		cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
		SetVerbose(HRVData, verbose)
	    }
	if (HRVData$Verbose) {
		cat("** Loading beats positions for record:", RecordName,"**\n")
		cat("   Path:",RecordPath,"\n")
	}

	setwd(RecordPath)

	# Extracts time and date information from file
	headerinfo=scan(RecordName,what=character(0),sep="=",skip=4,nlines=2,quiet=TRUE)	
	
	regexpdate="[0-9]{4}[0-9]{2}[0-9]{2}"
	if (length(headerinfo[regexpr(regexpdate,headerinfo)==1])) {
		dateinfo=headerinfo[regexpr(regexpdate,headerinfo)==1]
	} else {
		dateinfo="19000101"
	}
	
	regexptime="[0-9]{2}:[0-9]{2}:[0-9]{2}.[0-9]{1}"
	if (length(headerinfo[regexpr(regexptime,headerinfo)==1])) {
		timeinfo=headerinfo[regexpr(regexptime,headerinfo)==1]
	} else {
		timeinfo="00:00:00.0"
	}
	
	datetimeinfo = paste(dateinfo,timeinfo)
	datetimeaux = strptime(datetimeinfo,"%Y%m%d %H:%M:%S")
	
	if (HRVData$Verbose) {
		cat("   Date: ",sprintf("%02d",datetimeaux$mday),"/",
			sprintf("%02d",1+datetimeaux$mon),"/",
			1900+datetimeaux$year,"\n",sep="")
		cat("   Time: ",sprintf("%02d",datetimeaux$hour),":",
			sprintf("%02d",datetimeaux$min),":",
			sprintf("%02.01f",datetimeaux$sec),"\n",sep="")
	}
	HRVData$datetime=datetimeaux

	aux=scan(RecordName,what=character(0),sep="=",skip=39,quiet=TRUE)
	HRVData$Beat$RR=as.numeric(aux[-(1:which(aux=="[HRData]"))])

	HRVData$Beat$Time=cumsum(HRVData$Beat$RR)/1000
	if (HRVData$Verbose) {
        	cat("   Number of beats:", length(HRVData$Beat$Time), "\n")
    	}	

  return(HRVData)
}
