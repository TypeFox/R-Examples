LoadHeaderWFDB <-
function(HRVData, RecordName, RecordPath=".", verbose=NULL) {
#------------------------------------
# Loads header info from an wfdb file
#------------------------------------
#	RecordName -> record containing beat positions
#	RecordPath -> path

	if (!is.null(verbose)) {
		cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
		SetVerbose(HRVData,verbose)
	}
	
	dir=getwd()
  on.exit(setwd(dir))

	if (HRVData$Verbose) {
		cat("   Path:",RecordPath,"\n")
	}
	setwd(RecordPath)

	# Extracts time and date information from wfdb header
	headerfile=paste(RecordName,".hea",sep="")
	if (HRVData$Verbose) {
		cat("   Opening header file:",headerfile,"\n")
	}
	headerinfo=scan(headerfile,what=character(0),nlines=1,quiet=TRUE)	
	
	regexptime="[[:digit:]]{2}:[[:digit:]]{2}:[[:digit:]]{2}"
	if (length(headerinfo[regexpr(regexptime,headerinfo)==1])) {
		timeinfo=headerinfo[regexpr(regexptime,headerinfo)==1]
		if (HRVData$Verbose) {
			cat("      Time information in header:",timeinfo,"\n")
		}
	} else {
		timeinfo="00:00:00"
		if (HRVData$Verbose) {
			cat("      No time information in header:",timeinfo,"\n")
		}
	}
	
	regexpdate="[[:digit:]]{2}/[[:digit:]]{2}/[[:digit:]]{4}"
	if (length(headerinfo[regexpr(regexpdate,headerinfo)==1])) {
		dateinfo=headerinfo[regexpr(regexpdate,headerinfo)==1]
		if (HRVData$Verbose) {
			cat("      Date information in header:",dateinfo,"\n")
		}
	} else {
		dateinfo="01/01/1900"
		if (HRVData$Verbose) {
			cat("      No date information in header:",dateinfo,"\n")
		}
	}
	
	datetimeinfo = paste(dateinfo,timeinfo)
	datetimeaux = strptime(datetimeinfo,"%d/%m/%Y %H:%M:%S")
	
	if (HRVData$Verbose) {
		cat("   Date: ",sprintf("%02d",datetimeaux$mday),"/",
			sprintf("%02d",1+datetimeaux$mon),"/",
			1900+datetimeaux$year,"\n",sep="")
		cat("   Time: ",sprintf("%02d",datetimeaux$hour),":",
			sprintf("%02d",datetimeaux$min),":",
			sprintf("%02d",datetimeaux$sec),"\n",sep="")
	}
	HRVData$datetime=datetimeaux

	return(HRVData)

}

