LoadBeatAscii <- function (HRVData, RecordName, RecordPath=".", scale = 1, datetime = "1/1/1900 0:0:0", verbose = NULL) {
#-------------------------------
# Loads beats from an ASCII file
#-------------------------------
#	RecordName -> file containing values
#	RecordPath -> path
#-------------------------------

    if (!is.null(verbose)) {
        cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
        SetVerbose(HRVData, verbose)
    }
    if (HRVData$Verbose) {
        cat("** Loading beats positions for record:", RecordName,"**\n")
    }

    dir = getwd()
    on.exit(setwd(dir))
    if (HRVData$Verbose) {
       	cat("   Path:", RecordPath, "\n")
	cat("   Scale:", scale, "\n")
    }
    setwd(RecordPath)

    x = read.table(RecordName)
    beatsaux = x$V1
    beats = beatsaux[!duplicated(beatsaux)]
    if (length(beatsaux) != length(beats)) {
        if (HRVData$Verbose) {
            cat("   Removed", length(beatsaux) - length(beats),"duplicated beats\n")
        }
    }
    datetimeaux = strptime(datetime, "%d/%m/%Y %H:%M:%S")
    if (is.na(datetimeaux)) {
        cat("   --- ERROR: Date/time format is dd/mm/yyyy HH:MM:SS ---\n")
        return(HRVData)
    }
    if (HRVData$Verbose) {
        cat("   Date: ", sprintf("%02d", datetimeaux$mday), "/", 
            sprintf("%02d", 1 + datetimeaux$mon), "/", 1900 + 
                datetimeaux$year, "\n", sep = "")
        cat("   Time: ", sprintf("%02d", datetimeaux$hour), ":", 
            sprintf("%02d", datetimeaux$min), ":", sprintf("%02d", 
                datetimeaux$sec), "\n", sep = "")
    }
    HRVData$datetime = datetimeaux
    HRVData$Beat = data.frame(Time = beats * scale)
    if (HRVData$Verbose) {
        cat("   Number of beats:", length(HRVData$Beat$Time), "\n")
    }

   return(HRVData)
}
