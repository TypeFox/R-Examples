LoadApneaWFDB <-
function(HRVData, RecordName, RecordPath=".", Tag="APNEA", verbose=NULL) {
#--------------------------------------- 
# Loads apnea episodes from an wfdb file
#	Uses rdann from wfdbtools
#---------------------------------------
#	RecordName -> record containing beat positions
#	RecordPath -> path
#  Tag -> tag to include in episodes

	if (!is.null(verbose)) {
		cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
		SetVerbose(HRVData,verbose)
	}
	
	if (HRVData$Verbose) {
		cat("** Loading apnea episodes for record:",RecordName,"**\n")
	}
	
	dir=getwd()
  on.exit(setwd(dir))

	if (HRVData$Verbose) {
		cat("   Path:",RecordPath,"\n")
	}
	setwd(RecordPath)

   # Reads header, verbose=FALSE
   if (is.null(HRVData$datetime)) {
      if (HRVData$Verbose) {
         cat("   Reading header info for:",RecordName,"\n")
      }
      HRVData = LoadHeaderWFDB(HRVData,RecordName,RecordPath)
   } else {
      if (HRVData$Verbose) {
         cat("   Header info already present for:",RecordName,"\n")
      }
   }

   # Calls rdann to read apnea annotations
	command=paste("rdann -r",RecordName,"-a apn")
	if (HRVData$Verbose) {
		cat("   Command:",command,"\n")
	}
	x1=system(command,intern=TRUE)
   xlabels=substring(x1,27,27)
   xtimes=seq(from=60,to=60*length(xlabels),length.out=length(xlabels))

   if (HRVData$Verbose) {
      cat("   Number of labels:",length(xlabels),"\n")
   }

   index=c(TRUE,xlabels[2:length(xlabels)]!=xlabels[1:(length(xlabels)-1)])

   ylabels=xlabels[index]
   ytimes=xtimes[index]
   # Detects changes between labels "A" and "N"

   if (tail(ylabels,1)=="A") {
      ylabels=c(ylabels,"N")
      ytimes=c(ytimes,60*length(xlabels)+30)
   } # If the last point is "A", an "N" is added

   if (head(ylabels,1)=="N") {
      l=length(ylabels)
      ylabels=ylabels[2:l]
      ytimes=ytimes[2:l]
   } # If the first point is "N", it is removed

	
	 indexInit=seq(from=1,to=length(ytimes)-1,by=2) # Odd elements
   indexEnd=seq(from=2,to=length(ytimes),by=2) # Even elements

   HRVData=AddEpisodes(HRVData,
      InitTimes=ytimes[indexInit]-30,
      Tags=Tag,
      Durations=ytimes[indexEnd]-ytimes[indexInit],
      Values=0
   )

   return(HRVData)
}

