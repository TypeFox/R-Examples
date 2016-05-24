AddEpisodes <-
function(HRVData, InitTimes,Tags, Durations ,Values, verbose=NULL) {	
#------------------
# Adds new episodes 
#------------------
#	InitTimes -> Vector containing initial times in seconds
#	Tags -> Vector containing types of episodes
#	Durations -> Vector containing durations in seconds
#	Values -> Vector containing numerical values for episodes

	if (!is.null(verbose)) {
		cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
		SetVerbose(HRVData,verbose)
	}
	
   	if (HRVData$Verbose) {
		cat("** Adding new episodes **\n")
	}


   	NewEpisodes=data.frame(InitTime=InitTimes,Type=Tags,Duration=Durations,Value=Values)
   	if (HRVData$Verbose) {
      	cat("   Added",length(NewEpisodes$InitTime),"episodes from file\n")
   	}

   	HRVData$Episodes=rbind(HRVData$Episodes,NewEpisodes)
   	HRVData$Episodes=HRVData$Episodes[order(HRVData$Episodes$InitTime),]  # Sort episodes by InitTime
   	HRVData$Episodes=HRVData$Episodes[!duplicated(HRVData$Episodes),]  # Remove duplicated episodes

   	if (HRVData$Verbose) {
      	cat("   Number of episodes:",length(HRVData$Episodes$InitTime),"\n")
   	}

   	return(HRVData)
}

