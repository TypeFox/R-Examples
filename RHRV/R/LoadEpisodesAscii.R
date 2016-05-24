LoadEpisodesAscii <-
function(HRVData, FileName, RecordPath=".", Tag="", InitTime="0:0:0", verbose=NULL, header = TRUE) {	
#-------------------------------
# Loads episodes from ascii file
#-------------------------------
#	FileName -> file containing episodes
# RecordPath -> path to the file containing episodes
#   Tag -> specifies type of episodes
#	InitTime -> time (HH:MM:SS) absolute time of beginning of the record (subtracted from time of episodes)

#  	Example of file containing episodes:

#  	Init_Time	Resp_Events	Durat	SaO2
#  	00:33:00        GEN_HYPO	120.0	82.9
#  	01:30:00        OBS_APNEA	60.0	81.0
#  	...

#  	First line of file is discarded if the header argument is set to TRUE
#  	Duration in seconds
	
  dir = getwd()
  on.exit(setwd(dir))
  setwd(RecordPath)

	if (!is.null(verbose)) {
		cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
		SetVerbose(HRVData,verbose)
	}
	
	if (HRVData$Verbose) {
		cat("** Loading episodes file:",FileName,"**\n")
		cat("   Path:", RecordPath, "\n")
	}
  if (header){
    # skip first line
    x=read.table(FileName,skip=1)  
  }else{
    # Do not skip the first line,
    x=read.table(FileName,skip=0)  
  }
	

	if (HRVData$Verbose) {
   	if (Tag=="") {
      	cat("   No tag specified\n")
   	} else {
      	cat("   Tag:",Tag,"\n")
   	}
	}

	# obtaining time
  options(digits.secs=3)
	timeaux = strptime(InitTime,"%H:%M:%OS")
	if (is.na(timeaux)) {
		cat("   --- ERROR: Time format is HH:MM:SS ---\n")
		return(HRVData)
	}	
	
	if (HRVData$Verbose) {
		cat("   Initial time: ",sprintf("%02d",timeaux$hour),":",
			sprintf("%02d",timeaux$min),":",
			sprintf("%02.3f",timeaux$sec),"\n",sep="")
	}

	# calculating time in seconds, considering the initial time for the register

	EpisodeTimeAbs=strptime(x$V1,"%H:%M:%OS")
	EpisodeTimeRel=difftime(EpisodeTimeAbs,timeaux, units="secs")
	x$V1=as.numeric(EpisodeTimeRel)

	if (Tag=="") {
   	y=x
	} else {
   	y=subset(x,x$V2==Tag)
	}
   	
	added=length(y$V1)

   if (HRVData$Verbose) {
      if (is.null(y$V4)) {
         cat("   Data does not include values associated to episodes\n")
      } else {
         cat("   Data includes values associated to episodes\n")
      }
   }

	if (added==0) {
   	if (HRVData$Verbose) {
      	cat("   No episode was loaded\n")
   	}
	} else {

      if (is.null(y$V4)) {
         HRVData$Episodes=rbind(HRVData$Episodes,data.frame(InitTime=y$V1,Type=y$V2,Duration=y$V3))
      } else {
   	  HRVData$Episodes=rbind(HRVData$Episodes,data.frame(InitTime=y$V1,Type=y$V2,Duration=y$V3,Value=y$V4))
      }

   	HRVData$Episodes=HRVData$Episodes[order(HRVData$Episodes$InitTime),]  # Sorts episodes by InitTime
   	HRVData$Episodes=HRVData$Episodes[!duplicated(HRVData$Episodes),]  # Removes duplicated episodes
   
   	if (HRVData$Verbose) {
      	cat("   Loaded",added,"episodes from file\n")
   	}
	}

	if (HRVData$Verbose) {
   	cat("   Number of episodes:",length(HRVData$Episodes$InitTime),"\n")
	}

	return(HRVData)
}

