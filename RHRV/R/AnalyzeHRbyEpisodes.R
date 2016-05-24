AnalyzeHRbyEpisodes <-
function(HRVData, Tag="", func, ..., verbose=NULL) {
# ----------------------------------------------
# Analyzes Heart Rate using Episodes information
# ----------------------------------------------
#  Tag -> specifies tag of episodes
#  func -> function to apply 
#  Returns a list with two objects result

#  Function func musts receive a vector and returns an object

  funcToApply =  match.fun(func)
  nameFunc = deparse(substitute(func))

	if (!is.null(verbose)) {
		cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
		SetVerbose(HRVData,verbose)
	}
	
	if (HRVData$Verbose) {
		cat("** Applying function to heart rate signal using episodic information **\n");
    cat("   Function: ",nameFunc,"()\n",sep="")
   }

   if (is.null(HRVData$Episodes)) {
      stop("  --- Episodes not present\n    --- Quitting now!! ---\n")
   }

	if (is.null(HRVData$HR)) { 
      stop("  --- Interpolated heart rate not present\n    --- Quitting now!! ---\n")
	}

	if (HRVData$Verbose) {
      if (Tag=="") {
		   cat("   No tag was specified\n")
      } else {
		   cat("   Using episodes with tag:",Tag,"\n")
      }
	}

   vectors=SplitHRbyEpisodes(HRVData,Tag=Tag)

   resultIn=funcToApply(vectors$InEpisodes,...)
   resultOut=funcToApply(vectors$OutEpisodes,...)
   
   result=list(resultIn=resultIn,resultOut=resultOut)

   return(result)



}

