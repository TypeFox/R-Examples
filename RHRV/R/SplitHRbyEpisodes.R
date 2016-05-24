SplitHRbyEpisodes <-
function(HRVData, Tag="", verbose=NULL) {
# -------------------------------------------------
# Splits Heart Rate Data using Episodes information
# -------------------------------------------------
#  Tag -> specifies tag of episodes
#  Returns a list with two vectors: InEpisodes and OutEpisodes

	if (!is.null(verbose)) {
		cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
		SetVerbose(HRVData,verbose)
	}
	
	if (HRVData$Verbose) {
		cat("** Splitting heart rate signal using episodes **\n");
   }

   if (is.null(HRVData$Episodes)) {
      stop("  --- Episodes not present\n    --- Quitting now!! ---\n")
   }

	if (is.null(HRVData$HR)) { 
      stop("  --- Interpolated heart rate not present\n    --- Quitting now!! ---\n")
	}

	if (HRVData$Verbose) {
      if (Tag=="") {
		   cat("   No tag was specified\n");
      } else {
		   cat("   Using episodes with tag:",Tag,"\n");
      }
	}

   # Select episodes to split signal
   if (Tag=="") {
      ActiveEpisodes=HRVData$Episodes
   } else {
      ActiveEpisodes=subset(HRVData$Episodes,HRVData$Episodes$Type==Tag)
   }

   if (HRVData$Verbose) {
      cat("   Number of episodes:",length(ActiveEpisodes$InitTime),"\n")
   }

   Beg=ActiveEpisodes$InitTime
   End=ActiveEpisodes$InitTime+ActiveEpisodes$Duration

	npoints = length(HRVData$HR)
	first = head(HRVData$Beat$Time,1)
	last = tail(HRVData$Beat$Time,1)
	x=seq(first,last,length.out=npoints)

   # Auxiliary signal used to mark points inside episodes
   Aux=rep(0,times=npoints)
   for (i in 1:length(Beg)) {
      Aux[x>=Beg[i] & x<=End[i]] = 1
   }

   l=list(InEpisodes=HRVData$HR[Aux==1],OutEpisodes=HRVData$HR[Aux==0])

   if (HRVData$Verbose) {
      cat("   Inside episodes:",length(l$InEpisodes),"points\n")
      cat("   Outside episodes:",length(l$OutEpisodes),"points\n")
   }

   return(l)
}

