SplitPowerBandByEpisodes <-
function(HRVData, indexFreqAnalysis = length(HRVData$FreqAnalysis), Tag="", verbose=NULL) {
# ------------------------------------------------
# Splits Power Per Band using Episodes information
# ------------------------------------------------
#  Tag -> specifies tag of episodes
#  Returns a list with two lists: InEpisodes and OutEpisodes
#    Both lists include ULF, VLF, LF and HF bands

	if (!is.null(verbose)) {
		cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
		SetVerbose(HRVData,verbose)
	}
	
   if (HRVData$Verbose) {
		cat("** Splitting power bands using episodes**\n")
	}

   if (is.null(HRVData$Episodes)) {
      stop("  --- Episodes not present\n    --- Quitting now!! ---\n")
   }

	if (is.null(HRVData$FreqAnalysis[[indexFreqAnalysis]]$ULF)) {
      stop("  --- Power per band not present\n    --- Quitting now!! ---\n")
	}

	if (HRVData$Verbose) {
      if (Tag=="") {
		   cat("   No tag was specified\n");
      } else {
		   cat("   Using episodes with tag:",Tag,"\n");
      }
	}

   # Select episodes to split bands
   if (Tag=="") {
      ActiveEpisodes=HRVData$Episodes
   } else {
      ActiveEpisodes=subset(HRVData$Episodes,HRVData$Episodes$Type==Tag)
   }

   if (HRVData$Verbose) {
      cat("   Number of episodes:",length(ActiveEpisodes$InitTime),"\n")
   }

   lframes=length(HRVData$FreqAnalysis[[indexFreqAnalysis]]$HRV)
   # lframes is the number of frames

   EpisodesLeft=ActiveEpisodes$InitTime # Beg of episodes (seconds)
   EpisodesLeftFrame=EpisodesLeft*lframes/(tail(HRVData$Beat$Time,1)-head(HRVData$Beat$Time,1)) # Beg of episodes (frames)
   EpisodesRight=ActiveEpisodes$InitTime+ActiveEpisodes$Duration # Beg of episodes (seconds)
   EpisodesRightFrame=EpisodesRight*lframes/(tail(HRVData$Beat$Time,1)-head(HRVData$Beat$Time,1)) # Beg of episodes (frames)

   index=c()
   for (i in 1:length(EpisodesLeft)) {
      index=c(index,EpisodesLeftFrame[i]:EpisodesRightFrame[i])
   }

   l=list()

   l$InEpisodes=list(ULF=HRVData$FreqAnalysis[[indexFreqAnalysis]]$ULF[index],
      VLF=HRVData$FreqAnalysis[[indexFreqAnalysis]]$VLF[index],
      LF=HRVData$FreqAnalysis[[indexFreqAnalysis]]$LF[index],
      HF=HRVData$FreqAnalysis[[indexFreqAnalysis]]$HF[index]
   )

   l$OutEpisodes=list(ULF=HRVData$FreqAnalysis[[indexFreqAnalysis]]$ULF[-index],
      VLF=HRVData$FreqAnalysis[[indexFreqAnalysis]]$VLF[-index],
      LF=HRVData$FreqAnalysis[[indexFreqAnalysis]]$LF[-index],
      HF=HRVData$FreqAnalysis[[indexFreqAnalysis]]$HF[-index]
   )

   if (HRVData$Verbose) {
      cat("   No. of frames:",lframes,"\n")
      cat("   No. of frames in episodes:",length(l$InEpisodes$ULF),"\n")
      cat("   No. of frames outside episodes:",length(l$OutEpisodes$ULF),"\n")
   }

   return(l)

}

