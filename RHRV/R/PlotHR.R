PlotHR <-
function(HRVData, Tag=NULL, verbose=NULL) {
# -----------------------------
# Plots interpolated Heart Rate
# -----------------------------
#	Tag -> Tags of episodes to include in the plot
#    "all" includes all types


	if (!is.null(verbose)) {
		cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
		SetVerbose(HRVData,verbose)
	}
	
	if (HRVData$Verbose) {
		cat("** Plotting interpolated instantaneous heart rate **\n");
	}
   
   if (!is.null(Tag) & is.null(HRVData$Episodes)) {
      stop("  --- Episodes not present ---\n    --- Quitting now!! ---\n")
   }

	if (is.null(HRVData$HR)) { 
      stop("  --- Interpolated heart rate not present ---\n    --- Quitting now!! ---\n")
	}
	
	npoints = length(HRVData$HR)
	if (HRVData$Verbose) {
		cat("   Number of points:",npoints,"\n");
	}

	
	first = head(HRVData$Beat$Time,1)
	last = tail(HRVData$Beat$Time,1)
	x=seq(first,last,length.out=npoints)

	HRMin=min(HRVData$HR)
	HRMax=max(HRVData$HR)
	HRDiff=HRMax-HRMin

	if (!is.null(Tag)) {
		if (Tag[1]=="all") {
			Tag=levels(HRVData$Episodes$Type)
		}

		if (HRVData$Verbose) {
			cat("   Episodes in plot:",Tag,"\n")
		}
	}


	plot(x,HRVData$HR,type="l",xlab="time (sec.)",ylab="HR (beats/min.)",ylim=c(HRMin-0.1*HRDiff,HRMax))
	grid()
	
	if (!is.null(Tag)) {
		# Data for representing episodes
		EpisodesAuxLeft=HRVData$Episodes$InitTime[HRVData$Episodes$Type %in% Tag]
		EpisodesAuxBottom=c(HRMin-0.09*HRDiff,HRMin-0.04*HRDiff)
		EpisodesAuxRight=HRVData$Episodes$InitTime[HRVData$Episodes$Type %in% Tag] + 
			HRVData$Episodes$Duration[HRVData$Episodes$Type %in% Tag]
		EpisodesAuxTop=c(HRMin-0.07*HRDiff,HRMin-0.02*HRDiff)
		EpisodesAuxType=HRVData$Episodes$Type[HRVData$Episodes$Type %in% Tag]

		Pal=rainbow(length(Tag))
		Bor=Pal[match(EpisodesAuxType,Tag)]

		cat("   No of episodes:",length(EpisodesAuxLeft),"\n")
		cat("   No of classes of episodes:",length(Pal),"\n")

		rect(EpisodesAuxLeft,EpisodesAuxBottom,EpisodesAuxRight,EpisodesAuxTop,border=Bor,col=Bor)
		for (i in 1:length(EpisodesAuxLeft)) {

			lines(rep(EpisodesAuxLeft[i],times=2),c(HRMin-0.1*HRDiff,HRMax),lty=2,col=Bor[i])
			lines(rep(EpisodesAuxRight[i],times=2),c(HRMin-0.1*HRDiff,HRMax),lty=2,col=Bor[i])
		}

		legend("topright",inset=0.01,legend=Tag,fill=Pal,cex=0.6,horiz=FALSE,bg='white')
	}

	title(main="Interpolated instantaneous heart rate")

}

