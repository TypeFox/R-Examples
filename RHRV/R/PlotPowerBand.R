PlotPowerBand <-
function(HRVData, indexFreqAnalysis = length(HRVData$FreqAnalysis),  normalized=FALSE, hr=FALSE, ymax=NULL, ymaxratio=NULL, ymaxnorm=1, Tag=NULL, verbose=NULL) {
# --------------------
# Plots power per band
# --------------------
#  indexFreqAnalysis: index of an existing frequency analysis to use
# 	normalized: plots normalized powers if TRUE
# 	hr: plots heart rate signal if TRUE
# 	ymax: maximum value for y axis (unnormalized plots)
# 	ymaxratio: maximum value for y axis in LF/HF band (normalized and unnormalized plots)
#	Tag -> Tags of episodes to include in the plot
#    "all" includes all types

	if (!is.null(verbose)) {
		cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
		SetVerbose(HRVData,verbose)
	}
	
	if (HRVData$Verbose) {
		cat("** Plotting power per band **\n")
	}

   if ((length(HRVData$FreqAnalysis) < indexFreqAnalysis) || (indexFreqAnalysis<1) ) {
		stop("  --- Frequency analysis no.",indexFreqAnalysis,"not present!! ---\n    --- Quitting now!! ---\n")
   }

   if (!is.null(Tag) & is.null(HRVData$Episodes)) {
		stop("  --- Episodes not present!! ---\n    --- Quitting now!! ---\n")
   }
	
	if (is.null(HRVData$FreqAnalysis[[indexFreqAnalysis]]$ULF)) {
		cat("   --- ERROR: Power per band not present!! ---\n")
		return(HRVData)
	}

	

	# if (is.null(ymax)) {
	# 	ymax = max(max(HRVData$FreqAnalysis[[indexFreqAnalysis]]$ULF) , max(HRVData$FreqAnalysis[[indexFreqAnalysis]]$VLF) , max(HRVData$FreqAnalysis[[indexFreqAnalysis]]$LF) , max(HRVData$FreqAnalysis[[indexFreqAnalysis]]$HF) )
	# }
	
	 
	# normalization
	 if (normalized) {
		HRVData$FreqAnalysis[[indexFreqAnalysis]]$ULF=(HRVData$FreqAnalysis[[indexFreqAnalysis]]$ULF-min(HRVData$FreqAnalysis[[indexFreqAnalysis]]$ULF))/(max(HRVData$FreqAnalysis[[indexFreqAnalysis]]$ULF)-min(HRVData$FreqAnalysis[[indexFreqAnalysis]]$ULF))
		HRVData$FreqAnalysis[[indexFreqAnalysis]]$VLF=(HRVData$FreqAnalysis[[indexFreqAnalysis]]$VLF-min(HRVData$FreqAnalysis[[indexFreqAnalysis]]$VLF))/(max(HRVData$FreqAnalysis[[indexFreqAnalysis]]$VLF)-min(HRVData$FreqAnalysis[[indexFreqAnalysis]]$VLF))
		HRVData$FreqAnalysis[[indexFreqAnalysis]]$LF=(HRVData$FreqAnalysis[[indexFreqAnalysis]]$LF-min(HRVData$FreqAnalysis[[indexFreqAnalysis]]$LF))/(max(HRVData$FreqAnalysis[[indexFreqAnalysis]]$LF)-min(HRVData$FreqAnalysis[[indexFreqAnalysis]]$LF))
		HRVData$FreqAnalysis[[indexFreqAnalysis]]$HF=(HRVData$FreqAnalysis[[indexFreqAnalysis]]$HF-min(HRVData$FreqAnalysis[[indexFreqAnalysis]]$HF))/(max(HRVData$FreqAnalysis[[indexFreqAnalysis]]$HF)-min(HRVData$FreqAnalysis[[indexFreqAnalysis]]$HF))
		# HRVData$FreqAnalysis[[indexFreqAnalysis]]$LFHF=(HRVData$FreqAnalysis[[indexFreqAnalysis]]$LFHF-min(HRVData$FreqAnalysis[[indexFreqAnalysis]]$LFHF))/(max(HRVData$FreqAnalysis[[indexFreqAnalysis]]$LFHF)-min(HRVData$FreqAnalysis[[indexFreqAnalysis]]$LFHF))
		if (HRVData$Verbose) {
			cat("   Power per band normalized\n")
 		}
	}
			
	
	if (hr)
		numfilas=6
	else
		numfilas=5


	# lframes is the number of frames for plotting power per band
	lframes=length(HRVData$FreqAnalysis[[indexFreqAnalysis]]$HRV)

   # Episodes
   if (!is.null(Tag)) {

       if (Tag[1]=="all") {
         Tag=levels(HRVData$Episodes$Type)
      }

      if (HRVData$Verbose) {
         cat("   Episodes in plot:",Tag,"\n")
      }

      # Data for representing episodes
      EpisodesAuxLeft=HRVData$Episodes$InitTime[HRVData$Episodes$Type %in% Tag]
      EpisodesAuxLeftFrame=EpisodesAuxLeft*lframes/(tail(HRVData$Beat$Time,1)-head(HRVData$Beat$Time,1)) # Beg of episodes (frames)
      EpisodesAuxRight=HRVData$Episodes$InitTime[HRVData$Episodes$Type %in% Tag] + 
         HRVData$Episodes$Duration[HRVData$Episodes$Type %in% Tag]
      EpisodesAuxRightFrame=EpisodesAuxRight*lframes/(tail(HRVData$Beat$Time,1)-head(HRVData$Beat$Time,1)) # Beg of episodes (frames)
      EpisodesAuxType=HRVData$Episodes$Type[HRVData$Episodes$Type %in% Tag]
      if (HRVData$Verbose) {
      	cat("   No of episodes:",length(EpisodesAuxLeft),"\n")
      }
      
      Pal=rainbow(length(Tag))
      Bor=Pal[match(EpisodesAuxType,Tag)]


      EpisodesLeft=HRVData$Episodes$InitTime # Beg of episodes (seconds)
      EpisodesLeftFrame=EpisodesLeft*lframes/(tail(HRVData$Beat$Time,1)-head(HRVData$Beat$Time,1)) # Beg of episodes (frames)
      EpisodesRight=HRVData$Episodes$InitTime+HRVData$Episodes$Duration # Beg of episodes (seconds)
      EpisodesRightFrame=EpisodesRight*lframes/(tail(HRVData$Beat$Time,1)-head(HRVData$Beat$Time,1)) # Beg of episodes (frames)
   }
  previousPar = par()
	par(mfrow=c(numfilas,1),omi=c(0,0,0,0),mai=c(0,0,0,0),mar=c(2,4,1,1),oma=c(1,0,2,0),mgp=c(1.5,.5,0))

# ---------- LF/HF ----------
	if (is.null(ymaxratio)) {
		ymaxratio = max(HRVData$FreqAnalysis[[indexFreqAnalysis]]$LFHF)
	}

	mfg=c(1,1,numfilas,1)
	plot(seq(from=0,to=lframes,length.out=length(HRVData$FreqAnalysis[[indexFreqAnalysis]]$HRV)),HRVData$FreqAnalysis[[indexFreqAnalysis]]$LFHF,type='l',xlab="",ylab="LF/HF",ylim=c(0,ymaxratio*1.1))
	if (!is.null(Tag)) {
		EpisodesAuxTop=c(ymaxratio*1.09,ymaxratio*1.04)
		EpisodesAuxBottom=c(ymaxratio*1.06,ymaxratio*1.01)
		rect(EpisodesAuxLeftFrame,EpisodesAuxBottom,EpisodesAuxRightFrame,EpisodesAuxTop,border=Bor,col=Bor)

		for (i in 1:length(EpisodesAuxLeftFrame)) {
			lines(rep(EpisodesAuxLeftFrame[i],times=2),c(0,ymaxratio*1.1),lty=2,col=Bor[i])
			lines(rep(EpisodesAuxRightFrame[i],times=2),c(0,ymaxratio*1.1),lty=2,col=Bor[i])
		}


		par(xpd=NA) 
		legend(lframes/2,ymaxratio,legend=Tag,fill=Pal,cex=0.9,ncol=length(Tag),xjust=0.5,yjust=-0.2,bty="n")

	}
	if (HRVData$Verbose) {
		cat("   Plotted LF/HF\n")
	}

# ---------- ULF ----------
	if (normalized==TRUE) {
		ymaxv=c(0,ymaxnorm)
	} else if (!is.null(ymax)) {
		ymaxv=c(0,ymax)
	} else {
		ymaxv = c(0,max(HRVData$FreqAnalysis[[indexFreqAnalysis]]$ULF))
	}

	mfg=c(1,2,numfilas,1)
	plot(seq(from=0,to=lframes,length.out=length(HRVData$FreqAnalysis[[indexFreqAnalysis]]$HRV)),
			HRVData$FreqAnalysis[[indexFreqAnalysis]]$ULF,type='l',xlab="",ylab="ULF",ylim=ymaxv)
	if (!is.null(Tag)) {
		for (i in 1:length(EpisodesAuxLeftFrame)) {
			lines(rep(EpisodesAuxLeftFrame[i],times=2),c(ymaxv[1],ymaxv[2]),lty=2,col=Bor[i])
			lines(rep(EpisodesAuxRightFrame[i],times=2),c(ymaxv[1],ymaxv[2]),lty=2,col=Bor[i])
		}
	}
	if (HRVData$Verbose) {
		cat("   Plotted ULF\n")
	}

# ---------- VLF ----------
	if (normalized==TRUE) {
		ymaxv=c(0,ymaxnorm)
	} else if (!is.null(ymax)) {
		ymaxv=c(0,ymax)
	} else {
		ymaxv = c(0,max(HRVData$FreqAnalysis[[indexFreqAnalysis]]$VLF))
	}

	mfg=c(1,3,numfilas,1)
	plot(seq(from=0,to=lframes,length.out=length(HRVData$FreqAnalysis[[indexFreqAnalysis]]$HRV)),
			HRVData$FreqAnalysis[[indexFreqAnalysis]]$VLF,type='l',xlab="",ylab="VLF",ylim=ymaxv)
	if (!is.null(Tag)) {
		for (i in 1:length(EpisodesAuxLeftFrame)) {
			lines(rep(EpisodesAuxLeftFrame[i],times=2),c(ymaxv[1],ymaxv[2]),lty=2,col=Bor[i])
			lines(rep(EpisodesAuxRightFrame[i],times=2),c(ymaxv[1],ymaxv[2]),lty=2,col=Bor[i])
		}
	}
	if (HRVData$Verbose) {
		cat("   Plotted VLF\n")
	}

# ---------- LF ----------
	if (normalized==TRUE) {
		ymaxv=c(0,ymaxnorm)
	} else if (!is.null(ymax)) {
		ymaxv=c(0,ymax)
	} else {
		ymaxv = c(0,max(HRVData$FreqAnalysis[[indexFreqAnalysis]]$LF))
	}

	mfg=c(1,4,numfilas,1)
	plot(seq(from=0,to=lframes,length.out=length(HRVData$FreqAnalysis[[indexFreqAnalysis]]$HRV)),
			HRVData$FreqAnalysis[[indexFreqAnalysis]]$LF,type='l',xlab="",ylab="LF",ylim=ymaxv)
	if (!is.null(Tag)) {
		for (i in 1:length(EpisodesAuxLeftFrame)) {
			lines(rep(EpisodesAuxLeftFrame[i],times=2),c(ymaxv[1],ymaxv[2]),lty=2,col=Bor[i])
			lines(rep(EpisodesAuxRightFrame[i],times=2),c(ymaxv[1],ymaxv[2]),lty=2,col=Bor[i])
		}
	}
	if (HRVData$Verbose) {
		cat("   Plotted LF\n")
	}

# ---------- HF ----------
	if (normalized==TRUE) {
		ymaxv=c(0,ymaxnorm)
	} else if (!is.null(ymax)) {
		ymaxv=c(0,ymax)
	} else {
		ymaxv = c(0,max(HRVData$FreqAnalysis[[indexFreqAnalysis]]$HF))
	}

	mfg=c(1,5,numfilas,1)
	texto4="No. of frames"
	plot(seq(from=0,to=lframes,length.out=length(HRVData$FreqAnalysis[[indexFreqAnalysis]]$HRV)),
			HRVData$FreqAnalysis[[indexFreqAnalysis]]$HF,type='l',xlab=texto4,ylab="HF",ylim=ymaxv)
	if (!is.null(Tag)) {
		for (i in 1:length(EpisodesAuxLeftFrame)) {
			lines(rep(EpisodesAuxLeftFrame[i],times=2),c(ymaxv[1],ymaxv[2]),lty=2,col=Bor[i])
			lines(rep(EpisodesAuxRightFrame[i],times=2),c(ymaxv[1],ymaxv[2]),lty=2,col=Bor[i])
		}
	}
	if (HRVData$Verbose) {
		cat("   Plotted HF\n")
	} 

# ---------- HR ----------
	if (numfilas==6) {
		mfg=c(1,6,numfilas,1)
		# lsecs is the duration of the record in seconds for plotting heart rate signal
		lsecs=tail(HRVData$Beat$Time,1)-head(HRVData$Beat$Time,1)
		plot(seq(from=0,to=lsecs,length.out=length(HRVData$HR)),
				 HRVData$HR,type='l',xlab="Time (sec.)",ylab="HR (bps)")
		if (HRVData$Verbose) {
			cat("   Plotted HRV\n")
		}

		if (!is.null(Tag)) {
		for (i in 1:length(EpisodesAuxLeftFrame)) {
			lines(rep(EpisodesAuxLeft[i],times=2),c(min(HRVData$HR),max(HRVData$HR)),lty=2,col=Bor[i])
			lines(rep(EpisodesAuxRight[i],times=2),c(min(HRVData$HR),max(HRVData$HR)),lty=2,col=Bor[i])
		}
			#rect(EpisodesAuxLeft,rep(min(HRVData$HR),times=length(EpisodesAuxLeft)),EpisodesAuxRight,rep(max(HRVData$HR),times=length(EpisodesAuxLeft)),border=Bor)
		}
# --------------------
	}

	if (HRVData$Verbose & !is.null(Tag)) {
		cat("   Episodes plotted\n")
	}


	if ((normalized == TRUE) && (numfilas == 6)) 
		title(main="Normalized power bands of HRV and Heart Rate Signal",outer=TRUE)
	else if ((normalized == FALSE) && (numfilas == 6))
		title(main="Power bands of HRV and Heart Rate Signal",outer=TRUE)
	else if ((normalized == TRUE) && (numfilas != 6))
		title(main="Normalized power bands of HRV",outer=TRUE)
	else
		title(main="Power bands of HRV",outer=TRUE)
	
	if (HRVData$Verbose) {
		cat("   Power per band plotted\n")
	}	 

  #restore previous graphical parameters
  par(mfrow=previousPar$mfrow, omi=previousPar$omi, mai=previousPar$mai,
      mar=previousPar$mar,oma=previousPar$oma,mgp=previousPar$mgp)

}

