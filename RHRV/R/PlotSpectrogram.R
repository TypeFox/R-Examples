PlotSpectrogram <-
function(HRVData, size, shift, sizesp=NULL, freqRange =NULL, scale="linear", verbose=NULL) {
# -----------------
# Plots spectrogram
# -----------------
#    size, disp: size and displacement of window (sec.)
#    sizesp: seconds for calculating spectrogram (zero padding)
#	   scale: linear or logarithmic

	if (!is.null(verbose)) {
		cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
		SetVerbose(HRVData,verbose)
	}
	
	if (HRVData$Verbose) {
    	cat("** Plotting spectrogram **\n")
	}
	
	specgr=CalculateSpectrogram(HRVData,size,shift,sizesp)
  
	if(scale=="logaritmic"){
    specgr=log(specgr)
	}
	
  frequency = seq(from=0,to=HRVData$Freq_HR/2,length.out=ncol(specgr))
  time = seq(from=head(HRVData$Beat$Time,1),to=tail(HRVData$Beat$Time,1),length.out = nrow(specgr))
 
  if (is.null(freqRange)){ freqRange = range(frequency)}
  indx = which(frequency >= freqRange[[1]] & frequency <= freqRange[[2]])
  filled.contour(time,
		frequency[indx],
		specgr[,indx],
		xlab="Time (sec.)", ylab="Frequency (Hz.)", main="Spectrogram of the HR  series",
    color.palette=topo.colors
		#col=gray((256:0)/256)
	)
	if (HRVData$Verbose) {
		cat("   Spectrogram plotted\n")
	}
  return(specgr)
}

