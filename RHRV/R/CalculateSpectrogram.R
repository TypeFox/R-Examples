CalculateSpectrogram <-
function(HRVData, size, shift, sizesp=1024, verbose=NULL) {
# ---------------------------------------------------------------
# Calculates the spectrogram of an interpolated heart rate signal
# ---------------------------------------------------------------
#    size, disp: size and displacement of window (sec.)
#    sizesp: points for calculating spectrogram (zero padding)
# Used by PlotSpectrogram and CalculatePowerPerBand

	if (!is.null(verbose)) {
		cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
		SetVerbose(HRVData,verbose)
	}
	
	if (HRVData$Verbose) {
		cat("   Calculating spectrogram\n")
	}
	
	shift=shift*HRVData$Freq_HR
	size=floor(size*HRVData$Freq_HR)
	if (HRVData$Verbose) {
		cat("      Window: ",size," samples (shift ",shift," samples)\n",sep="")
	}
	
	
  if (is.null(sizesp)){
    sizesp = 2^ceiling(log2(size))
  }
  if (sizesp <= size){
    sizesp = size
  }

	sizezp=sizesp-size
	if (HRVData$Verbose) {
		cat("      Window size for calculation: ",sizesp," samples (zero padding: ",sizezp," samples)\n",sep="")
	}
	
  signal=1000.0/(HRVData$HR/60.0)  # msec.
 
	if (HRVData$Verbose) {
	
		cat("      Signal size: ",length(signal)," samples\n",sep="")
	}
	hamming=0.54-0.46*cos(2*pi*(0:(size-1))/(size-1))
  hammingfactor=1.586
	
	# Calculates the number of windows
	nw=1
	begnw=1
	repeat {
    	begnw=begnw+shift
    	if ((begnw+size-1)>=length(signal)) {
      		break
    	}
		nw=nw+1
	}
	
	if (HRVData$Verbose) {
		cat("      Windowing signal... ",nw," windows \n",sep="")
	}
	
	zp=matrix(nrow=nw,ncol=sizesp)
	for (i in 1:nw) {
		beg=1+(shift*(i-1))
		windowedSignal = signal[beg:(beg + size - 1)]
		windowedSignal = windowedSignal - mean(windowedSignal)
    zp[i, ] = c(windowedSignal * hamming, rep(0, length = sizezp))
	}
	
	z=matrix(nrow=nw,ncol=floor(sizesp/2))
	for (i in 1:nw) {
    	f= hammingfactor * abs(fft(zp[i,]))^2
      f = f/(2*(sizesp^2))
    	f=f[1:(length(f)/2)]
    	z[i,]=f
	}
	
	if (HRVData$Verbose) {
		cat("      Spectrogram calculated\n")
	}

	return(z)
}

