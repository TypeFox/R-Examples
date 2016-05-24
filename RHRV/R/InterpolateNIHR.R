InterpolateNIHR <-
function(HRVData, freqhr=4, method="linear", verbose=NULL) {
# -------------------------------------
# Interpolates instantaneous heart rate
# -------------------------------------
#	freqhr -> frequency for interpolating heart rate
#	method -> "linear" interpolation or "spline" monotone interpolation

	if (!is.null(verbose)) {
		cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
		SetVerbose(HRVData,verbose)
	}
	
	if (HRVData$Verbose) {
		cat("** Interpolating instantaneous heart rate **\n");
		cat("   Frequency: ",freqhr,"Hz\n",sep="")
	}
	
	HRVData$Freq_HR=freqhr

	first = head(HRVData$Beat$Time,1)
	last = tail(HRVData$Beat$Time,1)
 	npoints=as.integer((last-first)*HRVData$Freq_HR+1)

	if (HRVData$Verbose) {
		cat("   Number of beats:",length(HRVData$Beat$niHR),"\n")
		cat("   Number of points:",npoints,"\n");
	}

  	if (method=="linear")
		fun=approxfun(HRVData$Beat$Time,HRVData$Beat$niHR,method="linear",ties="ordered")
	else 
		fun=splinefun(HRVData$Beat$Time,HRVData$Beat$niHR,method="monoH.FC",ties="ordered")

	vectorxint=seq(first,last,1/HRVData$Freq_HR)
    
	HRVData$HR=fun(vectorxint)

	# limit indicates the maximum width in seconds of an interval without beats
	limit = 30 # seconds
	# beg and end are the beginning and end of the interval without beats
	begindex = which(diff(HRVData$Beat$Time)>limit)
	beg = HRVData$Beat$Time[begindex]
	end = HRVData$Beat$Time[begindex+1]

	# the value of HR in these intervals is set to zero
	if (length(begindex) > 0) {
		for (i in 1:length(beg)) {
			HRVData$HR[vectorxint>beg[i] & vectorxint<end[i]] = 0
			cat("  --- Warning: interval without beats detected!! ---\n")
		}
	}

	return(HRVData)
}

