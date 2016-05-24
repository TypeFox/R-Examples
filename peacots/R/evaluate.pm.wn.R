evaluate.pm.wn <-
function(times, signal, minPeakFreq=0, minFitFreq=0){
	# Calculate Lomb-Scargle periodogram from time series
	spectrum 	= LombScarglePeriodogram(times, signal);
	frequencies	= spectrum$frequencies;
	powers		= spectrum$powers;
	frequencies = frequencies[-1]; 		#discard zero-frequency mode
	powers	 	= powers[-1]; 			#discard zero-frequency mode
	if(tail(frequencies,n=1)<minPeakFreq) return(list(error=TRUE, errorMessage="All periodogram frequencies are below minPeakFreq"));
	if(tail(frequencies,n=1)<minFitFreq) return(list(error=TRUE, errorMessage="All periodogram frequencies are below minFitFreq"));
	minPeakMode	= which(frequencies>=minPeakFreq)[1];	#minimum mode considered for peak search. Does not influence ML fitting.
	minFitMode	= which(frequencies>=minFitFreq)[1];	#minimum mode considered for maximum-likelihood fitting. Does not influence peak search.
	
	# Detect spectral peak
	peakMode 	= (minPeakMode-1) + which.max(powers[minPeakMode:length(powers)]);											
	peakPower 	= powers[peakMode];
	peakFreq 	= frequencies[peakMode];	
	
	#Estimate white noise power
	#Linear least squares is equivalent to taking arithmetic mean, and the same as estimating the variance directly from the time series
	powerEstimate = mean(powers[minFitMode:length(powers)]);	#mean power in signal (proportional to signal variance)
	if(powerEstimate<1E-20) return(list(error=TRUE, errorMessage="Periodogram nearly zero"));
	significanceWN = 1 - (1 - exp(-peakPower/powerEstimate)) ^ (length(powers)-minPeakMode+1);
	
	#Return all results
	return(list(error			= FALSE,
				errorMessage	= "",
				frequencies		= frequencies, 
				periodogram 	= powers, 
				peakMode		= peakMode, 
				powerEstimate	= powerEstimate,
				minPeakMode		= minPeakMode, 
				minFitMode		= minFitMode,
				RSS				= sum((powers-powerEstimate)^2),
				P				= significanceWN));
}
