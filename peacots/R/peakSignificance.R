peakSignificance <-
function(power_o, lambda, power_e, times, frequencies, peakFreq, peakPower, accuracy){
	
	#Number of trials for estimation of false alarm probability P.
	#The variance of the estimator will be P(1-P)/trials and at most 0.25/trials
	trials 		= max(10,ceiling(0.25/accuracy^2));
	
	#Preliminary preparations
	series_size		= length(times);
	times			= sort(times);
	time_step		= (times[series_size] - times[1])/(series_size-1);
	EpeakPower 		= ps_ouss(peakFreq, power_o=power_o, lambda=lambda, power_e=power_e, time_step=time_step, series_size=series_size);
	Epowers 		= ps_ouss(frequencies, power_o=power_o, lambda=lambda, power_e=power_e, time_step=time_step, series_size=series_size);
	countPositives	= 0;
#	countPositivesLocal = 0;
	
	
	# DEPRECATED: ONLY VALID IN THE INFINITE TIME SERIES LIMIT
	#Generate all exponentially distributed variables corresponding to power estimates in a hypothetical periodogram
	NF = length(frequencies);
	rexps = sapply(1:NF, FUN=function(m){ if(Epowers[m]==0){ return(rep(0,trials)); }else{ return(rexp(n=trials,rate=1/Epowers[m])); } } );
	for(m in 1:length(frequencies)){
		if(Epowers[m]==0){
			rexps[,m] = rep(0, trials);
		}else{
			rexps[,m] = rexp(n=trials, rate=1/Epowers[m]);
		}
	}

	# DEPRECATED: ONLY VALID IN THE INFINITE TIME SERIES LIMIT
	#Evaluate trial periodograms, keeping track of how many ended up having peaks at least as extreme as the case given
	for(n in 1:trials){
		m = which.max(rexps[n,]);		# detect peak in current trial
		r_peakFreq 	= frequencies[m];
		r_peakPower = rexps[n,m];
		if(r_peakPower^2/Epowers[m] >= peakPower^2/EpeakPower){
			countPositives = countPositives + 1;
		}
#		if(r_peakPower/Epowers[m] >= peakPower/EpeakPower){
#			countPositivesLocal = countPositivesLocal + 1;
#		}
	}


	# Alternative:
	# generate OUSS series with given parameters, calculate periodograms and test for the presence of 'more extreme' peaks
#	for(n in 1:trials){
#		spectrum 	= LombScarglePeriodogram(times, generate_ouss(times,0,power_o=power_o,lambda=lambda,power_e=power_e),frequencies);
#		m 			= which.max(spectrum$powers);	# detect peak in current trial periodogram
#		r_peakFreq 	= frequencies[m];
#		r_peakPower = (spectrum$powers)[m];
#		if((r_peakPower^2/Epowers[m] >= peakPower^2/EpeakPower){
#			countPositives = countPositives + 1;
#		}
#		if(r_peakPower/Epowers[m] >= peakPower/EpeakPower){
#			countPositivesLocal = countPositivesLocal + 1;
#		}
#	}	

	#Estimated false alarm probability is fraction of random peridograms with "more extreme" peaks than case given 
#	return(list(P=countPositives/trials, 
#				Plocal=countPositivesLocal/trials));
	return(countPositives/trials)
}
