significanceOfLocalPeak <-
function(power_o, lambda, power_e, time_step, series_size, Nfreq, peakFreq, peakPower){
	Epower 	= ps_ouss(freq=peakFreq, power_o=power_o, lambda=lambda, power_e=power_e, time_step=time_step, series_size=series_size);
	ratio 	= peakPower/Epower;
	if(is.nan(ratio)) return(NaN);
	return(1 - (1-exp(-ratio))^Nfreq);
}
