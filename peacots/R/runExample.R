runExample <-
function(){
	# Parameters
	N 			= 100; 	# length of time series
	duration 	= 20;	# duration of time series
	lambda		= 1;	# Resilience (inverse correlation time) of OU process
	sigma		= 0.8;	# Standard deviation of OU process
	epsilon		= 0.5;	# Standard deviation of measurement errors
	cycleT		= 1;	# Cycle period for the cyclic case
	cycleA		= 0.6;	# Amplitude of cyclic equilibrium
	times 		= seq(0,N-1) * duration/(N-1);
	
	old.par <- par(mfrow=c(3, 2))	


	
	# Example 1 - Non cyclic time series
	
	cat(sprintf("Example 01: Generating random non-cyclic time series..\n"));
		
	# generate non-cyclic time series
	signal 	= generate_ouss(times, mu=0, sigma=sigma, lambda=lambda, epsilon=epsilon);

	#plot time series
	plot(ts(times), ts(signal), xy.label=FALSE, type="l", ylab="signal", xlab="time", main="Time series (non-cyclic)", cex=0.8, cex.main=0.9)

	cat(sprintf("            Evaluating periodogram..\n"));

	# find peak and estimate statistical significance
	report 	= peacots::evaluate.pm(times=times, signal=signal, minPeakFreq=0, minFitFreq=0, accuracy=0.01, startRadius=2, verbose=FALSE);

	#plot periodogram
	plot(ts(report$frequencies), ts(report$periodogram), xy.label=FALSE, type="l", ylab="power", xlab="frequency", main=sprintf("peacots periodogram analysis\n(peak freq=%.3g, P=%.2g, Plocal=%.2g)",report$frequencies[report$peakMode],report$P,report$Plocal), col="black", cex=0.8, cex.main=0.9);
	
	#plot fitted OUSS periodogram
	lines(report$frequencies[report$minFitMode:length(report$frequencies)], report$fittedPS[report$minFitMode:length(report$fittedPS)], col="red");
	
	#plot legend
	legend((0.6*report$frequencies[1]+0.4*tail(report$frequencies,1)), (0.85*max(report$periodogram)), c("periodogram", "fitted OUSS power spectrum"), lty=c(1,1), col=c("black", "red"), bty="n", cex=0.8)
	
	
	
	# Example 2a - Cyclic time series
	# In this example we use low-frequency trimming to avoid the low-frequency maximum

	cat(sprintf("Example 02: Generating cyclic time series..\n"));
		
	# generate cyclic time series
	signal 	= cycleA * cos(2*pi*times/cycleT) + generate_ouss(times, mu=0, sigma=sigma, lambda=lambda, epsilon=epsilon);
		
	#plot time series
	plot(ts(times), ts(signal), xy.label=FALSE, type="l", ylab="signal", xlab="time", main="Time series (cyclic)", cex=0.8, cex.main=0.9)

	cat(sprintf("            Evaluating periodogram..\n"));

	# find peak and estimate statistical significance
	# ignore frequencies lower than a pre-defined threshold to avoid masking by low-frequency maxima
	minFreq = 0.5/cycleT;
	report 	= peacots::evaluate.pm(times=times, signal=signal, minPeakFreq=minFreq, minFitFreq=minFreq, startRadius=2, accuracy=0.01, verbose=FALSE);
	
	#plot periodogram
	plot(ts(report$frequencies), ts(report$periodogram), xy.label=FALSE, type="l", ylab="power", xlab="frequency", main=sprintf("peacots periodogram analysis\nusing low-frequency trimming at threshold %.2g\n(peak freq=%.3g, P=%.2g, Plocal=%.2g)",minFreq,report$frequencies[report$peakMode],report$P,report$Plocal), col="black", cex=0.8, cex.main=0.9);
	
	#plot fitted OUSS periodogram
	lines(report$frequencies[report$minFitMode:length(report$frequencies)], report$fittedPS[report$minFitMode:length(report$fittedPS)], col="red");
	
	#plot legend
	legend((0.6*report$frequencies[1]+0.4*tail(report$frequencies,1)), (0.85*max(report$periodogram)), c("periodogram", "fitted OUSS power spectrum"), lty=c(1,1), col=c("black", "red"), bty="n", cex=0.8)
	
	
	
	
	
	# Example 2b - Cyclic time series
	# In this example we don't use low-frequency trimming
	# Instead, we pick the periodogram peak of interest and evaluate its local significance

	cat(sprintf("Example 03: Re-evaluating periodogram of previous example..\n"));
	
	# calculate periodogram and fit OUSS model
	report 		= peacots::evaluate.pm(times=times, signal=signal, minPeakFreq=0, minFitFreq=0, startRadius=2, accuracy=0.01, verbose=FALSE);
	
	# find which periodogram mode approximately corresponds to the frequency we are interested in
	cycleMode 	= which(report$frequencies>=0.99/cycleT)[1]; 
	
	# calculate P-value for local peak
	Pvalue 		= significanceOfLocalPeak(	power_o		= report$power_o, 
											lambda		= report$lambda, 
											power_e		= report$power_e, 
											time_step	= report$time_step,
											series_size = length(times),
											Nfreq		= length(report$frequencies), 
											peakFreq	= report$frequencies[cycleMode], 
											peakPower	= report$periodogram[cycleMode]);

	#plot time series
	plot(ts(times), ts(signal), xy.label=FALSE, type="l", ylab="signal", xlab="time", main="Time series (cyclic)", cex=0.8, cex.main=0.9)

	#plot periodogram
	plot(ts(report$frequencies), ts(report$periodogram), xy.label=FALSE, type="l", ylab="power", xlab="frequency", main=sprintf("peacots periodogram analysis\nfocusing on peak at freq=%.3g\nPlocal=%.2g",report$frequencies[cycleMode],Pvalue), col="black", cex=0.8, cex.main=0.9);
	
	#plot fitted OUSS periodogram
	lines(report$frequencies[report$minFitMode:length(report$frequencies)], report$fittedPS[report$minFitMode:length(report$fittedPS)], col="red");
	
	#plot legend
	legend((0.6*report$frequencies[1]+0.4*tail(report$frequencies,1)), (0.85*max(report$periodogram)), c("periodogram", "fitted OUSS"), lty=c(1,1), col=c("black", "red"), bty="n", cex=0.8)

	par(old.par)
}
