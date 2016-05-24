generate_ouss <-
function(	times,		# time points in non-decreasing order
			mu,			# deterministic equilibrium
			power_o,
			sigma,		# standard deviation of OU process (i.e. of its stationary distribution)
			lambda,		# relaxation rate of OU process
			power_e,
			epsilon){	# standard deviation of Gaussian measurement error. Set this to 0 for a classical OU process.
	N = length(times)
	if(N==0) return(c());
	T = times[N]-times[1];
	average_time_step = T/max(1,N-1);	
	if(missing(sigma)){
		# requires power_o
		if(missing(power_o)){ stop("Missing power_o or sigma"); }
		sigma = power_o_to_sigma(power_o, lambda, average_time_step);
	}	
	if(missing(epsilon)){
		# requires power_e
		if(missing(power_e)){ stop("Missing power_e or epsilon"); }
		epsilon = sqrt(power_e/average_time_step);
	}
	
	
	uncorrelated_signal = rnorm(N,mean=0,sd=sigma);
    signal		= rep(0,N);
    signal[1] 	= uncorrelated_signal[1]; # random starting point
    dtimes		= diff(times);
    rhos		= exp(-lambda*dtimes); # correlation between consecutive time steps
    increments 	= sqrt(1-rhos^2)*uncorrelated_signal[2:N];
    
    # simulate OU time series through correlated draws
    for(k in 2:N){
        signal[k] = rhos[k-1]*signal[k-1] + increments[k-1];
	}
	
	# shift uniformly and add Gaussian measurement errors to obtain an OU state space model
	signal = signal + mu + rnorm(N,mean=0,sd=epsilon);
		
	return(signal);
}
