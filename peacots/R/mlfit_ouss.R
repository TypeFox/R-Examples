mlfit_ouss <-
function(frequencies, periodogram, time_step, series_size, startRadius=10, iterations=200){
	power_e_scales		= mean(tail(periodogram,n=min(5,length(periodogram))));		# first guess for power_e
	powerInt 			= aux_trapezoid(frequencies, periodogram); 	# total (integrated) power in signal
	power_o_scales 		= mean(periodogram[1:min(5,length(periodogram))]);			# first guess for power_o
	lambda_scales 		= 4*powerInt/power_o_scales;	# first guess for lambda	
	
	# Define a grid of starting points for the fitting
	startRadius			= abs(startRadius);
	suspected_power_o 	= (4^(-startRadius:startRadius)) * power_o_scales;	# try out a range of start values for the fit
	suspected_lambda 	= (4^(-startRadius:startRadius)) * lambda_scales; 	# try out a range of start values for the fit
	suspected_power_e	= c(power_e_scales);	# power_e is a uniform translation parameter and so very easy to deal with by the optimizer
	
	# Go through all these starting points and attempt fit. Will choose the best fit at the end.
	ML_value = c(); ML_power_o =c(); ML_lambda =c(); ML_power_e =c();
	for(start_power_o in suspected_power_o){
		for(start_lambda in suspected_lambda){
			for(start_power_e in suspected_power_e){
				#try to fit theoretical spectrum of OU process to periodogram. This might fail due to numerical reasons.
				#	fitting parameters are params[1]^2=power_o, abs(params[2]) = lambda, params[3]^2 = power_e
				#	these transformations are used to keep the parameters unbounded (i.e. use unconstrained optimization)
				start_params = c(start_power_o=sqrt(start_power_o), lambda=start_lambda, power_e=sqrt(start_power_e));
				ML = try(stats::optim(	par=start_params, fn=nll_ps_ouss, gr=nll_ps_ouss_gradient,
										method="BFGS",
										control=list(	fnscale=+1,
														parscale=c(sqrt(power_o_scales),lambda_scales,sqrt(power_e_scales)), 
														maxit=iterations, reltol=1e-6),
										frequencies=frequencies, powers=periodogram, time_step=time_step, series_size=series_size), 
							silent=TRUE);
				if(!("try-error" %in% class(ML)) && (ML$convergence==0)){
					ML_value 	= c(ML_value, ML$value);	# note that ML_value is actually the negated maximum log-likelihood
					fits		= ML$par;
					ML_value 	= c(ML_value, ML$objective);
					ML_power_o	= c(ML_power_o, (fits[1])^2);
					ML_lambda 	= c(ML_lambda, abs(fits[2]));
					ML_power_e	= c(ML_power_e, (fits[3])^2);
				}else{
					next;
				}
			}
		}
	}
	if(length(ML_value)==0){
		# none of the fits worked
		return(list(error=TRUE, errorMessage="Maximum-likelihood fit failed"));
	}else{
		# pick best fit
		bestFit 		= which.min(ML_value);	
		fitted_LL		= -ML_value[bestFit];
		fitted_power_o	= ML_power_o[bestFit];
		fitted_lambda	= ML_lambda[bestFit];
		fitted_power_e	= ML_power_e[bestFit];
	}
	
	return(list(error			=FALSE,
				errorMessage	= "",
				power_o			= fitted_power_o, 
				lambda			= fitted_lambda,
				power_e			= fitted_power_e,
				MLL				= fitted_LL));
}
