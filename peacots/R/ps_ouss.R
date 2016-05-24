ps_ouss = function(freq, power_o, sigma, rho, lambda, power_e, epsilon, time_step, series_size){
	if(missing(freq)){ stop("Missing freq"); }
	if(missing(time_step)){ stop("Missing time_step"); }
	if(missing(rho)){
		# requires lambda
		if(missing(lambda)){ stop("Missing rho or lambda"); }
		rho = exp(-lambda*time_step);
	}else{
		lambda = -log(rho)/time_step;
	}
	N = series_size;
	T = time_step*(N-1);
	if(missing(sigma)){
		# requires power_o
		if(missing(power_o)){ stop("Missing power_o or sigma"); }
		sigma = power_o_to_sigma(power_o, lambda, time_step);
	}	
	if(missing(power_e)){
		# requires epsilon
		if(missing(epsilon)){ stop("Missing power_e or epsilon"); }
		power_e = time_step * epsilon^2;
	}
	r = rho*exp(complex(real=0,imaginary=-2*pi*freq*time_step));
	q = rho*exp(complex(real=0,imaginary=+2*pi*freq*time_step));	
	return(power_e + (T/N^2) * sigma^2 * (Re(N * (1/(1-r) + q/(1-q))) - 2*Re(r*(1-r^N)/(1-r)^2)));
}
