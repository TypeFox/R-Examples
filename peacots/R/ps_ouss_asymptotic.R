ps_ouss_asymptotic <-
function(freq, power_o, sigma, rho, lambda, power_e, epsilon, time_step){
	if(missing(freq)){ stop("Missing freq"); }
	if(missing(time_step)){ stop("Missing time_step"); }
	if(missing(rho)){
		# requires lambda
		if(missing(lambda)){ stop("Missing rho or lambda"); }
		rho = exp(-lambda*time_step);
	}else{
		lambda = -log(rho)/time_step;
	}
	if(missing(power_o)){
		# requires sigma
		if(missing(sigma)){ stop("Missing power_o or sigma"); }
		power_o = sigma_to_power_o(sigma, lambda, time_step);
	}
	if(missing(power_e)){
		# requires epsilon
		if(missing(epsilon)){ stop("Missing power_e or epsilon"); }
		power_e = time_step * epsilon^2;
	}
	return(power_o*(1-rho)^2/(1+rho^2-2*rho*cos(2*pi*freq*time_step)) + power_e);
}
