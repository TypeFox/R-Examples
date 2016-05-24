ps_ouss_gradient_wrt_power_o <-
function(freq, power_o, lambda, power_e, time_step, series_size){
	rho		= exp(-lambda*time_step);
	N 		= series_size;
	T		= time_step*(N-1);
	r 		= rho*exp(complex(real=0,imaginary=-2*pi*freq*time_step));
	q 		= rho*exp(complex(real=0,imaginary=+2*pi*freq*time_step));
	return(Re(T*(1-rho)/(N^2*time_step*(1+rho)) * (N/(1-r) + N*q/(1-q) - 2*Re(r*(1-r^N)/(1-r)^2))));
}
