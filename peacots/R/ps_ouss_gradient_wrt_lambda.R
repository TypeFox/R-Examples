ps_ouss_gradient_wrt_lambda <-
function(freq, power_o, lambda, power_e, time_step, series_size){
	rho = exp(-lambda*time_step);
	A = exp(complex(real=0,imaginary=+2*pi*freq*time_step));
	B = Conj(A);
	N = series_size;
	DOUdrho = (power_o*T*((2*N + A*((-1 + (rho*A)^N)*(-1 + rho*(2 + rho + A*(-1 + (-2 + rho)*rho))) + 
				N*(-1 - (A*rho)^N + rho*(-4 + A + rho + 2*A*rho - A*rho^2 + 
				(A*rho)^N *(A + rho - A*rho^2)))))/(-1 + A*rho)^3 + 
				(B*((-1 + N)*(-1 + rho*(2 + rho)) + 2*N*rho^3*B^2 + 
				(-1 + rho*(2 + rho) + N*(-1 + rho^2))*(rho*B)^N + 
				rho*B*((-1 + (-2 + rho)*rho)*(-1 + (rho*B)^N) + 
				N*(1 + (rho*B)^N - rho*(4 + rho + rho*(rho*B)^N)))))/
				(-1 + rho*B)^3))/(time_step*N^2*(1 + rho)^2);
	return(Re(DOUdrho * rho * (-time_step)));
}
