nll_ps_ouss_gradient <-
function(params, frequencies, powers, time_step, series_size){
	PSexpected 	= ps_ouss(frequencies, power_o=params[1]^2, lambda=abs(params[2]), power_e=params[3]^2, time_step=time_step, series_size=series_size);
	alpha 		= -1/PSexpected + powers/PSexpected^2;
	LL_grad		= c(0,0,0);
	LL_grad[1] 	= alpha %*% (2*params[1] 		* ps_ouss_gradient_wrt_power_o(frequencies, power_o=params[1]^2, lambda=abs(params[2]), power_e=params[3]^2, time_step=time_step, series_size=series_size));
	LL_grad[2]	= alpha %*% (sign(params[2]) 	* ps_ouss_gradient_wrt_lambda(frequencies, power_o=params[1]^2, lambda=abs(params[2]), power_e=params[3]^2, time_step=time_step, series_size=series_size));
	LL_grad[3]	= alpha %*% (2*params[3]		* ps_ouss_gradient_wrt_power_e(frequencies, power_o=params[1]^2, lambda=abs(params[2]), power_e=params[3]^2, time_step=time_step, series_size=series_size));
	return(-LL_grad);
}
