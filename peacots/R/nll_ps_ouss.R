nll_ps_ouss <-
function(params, frequencies, powers, time_step, series_size){
	PSexpected = ps_ouss(freq=frequencies, power_o=params[1]^2, lambda=abs(params[2]), power_e=params[3]^2, time_step=time_step, series_size=series_size);
	LL = sum(- log(PSexpected) - powers/PSexpected);
	return(-LL);
}
