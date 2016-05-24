power_o_to_sigma <-
function(power_o, lambda, time_step){
	return(sqrt(power_o/(time_step *(1+exp(-time_step*lambda))/(1-exp(-time_step*lambda)))));
}
