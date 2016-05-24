sigma_to_power_o <-
function(sigma, lambda, time_step){
	return(sigma^2 * time_step *(1+exp(-time_step*lambda))/(1-exp(-time_step*lambda)));
}
