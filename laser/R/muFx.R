`muFx` <-
function(mu0, z, time){
	mu <- mu0*(1-exp(-z*time));
	return(mu);	
}

