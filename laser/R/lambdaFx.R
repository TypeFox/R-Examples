`lambdaFx` <-
function(lam0, k, time){
	lambda <- lam0 * exp(-k*time);	
	return(lambda);
}

