`pTt.SPVAR` <-
function(lam0, k, mu, t0, Tmax){
	if (mu > lambdaFx(lam0, k, Tmax))
		stop('Error - mu exceeds lambda\n')
	
	fx1 <- function(x){
		(exp((mu*(x-t0))  + (lam0/k)*(exp(-k*x) - exp(-k*t0))));
	}
	temp_fx <- integrate(fx1, t0, Tmax, stop.on.error = FALSE);
	res <- (1/(1+mu*temp_fx$value));
	return(res);
}

