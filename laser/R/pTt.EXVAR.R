`pTt.EXVAR` <-
function(mu0, z, lam, t0, Tmax){
	if (muFx(mu0, z, Tmax) > lam)
		stop('Error - mu exceeds lambda\n');	fx1 <- function(x){
		 ( (1-exp(-z*x)) * exp(mu0*x + (mu0/z)*exp(-z*x) - lam*x - mu0*t0 - (mu0/z)*exp(-z*x) +lam*t0) );
	}
	temp_fx <- integrate(fx1, t0, Tmax, stop.on.error = FALSE);
	res <- 1/(1+mu0*temp_fx$value); 
	return(res);
}

