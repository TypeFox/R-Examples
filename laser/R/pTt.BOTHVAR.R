`pTt.BOTHVAR` <-
function(lam0, k, mu0, z, t0, Tmax){
	if (muFx(mu0, z, Tmax) > lambdaFx(lam0, k, Tmax))
		stop('Error - mu exceeds lambda\n');
	FXrho <- function(x){
		( (mu0*x) + (mu0/z)*exp(-z*x) + (lam0/k)*exp(-k*x) );
	}
	fx1 <- function(x){
		( (1-exp(-z*x)) * exp(FXrho(x) - FXrho(t0)) );
	}
	temp_fx <- integrate(fx1, t0, Tmax, stop.on.error = FALSE);
	res <- (1/(1+mu0*temp_fx$value));
	return(res);	
}

