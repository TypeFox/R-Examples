`fitBOTHVAR` <-
function(bt, init = c(0.3, .5, .1, .5))
{
	
	STIMES <- getSpeciationtimes(bt);
	TMAX <- max(bt);
	
	getLam0.BOTHVAR <- function(x)
	{
		k <- x[2];
		mu0 <- x[3];
		z <- x[4];
		lam0 <- (x[1] + mu0*(1-exp(-z*TMAX)))*exp(TMAX*k);
		return(lam0);
	}
	
	optimLH.BOTHVAR <- function(init) #c(r, k, mu0, z)
	{
		k <- init[2];
		mu0 <- init[3];
		z <- init[4];
		lam0 <- (init[1] + mu0*(1-exp(-z*TMAX)))*exp(TMAX*k);
		LH <- getLikelihood.BOTHVAR(lam0, k, mu0, z, STIMES, TMAX);
		return(-LH);
	}
	
	low <- c(0.001, 0.001, 0.001, .001);
	up <- c(Inf, Inf, Inf, Inf);
	
	init <- c(0.3, .5, .1, .5);
	temp <- optim(init, optimLH.BOTHVAR, method='L-BFGS-B', lower= low, upper=up );
	
	
	res <- list(model = 'BOTHVAR', LH = -temp$value, aic=(2*temp$value)+8, lam0 = getLam0.BOTHVAR(temp$par), k = temp$par[2], mu0=temp$par[3], z = temp$par[4]);
		
	return(res);
}

