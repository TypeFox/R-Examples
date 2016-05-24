`fitSPVAR` <-
function(bt, init = c(2, .2, 0.1))
{
	
	STIMES <- getSpeciationtimes(bt);
	TMAX <- max(bt);
	
	getLam0.SPVAR <- function(x)
	{
		k <- x[2];
		mu0 <- x[3];
		lam0 <- (x[1] + mu0) * exp(TMAX* k); 	
		return(lam0);
	}

	#the likelihood function, constrained version
	optimLH.SPVAR <- function(init) #c(r, k, mu0)
	{
		k <- init[2];
		mu0 <- init[3];
		lam0 <- (init[1] + mu0) * exp(TMAX* k); 	
		LH <- getLikelihood.SPVAR(lam0, k, mu0, STIMES, TMAX);
		return(-LH);
	}	
	
	low <- c(0.001, 0.001, 0.001);
	up <- c(Inf, Inf, Inf);
	
	temp <- optim(init, optimLH.SPVAR, method='L-BFGS-B', lower= low, upper=up );

	res <- list(model = 'SPVAR', LH = -temp$value, aic=(2*temp$value)+6, lam0 = getLam0.SPVAR(temp$par), k=temp$par[2], mu0 = temp$par[3]);
		
	return(res);
}

