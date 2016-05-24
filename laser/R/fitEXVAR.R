`fitEXVAR` <-
function(bt, init = c(0.3, .01, 1))
{
	
	STIMES <- getSpeciationtimes(bt);
	TMAX <- max(bt);
	
	
	optimLH.EXVAR <- function(init) #c(r, mu0, z)
	{
		mu0 <- init[2];
		z <- init[3];
		lam0 <- init[1] + mu0*(1-exp(-z*TMAX));
		LH <- getLikelihood.EXVAR(mu0, z, lam0, STIMES, TMAX);
		return(-LH);
	}

	getLam0.EXVAR <- function(x)
	{
		mu0 <- x[2];
		z <- x[3];
		lam0 <- x[1] + mu0*(1-exp(-z*TMAX));  	
		return(lam0);
	}

		
	low <- c(0.001, 0.001, 0.001);
	up <- c(Inf, Inf, Inf);
	temp <- optim(init, optimLH.EXVAR, method='L-BFGS-B', lower= low, upper=up );

	res <- list(model = 'EXVAR', LH = -temp$value, aic=(2*temp$value)+6, lam0 = getLam0.EXVAR(temp$par), mu0=temp$par[2], z = temp$par[3]);
		
	return(res);
}

