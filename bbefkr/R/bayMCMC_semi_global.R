bayMCMC_semi_global <-
function(data_x, data_y, data_xnew, Xvar, Xvarpred, warm=1000, M=1000, mutprob=0.44, errorprob=0.44, mutsizp=1.0, errorsizp=1.0, 
	prior_alpha=1.0, prior_beta=0.05, err_int = c(-10,10), err_ngrid=10001, num_batch=20, step=10, alpha=0.95, ...)
{	
	data_y <- as.vector(data_y)	
	 if (is.vector(data_xnew)) 
        		data_xnew <- as.matrix(t(data_xnew))
	    testfordim <- sum(dim(data_x) == dim(data_xnew)) == 2
	    twodatasets <- TRUE
	    if (testfordim) 
	        twodatasets <- sum(data_x == data_xnew) != prod(dim(data_x))
	
	SPECURVES1 = data_x
	Specresp1 = data_y

	########################
	# negative log posterior
	########################

	cost = function(xp)
	{
		n = length(Specresp1)		
		Wmat = funopare.kernel(Specresp1, SPECURVES1, SPECURVES1, bandwidth = exp(xp[1]), ...)$NWweit
		if(is.null(Wmat))
		{	
			result = -100000
		}
		else
		{
			tildacurves = (diag(n) - Wmat) %*% Xvar
			tilday = (diag(n) - Wmat) %*% as.matrix(Specresp1)
			betahat = solve(t(tildacurves)%*%tildacurves)%*%t(tildacurves)%*%tilday
		
			resi = Specresp1 - Xvar%*%betahat
			NW = funopare.kernel(resi, SPECURVES1, SPECURVES1, bandwidth = exp(xp[1]), ...)$Estimated.values
			yhat = NW + Xvar%*%betahat	
			resid = as.vector(Specresp1 - yhat)
			epsilon = as.vector(scale(resid))
			std = sd(resid)
			cont = (2.0*pi)^(-0.5)
		
			b = exp(xp[2])
			logf = vector(,length(resid))
			for(i in 1:length(resid))
			{
				temp = epsilon[i] - epsilon[-i]
				res = sum(cont*exp(-0.5*((temp/b)^2))/b)
				logf[i] = log(res/length(temp)/std)
			}
			sumlogf = sum(logf)
			# log Jacobi and log prior
			priorJacobi = vector(,2)
			for(i in 1:2)
			{
				priorJacobi[i] = xp[i] + logpriorh2((exp(xp[i]))^2, prior_alpha=prior_alpha, prior_beta=prior_beta)
			}
			result = sumlogf + sum(priorJacobi)
		}
		return(-result)
	}
			
	gibbs_mean = function(xp, k, mutsizp)
	{
		fx = xp[3]
		dv = rnorm(1) * mutsizp

		xp[1] = xp[1] + dv
		fy = cost(xp)
		rvalue = fx - fy
		if (class(rvalue) != "numeric") 
		{
			accept = 0
		}
		else
		{
			if(fx > fy) 
			{
				accept = 1
			}
			else
			{
				un = runif(1)
				if(un < exp(rvalue)) accept = 1
				else accept = 0
			}
		}
		accept_mean=0
		mutsizc = mutsizp/(mutprob * (1.0 - mutprob))
		if(accept == 1)
		{
			accept_mean = accept_mean + 1
			xp[3] = fy
			mutsizp = mutsizp + mutsizc * (1.0 - mutprob)/k
		}
		else
		{
			xp[1] = xp[1] - dv
			mutsizp = mutsizp - mutsizc * mutprob/k
		}
		return(list(xpnew = xp, mutsizpnew = mutsizp, mutsizcnew = mutsizc, acceptnw = accept_mean))
	}

	###################################################
	# sampling b used in the kernel-form error density
	###################################################

	gibbs_erro = function(xp, k, errorsizp)
	{
		fx = xp[3]
		dv = rnorm(1) * errorsizp

		xp[2] = xp[2] + dv
		fy = cost(xp)
		rvalue = fx - fy
		if (class(rvalue) != "numeric") {
			accept = 0
		}			
		else
		{
			if(fx > fy)
			{
				accept = 1
			}
			else
			{
				un = runif(1)
				if(un < exp(rvalue)) accept = 1
				else accept = 0
			}
		}
		accept_erro = 0
		errorsizc = errorsizp/(errorprob * (1.0 - errorprob))
		if(accept == 1)
		{
			accept_erro = accept_erro + 1
			xp[3] = fy
			errorsizp = errorsizp + errorsizc * (1.0 - errorprob)/k
		}
		else
		{
			xp[2] = xp[2] - dv
			errorsizp = errorsizp - errorsizc * errorprob/k
		}
		return(list(xperronew = xp, errorsizpnew = errorsizp, errorsizcnew = errorsizc, accepterro = accept_erro))  	
	}

	## warm-up stage		
	# initial values		
	ini_val = runif(2,min=1,max=3)
	xp = c(ini_val, cost(xp = ini_val))
	acceptnw = accepterro = vector(,warm)
	xpwarm = matrix(,warm,3)
	# burn-in period
	for(k in 1:warm)
	{
			dum = gibbs_mean(xp, k, mutsizp)
			xp = dum$xpnew
			acceptnw[k] = dum$acceptnw
			mutsizp = dum$mutsizpnew
			
			dum2 = gibbs_erro(xp, k, errorsizp)
			xp = dum2$xperronew
			accepterro[k] = dum2$accepterro
			errorsizp = dum2$errorsizpnew
			xpwarm[k,] = xp		
	}
	# MCMC recording
	acceptnwMCMC = accepterroMCMC  = vector(,M)
	xpM = xpMsquare = matrix(,M,3)
	cpost = matrix(, M/step, 3)
	for(k in 1:M)
	{
			dumMCMC = gibbs_mean(xp, k+warm, mutsizp)
			xp = dumMCMC$xpnew
			acceptnwMCMC[k] = dumMCMC$acceptnw
			mutsizp = dumMCMC$mutsizpnew
			
			dum2MCMC = gibbs_erro(xp, k+warm, errorsizp)
			xp = dum2MCMC$xperronew
			accepterroMCMC[k] = dum2MCMC$accepterro
			errorsizp = dum2MCMC$errorsizpnew
			
			xpM[k,] = xp	
			xpMsquare[k,] = exp(xp)^2
			index = ceiling(k/step)
			cpost[index,] = exp(xp)^2	
	}
	# ergodic average
	xpfinalres = colMeans(xpM)
	# obtaining the bandwidth of regression and residuals,
	kernelestfinal = funopare.kernel(Specresp1, SPECURVES1, SPECURVES1, bandwidth = exp(xpfinalres[1]), ...)
	Wmatkernelestfinal = kernelestfinal$NWweit							
	n = dim(Wmatkernelestfinal)[1]
		
	tildacurves = (diag(n) - Wmatkernelestfinal) %*% Xvar
	tilday = (diag(n) - Wmatkernelestfinal) %*% as.matrix(Specresp1)
	betahat = solve(t(tildacurves)%*%tildacurves)%*%t(tildacurves)%*%tilday
		
	resi = Specresp1 - Xvar%*%betahat
	NW = funopare.kernel(Response = resi, CURVES = SPECURVES1, PRED = data_xnew, bandwidth = exp(xpfinalres[1]), ...)
	regressionfunctionestimate = (Xvar%*%betahat + NW$Estimated.values)	
	residfinal = Specresp1 - regressionfunctionestimate	
	sif_value = SIF(exp(xpM[,1:(ncol(xpM)-1)]), M, num_batch)
    log_likelihood_Chib = loglikelihood_global_admkr(exp(xpfinalres[1:2]), residfinal)
    log_prior_Chib = logpriors_admkr(exp(xpfinalres[1:2])^2, prior_alpha=prior_alpha, prior_beta=prior_beta)
    log_density_Chib = logdensity_admkr(colMeans(xpMsquare[,1:2]), cpost[,1:2])
    mlikeres = log_likelihood_Chib + log_prior_Chib - log_density_Chib

	y = seq(err_int[1], err_int[2], by = diff(err_int)/(err_ngrid-1))
	fore.den.mkr = fore.cdf.mkr = vector(,length(y))
	for(i in 1:(length(y)))
	{
		eps = y[i]
		fore.den.mkr[i] = error.den(exp(xpfinalres[2]), eps, residfinal)
		fore.cdf.mkr[i] = error.cdf(exp(xpfinalres[2]), eps, residfinal)
	}
	if (twodatasets) 
	{
		NWpred = NW$Predicted.values + Xvarpred%*%betahat					
        lb = NWpred + y[which.min(abs(fore.cdf.mkr - (1-alpha)/2))]
        ub = NWpred + y[which.min(abs(fore.cdf.mkr - (1+alpha)/2))]
        PI = cbind(lb, ub)
		return(list(xpfinalres = exp(xpfinalres[1:2]), mhat = regressionfunctionestimate, betahat = betahat, sif_value = sif_value, mlikeres = mlikeres,
			acceptnwMCMC = mean(acceptnwMCMC), accepterroMCMC = mean(accepterroMCMC), 
			fore.den.mkr = fore.den.mkr, fore.cdf.mkr = fore.cdf.mkr, pointforecast = NWpred, PI = PI))
	}
	else
	{
		return(list(xpfinalres = exp(xpfinalres[1:2]), mhat = regressionfunctionestimate, betahat = betahat, sif_value = sif_value, mlikeres = mlikeres, log_likelihood_Chib = log_likelihood_Chib, log_prior_Chib = log_prior_Chib, 
            log_density_Chib = log_density_Chib, acceptnwMCMC = mean(acceptnwMCMC), accepterroMCMC = mean(accepterroMCMC), 
			fore.den.mkr = fore.den.mkr, fore.cdf.mkr = fore.cdf.mkr))
	}
}

