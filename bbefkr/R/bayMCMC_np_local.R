bayMCMC_np_local <-
function(data_x, data_y, data_xnew, warm=1000, M=1000, mutprob=0.44, errorprob=0.44, epsilonprob=0.44, mutsizp=1.0, errorsizp=1.0, epsilonsizp=1.0, prior_alpha=1.0, prior_beta=0.05, err_int = c(-10,10), err_ngrid=10001, num_batch=20, step=10, alpha=0.95, ...)
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
	
	# negative log posterior
	cost = function(xp)
	{
		kernelest = funopare.kernel(Specresp1, SPECURVES1, SPECURVES1, bandwidth = exp(xp[1]), ...)
		if(kernelest$Mse == 0)
		{
			result = -1000000
		}
		else
		{
			resid = Specresp1 - kernelest$Estimated.values
			epsilon = scale(resid)
			std = sd(resid)
			cont = (2.0*pi)^(-0.5)
			b = exp(xp[2])
			badj = exp(xp[3])/(1.0+exp(xp[3]))
			logf = vector(,length(resid))
			for(i in 1:length(resid))
			{
				temp = epsilon[i] - epsilon[-i]
				res = sum(cont*exp(-0.5*((temp/(b*(1+badj*abs(epsilon[-i]))))^2))/(b*(1+badj*abs(epsilon[-i]))))
				logf[i] = log(res/length(temp)/std)
			}
			sumlogf = sum(logf)
			priorJacobi = vector(,2)
			for(i in 1:2)
			{
				priorJacobi[i] = xp[i] + logpriorh2((exp(xp[i]))^2, prior_alpha=prior_alpha, prior_beta=prior_beta)
			}
			result = sumlogf + sum(priorJacobi) + log(badj) + log(1.0 - badj)
		}
		return(-result)
	}
	
	# sampling bandwidth h used in the functional NW estimator
	gibbs_mean = function(xp, k, mutsizp)
	{
		fx = xp[4]
		dv = rnorm(1)*mutsizp
		xp[1] = xp[1] + dv
		fy = cost(xp)
		rvalue = fx - fy
		if(is.nan(rvalue)) accept = 0
		else
		{
			if(fx > fy) accept = 1
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
			xp[4] = fy
			mutsizp = mutsizp + mutsizc*(1.0 -mutprob)/k
		}	
		else
		{
			xp[1] = xp[1] - dv
			mutsizp = mutsizp - mutsizc*mutprob/k
		}		
		return(list(xpnew = xp, mutsizpnew = mutsizp, mutsizcnew = mutsizc, acceptnw = accept_mean)) 
	}
    
	# sampling bandwidth b used in the kernel-form error density
	gibbs_erro = function(xp, k, errorsizp)
	{
		fx = xp[4]
		dv = rnorm(1)*errorsizp
		xp[2] = xp[2] + dv
		fy = cost(xp)
		rvalue = fx - fy
		if(is.nan(rvalue)) accept = 0
		else
		{
			if(fx > fy) accept = 1
			else
			{
				un = runif(1)
				if(un < exp(rvalue)) accept = 1
				else accept = 0
			}
		}
		accept_erro=0
		errorsizc = errorsizp/(errorprob * (1.0 - errorprob))
		if(accept == 1)
		{
			accept_erro = accept_erro + 1
			xp[4] = fy
			errorsizp = errorsizp + errorsizc*(1-errorprob)/k
		}
		else	
		{
			xp[2] = xp[2] - dv
			errorsizp = errorsizp - errorsizc*errorprob/k
		}	
		return(list(xperronew = xp, errorsizpnew = errorsizp, errorsizcnew = errorsizc, accepterro = accept_erro))
	}
    
	# sampling localized bandwidth adjustment tau_{\varepsilon}
	gibbs_epsilon = function(xp, k, epsilonsizp)
	{
		fx = xp[4]
		dv = rnorm(1)*epsilonsizp		
		xp[3] = xp[3] + dv
		fy = cost(xp)
		rvalue = fx - fy
		if(is.nan(rvalue)) accept = 0
		else
		{
			if(fx > fy) accept = 1
			else
			{
				un = runif(1)
				if(un < exp(rvalue)) accept = 1
				else accept = 0
			}
		}
		accept_epsilon = 0
		epsilonsizc = epsilonsizp/(epsilonprob * (1.0 - epsilonprob))
		if(accept == 1)
		{
			accept_epsilon = accept_epsilon + 1
			xp[4] = fy
			epsilonsizp = epsilonsizp + epsilonsizc*(1-epsilonprob)/k
		}
		else
		{
			xp[3] = xp[3] - dv
			epsilonsizp = epsilonsizp - epsilonsizc*epsilonprob/k
		}
		return(list(xpepsilonnew = xp, epsilonsizpnew = epsilonsizp, epsilonsizcnew = epsilonsizc, acceptepsilon = accept_epsilon))
	}
	## warm-up stage		
	# initial values		
	ini_val = runif(3, min=1, max=3)
	xp = c(ini_val, cost(xp = ini_val))
	acceptnw = accepterro = acceptepsilon = vector(,warm)
	xpwarm = matrix(,warm,4)
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
			
			dum3 = gibbs_epsilon(xp, k, epsilonsizp)
			xp = dum3$xpepsilonnew
			acceptepsilon[k] = dum3$acceptepsilon
			epsilonsizp = dum3$epsilonsizpnew
			xpwarm[k,] = xp		
	}
	# MCMC recording
	acceptnwMCMC = accepterroMCMC = acceptepsilonMCMC = vector(,M)
	xpM = xpMsquare = matrix(,M,4)
	cpost = matrix(, M/step, 4)
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
			
			dum3MCMC = gibbs_epsilon(xp, k+warm, epsilonsizp)
			xp = dum3MCMC$xpepsilonnew
			acceptepsilonMCMC[k] = dum3MCMC$acceptepsilon
			epsilonsizp = dum3MCMC$epsilonsizpnew
			xpM[k,] = xp	
			xpMsquare[k,] = exp(xp)^2
			index = ceiling(k/step)
			cpost[index,] = exp(xp)^2	
	}
	# ergodic average
	xpfinalres = colMeans(xpM)
	# obtaining the bandwidth of regression and residuals,
	kernelestfinal = funopare.kernel(Specresp1, SPECURVES1, data_xnew, bandwidth = exp(xpfinalres[1]), ...)
	residfinal = Specresp1 - kernelestfinal$Estimated.values
	sif_value = SIF(exp(xpM[,1:(ncol(xpM)-1)]), M, num_batch)
    log_likelihood_Chib = loglikelihood_local_admkr(exp(xpfinalres[1:3]), residfinal)
    log_prior_Chib = logpriors_admkr(exp(xpfinalres[1:3])^2, prior_alpha=prior_alpha, prior_beta=prior_beta)
    log_density_Chib = logdensity_admkr(colMeans(xpMsquare[,1:3]), cpost[,1:3])
    mlikeres =  log_likelihood_Chib + log_prior_Chib - log_density_Chib

	# approximate ISE	
	y = seq(err_int[1], err_int[2], by = diff(err_int)/(err_ngrid-1))
	fore.den.mkr = fore.cdf.mkr = vector(,length(y))
	for(i in 1:(length(y)))
	{
		eps = y[i]
		fore.den.mkr[i] = error.denadj(exp(xpfinalres[2]), exp(xpfinalres[3]), eps, residfinal)
		fore.cdf.mkr[i] = error.cdfadj(exp(xpfinalres[2]), exp(xpfinalres[3]), eps, residfinal)
	}
	if(twodatasets)
	{
		pointforecast = kernelestfinal$Predicted.values
        lb = pointforecast + y[which.min(abs(fore.cdf.mkr - (1-alpha)/2))]
        ub = pointforecast + y[which.min(abs(fore.cdf.mkr - (1+alpha)/2))]
        PI = cbind(lb, ub)
		return(list(xpfinalres = c(exp(xpfinalres[1:2]), exp(xpfinalres[3])/(1.0+exp(xpfinalres[3]))), mhat = kernelestfinal$Estimated.values, sif_value = sif_value, 
			mlikeres = mlikeres, acceptnwMCMC = mean(acceptnwMCMC), accepterroMCMC = mean(accepterroMCMC), 
			acceptepsilonMCMC = mean(acceptepsilonMCMC), fore.den.mkr = fore.den.mkr, fore.cdf.mkr = fore.cdf.mkr,
			pointforecast = pointforecast, PI = PI))
	}
	else
	{
		return(list(xpfinalres = c(exp(xpfinalres[1:2]), exp(xpfinalres[3])/(1.0+exp(xpfinalres[3]))), mhat = kernelestfinal$Estimated.values, sif_value = sif_value, 
			mlikeres = mlikeres, log_likelihood_Chib = log_likelihood_Chib, log_prior_Chib = log_prior_Chib, log_density_Chib = log_density_Chib,
			acceptnwMCMC = mean(acceptnwMCMC), accepterroMCMC = mean(accepterroMCMC), 
			acceptepsilonMCMC = mean(acceptepsilonMCMC),
			fore.den.mkr = fore.den.mkr, fore.cdf.mkr = fore.cdf.mkr))
	}
}
