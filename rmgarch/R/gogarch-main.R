#################################################################################
##
##   R package rmgarch by Alexios Ghalanos Copyright (C) 2008-2013.
##   This file is part of the R package rmgarch.
##
##   The R package rmgarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rmgarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################


#-----------------------------------------------------------------------------------------------

.gogarchfit = function(spec, data, out.sample = 0, solver = "solnp",
		fit.control = list(stationarity = 1), solver.control = list(), 
		cluster = NULL, VAR.fit = NULL, ARcoef = NULL, ...)
{
	tic = Sys.time()
	model = spec@model
	umodel = spec@umodel
	ufit.control = list()
	if(is.null(fit.control$stationarity)){
		ufit.control$stationarity = TRUE 
	} else {
		ufit.control$stationarity = fit.control$stationarity
		fit.control$stationarity = NULL
	}
	if(is.null(fit.control$scale)){
		ufit.control$scale = TRUE 
	} else{
		ufit.control$scale = fit.control$scale
		fit.control$scale = NULL
	}
	#-----------------------------------------------------------------------------------
	# Data Extraction
	m = dim(data)[2]
	if( is.null( colnames(data) ) ) cnames = paste("Asset_", 1:m, sep = "") else cnames = colnames(data)
	
	xdata = .extractmdata(data)
	if(!is.numeric(out.sample)) 
		stop("\ngogarchfit-->error: out.sample must be numeric\n")
	if(as.numeric(out.sample) < 0) 
		stop("\ngogarchfit-->error: out.sample must be positive\n")
	n.start = round(out.sample, 0)
	n = dim(xdata$data)[1]
	if( (n-n.start) < 100) 
		stop("\ngogarchfit-->error: function requires at least 100 data\n points to run\n")
	data   = xdata$data
	index  = xdata$index
	period = xdata$period

	model$modeldata$data = data
	model$modeldata$index = index
	model$modeldata$period = period
	T = model$modeldata$T = n - n.start
	model$modeldata$n.start = n.start
	model$modeldata$asset.names = cnames
	fit = list()
	varcoef = arcoef = NULL
	#-----------------------------------------------------------------------------------
	# Conditional Mean Equation
	if( spec@model$modelinc[3]>0 ){
		tmp = mvmean.varfit(model, data, VAR.fit, T, out.sample, cluster = cluster, ggn = 2)
		model = tmp$model
		zdata = tmp$zdata
		mu = tmp$mu
		varcoef = tmp$varcoef
		p = tmp$p
		N = tmp$N		
	} else{
		tmp = mvmean.arfit(model, data, ARcoef, T, out.sample, solver, solver.control, fit.control, 
				cluster = cluster)
		model = tmp$model
		zdata = tmp$zdata
		mu = tmp$mu
		arcoef = tmp$arcoef
		p = tmp$p
		N = tmp$N
	}
	T = dim(zdata)[1] - out.sample
	
	#-----------------------------------------------------------------------------------
	# ICA estimation
	# 	iidret = rmgarch:::.makeiid(zdata[1:T, ], method = spec@model$ica, ica.fix = spec@model$ica.fix, gfun = "tanh", maxiter1 = 100000, e)
	iidret = .makeiid(zdata[1:T, ], method = spec@model$ica, ica.fix = spec@model$ica.fix, ...)
	A = iidret$A
	Y = iidret$Y
	W = iidret$W
	K = iidret$K
	Kinv = iidret$Kinv
	U = iidret$U
	
	#-----------------------------------------------------------------------------------
	# Factor Dynamics
	xspec = ugarchspec(mean.model = list(include.mean = FALSE, armaOrder=c(0,0)),
			variance.model = list(model = umodel$vmodel, garchOrder = umodel$garchOrder, 
					submodel = umodel$vsubmodel, variance.targeting = umodel$variance.targeting), 
			distribution.model = umodel$distribution)
	# n.comp defines speclist length
	# A (mixing matrix) = n.comp (factors) x n.assets
	speclist = multispec( replicate(n = NCOL(A), xspec) )
	# Expand the Signals if needed to the full length (to include out.sample) without
	# any lookeahead bias.
	Y = zdata %*% t(W)
	
	fitlist = multifit(multispec = speclist, data = Y, out.sample = out.sample, solver = solver, 
			solver.control = solver.control, fit.control = fit.control, cluster = cluster)
	
	garchcoef = coef(fitlist)
	converge = sapply(fitlist@fit, FUN = function(x) x@fit$convergence)
	if(any(converge == 1)){
		# it's useless returning any more info since we have an ICA factor matrix for which the
		# user can do very little other than probably run it again hoping that the new ICA matrix
		# will return a 'nicer' set (perhaps return the random number seed used to initialize ICA).
		stop("\ngogarchfit-->warning: convergence problem in univariate fit...", call. = FALSE)
	}
	#-----------------------------------------------------------------------------------
	# Multivariate Affine Distribution
	tmp = list()
	xtmp = .gogarchllh.independent(W, fitlist)
	tmp$factor.sigmas = sapply(fitlist@fit, FUN = function(x) x@fit$sigma, simplify = TRUE)
	#tmp$factor.resids = resids
	tmp$mu = mu
	tmp$llh = xtmp$llh
	tmp$factor.likelihoods = xtmp$likelihoods
	distr = model$modeldesc$distribution
	# Return List
	tmp$ufit = fitlist
	tmp$arcoef = arcoef
	tmp$varcoef = varcoef
	tmp$garchcoef = garchcoef
	tmp$Y = Y
	tmp$W = W
	tmp$A = A
	tmp$K = K
	tmp$Kinv = Kinv
	tmp$U = U
	tmp$E = iidret$E
	tmp$D = iidret$D
	tmp$residuals = zdata
	model$umodel = umodel
	tmp$timer = Sys.time() - tic
	
	ans = new("goGARCHfit",
			mfit = tmp,
			model = model)
	return(ans)
}

.gogarchfilter = function(fit, data, out.sample, n.old, cluster = NULL, ...)
{
	tic = Sys.time()
	#-----------------------------------------------------------------------------------
	# Data Extraction
	if( is.null( colnames(data) ) ) cnames = paste("Asset_", 1:m, sep = "") else cnames = colnames(data)	
	xdata = .extractmdata(data)
	if(!is.numeric(out.sample)) 
		stop("\ngogarchfit-->error: out.sample must be numeric\n")
	if(as.numeric(out.sample) < 0) 
		stop("\ngogarchfit-->error: out.sample must be positive\n")
	n.start = round(out.sample, 0)
	n = dim(xdata$data)[1]
	data   = xdata$data
	index  = xdata$index
	period = xdata$period
	model = fit@model
	model$modeldata$data = data
	model$modeldata$index = index
	model$modeldata$period = period
	T = model$modeldata$T = n - n.start
	model$modeldata$n.start = n.start
	model$modeldata$asset.names = cnames
	umodel = model$umodel
	varcoef = arcoef = NULL
	#-----------------------------------------------------------------------------------
	# Conditional Mean Equation
	if( model$modelinc[3]>0 ){
		tmp = mvmean.varfilter(model, data, varcoef = fit@mfit$varcoef, T, out.sample, ggn = 2)
		zdata = tmp$zdata
		mu = tmp$mu
		p = tmp$p
		N = tmp$N
	} else{
		tmp = mvmean.arfilter(model, data, arcoef = fit@mfit$arcoef, T, out.sample, 
				cluster = cluster)
		zdata = tmp$zdata
		mu = tmp$mu
		p = tmp$p
		N = tmp$N
	}
	#-----------------------------------------------------------------------------------
	# ICA
	A = fit@mfit$A
	W = fit@mfit$W
	K = fit@mfit$K
	Kinv = fit@mfit$Kinv
	U = fit@mfit$U
	# should we update the covariance matrix with new information?
	# --> This would change the whitening matrix and hence A, W, and the factors...
	Y = zdata %*% t(W)
	#-----------------------------------------------------------------------------------
	# Factor Dynamics
	m = NCOL(A)
	specx = vector(mode = "list", length = m)
	for(i in 1:m) specx[[i]] = ugarchspec(mean.model = list(include.mean = FALSE, armaOrder=c(0,0)),
				variance.model = list(model = umodel$vmodel, garchOrder = umodel$garchOrder, 
						submodel = umodel$vsubmodel, variance.targeting = FALSE), 
				distribution.model = umodel$distribution, fixed.pars = as.list(fit@mfit$garchcoef[,i]))
	mspec = multispec( specx )
	filterlist = multifilter(multifitORspec = mspec, data = Y, out.sample = out.sample, n.old = NULL, 
			cluster = cluster)
	#-----------------------------------------------------------------------------------
	# Multivariate Affine Distribution
	tmp = list()
	tmp$factor.sigmas = sapply(filterlist@filter, FUN = function(x) x@filter$sigma, simplify = TRUE)
	#tmp$factor.resids = resids
	tmp$mu = mu
	tmp$factor.likelihoods = likelihood(filterlist)
	distr = model$modeldesc$distribution
	
	xtmp = .gogarchllh.independent(W, filterlist)
	tmp$llh = xtmp$llh
	
	#-----------------------------------------------------------------------------------
	# Return List
	tmp$ufilter = filterlist
	tmp$arcoef = arcoef
	tmp$varcoef = varcoef
	tmp$garchcoef = fit@mfit$garchcoef
	tmp$Y = Y
	tmp$W = W
	tmp$A = A
	tmp$K = K
	tmp$Kinv = Kinv
	tmp$U = U
	tmp$E = fit@mfit$E
	tmp$D = fit@mfit$D
	tmp$residuals = zdata
	tmp$timer = Sys.time() - tic
	model$umodel = fit@model$umodel
	ans = new("goGARCHfilter",
			mfilter = tmp,
			model = model)
	return(ans)
}

.gogarchforecast = function(fit, n.ahead, n.roll, external.forecasts, cluster = NULL, ...)
{
	tic = Sys.time()
	ns = fit@model$modeldata$n.start
	if( n.roll > ns )
		stop("\ngogarchforecast-->error: n.roll must not be greater than out.sample!")
	if(n.roll>0 && n.ahead>1) 
		stop("\ngogarchforecast-->error: n.ahead must be equal to 1 when using n.roll\n")
	
	#-----------------------------------------------------------------------------------
	# Data Extraction
	model = fit@model
	# checks for external forecasts
	tf = n.ahead + n.roll
	if( !is.null( external.forecasts$mregfor ) ){
		mregfor = external.forecasts$mregfor
		if( !is.matrix(mregfor) ) stop("\nmregfor must be a matrix.")
		if( dim(mregfor)[1] < tf ) stop("\nmregfor must have at least n.ahead + n.roll observations to be used")
		mregfor = mregfor[1:tf, , drop = FALSE]
	} else{
		mregfor = NULL
	}
	#-----------------------------------------------------------------------------------
	# Conditional Mean Equation
	if( model$modelinc[3] > 0 ){
		Mu = mvmean.varforecast(model, mregfor, fit@mfit$varcoef, n.ahead, n.roll, ns, ggn = 2)
	} else{
		Mu = mvmean.arforecast(model, mregfor, fit@mfit$arcoef, n.ahead, n.roll, 
				ns, cluster = cluster)
	}
	#-----------------------------------------------------------------------------------
	# ICA & Factor Dynamics
	A = fit@mfit$A
	m = NCOL(A)
	
	specx = vector(mode = "list", length = m)
	for(i in 1:m) specx[[i]] = ugarchspec(mean.model = list(include.mean = FALSE, armaOrder=c(0,0)),
				variance.model = list(model = model$umodel$vmodel, garchOrder = model$umodel$garchOrder, 
						submodel = model$umodel$vsubmodel, variance.targeting = FALSE), 
				distribution.model = model$umodel$distribution, fixed.pars = as.list(fit@mfit$garchcoef[,i]))
	mspec = multispec( specx )
	forclist = multiforecast(multifitORspec = mspec, data = fit@mfit$Y, 
			n.ahead = n.ahead, out.sample = ns, n.roll = n.roll, cluster = cluster, ...)
	#-----------------------------------------------------------------------------------
	# Multivariate Affine Distribution
	tmp = list()
	xS = sigma(forclist)
	fS = array(NA, dim = c(n.ahead, m, n.roll+1))
	for(i in 1:(n.roll+1)){
		fS[,,i] = apply(xS, 3, function(x) x[,i])
	}
	tmp$factor.sigmas = fS
	tmp$mu = Mu
	#-----------------------------------------------------------------------------------
	# Return List
	tmp$A = A
	tmp$Y = fit@mfit$Y
	tmp$W = fit@mfit$W
	tmp$K = fit@mfit$K
	tmp$Kinv = fit@mfit$Kinv
	tmp$U = fit@mfit$U
	#tmp$E = fit@mfit$E
	#tmp$D = fit@mfit$D
	tmp$arcoef = fit@mfit$arcoef
	tmp$varcoef = fit@mfit$varcoef
	tmp$garchcoef = fit@mfit$garchcoef
	model$n.ahead = n.ahead
	model$n.roll = n.roll
	tmp$timer = Sys.time() - tic
	
	ans = new("goGARCHforecast",
			mforecast = tmp,
			model = model)
	return(ans)
}

.gogarchsim = function(fit, n.sim, n.start, m.sim, startMethod, prereturns, 
		preresiduals, presigma, mexsimdata, rseed, cluster = NULL, ...)
{
	tic = Sys.time()
	T = fit@model$modeldata$T
	Data = fit@model$modeldata$data[1:T, ]
	A = fit@mfit$A
	m = NCOL(A)
	simlist = vector(mode = "list", length = m)
	model = fit@model
	p = sum(model$modelinc[2:3])
	idx = 1:T
	#mo = fit@ufit@fit[[1]]@model$maxOrder
	if(is.null(startMethod[1])) startMethod = "unconditional" else startMethod = startMethod[1]
	if( !is.null(rseed) ){
		if( !is.matrix(rseed) ){
			xseed = matrix(NA, ncol = m, nrow = m.sim)
			if( length(rseed) != 1 ) 
				stop("\ngogarchsim-->error: if not matrix, rseed must be a single value!\n")
			xseed = matrix(rseed+seq(1,m.sim*m,by=1), ncol = m, nrow = m.sim, byrow=TRUE)
		} else{
			if(dim(rseed)[2] != m) 
				stop("\ngogarchsim-->error: rseed must have N (= n.assets) columns!\n")
			if(dim(rseed)[1] != m.sim) 
				stop("\ngogarchsim-->error: rseed must have m.sim rows!\n")
			xseed = rseed
		}
	} else{
		xseed = matrix( as.integer( runif( m*m.sim, 0, as.integer(Sys.time()) ) ), ncol = m, nrow = m.sim, byrow = TRUE)
	}
	mx = max(c(model$umodel$garchOrder), p)
	if(startMethod == "sample"){
		if(!is.na(presigma[1])){
			presigma = as.matrix(presigma)
			if(dim(presigma)[1]<mx) 
				stop(paste("\ngogarchsim-->error: presigma must have row dimension ", mx, sep = ""))
			presig = tail(presigma, max(model$umodel$garchOrder))
		} else{
			presig = tail(fit@mfit$factor.sigmas, max(model$umodel$garchOrder))
		}
		sm = TRUE
	} else{
		presig = NA
		sm = FALSE
	}
	if(startMethod == "sample"){
		if(!is.na(preresiduals[1])){
			preresiduals = as.matrix(preresiduals)
			if(dim(preresiduals)[1]<mx) 
				stop(paste("\ngogarchsim-->error: preresiduals must have row dimension ", mx, sep = ""))
			factor.preres = tail(preresiduals, max(model$umodel$garchOrder))
			preres  = tail(preresiduals, mx) %*% t(A)
		} else{
			factor.preres = tail(fit@mfit$Y[idx, ], max(model$umodel$garchOrder))
			preres = tail(fit@mfit$residuals[idx, ], mx)
		}
		sm = TRUE
	} else{
		factor.preres = NA
		sm = FALSE
		preres = tail(fit@mfit$residuals[idx, ], mx)
	}
	#-----------------------------------------------------------------------------------
	# Factor Dynamics
	if( !is.null(cluster) ){
		specx = vector(mode = "list", length = m)
		mo = max(model$umodel$garchOrder)
		for(i in 1:m) specx[[i]] = ugarchspec(mean.model = list(include.mean = FALSE, armaOrder=c(0,0)),
					variance.model = list(model = model$umodel$vmodel, garchOrder = model$umodel$garchOrder, 
							submodel = model$umodel$vsubmodel, variance.targeting = FALSE), 
					distribution.model = model$umodel$distribution, fixed.pars = as.list(fit@mfit$garchcoef[,i]))
		clusterEvalQ(cluster, loadNamespace("rugarch"))
		clusterExport(cluster, c("presig", "factor.preres", "sm", "specx", 
						"n.sim", "n.start", "m.sim", "startMethod", "xseed"), 
				envir = environment())
		simlist = parLapply(cluster, as.list(1:m), fun = function(i){
					rugarch::ugarchpath(specx[[i]], n.sim = n.sim + n.start, n.start = 0, 
							m.sim = m.sim, rseed = xseed[,i],
							presigma = if(sm) presig[,i] else NA, 
							preresiduals  = if(sm) factor.preres[,i] else NA)
				})
	} else{
		specx = vector(mode = "list", length = m)
		for(i in 1:m){
			specx = ugarchspec(
					mean.model = list(include.mean = FALSE, armaOrder=c(0,0)),
					variance.model = list(model = model$umodel$vmodel, garchOrder = model$umodel$garchOrder, 
							submodel = model$umodel$vsubmodel, variance.targeting = FALSE), 
					distribution.model = model$umodel$distribution, fixed.pars = as.list(fit@mfit$garchcoef[,i]))
			simlist[[i]] = ugarchpath(specx, n.sim = n.sim + n.start, 
								n.start = 0, m.sim = m.sim, rseed = xseed[,i],
								presigma = if(sm) presig[,i] else NA, 
								preresiduals  = if(sm) factor.preres[,i] else NA)
		}
	}
	
	#-----------------------------------------------------------------------------------
	# ICA & Multivariate Affine Distribution
	factor.sigmas = vector(mode = "list", length = m.sim)
	res = factor.res = vector(mode = "list", length = m.sim)
	for(i in 1:m.sim){
		factor.sigmas[[i]] = matrix(sapply(simlist, FUN = function(x) tail(x@path$sigmaSim[,i], n.sim + n.start)), ncol = m)
		factor.res[[i]] = matrix(sapply(simlist, FUN = function(x) tail(x@path$residSim[,i], n.sim + n.start)), ncol = m)
		res[[i]] = matrix(t(apply(factor.res[[i]], 1, FUN = function(x) x %*% t(A))), ncol = NROW(A))		
		factor.res[[i]] = tail(factor.res[[i]], n.sim)
		factor.sigmas[[i]] = tail(factor.sigmas[[i]], n.sim)
	}
	#-----------------------------------------------------------------------------------
	# Conditional Mean Equation
	if( model$modelinc[3] > 0 ){
		seriesSim = mvmean.varsim(model, Data, res, mexsimdata, prereturns, m.sim, n.sim, 
				n.start, startMethod = startMethod, cluster =  cluster, ggn = 2)
	} else{
		seriesSim = mvmean.arsim(model, Data, res, arcoef = fit@mfit$arcoef, 
				 mxn = fit@mfit$mxn, mexdata = fit@modeldata$mexdata, mexsimdata, 
				 prereturns, preres, m.sim, n.sim, n.start, startMethod = startMethod, 
				 cluster =  cluster)
	}
	#-----------------------------------------------------------------------------------
	# Return List
	tmp = list()
	tmp$distribution = fit@mfit$distr$model
	tmp$A = A
	tmp$Y = fit@mfit$Y
	tmp$W = fit@mfit$W
	tmp$K = fit@mfit$K
	tmp$Kinv = fit@mfit$Kinv
	tmp$U = fit@mfit$U
	tmp$arcoef = fit@mfit$arcoef
	tmp$varcoef = fit@mfit$varcoef
	tmp$garchcoef = fit@mfit$garchcoef
	tmp$factor.sigmaSim = factor.sigmas
	tmp$factor.residSim = factor.res
	tmp$residSim = res
	tmp$seriesSim = seriesSim
	tmp$rseed = xseed
	tmp$n.sim = n.sim
	tmp$m.sim = m.sim
	ans = new("goGARCHsim",
			msim = tmp,
			model = model)
	return(ans)
}


#-----------------------------------------------------------------------------------------------
# Rolling
.rollgogarch = function(spec, data, n.ahead = 1, forecast.length = 50, 
		n.start = NULL, refit.every = 25, refit.window = c("recursive", "moving"), 
		window.size = NULL, solver = "solnp", solver.control = list(), 
		fit.control = list(), rseed = NULL, cluster = NULL, save.fit = FALSE, 
		save.wdir = NULL, ...)
{
	model = spec@model
	if(is.null(solver.control$trace)) trace = 0 else trace = solver.control$trace
	if(is.null(fit.control$stationarity)) fit.control$stationarity = TRUE
	if(is.null(fit.control$fixed.se)) fit.control$fixed.se = FALSE
	if(is.null(fit.control$scale)) fit.control$scale = FALSE
	mm = match(names(fit.control), c("stationarity", "fixed.se", "scale"))
	if(any(is.na(mm))){
		idx = which(is.na(mm))
		enx = NULL
		for(i in 1:length(idx)) enx = c(enx, names(fit.control)[idx[i]])
		warning(paste(c("unidentified option(s) in fit.control:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
	}
	asset.names = colnames(data)
	xdata = .extractmdata(data)
	data   = xdata$data
	index  = xdata$index
	period = xdata$period
	if(is.null(fit.control$stationarity)) fit.control$stationarity = 1
	if(is.null(fit.control$fixed.se)) fit.control$fixed.se = 0
	T = NROW(data)
	if(n.ahead>1) 
		stop("\ngogarchroll:--> n.ahead>1 not supported...try again.")
	if(is.null(n.start)){
		if(is.null(forecast.length)) 
			stop("\ngogarchroll:--> forecast.length amd n.start are both NULL....try again.")
		n.start = T - forecast.length
	} else{
		forecast.length = T - n.start
	}
	if(T<=n.start) 
		stop("\ngogarchroll:--> start cannot be greater than length of data")
	# the ending points of the estimation window
	s = seq(n.start+refit.every, T, by = refit.every)
	m = length(s)
	# the rolling forecast length
	out.sample = rep(refit.every, m)
	# adjustment to include all the datapoints from the end
	if(s[m]<T){
		s = c(s,T)
		m = length(s)
		out.sample = c(out.sample, s[m]-s[m-1])
	}
	if(refit.window == "recursive"){
		rollind = lapply(1:m, FUN = function(i) 1:s[i])
	} else{
		if(!is.null(window.size)){
			if(window.size<100) stop("\ngogarchroll:--> window size must be greater than 100.")
			rollind = lapply(1:m, FUN = function(i) max(1, (s[i]-window.size-out.sample[i])):s[i])
		} else{
			rollind = lapply(1:m, FUN = function(i) (1+(i-1)*refit.every):s[i])
		}
	}
		
	if(is.null(rseed)) rseed = as.numeric(Sys.time()) else rseed = as.integer(rseed)
	WD <- getwd()
	if(is.null(save.wdir)){
		if (!is.null(WD)) setwd(WD)
	} else{
		ND = save.wdir
		if (!is.null(ND)) setwd(ND)
	}
	# steps:
	# 1. Fit using out.sample (gofit)
	# 2. Forecast with rolling (goforecast)
	# 3. Convolution (convolution)
	# 4. VaR (makeQ)
	# Setup Indexed Data for Fitting
	forc = icaA = vector(mode = "list", length = m)
	for(i in 1:m){
		if(i == 1){
			set.seed(rseed)
			gofit = gogarchfit(spec, data = data[rollind[[i]],], out.sample = out.sample[i], 
					solver = solver, fit.control = fit.control, 
					solver.control = solver.control, cluster = cluster, ...)
			A = gofit@mfit$A
			icaA[[i]] = A
			forc[[i]] = gogarchforecast(gofit, n.ahead = 1, n.roll = out.sample[i]-1, cluster = cluster)
			spec@model$ica.fix = list(A = NULL, K = NULL)
		} else{
			gofit = gogarchfit(spec, data = data[rollind[[i]],], out.sample = out.sample[i], 
					solver = solver, fit.control = fit.control, 
					solver.control = solver.control, cluster = cluster, 
					A.init = A, ...)
			A = gofit@mfit$A
			icaA[[i]] = A
			forc[[i]] = gogarchforecast(gofit, n.ahead = 1, n.roll = out.sample[i]-1, cluster = cluster)
		}
		if(save.fit){
			eval(parse(text = paste("ggroll_",i,"=gofit",sep = "")))
			eval(parse(text = paste("save(ggroll_",i,",file='ggroll_",i,".rda')",sep = "")))
		}
	}
	model$n.start = n.start
	model$n.refits = m
	model$refit.every = refit.every
	model$refit.window = refit.window
	model$window.size = window.size
	model$forecast.length = forecast.length
	model$n.start = n.start
	model$asset.names = asset.names
	model$umodel = spec@umodel
	model$rollind = rollind
	model$out.sample = out.sample
	model$A = icaA
	model$index = index
	ans = new("goGARCHroll",
			forecast = forc,
			model = model)
	# return to the current directory
	setwd(WD)
	return(ans)
}

# The likelihood of the affine linear model with independent margins is simply the
# sum of the univariate/independent log-likelihoods plus a term for the mixing matrix A
# It is not yet clear what to do when W is not square (PCA based dimensionality 
# reduced system).
.gogarchllh.independent = function(W, fit)
{
	# from the density transformation theorem (see for example Schmidt Technical Report P.14)
	if(is(fit, "uGARCHmultifit")){
		lik = sapply(fit@fit, FUN = function(x) -x@fit$log.likelihoods)
	} else{
		lik = sapply(fit@filter, FUN = function(x) -x@filter$log.likelihoods)
	}
	n = dim(lik)[1]
	if(NCOL(W) == NROW(W)){
		dt = log(abs(det(W)))
		likelihoods = apply(lik, 1, FUN = function(x) sum(x))
		llh = n*dt + sum(likelihoods)
	} else{
		likelihoods = apply(lik, 1, FUN = function(x) sum(x))
		llh = NA
	}
	return(list(llh = llh, likelihoods = likelihoods))
}



# make changes to subdivisions and rel.tol
#.safeintegrate = function(fun, lower, upper, ...)
#{
#	ans = try(integrate(fun, lower, upper, ...))
#	if(inherits(ans, "try-error")){
#		ans = integrate(fun, lower, upper, ...)
#	}
#	return(ans)
#}

# quadrature based numerical moments from the convoluted density
.dcintmoments = function(z, y, sm1 = NULL, fix.sm1 = TRUE, subdivisions = 200, 
		rel.tol = .Machine$double.eps^0.5, n.approx = 100, error.reltol = 0.05, ...)
{
	# convoluted density arising from numerical inversion of characteristic function of multivariate affine distribution
	dmad = .makeDNew(x = z, dx = y, h = diff(z[1:10])[1], standM = "sum")
	if(!is.null(sm1) && fix.sm1){
		m1 = sm1
	} else{
		fm1 = function(x) x * dmad(x)
			m1 = integrate(fm1, -Inf, Inf, subdivisions = 150, rel.tol=.Machine$double.eps^0.25)$value
			if(!is.null(sm1)){
				if( ((m1-sm1)/m1)>error.reltol || ((m1-sm1)/m1)<(-error.reltol) ) {
					m1 = try(integrate(fm1, -Inf, Inf, subdivisions = 200, rel.tol=.Machine$double.eps^0.5)$value, silent = TRUE)
					if(inherits(m1, "try-error")) {
						m1 = sm1
					} else{
						if( ((m1-sm1)/m1)>error.reltol || ((m1-sm1)/m1)<(-error.reltol) ) m1 = sm1
					}
				}
			}
	}
	fm2 = function(x) ((x-m1)^2) * dmad(x)
	m2 = sqrt(integrate(fm2, -Inf, Inf, subdivisions = subdivisions, rel.tol = rel.tol, stop.on.error = FALSE)$value)
	fm3 =  function(x) ((x-m1)^3) * dmad(x)
	m3 = integrate(fm3, -Inf, Inf, subdivisions = subdivisions, rel.tol = rel.tol, stop.on.error = FALSE)$value
	fm4 =  function(x) ((x-m1)^4) * dmad(x)
	m4 = integrate(fm4, -Inf, Inf, subdivisions = subdivisions, rel.tol = rel.tol, stop.on.error = FALSE)$value
	return(c(mean = m1, sd = m2, m3 = m3, m4 = m4))
}

.dfun = function(z, y, fn, lower = -Inf, upper = Inf, subdivisions = 200, rel.tol = .Machine$double.eps^0.5)
{
	# convoluted density arising from numerical inversion of characteristic function of multivariate affine distribution
	dmad = .makeDNew(x = z, dx = y, h = diff(z[1:10])[1], standM = "sum")
	ff = eval(parse(text = paste("function(x) ", fn, "*dmad(x)", sep = "")))
	ans = integrate(ff, lower, upper, subdivisions = subdivisions, rel.tol = rel.tol)$value
	return(ans)
}


.simulate.constant = function(resid, Mu, n, k)
{
	ans = matrix(NA, ncol = k, nrow = n)
	for(i in 1:n){
		ans[i,] = Mu + resid[i,]
	}
	return(ans)
}
