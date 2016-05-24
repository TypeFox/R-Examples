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

#----------------------------------------------------------------------------------
# spec method
#----------------------------------------------------------------------------------
dccspec = function(uspec, VAR = FALSE, robust = FALSE, lag = 1, lag.max = NULL, 
				lag.criterion = c("AIC", "HQ", "SC", "FPE"), external.regressors = NULL, 
				robust.control = list("gamma" = 0.25, "delta" = 0.01, "nc" = 10, "ns" = 500), 
		dccOrder = c(1,1), model = c("DCC", "aDCC", "FDCC"), groups = rep(1, length(uspec@spec)), 
		distribution = c("mvnorm", "mvt", "mvlaplace"), start.pars = list(), fixed.pars = list())
{
	UseMethod("dccspec")
}

.xdccspec = function(uspec, VAR = FALSE, robust = FALSE, lag = 1, lag.max = NULL, 
		lag.criterion = c("AIC", "HQ", "SC", "FPE"), external.regressors = NULL, 
		robust.control = list("gamma" = 0.25, "delta" = 0.01, "nc" = 10, "ns" = 500), 
		dccOrder = c(1,1), model = c("DCC", "aDCC", "FDCC"), groups = rep(1, length(uspec@spec)), 
		distribution = c("mvnorm", "mvt", "mvlaplace"), 
		start.pars = list(), fixed.pars = list())
{
	if(tolower(model[1]) == "fdcc"){
		spec = .fdccspec(uspec = uspec, VAR = VAR, robust = robust, lag = lag, 
				lag.max = lag.max, lag.criterion = lag.criterion[1], 
				external.regressors = external.regressors, 
				robust.control = robust.control, fdccOrder = dccOrder, 
				fdccindex = groups, distribution = distribution[1], 
				start.pars = start.pars, fixed.pars = fixed.pars)
	} else{
		if(tolower(model[1]) == "adcc") asymmetric = TRUE else asymmetric = FALSE
		spec = .dccspec(uspec = uspec, VAR = VAR, robust = robust, lag = lag, 
				lag.max = lag.max, lag.criterion = lag.criterion[1], 
				external.regressors = external.regressors, 
				robust.control = robust.control, dccOrder = dccOrder, 
				asymmetric = asymmetric, distribution = distribution[1], 
				start.pars = start.pars, fixed.pars = fixed.pars)
	}
	return(spec)
}
	
setMethod(f = "dccspec", signature(uspec = "uGARCHmultispec"), definition = .xdccspec)

.setfixeddcc = function(object, value){
	model = object@model
	umodel = object@umodel
	ipars = model$pars
	pars = unlist(value)
	names(pars) = parnames = tolower(names(pars))
	# included parameters in model
	modelnames = rownames(ipars[which(ipars[,4]==1 | ipars[,2]==1), ,drop=FALSE])
	inc = NULL
	for(i in seq_along(parnames)){
		if(is.na(match(parnames[i], modelnames))){
			warning( (paste("Unrecognized Parameter in Fixed Values: ", parnames[i], "...Ignored", sep = "")))
		} else{
			inc = c(inc, i)
		}
	}
	fixed.pars = pars[inc]
	names(fixed.pars) = tolower(names(pars[inc]))
	
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, umodel$vt)
	# set parameter values
	if(object@model$DCC == "FDCC"){
		tmp = dccspec(uspec = mspec, VAR = ifelse(model$modelinc[1]>0, TRUE, FALSE), 
				robust = ifelse(!is.null(model$varmodel$robust), model$varmodel$robust, FALSE), 
				lag = model$modelinc[1], lag.max = model$varmodel$lag.max, lag.criterion = model$varmodel$lag.criterion, 
				external.regressors = if(model$modelinc[2]>0) model$modeldata$mexdata else NULL, 
				robust.control = if(!is.null(model$varmodel$robust.control)) model$varmodel$robust.control else list(), 
				dccOrder = model$modelinc[4:5], model="FDCC",  groups = model$fdccindex,
				distribution = model$modeldesc$distribution, 
				start.pars = if(is.null(model$start.pars)) model$start.pars else NULL, fixed.pars = as.list(fixed.pars))
	} else{
		tmp = dccspec(uspec = mspec, VAR = ifelse(model$modelinc[1]>0, TRUE, FALSE), 
				robust = ifelse(!is.null(model$varmodel$robust), model$varmodel$robust, FALSE), 
				lag = model$modelinc[1], lag.max = model$varmodel$lag.max, lag.criterion = model$varmodel$lag.criterion, 
				external.regressors = if(model$modelinc[2]>0) model$modeldata$mexdata else NULL, 
				robust.control = if(!is.null(model$varmodel$robust.control)) model$varmodel$robust.control else list(), 
				dccOrder = model$modelinc[3:4], model = ifelse(model$modelinc[5]>0, "aDCC", "DCC"), 
				distribution = model$modeldesc$distribution, 
				start.pars = if(is.null(model$start.pars)) model$start.pars else NULL, fixed.pars = as.list(fixed.pars))
	}
	return(tmp)
}

setReplaceMethod(f="setfixed", signature= c(object = "DCCspec", value = "vector"), definition = .setfixeddcc)


.setstartdcc = function(object, value){
	model = object@model
	umodel = object@umodel
	ipars = model$pars
	pars = unlist(value)
	names(pars) = parnames = tolower(names(pars))
	# included parameters in model
	modelnames = rownames(ipars[which(ipars[,4]==1 | ipars[,2]==1), ,drop=FALSE])
	inc = NULL
	for(i in seq_along(parnames)){
		if(is.na(match(parnames[i], modelnames))){
			warning( (paste("Unrecognized Parameter in Start Values: ", parnames[i], "...Ignored", sep = "")))
		} else{
			inc = c(inc, i)
		}
	}
	start.pars = pars[inc]
	names(start.pars) = tolower(names(pars[inc]))
	
	uspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, umodel$vt)
	
	# set parameter values
	if(object@model$DCC == "FDCC"){
		tmp =  dccspec(uspec, VAR = ifelse(model$modelinc[1]>0, TRUE, FALSE), 
				robust = ifelse(!is.null(model$varmodel$robust), model$varmodel$robust, FALSE), 
				lag = model$modelinc[1], lag.max = model$varmodel$lag.max, lag.criterion = model$varmodel$lag.criterion, 
				external.regressors = if(model$modelinc[2]>0) model$modeldata$mexdata else NULL, 
				robust.control = if(!is.null(model$varmodel$robust.control)) model$varmodel$robust.control else list(), 
				dccOrder = model$modelinc[4:5], model="FDCC",  groups = model$fdccindex,
				distribution = model$modeldesc$distribution, 
				start.pars = as.list(start.pars),  fixed.pars =  model$fixed.pars)
	} else{
		tmp =  dccspec(uspec, VAR = ifelse(model$modelinc[1]>0, TRUE, FALSE), 
				robust = ifelse(!is.null(model$varmodel$robust), model$varmodel$robust, FALSE), 
				lag = model$modelinc[1], lag.max = model$varmodel$lag.max, lag.criterion = model$varmodel$lag.criterion, 
				external.regressors = if(model$modelinc[2]>0) model$modeldata$mexdata else NULL, 
				robust.control = if(!is.null(model$varmodel$robust.control)) model$varmodel$robust.control else list(), 
				dccOrder = model$modelinc[3:4], model = ifelse(model$modelinc[5]>0, "aDCC", "DCC"), 
				distribution = model$modeldesc$distribution, 
				start.pars = as.list(start.pars),  fixed.pars =  model$fixed.pars)
	}
	return(tmp)
}

setReplaceMethod(f="setstart", signature= c(object = "DCCspec", value = "vector"), definition = .setstartdcc)
#----------------------------------------------------------------------------------
# fit method
#----------------------------------------------------------------------------------

dccfit = function(spec, data, out.sample = 0, solver = "solnp", solver.control = list(), 
		fit.control = list(eval.se = TRUE, stationarity = TRUE, scale = FALSE), 
		cluster = NULL, fit = NULL, VAR.fit = NULL, realizedVol = NULL, ...)
{
	UseMethod("dccfit")
}

.xdccfit = function(spec, data, out.sample = 0, solver = "solnp", solver.control = list(), 
		fit.control = list(eval.se = TRUE, stationarity = TRUE, scale = FALSE), 
		cluster = NULL, fit = NULL, VAR.fit = NULL, realizedVol = NULL, ...)
{
	if(spec@model$DCC == "FDCC"){
		ans = .fdccfit(spec = spec, data = data, out.sample = out.sample, 
				solver = solver, solver.control = solver.control, 
				fit.control = fit.control, cluster = cluster, fit = fit, 
				VAR.fit = VAR.fit, realizedVol = realizedVol, ...)
	} else{
		ans = .dccfit(spec = spec, data = data, out.sample = out.sample, 
				solver = solver, solver.control = solver.control, 
				fit.control = fit.control, cluster = cluster, fit = fit, 
				VAR.fit = VAR.fit, realizedVol = realizedVol, ...)
	}
	return(ans)
}
setMethod("dccfit", signature(spec = "DCCspec"), .xdccfit)

#----------------------------------------------------------------------------------
# filter method
#----------------------------------------------------------------------------------
dccfilter = function(spec, data, out.sample = 0, filter.control = list(n.old = NULL), 
		cluster = NULL, varcoef = NULL, realizedVol = NULL, ...)
{
	UseMethod("dccfilter")
}

.xdccfilter = function(spec, data, out.sample = 0, filter.control = list(n.old = NULL), 
		cluster = NULL, varcoef = NULL, realizedVol = NULL, ...)
{
	if(spec@model$DCC == "FDCC"){
		ans = .fdccfilter(spec = spec, data = data, out.sample = out.sample, 
				filter.control = filter.control, cluster = cluster, varcoef = varcoef, 
				realizedVol = realizedVol,...)
	} else{
		ans = .dccfilter(spec = spec, data = data, out.sample = out.sample, 
				filter.control = filter.control, cluster = cluster, varcoef = varcoef, 
				realizedVol = realizedVol,...)
	}
	return(ans)
}

setMethod("dccfilter", signature(spec = "DCCspec"), .xdccfilter)

#----------------------------------------------------------------------------------
# forecast method
#----------------------------------------------------------------------------------
dccforecast = function(fit, n.ahead = 1, n.roll = 0, 
		external.forecasts = list(mregfor = NULL, vregfor = NULL), cluster = NULL, ...)
{
	UseMethod("dccforecast")
}

.xdccforecast = function(fit, n.ahead = 1, n.roll = 0, 
		external.forecasts = list(mregfor = NULL, vregfor = NULL), cluster = NULL, ...)
{
	if(fit@model$DCC == "FDCC"){
		ans = .fdccforecast(fit, n.ahead = n.ahead, n.roll = n.roll, 
				external.forecasts = external.forecasts, cluster = cluster, ...)
	} else{
		ans = .dccforecast(fit, n.ahead = n.ahead, n.roll = n.roll, 
				external.forecasts = external.forecasts, cluster = cluster, ...)
	}
	return(ans)
}
setMethod("dccforecast", signature(fit = "DCCfit"), .xdccforecast)

#----------------------------------------------------------------------------------
# roll method
#----------------------------------------------------------------------------------

dccroll = function(spec, data, n.ahead = 1, forecast.length = 50, refit.every = 25, 
		n.start = NULL, refit.window = c("recursive", "moving"), window.size = NULL, 
		solver = "solnp", solver.control = list(), 
		fit.control = list(eval.se = TRUE, stationarity = TRUE, scale = FALSE), 
		cluster = NULL, save.fit = FALSE, save.wdir = NULL, 
		realizedVol = NULL, ...)
{
	UseMethod("dccroll")
}


setMethod("dccroll", signature(spec = "DCCspec"), .rolldcc)


#----------------------------------------------------------------------------------
# simulation method
#----------------------------------------------------------------------------------

dccsim = function(fitORspec, n.sim = 1000, n.start = 0, m.sim = 1, 
		startMethod = c("unconditional", "sample"), presigma = NULL, 
		preresiduals = NULL, prereturns = NULL, preQ = NULL, preZ = NULL, 
		Qbar = NULL, Nbar = NULL, rseed = NULL, mexsimdata = NULL, 
		vexsimdata = NULL, cluster = NULL, VAR.fit = NULL, 
		prerealized = NULL, ...)
{
	UseMethod("dccsim")
}

.xdccsim.fit = function(fitORspec, n.sim = 1000, n.start = 0, m.sim = 1, 
		startMethod = c("unconditional", "sample"), presigma = NULL, 
		preresiduals = NULL, prereturns = NULL, preQ = NULL, preZ = NULL, 
		Qbar = NULL, Nbar = NULL, rseed = NULL, mexsimdata = NULL, 
		vexsimdata = NULL, cluster = NULL, VAR.fit = NULL, prerealized = NULL, ...)
{
	if(fitORspec@model$DCC == "FDCC"){
		ans = .fdccsim.fit(fitORspec = fitORspec, n.sim = n.sim, n.start = n.start, 
				m.sim = m.sim, startMethod = startMethod, presigma = presigma, 
				preresiduals = preresiduals, prereturns = prereturns, 
				preQ = preQ, preZ = preZ, Qbar = Qbar, rseed = rseed, 
				mexsimdata = mexsimdata, vexsimdata = vexsimdata, 
				cluster = cluster, prerealized = prerealized, ...)
	} else{
		ans = .dccsim.fit(fitORspec = fitORspec, n.sim = n.sim, n.start = n.start, 
				m.sim = m.sim, startMethod = startMethod, presigma = presigma, 
				preresiduals = preresiduals, prereturns = prereturns, 
				preQ = preQ, preZ = preZ, Qbar = Qbar, Nbar = Nbar, 
				rseed = rseed, mexsimdata = mexsimdata, vexsimdata = vexsimdata, 
				cluster = cluster, prerealized = prerealized, ...)
	}
	return(ans)
}

.xdccsim.spec = function(fitORspec, n.sim = 1000, n.start = 0, m.sim = 1, 
		startMethod = c("unconditional", "sample"), presigma = NULL, 
		preresiduals = NULL, prereturns = NULL, preQ = NULL, preZ = NULL, 
		Qbar = NULL, Nbar = NULL, rseed = NULL, mexsimdata = NULL, 
		vexsimdata = NULL, cluster = NULL, VAR.fit = NULL, 
		prerealized = NULL, ...)
{
	if(fitORspec@model$DCC == "FDCC"){
		ans = .fdccsim.spec(fitORspec = fitORspec, n.sim = n.sim, n.start = n.start, 
				m.sim = m.sim, startMethod = startMethod, presigma = presigma, 
				preresiduals = preresiduals, prereturns = prereturns, 
				preQ = preQ, preZ = preZ, Qbar = Qbar, rseed = rseed, 
				mexsimdata = mexsimdata, vexsimdata = vexsimdata, 
				cluster = cluster, VAR.fit = VAR.fit, prerealized = prerealized, ...)
	} else{
		ans = .dccsim.spec(fitORspec = fitORspec, n.sim = n.sim, n.start = n.start, 
				m.sim = m.sim, startMethod = startMethod, presigma = presigma, 
				preresiduals = preresiduals, prereturns = prereturns, 
				preQ = preQ, preZ = preZ, Qbar = Qbar, Nbar = Nbar, 
				rseed = rseed, mexsimdata = mexsimdata, vexsimdata = vexsimdata, 
				cluster = cluster, VAR.fit = VAR.fit, prerealized = prerealized, ...)
	}
	return(ans)
}

setMethod("dccsim", signature(fitORspec = "DCCfit"), .xdccsim.fit)
setMethod("dccsim", signature(fitORspec = "DCCspec"), .xdccsim.spec)
#----------------------------------------------------------------------------------
# fitted
#----------------------------------------------------------------------------------
.fitted.dccfit = function(object)
{
	T = object@model$modeldata$T
	D = object@model$modeldata$index
	ans = xts(object@model$mu[1:T,], D[1:T])
	colnames(ans) = object@model$modeldata$asset.names
	return( ans )
}

setMethod("fitted", signature(object = "DCCfit"), .fitted.dccfit)
setMethod("fitted", signature(object = "DCCfilter"), .fitted.dccfit)

.fitted.dccforecast = function(object)
{
	T = object@model$modeldata$T
	m = NCOL(object@model$modeldata$data)
	n.roll = object@model$n.roll
	n.ahead = object@model$n.ahead
	T0 = object@model$modeldata$index[T:(T+n.roll)]
	ans = array(object@mforecast$mu, dim = c(n.ahead, m, n.roll+1), 
			dimnames = list(paste("T+", 1:n.ahead,sep=""), object@model$modeldata$asset.names, as.character(T0)))
	return( ans )
}

setMethod("fitted", signature(object = "DCCforecast"), .fitted.dccforecast)

.fitted.dccsim = function(object, sim = 1)
{
	n = length(object@msim$simR)
	m.sim = as.integer(sim)
	if( m.sim > n | m.sim < 1 ) stop("\rmgarch-->error: fitted (simulation) sim index out of bounds!")
	ans = object@msim$simX[[m.sim]]
	rownames(sim) = NULL
	colnames(ans) = object@model$modeldata$asset.names
	return( ans )
}

setMethod("fitted", signature(object = "DCCsim"), .fitted.dccsim)

.fitted.dccroll = function(object)
{
	n = length(object@mforecast)
	index = object@model$index
	m = NCOL(object@model$data)
	D = index[(object@model$n.start+1):(object@model$n.start+object@model$forecast.length)]
	M = matrix(unlist(sapply(object@mforecast, function(x) fitted(x))), ncol = m, byrow = TRUE)
	M = xts(M, D)
	colnames(M) = object@model$modeldata$asset.names
	return( M )
}

setMethod("fitted", signature(object = "DCCroll"), .fitted.dccroll)

#----------------------------------------------------------------------------------
# residuals
#----------------------------------------------------------------------------------
.residuals.dccfit = function(object)
{
	T = object@model$modeldata$T
	D = object@model$modeldata$index
	ans = xts(object@model$residuals[1:T,], D[1:T])
	colnames(ans) = object@model$modeldata$asset.names
	return( ans )
}

setMethod("residuals", signature(object = "DCCfit"), .residuals.dccfit)
setMethod("residuals", signature(object = "DCCfilter"), .residuals.dccfit)

#----------------------------------------------------------------------------------
# rcor
#----------------------------------------------------------------------------------
.rcor.dccfit = function(object, type = "R")
{
	T = object@model$modeldata$T
	D = object@model$modeldata$index
	nam = object@model$modeldata$asset.names
	if( type == "R"){
		if(class(object)[1]=="DCCfit") tmp = object@mfit$R else tmp = object@mfilter$R
		m = dim(tmp[[1]])[2]
		R = array(NA, dim = c(m, m, T))
		R[1:m,1:m, ] = sapply(tmp[1:T], FUN = function(x) x)
		dimnames(R)<-list(nam, nam, as.character(D[1:T]))
		return( R )
	} else{
		if(class(object)[1]=="DCCfit") tmp = object@mfit$Q else tmp = object@mfilter$Q
		m = dim(tmp[[1]])[2]
		Q = array(NA, dim = c(m, m, T))
		Q[1:m,1:m, ] = sapply(tmp[1:T], FUN = function(x) x)
		dimnames(Q)<-list(nam, nam, as.character(D[1:T]))
		return( Q )
	}
}

setMethod("rcor", signature(object = "DCCfit"), .rcor.dccfit)
setMethod("rcor", signature(object = "DCCfilter"), .rcor.dccfit)

.rcor.dccforecast = function(object, type = "R")
{
	n.roll = object@model$n.roll
	n.ahead = object@model$n.ahead
	T = object@model$modeldata$T
	D = as.character(object@model$modeldata$index[T:(T+n.roll)])
	nam = object@model$modeldata$asset.names
	if( type == "R"){
		tmp = object@mforecast$R
	} else{
		tmp = object@mforecast$Q
	}
	for(i in 1:(n.roll+1)){ dimnames(tmp[[i]]) = list(nam, nam, paste("T+",1:n.ahead,sep="")) }
	names(tmp) = D
	# attr(tmp, "T0") = D
	return( tmp )
}

setMethod("rcor", signature(object = "DCCforecast"), .rcor.dccforecast)

.rcor.dccsim = function(object, type = "R", sim = 1)
{
	n = object@model$m.sim
	m.sim = as.integer(sim)
	nam = object@model$modeldata$asset.names
	n.sim = object@model$n.sim
	if( m.sim > n | m.sim < 1 ) stop("\rmgarch-->error: rcor sim index out of bounds!")
	if( type == "R"){
		tmp = object@msim$simR[[m.sim]]
	} else{
		tmp = object@msim$simQ[[m.sim]]
	}
	dimnames(tmp) = list(nam, nam, 1:n.sim)
	return( tmp )
	
}

setMethod("rcor", signature(object = "DCCsim"), .rcor.dccsim)


.rcor.dccroll = function(object, type = "R")
{
	n = length(object@mforecast)
	index = object@model$index
	m = NCOL(object@model$data)
	fl = object@model$forecast.length
	nam = object@model$modeldata$asset.names
	D = index[(object@model$n.start+1):(object@model$n.start+fl)]
	Cf = array(unlist(sapply(object@mforecast, function(x) unlist(rcor(x, type=type)))), dim=c(m,m,fl),
			dimnames = list(nam, nam, as.character(D)))
	return( Cf )
}

setMethod("rcor", signature(object = "DCCroll"), .rcor.dccroll)
#----------------------------------------------------------------------------------
# rcov
#----------------------------------------------------------------------------------
.rcov.dccfit = function(object)
{
	T = object@model$modeldata$T
	D = object@model$modeldata$index
	nam = object@model$modeldata$asset.names
	if(class(object)[1]=="DCCfit") H = object@mfit$H[,,1:T] else H = object@mfilter$H[,,1:T]
	dimnames(H)<-list(nam, nam, as.character(D[1:T]))	
	return( H )
}

setMethod("rcov", signature(object = "DCCfit"), .rcov.dccfit)
setMethod("rcov", signature(object = "DCCfilter"), .rcov.dccfit)


.rcov.dccforecast = function(object)
{
	n.roll = object@model$n.roll
	n.ahead = object@model$n.ahead
	T = object@model$modeldata$T
	D = as.character(object@model$modeldata$index[T:(T+n.roll)])
	nam = object@model$modeldata$asset.names
	H = object@mforecast$H
	for(i in 1:(n.roll+1)){ dimnames(H[[i]]) = list(nam, nam, paste("T+",1:n.ahead,sep="")) }
	names(H) = D
	#attr(H, "T0") = D	
	return( H )
}

setMethod("rcov", signature(object = "DCCforecast"), .rcov.dccforecast)


.rcov.dccsim = function(object, sim = 1)
{
	n = object@model$m.sim
	nam = object@model$modeldata$asset.names
	n.sim = object@model$n.sim
	m.sim = as.integer(sim)
	if( m.sim > n | m.sim < 1 ) stop("\rmgarch-->error: rcov sim index out of bounds!")
	tmp = object@msim$simH[[m.sim]]
	dimnames(tmp) = list(nam, nam, 1:n.sim)
	return( tmp )
}

setMethod("rcov", signature(object = "DCCsim"), .rcov.dccsim)

.rcov.dccroll = function(object)
{
	n = length(object@mforecast)
	index = object@model$index
	m = NCOL(object@model$data)
	fl = object@model$forecast.length
	nam = object@model$modeldata$asset.names
	D = index[(object@model$n.start+1):(object@model$n.start+fl)]
	Cf = array(unlist(sapply(object@mforecast, function(x) unlist(rcov(x)))), dim=c(m,m,fl),
			dimnames = list(nam, nam, as.character(D)))
	return( Cf )
}

setMethod("rcov", signature(object = "DCCroll"), .rcov.dccroll)
#----------------------------------------------------------------------------------
# sigma
#----------------------------------------------------------------------------------
.sigma.dccfit = function(object)
{
	T = object@model$modeldata$T
	D = object@model$modeldata$index
	nam = object@model$modeldata$asset.names
	if(class(object)[1] == "DCCfit"){
		H = object@mfit$H[,,1:T,drop=FALSE]
		m = dim(H)[2]
		sig = sqrt(.Call("ArrayDiag", H, c(m,m,T), PACKAGE="rmgarch"))
	} else{
		H = object@mfilter$H[,,1:T,drop=FALSE]
		m = dim(H)[2]
		sig = sqrt(.Call("ArrayDiag", H, c(m,m,T), PACKAGE="rmgarch"))		
	}
	colnames(sig) = nam
	sig = xts(sig, D[1:T])
	return( sig )
}

setMethod("sigma", signature(object = "DCCfit"), .sigma.dccfit)
setMethod("sigma", signature(object = "DCCfilter"), .sigma.dccfit)

.sigma.dccforecast = function(object)
{
	T = object@model$modeldata$T
	D = object@model$modeldata$index
	n.ahead = object@model$n.ahead
	nam = object@model$modeldata$asset.names
	n.roll = object@model$n.roll
	m = dim(object@model$mpars)[2]-1
	sig = array(NA, dim = c(n.ahead, m, n.roll+1))
	for(i in 1:(n.roll+1)){
		sig[,,i] = sqrt(.Call("ArrayDiag", object@mforecast$H[[i]], c(m,m,n.ahead), PACKAGE="rmgarch"))
	}
	dimnames(sig) = list(paste("T+",1:n.ahead,sep=""), nam, as.character(D[T:(T+n.roll)]))
	return( sig )
}

setMethod("sigma", signature(object = "DCCforecast"), .sigma.dccforecast)

.sigma.dccsim = function(object, sim = 1)
{
	H = rcov(object, sim = sim)
	m = dim(H)[1]
	n = dim(H)[3]
	nam = object@model$modeldata$asset.names
	sig = sapply(1:m, FUN = function(i) sqrt(H[i,i, ]))
	colnames(sig) = nam
	rownames(sig) = NULL
	return( sig )
}

setMethod("sigma", signature(object = "DCCsim"), .sigma.dccsim)

.sigma.dccroll = function(object)
{
	n = length(object@mforecast)
	index = object@model$index
	m = NCOL(object@model$data)
	D = index[(object@model$n.start+1):(object@model$n.start+object@model$forecast.length)]
	S = matrix(unlist(sapply(object@mforecast, function(x) sigma(x))), ncol = m, byrow = TRUE)
	S = xts(S, D)
	colnames(S) = object@model$modeldata$asset.names
	return( S )
}

setMethod("sigma", signature(object = "DCCroll"), .sigma.dccroll)

#----------------------------------------------------------------------------------
# rskew
#----------------------------------------------------------------------------------
rskew = function(object, ...)
{
	UseMethod("rskew")
}

.skew.dccfit = function(object)
{
	sk = NA
	if( object@model$modelinc[7]>0 ){
		cf = object@model$mpars[,dim(object@model$mpars)[2]]
		nx = which( substr(names(cf), 1, 5) == "mskew" )
		sk = cf[nx]
	}
	return( sk )
}

setMethod("rskew", signature(object = "DCCfit"), .skew.dccfit)
setMethod("rskew", signature(object = "DCCfilter"), .skew.dccfit)
setMethod("rskew", signature(object = "DCCforecast"), .skew.dccfit)

.skew.dccroll = function(object)
{
	m = dim(object@model$modeldata$data)[1]
	n = length(object@mforecast)
	sk = rep(NA, n)
	if( object@model$modelinc[7]>0 ){
		sk = sapply(object@model$rollcoef, FUN = function(x) x["mskew",m+1])
		colnames(sk) = paste("roll-", 1:n, sep = "")
	}
	return( sk )
}

setMethod("rskew", signature(object = "DCCroll"), .skew.dccroll)
#----------------------------------------------------------------------------------
# rshape
#----------------------------------------------------------------------------------
rshape = function(object, ...)
{
	UseMethod("rshape")
}

.shape.dccfit = function(object)
{
	sh = NA
	if( object@model$modelinc[6]>0 ){
		cf = object@model$mpars[,dim(object@model$mpars)[2]]
		nx = which( substr(names(cf), 1, 6) == "mshape" )
		sh = cf[nx]
	}
	return( sh )
}

setMethod("rshape", signature(object = "DCCfit"), .shape.dccfit)
setMethod("rshape", signature(object = "DCCfilter"), .shape.dccfit)
setMethod("rshape", signature(object = "DCCforecast"), .shape.dccfit)

.shape.dccroll = function(object)
{
	m = dim(object@model$modeldata$data)[1]
	n = length(object@mforecast)
	sh = rep(NA, n)
	if( object@model$modelinc[6]>0 ){
		sh = sapply(object@model$rollcoef, FUN = function(x) x["mshape",m+1])
		colnames(sh) = paste("roll-", 1:n, sep = "")
	}
	return( sh )
}

setMethod("rshape", signature(object = "DCCroll"), .shape.dccroll)
#----------------------------------------------------------------------------------
# coef
#----------------------------------------------------------------------------------
.coef.dccfit = function(object, type = "all")
{
	mpars = object@model$mpars
	m = dim(mpars)[2]-1
	if( type == "all" ){
		if(class(object)[1]=="DCCfit") cf = object@mfit$coef else cf = object@mfilter$coef
	} else if( type == "dcc"){
		cf = object@model$mpars[which(object@model$eidx[,m+1]==1), m+1]
		if(class(object)[1]=="DCCfit"){
			names(cf) = object@mfit$dccnames
		} else{
			names(cf) = object@mfilter$dccnames
		}
	} else{
		cf = object@model$mpars[which(object@model$eidx[,1:m]==1)]
		if(class(object)[1]=="DCCfit"){
			names(cf) = object@mfit$garchnames
		} else{
			names(cf) = object@mfilter$garchnames
		}
	}
	return( cf )
}

setMethod("coef", signature(object = "DCCfit"), .coef.dccfit)
setMethod("coef", signature(object = "DCCfilter"), .coef.dccfit)

.coef.dccroll = function(object)
{
	model = object@model
	cnames = object@model$modeldata$asset.names
	m = dim(model$umodel$modelinc)[2]
	allnames = NULL
	midx = .fullinc(model$modelinc, model$umodel)
	
	for(i in 1:m){
		allnames = c(allnames, paste("[",cnames[i],"].", rownames(midx[midx[,i]==1,i, drop = FALSE]), sep = ""))
	}
	garchnames = allnames
	dccnames = rownames(midx[midx[,m+1]==1,m+1, drop = FALSE])
	allnames = c(garchnames, paste("[Joint]", rownames(midx[midx[,m+1]==1,m+1, drop = FALSE]), sep = ""))
	n = length(object@mforecast)
	rollmat = matrix(NA, ncol = n, nrow = length(allnames))
	for(i in 1:n){
		rollmat[,i] = model$rollcoef[[i]][midx==1]
	}
	colnames(rollmat) = paste("roll-", 1:n, sep = "")
	rownames(rollmat) = allnames
	return( rollmat )
}

setMethod("coef", signature(object = "DCCroll"), .coef.dccroll)
#----------------------------------------------------------------------------------
# likelihood
#----------------------------------------------------------------------------------
.likelihood.dccfit = function(object)
{
	switch(class(object)[1],
			DCCfit = object@mfit$llh,
			DCCfilter = object@mfilter$llh,
			DCCroll = {
				n = length(object@mforecast)
				cf = sapply( object@model$rolllik, FUN = function(x) x )
				names(cf) = paste("roll-", 1:n, sep = "")
				return(cf)
			})
}

setMethod("likelihood", signature(object = "DCCfit"), .likelihood.dccfit)
setMethod("likelihood", signature(object = "DCCfilter"), .likelihood.dccfit)
setMethod("likelihood", signature(object = "DCCroll"), .likelihood.dccfit)
#----------------------------------------------------------------------------------
# infocriteria
#----------------------------------------------------------------------------------
.dccinfocriteria.fit = function(object)
{
	n = object@model$modeldata$T
	m  = NCOL(object@model$modeldata$data)
	if(object@model$modelinc[1]>0){
		npvar = dim(object@model$varcoef)[1] * dim(object@model$varcoef)[2]
	} else{
		npvar = 0
	}
	# VAR + GARCH + DCC + Qbar
	estpars = NROW(object@mfit$matcoef)
	np = npvar + estpars + ( (m^2 - m)/2 )
	itest = rugarch:::.information.test(likelihood(object), nObs = n, nPars = np)
	itestm = matrix(0, ncol = 1, nrow = 4)
	itestm[1,1] = itest$AIC
	itestm[2,1] = itest$BIC
	itestm[3,1] = itest$SIC
	itestm[4,1] = itest$HQIC
	colnames(itestm) = ""
	rownames(itestm) = c("Akaike", "Bayes", "Shibata", "Hannan-Quinn")
	return(itestm)
}

setMethod("infocriteria", signature(object = "DCCfit"), .dccinfocriteria.fit)
#----------------------------------------------------------------------------------

dcc.symcheck = function(x, m, d = rep(1, m), tol = 1e-12)
{
	n1 = dim(x)[1]
	n2 = dim(x)[2]
	if( n1 != n2 ) stop("\nmatrix not square!")
	if( n1 != m )  stop("\nmatrix not of expected dimension!")
	if( max(abs(x - t(x))) > tol ) stop("\nmatrix is not symmetric!")
	if( !is.null( d ) ){
		if( any( diag(x) != d ) ) stop("\nwrong values on diagonal of matrix!")
	}
	return( 0 )
}

#----------------------------------------------------------------------------------
# show methods
#----------------------------------------------------------------------------------
setMethod("show",
		signature(object = "DCCspec"),
		function(object){
			m = dim(object@umodel$modelinc)[2]
			dccpars = sum(object@model$modelinc[3:7])
			garchpars = sum( object@umodel$modelinc[1:18,] )
			# VAR = mxm x lags + lags*constant + lags*mxreg
			varpars = object@model$modelinc[1]*(m*m + m + object@model$modelinc[2])
			cat(paste("\n*------------------------------*", sep = ""))
			cat(paste("\n*       DCC GARCH Spec         *", sep = ""))
			cat(paste("\n*------------------------------*", sep = ""))
			if(object@model$DCC=="FDCC"){
				cat("\nModel          : ", paste(object@model$DCC, "(", object@model$modelinc[4], ",", object@model$modelinc[5],")", sep=""))
				cat("\nNo.Groups      : ", object@model$modelinc[3])
			} else{
				cat("\nModel          : ", paste(object@model$DCC, "(", object@model$modelinc[3], ",", object@model$modelinc[4],")", sep=""))
			}
			cat("\nEstimation     : ", object@model$modeldesc$type)
			cat("\nDistribution   : ", object@model$modeldesc$distribution)
			cat("\nNo. Parameters : ", dccpars + garchpars + varpars + ( (m^2 - m)/2 ))
			cat("\nNo. Series     : ", m)
			cat("\n\n")
			invisible(object)
		})
		
# fit show
setMethod("show",
		signature(object = "DCCfit"),
		function(object){
			m = dim(object@model$modeldata$data)[2]
			T = object@model$modeldata$T
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\n*          DCC GARCH Fit          *", sep = ""))
			cat(paste("\n*---------------------------------*", sep = ""))	
			cat("\n\nDistribution         : ", object@model$modeldesc$distribution)
			if(object@model$DCC=="FDCC"){
				cat("\nModel                : ", paste(object@model$DCC, "(", object@model$modelinc[4], ",", object@model$modelinc[5],")", sep=""))
				cat("\nNo.Groups            : ", object@model$modelinc[3])
			} else{
				cat("\nModel                : ", paste(object@model$DCC, "(", object@model$modelinc[3], ",", object@model$modelinc[4],")", sep=""))
			}
			# Need to show parameters as Q + DCC + GARCH
			if(object@model$modelinc[1]>0){
				npvar = dim(object@model$varcoef)[1] * dim(object@model$varcoef)[2]
			} else{
				npvar = 0
			}
			NP = paste("[",npvar, "+", length(object@mfit$garchnames),"+",length(object@mfit$dccnames), "+",(m^2 - m)/2,"]", sep="") 
			cat("\nNo. Parameters       : ", npvar+NROW(object@mfit$matcoef) + ( (m^2 - m)/2 ))
			cat("\n[VAR GARCH DCC UncQ] :", NP)
			cat("\nNo. Series           : ", m)
			cat("\nNo. Obs.             : ", T)
			cat("\nLog-Likelihood       : ", object@mfit$llh)
			cat("\nAv.Log-Likelihood    : ", round(object@mfit$llh/T, 2), "\n")
			cat("\nOptimal Parameters")
			cat(paste("\n-----------------------------------\n", sep = ""))
			print(round(object@mfit$matcoef,6), digits = 5)
			itest = rugarch:::.information.test(object@mfit$llh, nObs = T, nPars =  npvar + (m^2 - m)/2 + length(object@mfit$matcoef[,1]))
			itestm = matrix(0, ncol = 1, nrow = 4)
			itestm[1,1] = itest$AIC
			itestm[2,1] = itest$BIC
			itestm[3,1] = itest$SIC
			itestm[4,1] = itest$HQIC
			colnames(itestm) = ""
			rownames(itestm) = c("Akaike", "Bayes", "Shibata", "Hannan-Quinn")
			cat("\nInformation Criteria")
			cat(paste("\n---------------------\n", sep = ""))
			print(itestm,digits=5)
			cat("\n")
			cat("\nElapsed time :", object@mfit$timer,"\n\n")
			invisible(object)
		})
		
# filter show
setMethod("show",
		signature(object = "DCCfilter"),
		function(object){
			m = dim(object@model$modeldata$data)[2]
			T = object@model$modeldata$T
			cat(paste("\n*------------------------------------*", sep = ""))
			cat(paste("\n*          DCC GARCH Filter          *", sep = ""))
			cat(paste("\n*------------------------------------*", sep = ""))	
			cat("\n\nDistribution         : ", object@model$modeldesc$distribution)
			if(object@model$DCC=="FDCC"){
				cat("\nModel                : ", paste(object@model$DCC, "(", object@model$modelinc[4], ",", object@model$modelinc[5],")", sep=""))
				cat("\nNo.Groups            : ", object@model$modelinc[3])
			} else{
				cat("\nModel                : ", paste(object@model$DCC, "(", object@model$modelinc[3], ",", object@model$modelinc[4],")", sep=""))
			}
			# Need to show parameters as Q + DCC + GARCH
			if(object@model$modelinc[1]>0){
				npvar = dim(object@model$varcoef)[1] * dim(object@model$varcoef)[2]
			} else{
				npvar = 0
			}
			NP = paste("[",npvar, "+", length(object@mfilter$garchnames),"+",length(object@mfilter$dccnames), "+",(m^2 - m)/2,"]", sep="") 			
			cat("\nNo. of Parameters    : ", length(object@mfilter$matcoef[,1]) + ( (m^2 - m)/2 ))
			cat("\n[VAR GARCH DCC UncQ] :", NP)
			cat("\nNo. of Series        : ", m)
			cat("\nNo. of Obs.          : ", T)
			cat("\nLog-Likelihood       : ", object@mfilter$llh)
			cat("\nAv.Log-Likelihood    : ", round(object@mfilter$llh/T, 2), "\n")
			cat("\nParameters")
			cat(paste("\n--------------------------------------\n", sep = ""))
			print(round(object@mfilter$matcoef[,1,drop=FALSE],6), digits = 5)
			cat("\n")
			cat("\nElapsed time :", object@mfilter$timer,"\n\n")
			invisible(object)
		})
# forecast show
setMethod("show",
		signature(object = "DCCforecast"),
		function(object){
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\n*       DCC GARCH Forecast        *", sep = ""))
			cat(paste("\n*---------------------------------*", sep = ""))
			cat("\n\nDistribution         : ", object@model$modeldesc$distribution)
			if(object@model$DCC=="FDCC"){
				cat("\nModel                : ", paste(object@model$DCC, "(", object@model$modelinc[4], ",", object@model$modelinc[5],")", sep=""))
				cat("\nNo.Groups            : ", object@model$modelinc[3])
			} else{
				cat("\nModel                : ", paste(object@model$DCC, "(", object@model$modelinc[3], ",", object@model$modelinc[4],")", sep=""))
			}
			n.ahead = object@model$n.ahead
			cat("\nHorizon              : ", n.ahead)
			cat("\nRoll Steps           : ", object@model$n.roll)
			cat("\n-----------------------------------")
			cat("\n\n0-roll forecast: \n")
			forc = object@mforecast$R[[1]]
			if( dim(forc)[3] > 5 ){
				cat("\nFirst 2 Correlation Forecasts\n")
				print(forc[, , 1:2], digits = 4)
				cat(paste(rep(".", dim(forc[,,1])[1], collapse = TRUE)))
				cat("\n")
				cat(paste(rep(".", dim(forc[,,1])[1], collapse = TRUE)))
				cat("\n")
				cat("\nLast 2 Correlation Forecasts\n")
				print(last(forc, 2), digits = 4)
			} else{
				print(forc, digits = 4)
			}
			cat("\n\n")
			invisible(object)
		})

setMethod("show",
		signature(object = "DCCsim"),
		function(object){
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\n*      DCC GARCH Simulation       *", sep = ""))
			cat(paste("\n*---------------------------------*", sep = ""))
			cat("\n\nDistribution\t\t:", object@msim$model$distribution)
			if(object@model$DCC=="FDCC"){
				cat("\nModel\t\t: ", paste(object@model$DCC, "(", object@model$modelinc[4], ",", object@model$modelinc[5],")", sep=""))
				cat("\nNo.Groups\t: ", object@model$modelinc[3])
			} else{
				cat("\nModel\t\t: ", paste(object@model$DCC, "(", object@model$modelinc[3], ",", object@model$modelinc[4],")", sep=""))
			}
			cat(paste("\nSimulation Horizon\t: ",  object@msim$model$n.sim, sep = ""))
			cat(paste("\nBurn In\t\t\t\t: ",  object@msim$model$n.start, sep = ""))
			cat(paste("\nNo. of Simulations\t: ",object@msim$model$m.sim, sep = ""))
			cat("\n\n")
			invisible(object)
		})

setMethod("show",
		signature(object = "DCCroll"),
		function(object){
			cat(paste("\n*---------------------------------------------------*", sep = ""))
			cat(paste("\n*                    DCC GARCH Roll                 *", sep = ""))
			cat(paste("\n*---------------------------------------------------*", sep = ""))
			cat("\n\nDistribution\t\t:", object@model$modeldesc$distribution)
			cat(paste("\nSimulation Horizon\t: ",  object@model$forecast.length, sep = ""))
			cat(paste("\nRefits\t\t\t\t: ",  object@model$n.refits, sep = ""))
			cat(paste("\nWindow\t\t\t\t: ",  object@model$refit.window, sep = ""))
			cat(paste("\nNo.Assets\t\t\t: ", length(object@model$modeldata$asset.names), sep = ""))
			cat("\n\nOptimal Parameters Across Rolls (First 2, Last 2)")
			cat(paste("\n---------------------------------------------------\n", sep = ""))
			tmp = coef(object)
			if(dim(tmp)[2]>4){
				n = dim(tmp)[2]
				print(round(cbind(tmp[,1:2], tmp[,(n-1):n]), 4))
			} else{
				print(round(tmp, 4))
			}
			cat("\n")
			invisible(object)
		})

#----------------------------------------------------------------------------------
# plot
#----------------------------------------------------------------------------------
setMethod(f = "plot", signature(x = "DCCfit", y = "missing"), .plotdccfit)
setMethod(f = "plot", signature(x = "DCCfilter", y = "missing"), .plotdccfit)
setMethod(f = "plot", signature(x = "DCCforecast", y = "missing"), .plotdccforecast)
#setMethod(f = "plot", signature(x = "DCCsim", y = "missing"), .plotdccsim)
setMethod(f = "plot", signature(x = "DCCroll", y = "missing"), .plotdccroll)


#----------------------------------------------------------------------------------
# First and Last object method
#----------------------------------------------------------------------------------
first = function(x, index = 1, ...)
{
	UseMethod("first")
}

.first.array = function(x, index = 1)
{
	T = dim(x)[3]
	if( index > T | index < 1 ) stop("\nindex out of bounds")
	x[, , 1:index, drop = FALSE]
}

setMethod("first", signature(x = "array"), .first.array)

last = function(x, index = 1, ...)
{
	UseMethod("last")
}

.last.array = function(x, index = 1)
{
	T = dim(x)[3]
	if( index > T | index < 1 ) stop("\nindex out of bounds")
	x[, , (T - index + 1):T, drop = FALSE]
}

setMethod("last", signature(x = "array"), .last.array)

#----------------------------------------------------------------------------------
# nisurface
#----------------------------------------------------------------------------------
.newsimpact.dcc = function(object, type = "cov", pair = c(1,2), plot = TRUE, plot.type = c("surface", "contour"))
{
	if(object@model$DCC=="FDCC"){
		ans = switch(type,
				cov = .newsimpact.fdcc.cov(object = object, pair = pair, plot = plot, type = plot.type),
				cor = .newsimpact.fdcc.cor(object = object, pair = pair, plot = plot, type = plot.type)
		)
	} else{
		ans = switch(type,
				cov = .newsimpact.dcc.cov(object = object, pair = pair, plot = plot, type = plot.type),
				cor = .newsimpact.dcc.cor(object = object, pair = pair, plot = plot, type = plot.type)
		)
	}
	return(ans)
}

setMethod("nisurface", signature(object = "DCCfit"), .newsimpact.dcc)
setMethod("nisurface", signature(object = "DCCfilter"), .newsimpact.dcc)

.newsimpact.dcc.cor = function(object, pair = c(1,2), plot = TRUE, type = c("surface", "contour")){

	if(is(object, "DCCfilter")){
		Z = object@mfilter$stdresid
		cnames = object@model$modeldata$asset.names
		m = dim(object@model$umodel$modelinc)[2]
		maxz = round(max(apply(object@mfilter$stdresid, 1, "max")) + 1, 0)
		minz = round(min(apply(object@mfilter$stdresid, 1, "min")) - 1, 0)
		zseq = seq(minz, maxz, length.out = 100)
		idx = object@model$pidx
		dcca = object@model$pars[idx["dcca",1], 1]
		dccg = object@model$pars[idx["dccg",1], 1]
		dccb = object@model$pars[idx["dccb",1], 1]
		U = object@mfilter$Qbar*(1 - dcca - dccb) - dccg*object@mfilter$Nbar
		Qbar = object@mfilter$Qbar
	} else{
		Z = object@mfit$stdresid
		cnames = object@model$modeldata$asset.names
		m = dim(object@model$umodel$modelinc)[2]
		maxz = round(max(apply(object@mfit$stdresid, 1, "max")) + 1, 0)
		minz = round(min(apply(object@mfit$stdresid, 1, "min")) - 1, 0)
		zseq = seq(minz, maxz, length.out = 100)
		idx = object@model$pidx
		dcca = object@model$pars[idx["dcca",1], 1]
		dccg = object@model$pars[idx["dccg",1], 1]
		dccb = object@model$pars[idx["dccb",1], 1]
		U = object@mfit$Qbar*(1 - dcca - dccb) - dccg*object@mfit$Nbar
		Qbar = object@mfit$Qbar
	}

	ni = matrix(0, 100, 100)
	for(i in 1:100){
		for(j in 1:100){
			z = za = matrix(0, ncol = m, nrow = 1)
			z[1,pair[1]] = zseq[i]
			z[1,pair[2]] = zseq[j]
			za = z * .asymI(z)
			tmp = U + dcca*t(z)%*%z + dccg*t(za)%*%za + dccb * Qbar
			tmp = tmp/(sqrt(diag(tmp)) %*% t(sqrt(diag(tmp))) )
			ni[i,j] = tmp[pair[1],pair[2]]
		}
	}
	type == type[1]
	if(plot){
		if(tolower(type[1]) == "surface"){
			x1 = shape::drapecol(ni, col = shape::femmecol(100), NAcol = "white")
			persp(  x = zseq,
					y = zseq,
					z = ni,  col = x1, theta = 45, phi = 25, expand = 0.5,
					ltheta = 120, shade = 0.75, ticktype = "detailed", xlab = "shock[z_1]",
					ylab = "shock[z_2]", zlab = "cor",
					cex.axis = 0.7,  cex.main = 0.8, main = paste("DCC News Impact Correlation Surface\n", 
							cnames[pair[1]], "-", cnames[pair[2]], sep = ""))
		} else{
			tcol <- terrain.colors(12)
			contour(x = zseq,
					y = zseq,
					z = ni,  col = tcol[2], lty = "solid", cex.axis = 0.7,  cex.main = 0.8,
					main = paste("DCC News Impact Correlation Contour\n", 
							cnames[pair[1]], "-", cnames[pair[2]], sep = ""))
		}
	}
	return(list(nisurface = ni, axis = zseq))
}


.newsimpact.dcc.cov = function(object, pair = c(1,2), plot = TRUE, type = c("surface", "contour")){
	# for covariance we need the unconditional variances of the univariate model
	
	umodel = object@model$umodel
	m = dim(umodel$modelinc)[2]
	fpars = lapply(1:m, FUN = function(i) object@model$mpars[object@model$midx[,i]==1,i])
	
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			fpars, NULL)
	D = rep(0, m)
	for(i in 1:m) D[i] = sqrt( as.numeric( uncvariance(mspec@spec[[i]])))
	

	if(is(object, "DCCfilter")){
		Z = object@mfilter$stdresid
		cnames = object@model$modeldata$asset.names
		maxz = round(max(apply(object@mfilter$stdresid, 1, "max")) + 1, 0)
		minz = round(min(apply(object@mfilter$stdresid, 1, "min")) - 1, 0)
		zseq = seq(minz, maxz, length.out = 100)
		idx = object@model$pidx
		dcca = object@model$pars[idx["dcca",1], 1]
		dccg = object@model$pars[idx["dccg",1], 1]
		dccb = object@model$pars[idx["dccb",1], 1]
		U = object@mfilter$Qbar*(1 - dcca - dccb) - dccg*object@mfilter$Nbar
		Qbar = object@mfilter$Qbar
		
	} else{
		Z = object@mfit$stdresid
		cnames = object@model$modeldata$asset.names
		maxz = round(max(apply(object@mfit$stdresid, 1, "max")) + 1, 0)
		minz = round(min(apply(object@mfit$stdresid, 1, "min")) - 1, 0)
		zseq = seq(minz, maxz, length.out = 100)
		idx = object@model$pidx
		dcca = object@model$pars[idx["dcca",1], 1]
		dccg = object@model$pars[idx["dccg",1], 1]
		dccb = object@model$pars[idx["dccb",1], 1]
		U = object@mfit$Qbar*(1 - dcca - dccb) - dccg*object@mfit$Nbar
		Qbar = object@mfit$Qbar
	}

	ni = matrix(0, 100, 100)
	for(i in 1:100){
		for(j in 1:100){
			z = za = matrix(0, ncol = m, nrow = 1)
			z[1,pair[1]] = zseq[i]
			z[1,pair[2]] = zseq[j]
			za = z * .asymI(z)
			tmp = U + dcca*t(z)%*%z + dccg*t(za)%*%za + dccb * Qbar
			tmp = tmp/(sqrt(diag(tmp)) %*% t(sqrt(diag(tmp))) )
			H = diag(D)%*%tmp%*%diag(D)
			ni[i,j] = H[pair[1],pair[2]]
		}
	}
	type == type[1]
	if(plot){
		if(tolower(type[1]) == "surface"){
			x1 = shape::drapecol(ni, col = shape::femmecol(100), NAcol = "white")
			persp(  x = zseq,
					y = zseq,
					z = ni,  col = x1, theta = 45, phi = 25, expand = 0.5,
					ltheta = 120, shade = 0.75, ticktype = "detailed", xlab = "shock[z_1]",
					ylab = "shock[z_2]", zlab = "cov",
					cex.axis = 0.7,  cex.main = 0.8, main = paste("DCC News Impact Covariance Surface\n", 
							cnames[pair[1]], "-", cnames[pair[2]], sep = ""))
		} else{
			tcol <- terrain.colors(12)
			contour(x = zseq,
					y = zseq,
					z = ni,  col = tcol[2], lty = "solid", cex.axis = 0.7,  cex.main = 0.8,
					main = paste("DCC News Impact Covariance Contour\n", 
							cnames[pair[1]], "-", cnames[pair[2]], sep = ""))
		}
	}
	return(list(nisurface = ni, axis = zseq))
}

.newsimpact.fdcc.cor = function(object, pair = c(1,2), plot = TRUE, type = c("surface", "contour"))
{
	tmp = getfdccpars(object@model$pars, object@model)
	fC = tmp$C
	fA = tmp$A
	fB = tmp$B
	if(is(object, "DCCfilter")){
		Z = object@mfilter$stdresid
		modelinc = object@model$modelinc
		cnames = object@model$modeldata$asset.names
		m = dim(object@model$umodel$modelinc)[2]
		maxz = round(max(apply(object@mfilter$stdresid, 1, "max")) + 1, 0)
		minz = round(min(apply(object@mfilter$stdresid, 1, "min")) - 1, 0)
		zseq = seq(minz, maxz, length.out = 100)
		Qbar = object@mfilter$Qbar
	} else{
		Z = object@mfit$stdresid
		cnames = object@model$modeldata$asset.names
		m = dim(object@model$umodel$modelinc)[2]
		maxz = round(max(apply(object@mfit$stdresid, 1, "max")) + 1, 0)
		minz = round(min(apply(object@mfit$stdresid, 1, "min")) - 1, 0)
		zseq = seq(minz, maxz, length.out = 100)
		Qbar = object@mfit$Qbar
	}
	U = (fC %*% t(fC)) * Qbar
	
	ni = matrix(0, 100, 100)
	for(i in 1:100){
		for(j in 1:100){
			z = za = matrix(0, ncol = m, nrow = 1)
			z[1,pair[1]] = zseq[i]
			z[1,pair[2]] = zseq[j]
			tmp = U + (fA[,1] %*% t(fA[,1]))* (t(z)%*%z) + (fB[,1]%*%t(fB[,1])) * Qbar
			tmp = tmp/(sqrt(diag(tmp)) %*% t(sqrt(diag(tmp))) )
			ni[i,j] = tmp[pair[1],pair[2]]
		}
	}
	type == type[1]
	if(plot){
		if(tolower(type[1]) == "surface"){
			x1 = shape::drapecol(ni, col = shape::femmecol(100), NAcol = "white")
			persp(  x = zseq,
					y = zseq,
					z = ni,  col = x1, theta = 45, phi = 25, expand = 0.5,
					ltheta = 120, shade = 0.75, ticktype = "detailed", xlab = "shock[z_1]",
					ylab = "shock[z_2]", zlab = "cor",
					cex.axis = 0.7,  cex.main = 0.8, main = paste("DCC News Impact Correlation Surface\n", 
							cnames[pair[1]], "-", cnames[pair[2]], sep = ""))
		} else{
			tcol <- terrain.colors(12)
			contour(x = zseq,
					y = zseq,
					z = ni,  col = tcol[2], lty = "solid", cex.axis = 0.7,  cex.main = 0.8,
					main = paste("DCC News Impact Correlation Contour\n", 
							cnames[pair[1]], "-", cnames[pair[2]], sep = ""))
		}
	}
	return(list(nisurface = ni, axis = zseq))
}

.newsimpact.fdcc.cov = function(object, pair = c(1,2), plot = TRUE, type = c("surface", "contour")){
	# for covariance we need the unconditional variances of the univariate model
	
	umodel = object@model$umodel
	m = dim(umodel$modelinc)[2]
	fpars = lapply(1:m, FUN = function(i) object@model$mpars[object@model$midx[,i]==1,i])
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			fpars, NULL)
	D = rep(0, m)
	for(i in 1:m) D[i] = sqrt( as.numeric( uncvariance(mspec@spec[[i]])))
	
	tmp = getfdccpars(object@model$pars, object@model)
	fC = tmp$C
	fA = tmp$A
	fB = tmp$B
	if(is(object, "DCCfilter")){
		Z = object@mfilter$stdresid
		cnames = object@model$modeldata$asset.names
		maxz = round(max(apply(object@mfilter$stdresid, 1, "max")) + 1, 0)
		minz = round(min(apply(object@mfilter$stdresid, 1, "min")) - 1, 0)
		zseq = seq(minz, maxz, length.out = 100)
		Qbar = object@mfilter$Qbar
		
	} else{
		Z = object@mfit$stdresid
		cnames = object@model$modeldata$asset.names
		maxz = round(max(apply(object@mfit$stdresid, 1, "max")) + 1, 0)
		minz = round(min(apply(object@mfit$stdresid, 1, "min")) - 1, 0)
		zseq = seq(minz, maxz, length.out = 100)
		Qbar = object@mfit$Qbar
	}
	U = (fC %*% t(fC)) * Qbar
	
	ni = matrix(0, 100, 100)
	for(i in 1:100){
		for(j in 1:100){
			z = za = matrix(0, ncol = m, nrow = 1)
			z[1,pair[1]] = zseq[i]
			z[1,pair[2]] = zseq[j]
			tmp = U + (fA[,1] %*% t(fA[,1]))* (t(z)%*%z) + (fB[,1]%*%t(fB[,1])) * Qbar
			tmp = tmp/(sqrt(diag(tmp)) %*% t(sqrt(diag(tmp))) )
			H = diag(D)%*%tmp%*%diag(D)
			ni[i,j] = H[pair[1],pair[2]]
		}
	}
	type == type[1]
	if(plot){
		if(tolower(type[1]) == "surface"){
			x1 = shape::drapecol(ni, col = shape::femmecol(100), NAcol = "white")
			persp(  x = zseq,
					y = zseq,
					z = ni,  col = x1, theta = 45, phi = 25, expand = 0.5,
					ltheta = 120, shade = 0.75, ticktype = "detailed", xlab = "shock[z_1]",
					ylab = "shock[z_2]", zlab = "cov",
					cex.axis = 0.7,  cex.main = 0.8, main = paste("DCC News Impact Covariance Surface\n", 
							cnames[pair[1]], "-", cnames[pair[2]], sep = ""))
		} else{
			tcol <- terrain.colors(12)
			contour(x = zseq,
					y = zseq,
					z = ni,  col = tcol[2], lty = "solid", cex.axis = 0.7,  cex.main = 0.8,
					main = paste("DCC News Impact Covariance Contour\n", 
							cnames[pair[1]], "-", cnames[pair[2]], sep = ""))
		}
	}
	return(list(nisurface = ni, axis = zseq))
}