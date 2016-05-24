#################################################################################
##
##   R package rugarch by Alexios Ghalanos Copyright (C) 2008-2015.
##   This file is part of the R package rugarch.
##
##   The R package rugarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rugarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

#---------------------------------------------------------------------------------
# SECTION sGARCH fit
#---------------------------------------------------------------------------------
.arfimafit = function(spec, data, out.sample = 0, solver = "solnp", solver.control = list(), 
		fit.control = list(fixed.se = 0, scale = 0), 
		numderiv.control = list(grad.eps=1e-4, grad.d=0.0001, grad.zero.tol=sqrt(.Machine$double.eps/7e-7),
				hess.eps=1e-4, hess.d=0.1, hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2))
{
	tic = Sys.time()
	
	default.numd = list(grad.eps=1e-4, grad.d=0.0001, grad.zero.tol=sqrt(.Machine$double.eps/7e-7),
			hess.eps=1e-4, hess.d=0.1, hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2)
	idx1 = na.omit(match(names(numderiv.control), names(default.numd)))
	idx2 = na.omit(match(names(default.numd), names(numderiv.control)))
	if(length(idx1)>0){
		default.numd[idx1] = numderiv.control[idx2]
	}
	numderiv.control = default.numd
	
	if(is.null(solver.control$trace)) trace = 0 else trace = solver.control$trace
	
	if(is.null(fit.control$stationarity)) fit.control$stationarity = TRUE
	if(is.null(fit.control$fixed.se)) fit.control$fixed.se = FALSE
	if(is.null(fit.control$scale)) fit.control$scale = FALSE
	# if we have external regressors in variance turn off scaling
	if(spec@model$modelinc[15] > 0) fit.control$scale = FALSE
	# if we have arch-in-mena turn off scaling
	if(spec@model$modelinc[5] > 0) fit.control$scale = FALSE
	
	# if there are fixed pars we do no allow scaling as there would be no way of mixing scaled
	# amd non scaled parameters	
	if(sum(spec@model$pars[,2]) > 0) fit.control$scale = FALSE
	xdata = .extractdata(data)
	if(!is.numeric(out.sample)) stop("\narfimafit-->error: out.sample must be numeric\n")
	if(as.numeric(out.sample)<0) stop("\narfimafit-->error: out.sample must be positive\n")
	n.start = round(out.sample,0)
	n = length(xdata$data)
	#if((n-n.start)<100) stop("\narfimafit-->error: function requires at least 100 data\n points to run\n")
	data = xdata$data[1:(n-n.start)]
	index = xdata$index[1:(n-n.start)]
	origdata = xdata$data
	origindex = xdata$index
	period = xdata$period
	# create a temporary environment to store values (deleted at end of function)
	garchenv = new.env(hash = TRUE)
	arglist = list()
	arglist$garchenv <- garchenv
	arglist$pmode = 0
	model = spec@model
	modelinc = model$modelinc
	pidx = model$pidx
	# expand the spec object and assign spec lists
	if(modelinc[6] > 0){
		mexdata = model$modeldata$mexdata[1:(n-n.start), , drop = FALSE]
	} else{
		mexdata = NULL
	}
	
	arglist$index = index
	arglist$trace = trace
	arglist$fit.control = fit.control
	m =  model$maxOrder
	model$modeldata$T = T = length(as.numeric(data))
	dist = model$modeldesc$distribution
	if(fit.control$scale) dscale = sd(data) else dscale = 1
	zdata = data/dscale
	arglist$data = zdata
	arglist$dscale = dscale
	arglist$model = model
	ipars = model$pars
	# Optimization Starting Parameters Vector & Bounds
	ipars = .arfimastart(ipars, arglist = arglist)
	arglist$ipars = ipars
	# we now split out any fixed parameters
	estidx = as.logical( ipars[,4] )
	arglist$estidx = estidx	
	npars = sum(estidx)
	if(any(ipars[,2]==1)){
		if(npars == 0){
			if(fit.control$fixed.se==0) {
				# if all parameters are fixed an no standard erros are to
				# be calculated then we return a ugarchfilter object
				warning("\narfimafit-->warning: all parameters fixed...returning ugarchfilter 
								object instead\n")
				return(arfimafilter(data = xts(origdata, origindex), spec = spec, out.sample = out.sample))
			} else{
				# if all parameters are fixed but we require standard errors, we
				# skip the solver
				use.solver = 0
				ipars[ipars[,2]==1, 4] = 1
				ipars[ipars[,2]==1, 2] = 0
				arglist$ipars = ipars
				estidx = as.logical( ipars[,4] )
				arglist$estidx = estidx	
			}
		} else{
			# with some parameters fixed we extract them (to be rejoined at end)
			# so that they do not enter the solver
			use.solver = 1
		}
	} else{
		use.solver = 1
	}
	
	# start counter
	assign("rugarch_llh", 1, envir = garchenv)

	if(fit.control$stationarity==1 && modelinc[2]>0){
		Ifn = function(pars, arglist){
			arx = 0.5
			modelinc = arglist$model$modelinc
			idx = arglist$model$pidx
			ipars = arglist$ipars
			estidx = arglist$estidx
			ipars[estidx, 1] = pars		
			if(modelinc[2] > 0 && modelinc[4] == 0) arx = ipars[idx["ar",1]:idx["ar",2],1]
			min(Mod(polyroot(c(1, -arx))))
		}
		ILB = 1.01
		IUB = 100
		if(solver == "solnp" | solver == "gosolnp" | solver == "hybrid") fit.control$stationarity = 0
	} else{
		Ifn = ILB = IUB = NULL
	}
	arglist$fit.control = fit.control
	# need to control for use of parallel environment which cannot handle the
	# calls between the locally stored environment
	if(use.solver){
		parscale = rep(1, length = npars)
		names(parscale) = rownames(ipars[estidx,])
		if(modelinc[1] > 0) parscale["mu"] = abs(mean(zdata))
		if(modelinc[7] > 0) parscale["sigma"] = sd(zdata)
		arglist$returnType = "llh"
		solution = .garchsolver(solver, pars = ipars[estidx, 1], fun = .arfimaLLH, 
				Ifn, ILB, IUB, gr = NULL, hessian = NULL, parscale = parscale, 
				control = solver.control, LB = ipars[estidx, 5], UB = ipars[estidx, 6], 
				ux = NULL, ci = NULL, mu = NULL, arglist = arglist)
		sol = solution$sol
		hess = solution$hess
		timer = Sys.time()-tic
		if(!is.null(sol$par)) ipars[estidx, 1] = sol$par else ipars[estidx, 1] = NA
		if(sum(ipars[,2]) == 0){
			if(modelinc[1] > 0) ipars[pidx["mu",1]:pidx["mu",2], 1] = ipars[pidx["mu",1]:pidx["mu",2], 1] * dscale
			if(modelinc[5] > 0){
				ipars[pidx["mxreg", 1]:pidx["mxreg", 2], 1] = ipars[pidx["mxreg", 1]:pidx["mxreg", 2], 1] * dscale
			}
			ipars[pidx["sigma",1], 1] = ipars[pidx["sigma",1],1] * dscale
		}
		arglist$ipars = ipars
		convergence = sol$convergence
	} else{
		hess = NULL
		timer = Sys.time()-tic
		convergence = 0
		sol = list()
		sol$message = "all parameters fixed"
	}
	fit = list()
	# check convergence else write message/return
	ipars2 = ipars
	if(convergence == 0){
		arglist$dscale = 1
		if(sum(ipars[,2]) > 0 && fit.control$fixed.se == 1){
			ipars[ipars[,2]==1, 4] = 1
			ipars[ipars[,2]==1, 2] = 0
			arglist$ipars = ipars
			estidx = as.logical( ipars[,4] )
			arglist$estidx = estidx	
		}
		arglist$data = data
		fit = .makearfimafitmodel(f = .arfimaLLH, T = T, m = m, timer = timer, 
				convergence = convergence, message = sol$message, hess, 
				arglist = arglist, numderiv.control = numderiv.control)
		model$modelinc[7] = modelinc[7]
		model$modeldata$data = origdata
		model$modeldata$index = origindex
		model$modeldata$period = period
		
		model$pars = ipars
		model$pars[, 1] = fit$ipars[,1]
		fit$ipars[, 4] = ipars2[, 4]
		fit$ipars[, 2] = ipars2[, 2]
	} else{
		fit$message = sol$message
		fit$convergence = 1
		model$modeldata$data = origdata
		model$modeldata$index = origindex
		model$modeldata$period = period
	}
	
	# make model list to return some usefule information which
	# will be called by other functions (show, plot, sim etc)
	model$n.start = n.start
	fit$solver = solution
	ans = new("ARFIMAfit",
			fit = fit,
			model = model)
	return(ans)
}

#---------------------------------------------------------------------------------
# SECTION sGARCH LLH
#---------------------------------------------------------------------------------
.arfimaLLH = function(pars, arglist)
{
	# prepare inputs
	# rejoin fixed and pars
	data = arglist$data
	returnType = arglist$returnType
	garchenv = arglist$garchenv
	# rejoin fixed and pars
	model = arglist$model
	estidx = arglist$estidx
	idx = model$pidx
	ipars = arglist$ipars	
	ipars[estidx, 1] = pars
	trace = arglist$trace
	T = length(data)
	
	fit.control = arglist$fit.control	
	m = model$maxOrder
	N = c(m,T)
	mexdata = model$modeldata$mexdata[1:T,, drop = FALSE]
	distribution = model$modeldesc$distribution
	modelinc = model$modelinc
	# distribution number
	# 1 = norm, 2=snorm, 3=std, 4=sstd, 5=ged, 6=sged, 7=nig
	dist = model$modeldesc$distno
	hm = 0
	if(modelinc[4]>0){
		rx = .arfimaxfilter(modelinc[1:21], ipars[,1], idx, mexdata = mexdata, h = 0, data = data, N = N)
		res = rx$res
		zrf = rx$zrf
		res[is.na(res) | !is.finite(res) | is.nan(res)] = 0
	} else{
		res = rep(0, T)
		zrf = 0
	}
	# unconditional sigma value
	if(fit.control$stationarity == 1 && modelinc[4] == 0 && modelinc[2] > 0){
		arx = ipars[idx["ar",1]:idx["ar",2],1]
		kappa = min(Mod(polyroot(c(1, -arx))))
		if(is.na(kappa) || (kappa < 1 || kappa > 100) ){
			if(arglist$pmode!=1){
				return(llh = get("rugarch_llh", garchenv) + 0.1*(abs(get("rugarch_llh", garchenv))))
			} else{
				return(llh = 1e10)
			}
		}
	}
	if(modelinc[6]>0) mexdata = as.double(as.vector(mexdata)) else mexdata = double(1)

	ans = try( .C("arfimafitC", model = as.integer(modelinc[1:21]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					x = as.double(data), res = as.double(res), mexdata = mexdata, 
					zrf = as.double(zrf), constm = double(T), condm = double(T), 
					m = as.integer(m), T = as.integer(T), z = double(T), 
					llh = double(1), LHT = double(T), PACKAGE = "rugarch"), silent = TRUE )
	
	if(inherits(ans, "try-error")){
		if(arglist$pmode!=1){
			assign("rugarch_csol", 1, envir = garchenv)
			assign("rugarch_filtermessage", ans, envir = garchenv)
			if( trace > 0 ) cat(paste("\narfimafit-->warning: ", get("rugarch_filtermessage", garchenv),"\n", sep=""))
			return(llh = (get("rugarch_llh", garchenv) + 0.1*(abs(get("rugarch_llh", garchenv)))))
		} else{
			return(llh = 1e10)
		}
	} else{
		if(arglist$pmode!=1){
			assign("rugarch_csol", 0, envir = garchenv)
		}	
	}
	z = ans$z
	epsx = ans$res
	llh = ans$llh
	
	if(is.finite(llh) && !is.na(llh) && !is.nan(llh)){
		if(arglist$pmode!=1) assign("rugarch_llh", llh, envir = garchenv)
	} else {
		if(arglist$pmode!=1) llh = (get("rugarch_llh", garchenv) + 0.1*(abs(get("rugarch_llh", garchenv)))) else llh = 1e10
	}
	
	# LHT = raw scores
	# LHT = -ans$LHT[(m):T]
	LHT = -ans$LHT
	ans = switch(returnType,
			llh = llh,
			LHT = LHT,
			all = list(llh = llh, res = epsx, z = z, LHT = LHT))
	return(ans)
}

#---------------------------------------------------------------------------------
# SECTION sGARCH filter
#---------------------------------------------------------------------------------
.arfimafilter = function(spec, data, out.sample = 0, n.old = NULL)
{
	# n.old is optional and indicates the length of the original dataseries (in
	# cases when this represents a dataseries augmented by newer data). The reason
	# for using this is so that the old and new datasets agree since the original
	# recursion uses the sum of the residuals to start the recursion and therefore
	# is influenced by new data. For a small augmentation the values converge after
	# x periods, but it is sometimes preferable to have this option so that there is
	# no forward looking information contaminating the study.
	xdata = .extractdata(data)
	data = xdata$data
	index = xdata$index
	period = xdata$period
	origdata = data
	origindex = index
	T = length(origdata)  - out.sample
	data = origdata[1:T]
	index = origindex[1:T]
	if(!is.null(n.old)) Nx = n.old else Nx = length(data)
	
	model = spec@model
	ipars = model$pars
	pars = unlist(model$fixed.pars)
	parnames = names(pars)
	modelnames = .checkallfixed(spec)
	if(is.na(all(match(modelnames, parnames), 1:length(modelnames)))) {
		cat("\narfimafilter-->error: parameters names do not match specification\n")
		cat("Expected Parameters are: ")
		cat(paste(modelnames))
		cat("\n")
		stop("Exiting", call. = FALSE)
	}
	
	# once more into the spec
	setfixed(spec)<-as.list(pars)
	model = spec@model
	model$modeldata$T = T
	ipars = model$pars
	idx = model$pidx
	modelinc = model$modelinc
	m = model$maxOrder
	N = c(m,T)
	mexdata = model$modeldata$mexdata[1:T, , drop = FALSE]
	distribution = model$modeldesc$distribution
	dist = model$modeldesc$distno
	
	
	if(modelinc[4]>0){
		rx = .arfimaxfilter(modelinc[1:21], ipars[,1], idx, mexdata = mexdata, h = 0, data = data, N = N)
		res = rx$res
		zrf = rx$zrf
		res[is.na(res) | !is.finite(res) | is.nan(res)] = 0
	} else{
		res = rep(0, T)
		zrf = 0
	}
	
	if(modelinc[6]>0) mexdata = as.double(as.vector(mexdata)) else mexdata = double(1)
	
	ans = try( .C("arfimafitC", model = as.integer(modelinc[1:21]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					x = as.double(data), res = as.double(res), mexdata = mexdata, 
					zrf = as.double(zrf), constm = double(T), condm = double(T), 
					m = as.integer(m), T = as.integer(T), z = double(T), 
					llh = double(1), LHT = double(T), PACKAGE = "rugarch"), silent = TRUE )
	
	if(inherits(ans, "try-error")) stop("\nugarchfilter-->error: problem in C code...exiting.")
	filter = list()
	filter$z = ans$z
	filter$residuals = ans$res
	filter$LLH = -ans$llh
	filter$log.likelihoods = ans$LHT
	filter$distribution = distribution
	filter$ipars = ipars
	model$modeldata$data = origdata
	model$modeldata$index = origindex
	model$modeldata$period = period
	model$n.start = out.sample
	
	sol = new("ARFIMAfilter",
			filter = filter,
			model = model)
	return(sol)
}

#---------------------------------------------------------------------------------
# SECTION sGARCH forecast
#---------------------------------------------------------------------------------
.arfimaforecast = function(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0, 
		external.forecasts = list(mregfor = NULL), ...)
{
	fit = fitORspec
	data   = fit@model$modeldata$data
	Nor    = length(as.numeric(data))
	index  = fit@model$modeldata$index
	period = fit@model$modeldata$period
	ns = fit@model$n.start
	N = Nor - ns
	model = fit@model
	ipars = fit@fit$ipars
	modelinc = model$modelinc
	idx = model$pidx
	if( n.roll > ns ) stop("\nugarchforecast-->error: n.roll must not be greater than out.sample!")
	pars = fit@fit$coef
	ipars = fit@fit$ipars
	
	# check if necessary the external regressor forecasts provided first
	xreg = .forcregressors(model = model, mregfor = external.forecasts$mregfor, vregfor = NULL, 
			n.ahead = n.ahead, N = Nor, out.sample = ns, n.roll = n.roll)
	mxf = xreg$mxf
	
	# filter data (check external regressor data - must equal length of origData)
	fcreq = ifelse(ns >= (n.ahead+n.roll), n.ahead+n.roll, ns)
	fspec = arfimaspec(
			mean.model = list(armaOrder = c(modelinc[2], modelinc[3]),
					include.mean = modelinc[1], arfima = modelinc[4], 
					external.regressors = mxf[1:(N + fcreq), , drop = FALSE]), 
			distribution.model = model$modeldesc$distribution, fixed.pars = as.list(pars))
	tmp =  xts(data[1:(N + fcreq)], index[1:(N + fcreq)])
	flt = .arfimafilter(data = tmp, spec = fspec, n.old = N)
	resfilter = flt@filter$residuals
	zfilter = flt@filter$z
	# forecast ARFIMA process
	seriesFor = matrix(NA, ncol = n.roll+1, nrow = n.ahead)
	colnames(seriesFor) = as.character(index[N:(N+n.roll)])
	rownames(seriesFor) = paste("T+", 1:n.ahead, sep="")
	for(i in 1:(n.roll+1)){
		np = N + i - 1
		if(modelinc[1] > 0){
			mu = rep(ipars[idx["mu",1]:idx["mu",2], 1], N+i+n.ahead-1)
		} else{
			mu = rep(0, N+i+n.ahead-1)
		}
		epsx = c(resfilter[1:(N+i-1)], rep(0, n.ahead))
		x = c(data[1:(N+i-1)], rep(0, n.ahead))
		z = c(zfilter[1:(N+i-1)], rep(0, n.ahead))
		# forecast of externals is provided outside the system
		mxfi = mxf[1:(N+i-1+n.ahead), , drop = FALSE]
		
		if(modelinc[4]>0){
			res = arfimaf(ipars, modelinc[1:21], idx, mu, mxfi, h=0, epsx, z, data = x, N = np, n.ahead)
		} else{
			res = armaf(ipars, modelinc[1:21], idx, mu, mxfi, h=0, epsx, z, data = x, N = np, n.ahead)
		}
		seriesFor[,i] = res[(np + 1):(np + n.ahead)]
	}
	fcst = list()
	fcst$n.ahead = n.ahead
	fcst$N = N+ns
	fcst$n.start = ns
	fcst$n.roll = n.roll
	fcst$seriesFor = seriesFor
	model$modeldata$residuals = flt@filter$residuals
	ans = new("ARFIMAforecast",
			forecast = fcst,
			model = model)
	return(ans)
}

#---------------------------------------------------------------------------------
# 2nd dispatch method for forecast
.arfimaforecast2 = function(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0, 
		external.forecasts = list(mregfor = NULL), ...)
{
	# first we filter the data to get the results:
	spec = fitORspec
	if(is.null(data)) stop("\nugarchforecast-->error: data must not be NULL when using a specification!")
	# we then extract the data/coefs etc
	xdata = .extractdata(data)
	Nor = length(as.numeric(xdata$data))
	data = xdata$data
	N = length(as.numeric(data))
	index = xdata$index
	period = xdata$period
	ns = out.sample
	if( n.roll > ns ) stop("\nugarchforecast-->error: n.roll must not be greater than out.sample!")
	N = Nor - ns
	
	model = spec@model
	ipars = model$pars
	pars = unlist(model$fixed.pars)
	parnames = names(pars)
	modelnames = .checkallfixed(spec)
	if(is.na(all(match(modelnames, parnames), 1:length(modelnames)))) {
		cat("\nugarchforecast-->error: parameters names do not match specification\n")
		cat("Expected Parameters are: ")
		cat(paste(modelnames))
		cat("\n")
		stop("Exiting", call. = FALSE)
	}
	# once more into the spec
	setfixed(spec)<-as.list(pars)
	model = spec@model
	idx = model$pidx
	ipars = model$pars
	modelinc = model$modelinc
	model$modeldata$data = data
	model$modeldata$index = index
	model$modeldata$period = period
	
	# check if necessary the external regressor forecasts provided first
	xreg = .forcregressors(model = model, mregfor = external.forecasts$mregfor, vregfor = NULL, 
			n.ahead = n.ahead, N = Nor, out.sample = ns, n.roll = n.roll)
	mxf = xreg$mxf
	vxf = xreg$vxf
	
	# filter data (check external regressor data - must equal length of origData)
	fcreq = ifelse(ns >= (n.ahead+n.roll), n.ahead+n.roll, ns)
	fspec = arfimaspec(
			mean.model = list(armaOrder = c(modelinc[2], modelinc[3]),
					include.mean = modelinc[1], arfima = modelinc[4], 
					external.regressors = mxf[1:(N + fcreq), , drop = FALSE]), 
			distribution.model = model$modeldesc$distribution, fixed.pars = as.list(pars))
	tmp =  xts(data[1:(N + fcreq)], index[1:(N + fcreq)])
	flt = .arfimafilter(data = tmp, spec = fspec, n.old = N)
	resfilter = flt@filter$residuals
	zfilter = flt@filter$z
	# forecast GARCH process
	seriesFor = matrix(NA, ncol = n.roll+1, nrow = n.ahead)
	colnames(seriesFor) = as.character(index[N:(N+n.roll)])
	rownames(seriesFor) = paste("T+", 1:n.ahead, sep="")

	
	for(i in 1:(n.roll+1)){
		np = N + i - 1
		if(modelinc[1] > 0){
			mu = rep(ipars[idx["mu",1]:idx["mu",2], 1], N+i+n.ahead-1)
		} else{
			mu = rep(0, N+i+n.ahead-1)
		}
		epsx = c(resfilter[1:(N+i-1)], rep(0, n.ahead))
		x = c(data[1:(N+i-1)], rep(0, n.ahead))
		z = c(zfilter[1:(N+i-1)], rep(0, n.ahead))
		# forecast of externals is provided outside the system
		mxfi = mxf[1:(N+i-1+n.ahead), , drop = FALSE]
		
		if(modelinc[4]>0){
			res = arfimaf(ipars, modelinc[1:21], idx, mu, mxfi, h=0, epsx, z, data = x, N = np, n.ahead)
		} else{
			res = armaf(ipars, modelinc[1:21], idx, mu, mxfi, h=0, epsx, z, data = x, N = np, n.ahead)
		}
		seriesFor[,i] = res[(np + 1):(np + n.ahead)]
	}
	fcst = list()
	fcst$n.ahead = n.ahead
	fcst$N = N+ns
	fcst$n.start = ns
	fcst$n.roll = n.roll
	fcst$seriesFor = seriesFor
	model$modeldata$residuals = flt@filter$residuals
	ans = new("ARFIMAforecast",
			forecast = fcst,
			model = model)
	return(ans)
	
}
#---------------------------------------------------------------------------------
# SECTION sGARCH simulate
#---------------------------------------------------------------------------------

.arfimasim = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, startMethod = 
				c("unconditional","sample"), prereturns = NA, 
		preresiduals = NA, rseed = NA, custom.dist = list(name = NA, distfit = NA, type = "z"), 
		mexsimdata = NULL)
{
	if(fit@model$modelinc[4]>0){
		if(n.start<fit@model$modelinc[3]){
			warning("\narfimasim-->warning: n.start>=MA order for arfima model...automatically setting.")
			n.start = fit@model$modelinc[3]
		}
	}
	# some checks
	if(is.na(rseed[1])){
		sseed = as.integer(runif(1,0,as.integer(Sys.time())))
	} else{
		if(length(rseed) != m.sim) sseed = as.integer(rseed[1]) else sseed = rseed[1:m.sim]
	}
	model = fit@model
	modelinc = model$modelinc
	# Enlarge Series:
	# need to allow for arfima case:
	if(modelinc[4]>0) n.start = max(modelinc[2]+modelinc[3], n.start)
	n = n.sim + n.start
	startMethod = startMethod[1]
	data = fit@model$modeldata$data
	N = length(as.numeric(data))
	data = data[1:(N - fit@model$n.start)]
	N = length(as.numeric(data))
	m = fit@model$maxOrder
	resids = fit@fit$residuals

	idx = model$pidx
	ipars = fit@fit$ipars
	# check if necessary the external regressor forecasts provided first
	xreg = .simregressors(model = model, mexsimdata = mexsimdata, vexsimdata = NULL, 
			N = N, n = n, m.sim = m.sim, m = m)	
	mexsim = xreg$mexsimlist
	Sigma = ipars[idx["sigma", 1], 1]
	if(N < n.start){
		startmethod[1] = "unconditional"
		warning("\narfimasim-->warning: n.start greater than length of data...using unconditional start method...\n")
	}
	
	# Random Samples from the Distribution
	# For the arfima method this is the residuals, NOT the standardized residuals
	if(length(sseed) == 1){
		zmatrix = data.frame(dist = model$modeldesc$distribution, lambda = ipars[idx["ghlambda",1], 1], 
				skew = ipars[idx["skew",1], 1], shape = ipars[idx["shape",1], 1], n = n * m.sim, seed = sseed[1])
		z = .custzdist(custom.dist, zmatrix, m.sim, n)
	} else{
		zmatrix = data.frame(dist = rep(model$modeldesc$distribution, m.sim), lambda = rep(ipars[idx["ghlambda",1], 1], m.sim), 
				skew = rep(ipars[idx["skew",1], 1], m.sim), shape = rep(ipars[idx["shape",1], 1], m.sim), 
				n = rep(n, m.sim), seed = sseed)
		z = .custzdist(custom.dist, zmatrix, m.sim, n)
	}
	if( is.matrix( custom.dist$distfit ) && custom.dist$type != "z") cres = TRUE else cres = FALSE
	
	if(startMethod == "unconditional"){
		z = rbind(matrix(0, nrow = m, ncol = m.sim), z)
	} else{
		z = rbind(matrix(tail(fit@fit$z, m), nrow = m, ncol = m.sim), z) 
	}
	
	if(!is.na(prereturns[1])){
		prereturns = as.vector(prereturns)
		if(length(prereturns)<m) stop(paste("\narfimasim-->error: prereturns must be of length ", m, sep=""))
	}
	if(!is.na(preresiduals[1])){
		preresiduals = as.vector(preresiduals)
		if(length(preresiduals)<m) stop(paste("\narfimasim-->error: preresiduals must be of length ", m, sep=""))
		preres = matrix(preresiduals, nrow = m)
	}
	if(is.na(prereturns[1])){
		if(startMethod[1] == "unconditional"){
			prereturns = as.numeric(rep(uncmean(fit), m))
		}
		else{
			prereturns = tail(data, m)
		}
	}
	
	# input vectors/matrices
	x = c(prereturns, rep(0, n))
	constm = matrix(ipars[idx["mu",1]:idx["mu",2], 1], ncol = m.sim, nrow = n + m)
	
	# outpus matrices
	seriesSim = matrix(0, ncol = m.sim, nrow = n.sim)
	residSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	
	for(i in 1:m.sim){
		if(any(is.na(preresiduals))){
			if(startMethod[1] == "sample"){
				preres 	= tail(resids, m)
				res 	= c(preres, if(cres) tail(z[(m+1):(n+m), i], n) else tail(z[(m+1):(n+m), i], n) * Sigma)
				residSim[,i] = tail(res, n.sim)
			} else{
				res 	= if(cres) as.numeric(z[,i]) else as.numeric(z[,i])*Sigma
				residSim[,i] = tail(res, n.sim)
			}
		}
		if(modelinc[6]>0){
			mxreg = matrix( ipars[idx["mxreg",1]:idx["mxreg",2], 1], ncol = modelinc[6] )
			constm[,i] = constm[,i] + mxreg %*%t( matrix( mexsim[[i]], ncol = modelinc[6] ) )
		}
		
		if(modelinc[4]>0){
			x = c(prereturns, rep(0, n +  modelinc[3]))
			res = c(if(cres) as.numeric(z[(m+1):(n+m),i]) else as.numeric(z[(m+1):(n+m),i]) * Sigma, if(modelinc[3]>0) rep(0,  modelinc[3]) else NULL)
			ans2 = .arfimaxsim(modelinc[1:21], ipars, idx, constm[1:n, i], res[1:(n+modelinc[3])], T = n)	
			seriesSim[,i] = head(ans2$series, n.sim)
		} else{
			ans2 = .armaxsim(modelinc[1:21], ipars, idx, constm[,i],  x, res, T = n + m, m)
			seriesSim[,i] = tail(ans2$x, n.sim)
		}
	}
	sim = list(seriesSim = seriesSim, residSim = residSim)
	sol = new("ARFIMAsim",
			simulation = sim,
			model = model,
			seed = as.integer(sseed))
	return(sol)
}


#---------------------------------------------------------------------------------
# SECTION sGARCH path
#---------------------------------------------------------------------------------

.arfimapath = function(spec, n.sim = 1000, n.start = 0, m.sim = 1, prereturns = NA, preresiduals = NA, 
		rseed = NA, custom.dist = list(name = NA, distfit = NA, type = "z"), mexsimdata = NULL)
{
	# some checks
	if(spec@model$modelinc[4]>0){
		if(n.start<spec@model$modelinc[3]){
			warning("\narfimapath-->warning: n.start>=MA order for arfima model...automatically setting.")
			n.start = spec@model$modelinc[3]
		}
	}
	if(is.na(rseed[1])){
		sseed = as.integer(runif(1,0,as.integer(Sys.time())))
	} else{
		if(length(rseed) != m.sim) sseed = as.integer(rseed[1]) else sseed = rseed[1:m.sim]
	}
	model = spec@model
	ipars = model$pars
	pars = unlist(model$fixed.pars)
	parnames = names(pars)
	modelnames = .checkallfixed(spec)
	if(is.na(all(match(modelnames, parnames), 1:length(modelnames)))) {
		cat("\narfimapath-->error: parameters names do not match specification\n")
		cat("Expected Parameters are: ")
		cat(paste(modelnames))
		cat("\n")
		stop("Exiting", call. = FALSE)
	}
	# once more into the spec
	setfixed(spec)<-as.list(pars)
	model = spec@model
	ipars = model$pars
	idx = model$pidx
	modelinc = model$modelinc
	Sigma = ipars[idx["sigma", 1], 1]
	# Enlarge Series:
	n = n.sim + n.start
	m = model$maxOrder
	N = 0
	if(modelinc[6]>0) {
		mexdata = matrix(model$modeldata$mexdata, ncol = modelinc[6])
		N = dim(mexdata)[1]
	} else { mexdata = NULL }
	distribution = model$modeldesc$distribution	
	# check if necessary the external regressor forecasts provided first
	xreg = .simregressors(model = model, mexsimdata = mexsimdata, vexsimdata = NULL, 
			N = N, n = n, m.sim = m.sim, m = m)	
	mexsim = xreg$mexsimlist
	
	if(length(sseed) == 1){
		zmatrix = data.frame(dist = model$modeldesc$distribution, lambda = ipars[idx["ghlambda",1], 1], 
				skew = ipars[idx["skew",1], 1], shape = ipars[idx["shape",1], 1], n = n * m.sim, seed = sseed[1])
		z = .custzdist(custom.dist, zmatrix, m.sim, n)
	} else{
		zmatrix = data.frame(dist = rep(model$modeldesc$distribution, m.sim), lambda = rep(ipars[idx["ghlambda",1], 1], m.sim), 
				skew = rep(ipars[idx["skew",1], 1], m.sim), shape = rep(ipars[idx["shape",1], 1], m.sim), 
				n = rep(n, m.sim), seed = sseed)
		z = .custzdist(custom.dist, zmatrix, m.sim, n)
	}
	z = rbind(matrix(0, nrow = m, ncol = m.sim), z)
	
	if( is.matrix( custom.dist$distfit ) && custom.dist$type != "z") cres = TRUE else cres = FALSE
	
		
	# create the presample information
	if(!is.na(prereturns[1])){
		prereturns = as.vector(prereturns)
		if(length(prereturns)<m) stop(paste("\narfimapath-->error: prereturns must be of length ", m, sep=""))
	}
	
	if(!is.na(preresiduals[1])){
		preresiduals = as.vector(preresiduals)
		if(length(preresiduals)<m) stop(paste("\narfimapath-->error: preresiduals must be of length ", m, sep=""))
		preres = matrix(preresiduals, nrow = m)
	}
	if(is.na(prereturns[1])){
		prereturns = as.numeric(rep(uncmean(spec), times = m))
	}
	
	
	# input vectors/matrices
	x = c(prereturns, rep(0, n))
	constm = matrix(ipars[idx["mu",1]:idx["mu",2],1], ncol = m.sim, nrow = n + m)

	# outpus matrices
	seriesSim = matrix(0, ncol = m.sim, nrow = n.sim)
	residSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	
	for(i in 1:m.sim){
		if(any(is.na(preresiduals))){
			preres = if(cres) as.numeric(z[1:m,i]) else as.numeric(z[1:m,i])*Sigma
		} 
		res = c(preres, if(cres) as.numeric(z[(m+1):(n+m),i]) else as.numeric(z[(m+1):(n+m),i]) * Sigma)
		residSim[,i] = tail(res, n.sim)
		
		if(modelinc[6]>0){
			mxreg = matrix( ipars[idx["mxreg",1]:idx["mxreg",2], 1], ncol = modelinc[6] )
			constm[,i] = constm[,i] + mxreg %*%t( matrix( mexsim[[i]], ncol = modelinc[6] ) )
		}
		if(modelinc[4]>0){
			x = c(prereturns, rep(0, n + modelinc[3]))
			res = c(if(cres) as.numeric(z[(m+1):(n+m),i]) else as.numeric(z[(m+1):(n+m),i]) * Sigma, if(modelinc[3]>0) rep(0, modelinc[3]) else NULL)
			ans2 = .arfimaxsim(modelinc[1:21], ipars, idx, constm[1:n, i], res[1:(n+modelinc[3])], T = n)
			seriesSim[,i] = head(ans2$series, n.sim)
		} else{
			ans2 = .armaxsim(modelinc[1:21], ipars, idx, constm[,i],  x, res, T = n + m, m)
			seriesSim[,i] = tail(ans2$x, n.sim)
		}
	}
	path = list(seriesSim = seriesSim, residSim = residSim)
	sol = new("ARFIMApath",
			path = path,
			model = model,
			seed = as.integer(sseed))
	return(sol)
}

.arfimaroll = function(spec, data, n.ahead = 1, forecast.length = 500, 
		n.start = NULL, refit.every = 25, refit.window = c("recursive", "moving"), 
		window.size = NULL, solver = "hybrid", fit.control = list(), solver.control = list(),
		calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), cluster = NULL,
		keep.coef = TRUE, ...)
{

	tic = Sys.time()
	if(is.null(solver.control$trace)) trace = 0 else trace = solver.control$trace
	if(is.null(fit.control$stationarity)) fit.control$stationarity = TRUE
	if(is.null(fit.control$fixed.se)) fit.control$fixed.se = FALSE
	if(is.null(fit.control$scale)) fit.control$scale = FALSE
	if(is.null(fit.control$rec.init)) fit.control$rec.init = 'all'
	mm = match(names(fit.control), c("stationarity", "fixed.se", "scale", "rec.init"))
	if(any(is.na(mm))){
		idx = which(is.na(mm))
		enx = NULL
		for(i in 1:length(idx)) enx = c(enx, names(fit.control)[idx[i]])
		warning(paste(c("unidentified option(s) in fit.control:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
	}
	datanames = names(data)
	xdata = .extractdata(data)
	data = xdata$data
	index = xdata$index
	period = xdata$period
	T = NROW(data)
	modelinc = spec@model$modelinc
	if( modelinc[6]>0 ){
		mexdata = spec@model$modeldata$mexdata
		mex = TRUE
	} else{
		mex = FALSE
		mexdata = NULL
	}
	if(n.ahead>1) stop("\narfimaroll:--> n.ahead>1 not supported...try again.")
	if(is.null(n.start)){
		if(is.null(forecast.length)) stop("\narfimaroll:--> forecast.length amd n.start are both NULL....try again.")
		n.start = T - forecast.length
	}
	if(T<=n.start) stop("\nugarchroll:--> start cannot be greater than length of data")
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
			if(window.size<100) stop("\narfimaroll:--> window size must be greater than 100.")
			rollind = lapply(1:m, FUN = function(i) max(1, (s[i]-window.size-out.sample[i])):s[i])
		} else{
			rollind = lapply(1:m, FUN = function(i) (1+(i-1)*refit.every):s[i])
		}
	}
	# distribution
	distribution = spec@model$modeldesc$distribution
	if(any(distribution==c("snorm","sstd","sged","jsu","nig","ghyp","ghst"))) skewed = TRUE else skewed = FALSE
	if(any(distribution==c("std","sstd","ged","sged","jsu","nig","ghyp","ghst"))) shaped = TRUE else shaped = FALSE
	if(any(distribution==c("ghyp"))) ghyp = TRUE else ghyp = FALSE
	
	if(!is.null(cluster)){
		clusterEvalQ(cl = cluster, library(rugarch))
		clusterExport(cluster, c("data", "index", "s","refit.every", 
						"keep.coef", "shaped", "skewed", "ghyp", 
						"rollind", "spec", "out.sample", "mex",
						"solver", "solver.control", "fit.control"), envir = environment())
		if(mex) clusterExport(cluster, c("mex"), envir = environment())
		tmp = parLapply(cluster, as.list(1:m), fun = function(i){
					if(mex) spec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
					fit = try(arfimafit(spec, data[rollind[[i]]], out.sample = out.sample[i], 
									solver = solver, solver.control = solver.control, 
									fit.control = fit.control), silent=TRUE)
					if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
						ans = list(y = NA, cf = NA, converge = FALSE)
					} else{
						if(mex) fmex = tail(mexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fmex = NULL
						f = arfimaforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1, 
								external.forecasts = list(mregfor = fmex))
						sig = as.numeric( rep(coef(fit)["sigma"], out.sample[i]) )
						ret = as.numeric( fitted(f) )
						if(shaped) shp = rep(coef(fit)["shape"], out.sample[i]) else shp = rep(0, out.sample[i])
						if(skewed) skw = rep(coef(fit)["skew"], out.sample[i]) else skw = rep(0, out.sample[i])
						if(ghyp) shpgig = rep(coef(fit)["ghlambda"], out.sample[i]) else shpgig = rep(0, out.sample[i])
						rlz = tail(data[rollind[[i]]], out.sample[i])
						# use xts for indexing the forecasts
						y = as.data.frame(cbind(ret, sig, skw, shp, shpgig, rlz))
						rownames(y) = tail(as.character(index[rollind[[i]]]), out.sample[i])						
						colnames(y) = c("Mu", "Sigma", "Skew", "Shape", "Shape(GIG)", "Realized")
						if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
						ans = list(y = y, cf = cf, converge = TRUE)
					}
					return(ans)})
	} else{
		tmp = lapply(as.list(1:m), FUN = function(i){
					if(mex) spec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
					fit = try(arfimafit(spec, data[rollind[[i]]], out.sample = out.sample[i], 
									solver = solver, solver.control = solver.control, 
									fit.control = fit.control), silent=TRUE)
					if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
						ans = list(y = NA, cf = NA, converge = FALSE)
					} else{
						if(mex) fmex = tail(mexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fmex = NULL
						f = arfimaforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1, 
								external.forecasts = list(mregfor = fmex))
						sig = as.numeric( rep(coef(fit)["sigma"], out.sample[i]) )
						ret = as.numeric( fitted(f) )
						if(shaped) shp = rep(coef(fit)["shape"], out.sample[i]) else shp = rep(0, out.sample[i])
						if(skewed) skw = rep(coef(fit)["skew"], out.sample[i]) else skw = rep(0, out.sample[i])
						if(ghyp) shpgig = rep(coef(fit)["ghlambda"], out.sample[i]) else shpgig = rep(0, out.sample[i])
						rlz = tail(data[rollind[[i]]], out.sample[i])
						# use xts for indexing the forecasts
						y = as.data.frame(cbind(ret, sig, skw, shp, shpgig, rlz))
						rownames(y) = tail(as.character(index[rollind[[i]]]), out.sample[i])
						colnames(y) = c("Mu", "Sigma", "Skew", "Shape", "Shape(GIG)", "Realized")
						if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
						ans = list(y = y, cf = cf, converge = TRUE)
					}
					return(ans)})
	}
	conv = sapply(tmp, FUN = function(x) x$converge)
	if(any(!conv)){
		warning("\nnon-converged estimation windows present...resubsmit object with different solver parameters...")
		noncidx = which(!conv)
		model = list()
		model$spec = spec
		model$data = data
		model$index = index
		model$period = period
		model$datanames = datanames
		model$n.ahead = n.ahead
		model$forecast.length = forecast.length 
		model$n.start = n.start
		model$n.refits = m
		model$refit.every = refit.every
		model$refit.window = refit.window
		model$window.size = window.size
		model$calculate.VaR = calculate.VaR
		model$VaR.alpha = VaR.alpha
		model$keep.coef = keep.coef
		model$noncidx = noncidx
		forecast = tmp
		toc = Sys.time()-tic
		model$elapsed = toc
		ans = new("ARFIMAroll",
				model = model,
				forecast = forecast)
		return(ans)
	} else{
		noncidx = NULL
		forc = tmp[[1]]$y
		for(i in 2:m){
			forc = rbind(forc, tmp[[i]]$y)
		}
		cf = vector(mode = "list", length = m)
		for(i in 1:m){
			cf[[i]]$index = index[tail(rollind[[i]],1) - out.sample[i]]
			cf[[i]]$coef = tmp[[i]]$cf
		}
		if(calculate.VaR){
			if(is.null(VaR.alpha)) VaR.alpha = c(0.01, 0.05)
			VaR.matrix = matrix(NA, ncol = length(VaR.alpha)+1, nrow = NROW(forc))
			for(i in 1:length(VaR.alpha)){
				VaR.matrix[,i] = qdist(VaR.alpha[i] , mu = forc[,1], sigma = forc[,2], 
						skew = forc[,3], shape = forc[,4], lambda = forc[,5], 
						distribution = spec@model$modeldesc$distribution)
			}
			VaR.matrix[,length(VaR.alpha)+1] = forc[,6]
			colnames(VaR.matrix) = c(paste("alpha(", round(VaR.alpha,2)*100, "%)",sep=""), "realized")
			VaR.matrix = as.data.frame(VaR.matrix)
			rownames(VaR.matrix) = rownames(forc)
		} else{
			VaR.matrix = NULL
		}
		model = list()
		model$spec = spec
		model$data = data
		model$index = index
		model$period = period
		model$n.ahead = n.ahead
		model$forecast.length = forecast.length 
		model$n.start = n.start
		model$refit.every = refit.every
		model$n.refits = m
		model$refit.window = refit.window
		model$window.size = window.size
		model$calculate.VaR = calculate.VaR
		model$VaR.alpha = VaR.alpha
		model$keep.coef = keep.coef
		model$noncidx = noncidx
		model$coef = cf
		forecast = list(VaR = VaR.matrix, density = forc)
	}
	toc = Sys.time()-tic
	model$elapsed = toc
	ans = new("ARFIMAroll",
			model = model,
			forecast = forecast)
	return( ans )
}

.resumeroll2 = function(object, solver = "hybrid", fit.control = list(), 
		solver.control = list(), cluster = NULL)
{
	if(!is.null(object@model$noncidx)){
		noncidx = object@model$noncidx
		tic = Sys.time()
		if(is.null(solver.control$trace)) trace = 0 else trace = solver.control$trace
		if(is.null(fit.control$stationarity)) fit.control$stationarity = TRUE
		if(is.null(fit.control$fixed.se)) fit.control$fixed.se = FALSE
		if(is.null(fit.control$scale)) fit.control$scale = FALSE
		if(is.null(fit.control$rec.init)) fit.control$rec.init = 'all'
		mm = match(names(fit.control), c("stationarity", "fixed.se", "scale", "rec.init"))
		if(any(is.na(mm))){
			idx = which(is.na(mm))
			enx = NULL
			for(i in 1:length(idx)) enx = c(enx, names(fit.control)[idx[i]])
			warning(paste(c("unidentified option(s) in fit.control:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
		}
		model = object@model
		keep.coef = model$keep.coef
		spec = model$spec
		datanames = model$datanames
		data = model$data
		index = model$index
		period= model$period
		T = NROW(data)
		modelinc = spec@model$modelinc
		calculate.VaR = model$calculate.VaR
		VaR.alpha = model$VaR.alpha
		if( modelinc[6]>0 ){
			mexdata = spec@model$modeldata$mexdata
			mex = TRUE
		} else{
			mex = FALSE
			mexdata = NULL
		}
		n.ahead = model$n.ahead
		n.start = model$n.start
		forecast.length = model$forecast.length
		refit.every = model$refit.every
		refit.window = model$refit.window
		window.size = model$window.size
		if(n.ahead>1) stop("\narfimaroll:--> n.ahead>1 not supported...try again.")
		if(is.null(n.start)){
			if(is.null(forecast.length)) stop("\narfimaroll:--> forecast.length amd n.start are both NULL....try again.")
			n.start = T - forecast.length
		}
		if(T<=n.start) stop("\narfimaroll:--> start cannot be greater than length of data")
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
				if(window.size<100) stop("\nugarchroll:--> window size must be greater than 100.")
				rollind = lapply(1:m, FUN = function(i) max(1, (s[i]-window.size-out.sample[i])):s[i])
			} else{
				rollind = lapply(1:m, FUN = function(i) (1+(i-1)*refit.every):s[i])
			}
		}
		# distribution
		distribution = spec@model$modeldesc$distribution
		if(any(distribution==c("snorm","sstd","sged","jsu","nig","ghyp","ghst"))) skewed = TRUE else skewed = FALSE
		if(any(distribution==c("std","sstd","ged","sged","jsu","nig","ghyp","ghst"))) shaped = TRUE else shaped = FALSE
		if(any(distribution==c("ghyp"))) ghyp = TRUE else ghyp = FALSE
		if(!is.null(cluster)){
			clusterEvalQ(cl = cluster, library(rugarch))
			clusterExport(cluster, c("data", "index","s","refit.every",
							"keep.coef", "shaped", "skewed", "ghyp", 
							"rollind", "spec", "out.sample", "mex", 
							"noncidx", "solver", "solver.control", "fit.control"),
					envir = environment())
			if(mex) clusterExport(cluster, c("mex"), envir = environment())
			tmp = parLapply(cluster, as.list(noncidx), fun = function(i){
						if(mex) spec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
						fit = try(arfimafit(spec, data[rollind[[i]]], out.sample = out.sample[i], 
										solver=solver, solver.control = solver.control, 
										fit.control = fit.control), silent=TRUE)
						if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
							ans = list(y = NA, cf = NA, converge = FALSE)
						} else{
							if(mex) fmex = tail(mexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fmex = NULL
							f = arfimaforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1, 
									external.forecasts = list(mregfor = fmex))
							sig = as.numeric( rep(coef(fit)["sigma"], out.sample[i]) )
							ret = as.numeric( fitted(f) )
							if(shaped) shp = rep(coef(fit)["shape"], out.sample[i]) else shp = rep(0, out.sample[i])
							if(skewed) skw = rep(coef(fit)["skew"], out.sample[i]) else skw = rep(0, out.sample[i])
							if(ghyp) shpgig = rep(coef(fit)["ghlambda"], out.sample[i]) else shpgig = rep(0, out.sample[i])
							rlz = tail(data[rollind[[i]]], out.sample[i])
							# use xts for indexing the forecasts
							y = as.data.frame(cbind(ret, sig, skw, shp, shpgig, rlz))
							rownames(y) = tail(as.character(index[rollind[[i]]]), out.sample[i])
							colnames(y) = c("Mu", "Sigma", "Skew", "Shape", "Shape(GIG)", "Realized")
							if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
							ans = list(y = y, cf = cf, converge = TRUE)
						}
						return(ans)})
		} else{
			tmp = lapply(as.list(noncidx), FUN = function(i){
						if(mex) spec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
						fit = try(arfimafit(spec, data[rollind[[i]]], out.sample = out.sample[i], 
										solver=solver, solver.control = solver.control, 
										fit.control = fit.control), silent=TRUE)
						if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
							ans = list(y = NA, cf = NA, converge = FALSE)
						} else{
							if(mex) fmex = tail(mexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fmex = NULL
							f = arfimaforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1, 
									external.forecasts = list(mregfor = fmex))
							sig = as.numeric( rep(coef(fit)["sigma"], out.sample[i]) )
							ret = as.numeric( fitted(f) )
							if(shaped) shp = rep(coef(fit)["shape"], out.sample[i]) else shp = rep(0, out.sample[i])
							if(skewed) skw = rep(coef(fit)["skew"], out.sample[i]) else skw = rep(0, out.sample[i])
							if(ghyp) shpgig = rep(coef(fit)["ghlambda"], out.sample[i]) else shpgig = rep(0, out.sample[i])
							rlz = tail(data[rollind[[i]]], out.sample[i])
							# use xts for indexing the forecasts
							y = as.data.frame(cbind(ret, sig, skw, shp, shpgig, rlz))
							rownames(y) = tail(as.character(index[rollind[[i]]]), out.sample[i])
							colnames(y) = c("Mu", "Sigma", "Skew", "Shape", "Shape(GIG)", "Realized")
							if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
							ans = list(y = y, cf = cf, converge = TRUE)
						}
						return(ans)})
		}
		forecast = object@forecast
		conv = sapply(tmp, FUN = function(x) x$converge)
		for(i in 1:length(noncidx)){
			if(conv[i]) forecast[[noncidx[i]]] = tmp[[i]]
		}
		if(any(!conv)){
			warning("\nnon-converged estimation windows present...resubsmit object with different solver parameters...")
			noncidx = noncidx[which(!conv)]
			model = list()
			model$spec = spec
			model$data = data
			model$index = index
			model$period = period
			model$n.ahead = n.ahead
			model$forecast.length = forecast.length 
			model$n.start = n.start
			model$refit.every = refit.every
			model$n.refits = m
			model$refit.window = refit.window
			model$window.size = window.size
			model$calculate.VaR = calculate.VaR
			model$VaR.alpha = VaR.alpha
			model$keep.coef = keep.coef
			model$noncidx = noncidx
			forecast = forecast
			toc = Sys.time()-tic
			model$elapsed = toc
			ans = new("ARFIMAroll",
					model = model,
					forecast = forecast)
			return( ans )			
		} else{
			noncidx = NULL
			forc = forecast[[1]]$y
			for(i in 2:m){
				forc = rbind(forc, forecast[[i]]$y)
			}
			cf = vector(mode = "list", length = m)
			for(i in 1:m){
				cf[[i]]$index = index[tail(rollind[[i]],1) - out.sample[i]]
				cf[[i]]$coef = forecast[[i]]$cf
			}
			if(calculate.VaR){
				if(is.null(VaR.alpha)) VaR.alpha = c(0.01, 0.05)
				VaR.matrix = matrix(NA, ncol = length(VaR.alpha)+1, nrow = NROW(forc))
				for(i in 1:length(VaR.alpha)){
					VaR.matrix[,i] = qdist(VaR.alpha[i], mu = forc[,1], sigma = forc[,2], 
							skew = forc[,3], shape = forc[,4], lambda = forc[,5], 
							distribution = spec@model$modeldesc$distribution)
				}
				VaR.matrix[,length(VaR.alpha)+1] = forc[,6]
				colnames(VaR.matrix) = c(paste("alpha(", round(VaR.alpha,2)*100, "%)",sep=""), "realized")
				VaR.matrix = as.data.frame(VaR.matrix)
				rownames(VaR.matrix) = rownames(forc)
			} else{
				VaR.matrix = NULL
			}
			model = list()
			model$spec = spec
			model$data = data
			model$index = index
			model$period = period
			model$n.ahead = n.ahead
			model$forecast.length = forecast.length 
			model$n.start = n.start
			model$refit.every = refit.every
			model$n.refits = m
			model$refit.window = refit.window
			model$window.size = window.size
			model$calculate.VaR = calculate.VaR
			model$VaR.alpha = VaR.alpha
			model$keep.coef = keep.coef
			model$noncidx = noncidx
			model$coef = cf
			forecast = list(VaR = VaR.matrix, density = forc)
		}
		toc = Sys.time()-tic
		model$elapsed = toc
		ans = new("ARFIMAroll",
				model = model,
				forecast = forecast)
	} else{
		# do nothing...all converged
		ans = object
	}
	return( ans )
}

.arfimadistribution = function(fitORspec, n.sim = 2000, n.start = 1, m.sim = 100, 
		recursive = FALSE, recursive.length = 6000, recursive.window = 1000, 
		prereturns = NA, preresiduals = NA, rseed = NA, 
		custom.dist = list(name = NA, distfit = NA, type = "z"), 
		mexsimdata = NULL, fit.control = list(), solver = "solnp", 
		solver.control = list(), cluster = NULL, ...)
{
	if(recursive){
		nwindows = 1 + round( (recursive.length - n.sim) / recursive.window )
		swindow = vector(mode = "list", length = nwindows)
		rwindow = vector(mode = "list", length = nwindows)
	} else{
		nwindows = 1
		swindow = vector(mode = "list", length = nwindows)
		rwindow = vector(mode = "list", length = nwindows)
		recursive.window = 0
	}
	if(is(fitORspec, "ARFIMAfit")){
		for(i in 1:nwindows){
			sim = arfimasim(fitORspec, n.sim = n.sim + (i-1)*recursive.window, 
					n.start = n.start, m.sim = m.sim, prereturns = prereturns, 
					preresiduals = preresiduals, rseed = rseed, custom.dist = custom.dist, 
					mexsimdata = mexsimdata)
			swindow[[i]]$path.df = as.data.frame(sim)
			swindow[[i]]$seed = sim@seed
		}
		xmodel = fitORspec@model
		fixpars = as.list(coef(fitORspec))
		truecoef = fitORspec@fit$robust.matcoef
		spec = getspec(fitORspec)
	}
	# simulate series paths
	if(is(fitORspec, "ARFIMAspec")){
		
		for(i in 1:nwindows){
			sim  = arfimapath(fitORspec, n.sim =  n.sim + (i-1)*recursive.window, 
					n.start = n.start, m.sim = m.sim, prereturns = prereturns, 
					preresiduals = preresiduals, rseed = rseed, custom.dist = custom.dist, 
					mexsimdata = mexsimdata)
			swindow[[i]]$path.df = fitted(sim)
			swindow[[i]]$seed = sim@seed
		}
		spec = fitORspec
		xmodel = spec@model
		setfixed(spec) <- list(NA)
		fixpars = fitORspec@model$fixed.pars
		truecoef = as.matrix(cbind(unlist(fitORspec@model$fixed.pars),rep(0,length(fixpars)),
						rep(10, length(fixpars)),rep(0,length(fixpars))))
	}
	fitlist = vector( mode = "list", length = m.sim )
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("spec", "swindow", "solver", "fit.control", "solver.control"), 
				envir = environment())
		for(i in 1:nwindows){
			nx = NCOL(swindow[[i]]$path.df)
			clusterExport(cluster, "i", envir = environment())
			rwindow[[i]]$fitlist = parLapply(cluster, as.list(1:nx), fun = function(j){
						.fitandextractarfima(spec, swindow[[i]]$path.df[,j], 
								out.sample = 0, solver = solver, fit.control = fit.control, 
								solver.control = solver.control)
					})
		}
	} else{
		for(i in 1:nwindows){
			rwindow[[i]]$fitlist = apply(swindow[[i]]$path.df, 2, FUN = function(x){
						.fitandextractarfima(spec, x, out.sample = 0, solver = solver, 
								fit.control = fit.control, solver.control = solver.control)
					})
		}
	}
	
	reslist = vector(mode = "list", length = nwindows)
	for(j in 1:nwindows){
		reslist[[j]]$simcoef = 	matrix(NA, ncol = length(fixpars), nrow = m.sim)
		reslist[[j]]$rmse = 	rep(NA, length = length(fixpars))
		reslist[[j]]$simcoefse = matrix(NA, ncol = length(fixpars), nrow = m.sim)
		reslist[[j]]$likelist = rep(NA, length = m.sim)
		reslist[[j]]$mlongrun = rep(NA, length = m.sim)
		reslist[[j]]$simmaxdata  = matrix(NA, ncol = 2, nrow = m.sim)
		reslist[[j]]$simmindata  = matrix(NA, ncol = 2, nrow = m.sim)
		reslist[[j]]$simmeandata  = matrix(NA, ncol = 2, nrow = m.sim)
		reslist[[j]]$simmomdata = matrix(NA, ncol = 2, nrow = m.sim)
		reslist[[j]]$convergence = 	rep(1, length = m.sim)
		reslist[[j]]$seeds = rep(1, length = m.sim)
		for(i in 1:m.sim){
			if(rwindow[[j]]$fitlist[[i]]$convergence!=0) next()
			reslist[[j]]$simcoef[i, ] = rwindow[[j]]$fitlist[[i]]$simcoef
			reslist[[j]]$simcoefse[i,] = rwindow[[j]]$fitlist[[i]]$simcoefse
			reslist[[j]]$likelist[i] = 	rwindow[[j]]$fitlist[[i]]$llh
			reslist[[j]]$mlongrun[i] = 	rwindow[[j]]$fitlist[[i]]$mlongrun
			reslist[[j]]$simmaxdata[i, ] = rwindow[[j]]$fitlist[[i]]$maxdata
			reslist[[j]]$simmindata[i, ] = rwindow[[j]]$fitlist[[i]]$mindata
			reslist[[j]]$simmeandata[i, ] = rwindow[[j]]$fitlist[[i]]$meandata
			reslist[[j]]$simmomdata[i, ] = rwindow[[j]]$fitlist[[i]]$momdata
			reslist[[j]]$convergence[i] = rwindow[[j]]$fitlist[[i]]$convergence
		}
		reslist[[j]]$seed = swindow[[j]]$seed
		reslist[[j]]$rmse = .rmse(reslist[[j]]$simcoef, unlist(fixpars))
	}
	reslist$details = list(n.sim = n.sim, n.start = n.start, m.sim = m.sim,  
			recursive = recursive, recursive.length = recursive.length, recursive.window = recursive.window,
			nwindows = nwindows)
	ans = new("ARFIMAdistribution",
			dist = reslist,
			truecoef = truecoef,
			model = xmodel)
	
	return(ans)
}

.fitandextractarfima = function(spec, x, out.sample = 0,  solver = "solnp", fit.control = list(), solver.control = list())
{
	dist = list()
	fit = .safefitarfima(spec, x, out.sample = 0, solver = solver, fit.control = fit.control, solver.control = solver.control)
	if( is.null(fit) || fit@fit$convergence == 1 || !is( fit, "ARFIMAfit" ) || any( is.na( coef( fit ) ) ) ){
		dist$convergence = 1
		return(dist)
	}
	dist$simcoef = coef(fit)
	dist$simcoefse = fit@fit$robust.matcoef[, 2]
	dist$llh = likelihood(fit)
	dist$mlongrun = uncmean(fit)
	tmp = cbind(fit@model$modeldata$data[1:fit@model$modeldata$T], as.numeric(fitted(fit)))
	dist$maxdata = apply(tmp, 2, "max")
	dist$mindata = apply(tmp, 2, "min")
	dist$meandata = apply(tmp, 2, "mean")
	dist$momdata = c(.kurtosis(tmp[,1]), .skewness(tmp[,1]))
	dist$convergence = fit@fit$convergence
	return(dist)
}

# autoarfima tests for best penalized insample fit based on information criteria
autoarfima = function(data, ar.max = 2, ma.max = 2, criterion = c("AIC", "BIC", "SIC", "HQIC"),
		method = c("partial", "full"), arfima = FALSE, include.mean = NULL, 
		distribution.model = "norm",
		cluster = NULL, external.regressors = NULL, solver = "solnp", 
		solver.control=list(), fit.control=list(), return.all = FALSE){
	m = tolower(method)
	ans = switch(m, 
			partial = .autoarfima1(Data = data, ar.max = ar.max, ma.max = ma.max, 
					criterion = criterion[1], arfima = arfima, include.mean = include.mean, 
					distribution.model = distribution.model, cluster = cluster,
					external.regressors = external.regressors, solver = solver, 
					solver.control=solver.control, fit.control=fit.control, 
					return.all = return.all),
			full = .autoarfima2(Data = data, ar.max = ar.max, ma.max = ma.max, 
					criterion = criterion[1], arfima = arfima, include.mean = include.mean, 
					distribution.model = distribution.model, cluster = cluster,
					external.regressors = external.regressors, solver = solver, 
					solver.control=solver.control, fit.control=fit.control, 
					return.all = return.all))
	return(ans)
}


.autoarfima1 = function(Data, ar.max = 2, ma.max = 2, criterion = c("AIC", "BIC", "SIC", "HQIC"),
		arfima = FALSE, include.mean = NULL, distribution.model = "norm", 
		cluster = NULL, external.regressors = NULL, solver = "solnp", 
		solver.control=list(), fit.control=list(), return.all = FALSE){
	# combinations
	ar = 0:ar.max
	ma = 0:ma.max
	if(is.null(include.mean)) im = c(0,1) else im = as.integer(include.mean)
	if(arfima) arf = c(0, 1) else arf = 0
	d = expand.grid(ar=ar, ma=ma, im = im, arf = arf)
	# eliminate the zero row
	check = apply(d, 1, "sum")
	if(any(check == 0)){
		idx=which(check==0)
		d = d[-idx,]
	}
	n = dim(d)[1]
	IC = match(criterion[1], c("AIC", "BIC", "SIC", "HQIC"))
	fitlist = vector(mode = "list", length = n)
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("d", "Data", "n", "solver", "solver.control", 
						"fit.control"), envir = environment())
		fitlist = parLapply(cluster, as.list(1:n), fun = function(i){
						spec = arfimaspec(
								mean.model = list(armaOrder = c(d[i,1], d[i,2]),
										include.mean =  as.logical(d[i,3]), 
										arfima = as.logical(d[i,4])),
								external.regressors = external.regressors, 
								distribution.model = distribution.model)
						fit = try(arfimafit(spec = spec, data = Data, solver = solver, 
										solver.control = solver.control, fit.control = fit.control), 
								silent = TRUE)
						if( !is(fit, "try-error") && fit@fit$convergence!=0 && solver == "solnp" ){
							fit = try(arfimafit(spec = spec, data = Data, 
											solver = "nlminb", solver.control = list(trace=0), 
											fit.control = fit.control), silent = TRUE)
						}
						return(fit)
					})
	} else{
		fitlist = lapply(as.list(1:n), FUN = function(i){
					spec = arfimaspec(
							mean.model = list(armaOrder = c(d[i,1], d[i,2]),
									include.mean =  as.logical(d[i,3]), 
									arfima = as.logical(d[i,4])),
							external.regressors = external.regressors, 
							distribution.model = distribution.model)
					fit = try(arfimafit(spec = spec, data = Data, solver = solver, 
									solver.control = solver.control, fit.control = fit.control), 
							silent = TRUE)
					if( !is(fit, "try-error") && fit@fit$convergence!=0 && solver == "solnp" ){
						fit = try(arfimafit(spec = spec, data = Data, 
										solver = "nlminb", solver.control = list(trace=0), 
										fit.control = fit.control), silent = TRUE)
					}
					return(fit)
				})
	}
	rankmat = matrix(NA, ncol = 6, nrow = n)
	colnames(rankmat) = c("AR", "MA", "Mean", "ARFIMA", criterion[1], "converged")
	rankmat[,1:4] = as.matrix(d[1:4])
	for(i in 1:n){
		if(inherits(fitlist[[i]], 'try-error') || fitlist[[i]]@fit$convergence!=0){
			rankmat[i,6] = FALSE
		} else{
			rankmat[i,5] = infocriteria(fitlist[[i]])[IC]
			rankmat[i,6] = TRUE
		}
	}
	rk = rankmat[order(rankmat[,5]),]
	rownames(rk) = 1:dim(rk)[1]
	i = which(rankmat[,5] == min(rankmat[,5], na.rm = TRUE))
	if(return.all){
		ans = list(fit = fitlist, rank.matrix = rankmat)
	} else{
		ans = list(fit = fitlist[[i]], rank.matrix = rk)	
	}
	return(ans)
}

.autoarfima2 = function(Data, ar.max = 2, ma.max = 2, criterion = c("AIC", "BIC", "SIC", "HQIC"),
		arfima = FALSE, include.mean = NULL, distribution.model = "norm", 
		cluster = NULL, external.regressors = NULL, solver = "solnp", 
		solver.control=list(), fit.control=list(), return.all = FALSE)
{
	# combinations
	arnames = paste("ar", 1:ar.max, sep = "")
	manames = paste("ma", 1:ma.max, sep = "")
	.str = NULL
	if(ar.max>0){
		for(i in 1:ar.max){
			.str = c(.str, paste(arnames[i],"=c(0,1),",sep=""))
		}
	}
	if(ma.max>0){
		for(i in 1:ma.max){
			.str = c(.str, paste(manames[i],"=c(0,1),",sep=""))
		}
	}
	if(is.null(include.mean)){
		.str = c(.str, "im = c(0,1)")
	} else{
		.str = c(.str, "im = as.integer(include.mean)")
	}
	if(is.null(arfima)){
		.str = c(.str, ",arf = c(0,1)")
	} else{
		.str = c(.str, ",arf = as.integer(arfima)")
	}
	str = c("d = expand.grid(", paste(.str), ')')
	xstr = paste(str, sep="", collapse="")
	eval(parse(text=xstr))
	# eliminate the zero row
	check = apply(d, 1, "sum")
	if(any(check == 0)){
		idx=which(check==0)
		d = d[-idx,]
	}
	sumar = apply(d[,1:ar.max,], 1, "sum")
	summa = apply(d[,(ar.max+1):(ma.max+ar.max),], 1, "sum")
	n = dim(d)[1]
	IC = match(criterion[1], c("AIC", "BIC", "SIC", "HQIC"))
	fitlist = vector(mode = "list", length = n)
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("d", "Data", "n", "solver", 
						"external.regressors",  "distribution.model",
						"ar.max", "ma.max", "solver.control", "fit.control"), 
				envir = environment())
		fitlist = parLapply(cluster, as.list(1:n), fun = function(i){
						if(ar.max>0){
							arr = d[i,1:ar.max]
							if(ma.max>0){
								mar = d[i,(ar.max+1):(ma.max+ar.max)]
							} else{
								mar = 0
							}
						} else{
							arr = 0
							if(ma.max>0){
								mar = d[i,1:ma.max]
							} else{
								mar = 0
							}
						}
						spec = .zarfimaspec( arOrder = arr, maOrder = mar, 
								include.mean = d[i,'im'], arfima = d[i,'arf'], 
								external.regressors = external.regressors, 
								distribution.model = distribution.model)
						fit = try(arfimafit(spec = spec, data = Data, solver = solver, 
										solver.control = solver.control, 
										fit.control = fit.control), silent = TRUE)
						return(fit)
					})
	} else{
		fitlist = lapply(as.list(1:n), FUN = function(i){
					if(ar.max>0){
						arr = d[i,1:ar.max]
						if(ma.max>0){
							mar = d[i,(ar.max+1):(ma.max+ar.max)]
						} else{
							mar = 0
						}
					} else{
						arr = 0
						if(ma.max>0){
							mar = d[i,1:ma.max]
						} else{
							mar = 0
						}
					}
					spec = .zarfimaspec( arOrder = arr, maOrder = mar, 
							include.mean = d[i,'im'], arfima = d[i,'arf'], 
							external.regressors = external.regressors, 
							distribution.model = distribution.model)
					fit = try(arfimafit(spec = spec, data = Data, solver = solver, 
									solver.control = solver.control, 
									fit.control = fit.control), silent = TRUE)
					return(fit)
				})	
	}
	m = dim(d)[2]
	rankmat = matrix(NA, ncol = m+2, nrow = n)
	colnames(rankmat) = c(colnames(d), criterion[1], "converged")
	rankmat[,1:m] = as.matrix(d)
	for(i in 1:n){
		if(inherits(fitlist[[i]], 'try-error') || fitlist[[i]]@fit$convergence!=0){
			rankmat[i,m+2] = 0
		} else{
			rankmat[i,m+1] = infocriteria(fitlist[[i]])[IC]
			rankmat[i,m+2] = 1
		}
	}
	rk = rankmat[order(rankmat[,m+1]),]
	rownames(rk) = 1:dim(rk)[1]
	i = which(rankmat[,m+1] == min(rankmat[,m+1], na.rm = TRUE))
	if(return.all){
		ans = list(fit = fitlist, rank.matrix = rankmat)
	} else{
	 	ans = list(fit = fitlist[[i]], rank.matrix = rk)	
	}
	return(ans)
}