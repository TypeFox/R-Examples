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
# SECTION eGARCH fit
#---------------------------------------------------------------------------------
.egarchfit = function(spec, data, out.sample = 0, solver = "solnp", solver.control = list(), 
		fit.control = list(stationarity = 1, fixed.se = 0, scale = 0, rec.init = 'all'), 
		numderiv.control = list(grad.eps=1e-4, grad.d=0.0001, grad.zero.tol=sqrt(.Machine$double.eps/7e-7),
				hess.eps=1e-4, hess.d=0.1, hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2))
{
	tic = Sys.time()
	if(is.null(solver.control$trace)) trace = 0 else trace = solver.control$trace
	
	if(is.null(fit.control$stationarity)) fit.control$stationarity = TRUE
	if(is.null(fit.control$fixed.se)) fit.control$fixed.se = FALSE
	# No scaling in the eGARCH model (because of the transformation)
	fit.control$scale = FALSE
	if(is.null(fit.control$rec.init)) fit.control$rec.init = 'all'
	mm = match(names(fit.control), c("stationarity", "fixed.se", "scale", "rec.init"))
	if(any(is.na(mm))){
		idx = which(is.na(mm))
		enx = NULL
		for(i in 1:length(idx)) enx = c(enx, names(fit.control)[idx[i]])
		warning(paste(c("unidentified option(s) in fit.control:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
	}
	
	# if there are fixed pars we do no allow scaling as there would be no way of mixing scaled
	# amd non scaled parameters	
	#if(sum(spec@model$pars[,2]) > 0) fit.control$scale = FALSE
	xdata = .extractdata(data)
	if(!is.numeric(out.sample)) stop("\nugarchfit-->error: out.sample must be numeric\n")
	if(as.numeric(out.sample)<0) stop("\nugarchfit-->error: out.sample must be positive\n")
	n.start = round(out.sample,0)
	n = length(xdata$data)
	if((n-n.start)<100) stop("\nugarchfit-->error: function requires at least 100 data\n points to run\n")
	data = xdata$data[1:(n-n.start)]
	index = xdata$index[1:(n-n.start)]
	origdata = xdata$data
	origindex = xdata$index
	period = xdata$period
	# arglist replaces the use of custom environments (resolves the problem of 
	# non-shared environment in parallel estimation, particularly windows)
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
	
	if(modelinc[15] > 0){
		vexdata = model$modeldata$vexdata[1:(n-n.start), ,drop = FALSE]
	} else{
		vexdata = NULL
	}
	arglist$index = index
	arglist$trace = trace
	m =  model$maxOrder
	model$modeldata$T = T = length(as.numeric(data))
	dist = model$modeldesc$distribution
	# no scaling so FALSE for sure
	if(fit.control$scale) dscale = sd(data) else dscale = 1
	zdata = data/dscale
	recinit = .checkrec(fit.control$rec.init, T)
	arglist$data = zdata
	arglist$recinit = recinit
	arglist$dscale = dscale
	arglist$model = model
	ipars = model$pars
	# Optimization Starting Parameters Vector & Bounds
	tmp = .garchstart(ipars, arglist)
	ipars = arglist$ipars = tmp$pars
	arglist$tmph  = tmp$tmph
	# we now split out any fixed parameters
	estidx = as.logical( ipars[,4] )
	arglist$estidx = estidx	
	arglist$fit.control = fit.control
	npars = sum(estidx)
	if(any(ipars[,2]==1)){
		if(npars == 0){
			if(fit.control$fixed.se==0) {
				# if all parameters are fixed an no standard erros are to
				# be calculated then we return a ugarchfilter object
				warning("\nugarchfit-->warning: all parameters fixed...returning ugarchfilter object instead\n")
				return(ugarchfilter(data = xts(origdata, origindex), spec = spec, out.sample = out.sample))
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
	# assisgn solver constraints (solnp directly else exterior type penalty
	# for other solvers)
	if(fit.control$stationarity == 1 && modelinc[15] == 0){
		cb = .garchconbounds()
		Ifn = .egarchcon
		ILB = 1.00000001
		IUB = 1000
		if(solver == "solnp" | solver == "gosolnp" | solver == "hybrid") fit.control$stationarity = 0
	} else{
		Ifn = NULL
		ILB = NULL
		IUB = NULL
	}
	# conditions controls the non-solnp solver penalty
	arglist$fit.control = fit.control
	if(use.solver){
		parscale = rep(1, length = npars)
		names(parscale) = rownames(ipars[estidx,])
		if(modelinc[1] > 0) parscale["mu"] = abs(mean(zdata))
		if(modelinc[7] > 0) parscale["omega"] = var(zdata)
		arglist$returnType = "llh"
		solution = .garchsolver(solver, pars = ipars[estidx, 1], fun = .egarchLLH, 
				Ifn, ILB, IUB, gr = NULL, hessian = NULL, parscale = parscale, 
				control = solver.control, LB = ipars[estidx, 5], 
				UB = ipars[estidx, 6], ux = NULL, ci = NULL, mu = NULL, 
				arglist)
		sol = solution$sol
		hess = solution$hess
		timer = Sys.time()-tic
		if(!is.null(sol$par)){
			ipars[estidx, 1] = sol$par
			if(modelinc[7]==0){
				# call it once more to get omega
				tmpx = .egarchLLH(sol$par, arglist)
				ipars[pidx["omega",1], 1] = get("omega", garchenv)
			}
			if(sum(ipars[,2]) == 0){
				if(modelinc[1] > 0) ipars[pidx["mu",1]:pidx["mu",2], 1] = ipars[pidx["mu",1]:pidx["mu",2], 1] * dscale
				if(modelinc[6] > 0){
					ipars[pidx["mxreg", 1]:pidx["mxreg", 2], 1] = ipars[pidx["mxreg", 1]:pidx["mxreg", 2], 1] * dscale
				}
				ipars[pidx["omega",1], 1] = ipars[pidx["omega",1],1] * dscale^2
			}
		} else{
			ipars[estidx, 1] = NA
		}

		arglist$ipars = ipars
		convergence = sol$convergence
		if(convergence != 0) warning("\nugarchfit-->warning: solver failer to converge.")
	} else{
		solution = NULL
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
		fit = .makefitmodel(garchmodel = "eGARCH", f = .egarchLLH, T = T, m = m, 
				timer = timer, convergence = convergence, message = sol$message, 
				hess, arglist = arglist, numderiv.control = numderiv.control)
		model$modelinc[7] = modelinc[7]
		model$modeldata$data = origdata
		model$modeldata$index = origindex
		model$modeldata$period = period
		model$pars[, 1] = fit$ipars[,1]
		model$pars[, 5:6] = ipars2[,5:6]
		fit$ipars[, 4] = ipars2[, 4]
		fit$ipars[, 2] = ipars2[, 2]
		fit$ipars[, 5:6] = ipars2[,5:6]
		# make sure omega is now included (for working with object post-estimation)
		fit$ipars["omega", 3] = 1
		model$pars["omega", 3] = 1
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
	ans = new("uGARCHfit",
			fit = fit,
			model = model)
	return(ans)
}

#---------------------------------------------------------------------------------
# SECTION eGARCH LLH
#---------------------------------------------------------------------------------
.egarchLLH = function(pars, arglist)
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
	dscale = arglist$dscale	
	recinit = arglist$recinit	
	fit.control = arglist$fit.control
	m = model$maxOrder
	N = c(m,T)
	mexdata = model$modeldata$mexdata[1:T,, drop = FALSE]
	vexdata = model$modeldata$vexdata[1:T,, drop = FALSE]
	distribution = model$modeldesc$distribution
	modelinc = model$modelinc
	# distribution number
	# 1 = norm, 2=snorm, 3=std, 4=sstd, 5=ged, 6=sged, 7=nig
	dist = model$modeldesc$distno
	hm = arglist$tmph
	rx = .arfimaxfilter(modelinc[1:21], ipars[,1], idx, mexdata = mexdata, h = hm, data = data, N = N)
	res = rx$res
	zrf = rx$zrf
	res[is.na(res) | !is.finite(res) | is.nan(res)] = 0
	# sgarch persistence value
	kappa = egarchKappa(ipars[idx["ghlambda",1],1], ipars[idx["shape",1],1], ipars[idx["skew",1],1], distribution)
	persist = sum(ipars[idx["beta",1]:idx["beta",2],1])
	# unconditional sigma value
	mvar = ifelse(recinit$type==1, mean(res[1:recinit$n]*res[1:recinit$n]), backcastv(res, T, recinit$n))
	if(modelinc[15]>0) {
		mv = sum(apply(matrix(vexdata, ncol = modelinc[15]), 2, "mean")*ipars[idx["vxreg",1]:idx["vxreg",2],1])
	} else{
		mv = 0
	}
	hEst = mvar
	if(modelinc[7]==0){
		mvar2 = ifelse(!is.na(modelinc[22]), modelinc[22]/dscale, mvar)
		ipars[idx["omega",1],1] = log(mvar2) * max(1 - persist, 0.001) - mv
		assign("omega", ipars[idx["omega",1],1], garchenv)
	}
	if(is.na(hEst) | !is.finite(hEst) | is.nan(hEst)) hEst = var(data)
	# if we have external regressors in variance equation we cannot have
	# stationarity checks in likelihood. Stationarity conditions not valid for
	# solnp sover which implements them internally as constraints.
	if(fit.control$stationarity == 1 && modelinc[15] == 0){
		if(!is.na(persist) && persist >= 1){
			if(arglist$pmode!=1){
				return(llh = get("rugarch_llh", garchenv) + 0.1*(abs(get("rugarch_llh", garchenv))))
			} else{
				return(llh = 1e10)
			}
		}
	}
	if(modelinc[6]>0) mexdata = as.double(as.vector(mexdata)) else mexdata = double(1)
	if(modelinc[15]>0) vexdata = as.double(as.vector(vexdata)) else vexdata = double(1)
	
	ans = try( .C("egarchfilterC", model = as.integer(modelinc[1:21]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					hEst = as.double(hEst), meanz = as.double(kappa), 
					x = as.double(data), res = as.double(res), e = double(T), 
					mexdata = mexdata, vexdata = vexdata, zrf = as.double(zrf),
					constm = double(T), condm = double(T), m = as.integer(m), 
					T = as.integer(T), h = double(T), z = double(T), llh = double(1), 
					LHT = double(T), PACKAGE = "rugarch"), silent = TRUE )
	
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
	h = ans$h
	epsx = ans$res
	llh = ans$llh
	
	if(is.finite(llh) && !is.na(llh) && !is.nan(llh)){
		if(arglist$pmode!=1) assign("rugarch_llh", llh, envir = garchenv)
	} else {
		if(arglist$pmode!=1) llh = (get("rugarch_llh", garchenv) + 0.1*(abs(get("rugarch_llh", garchenv)))) else llh = 1e10
	}
	
	# LHT = raw scores
	# ? -ans$LHT[(m+1):T]
	LHT = -ans$LHT
	ans = switch(returnType,
			llh = llh,
			LHT = LHT,
			all = list(llh = llh, h = h, epsx = epsx, z = z, kappa = kappa, 
					LHT = LHT, persistence = persist))
	return(ans)
}

#---------------------------------------------------------------------------------
# SECTION eGARCH filter
#---------------------------------------------------------------------------------
.egarchfilter = function(spec, data, out.sample = 0, n.old = NULL, rec.init = 'all')
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
	recinit = .checkrec(rec.init, Nx)
	model = spec@model
	model$modeldata$T = T
	ipars = model$pars
	pars = unlist(model$fixed.pars)
	parnames = names(pars)
	modelnames = .checkallfixed(spec)
	if(is.na(all(match(modelnames, parnames), 1:length(modelnames)))) {
		cat("\nugarchfilter-->error: parameters names do not match specification\n")
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
	vexdata = model$modeldata$vexdata[1:T, , drop = FALSE]
	distribution = model$modeldesc$distribution
	dist = model$modeldesc$distno
	
	kappa = egarchKappa(ipars[idx["ghlambda",1],1], ipars[idx["shape",1],1], ipars[idx["skew",1],1], distribution)
	persist = sum(ipars[idx["beta",1]:idx["beta",2],1])
	
	rx = .arfimaxfilter(modelinc[1:21], ipars[,1], idx, mexdata = mexdata, h = 0, data = data, N = N)
	res = rx$res
	zrf = rx$zrf
	
	if(!is.null(n.old)){
		rx2 = .arfimaxfilter(modelinc[1:21], ipars[,1], idx, mexdata = mexdata[1:Nx, , drop = FALSE], h = 0, data = origdata[1:Nx], N = c(m, Nx))
		res2 = rx2$res
		mvar = ifelse(recinit$type==1, mean(res2[1:recinit$n]*res2[1:recinit$n]), backcastv(res2, Nx, recinit$n))
	} else{
		mvar = ifelse(recinit$type==1, mean(res[1:recinit$n]*res[1:recinit$n]), backcastv(res, T, recinit$n))
	}
	hEst = mvar
	if(modelinc[15]>0) {
		mv = sum(apply(matrix(vexdata, ncol = modelinc[15]), 2, "mean")*ipars[idx["vxreg",1]:idx["vxreg",2],1])
	} else{
		mv = 0
	}
	if(modelinc[7]==0){
		mvar2 = ifelse(!is.na(modelinc[22]), modelinc[22], mvar)
		ipars[idx["omega",1],1] = log(mvar2) * (1 - persist) - mv
	}
	if(modelinc[6]>0) mexdata = as.double(as.vector(mexdata)) else mexdata = double(1)
	if(modelinc[15]>0) vexdata = as.double(as.vector(vexdata)) else vexdata = double(1)
	
	ans = try( .C("egarchfilterC", model = as.integer(modelinc[1:21]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					hEst = as.double(hEst), meanz = as.double(kappa), 
					x = as.double(data), res = as.double(res), e = double(T), 
					mexdata = mexdata, vexdata = vexdata, zrf = as.double(zrf),
					constm = double(T), condm = double(T), m = as.integer(m), 
					T = as.integer(T), h = double(T), z = double(T), llh = double(1), 
					LHT = double(T), PACKAGE = "rugarch"), silent = TRUE )
		
	filter = list()
	filter$z = ans$z
	filter$sigma = sqrt(ans$h)
	filter$residuals = ans$res
	filter$LLH = -ans$llh
	filter$log.likelihoods = ans$LHT
	filter$persistence = persist
	filter$distribution = distribution
	filter$ipars = ipars
	model$modeldata$data = origdata
	model$modeldata$index = origindex
	model$modeldata$period = period
	model$n.start = out.sample
	sol = new("uGARCHfilter",
			filter = filter,
			model = model)
	return(sol)
}

#---------------------------------------------------------------------------------
# SECTION eGARCH forecast
#---------------------------------------------------------------------------------
.egarchforecast = function(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0, 
		external.forecasts = list(mregfor = NULL, vregfor = NULL), ...)
{
	fit    = fitORspec
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
	distribution = model$modeldesc$distribution
	if( n.roll > ns ) stop("\nugarchforecast-->error: n.roll must not be greater than out.sample!")
	pars = fit@fit$coef
	ipars = fit@fit$ipars
	# check if necessary the external regressor forecasts provided first
	xreg = .forcregressors(model, external.forecasts$mregfor, external.forecasts$vregfor, n.ahead, Nor, out.sample = ns, n.roll)
	mxf = xreg$mxf
	vxf = xreg$vxf
	
	fcreq = ifelse(ns >= (n.ahead+n.roll), n.ahead+n.roll, ns)
	fspec = ugarchspec(variance.model = list(model = "eGARCH", 
					garchOrder = c(modelinc[8], modelinc[9]), submodel = NULL, 
					external.regressors = vxf[1:(N + fcreq), , drop = FALSE]), 
			mean.model = list(armaOrder = c(modelinc[2], modelinc[3]),
					include.mean = modelinc[1], 
					archm = ifelse(modelinc[5]>0,TRUE,FALSE), archpow = modelinc[5], arfima = modelinc[4], 
					external.regressors = mxf[1:(N + fcreq), , drop = FALSE], archex = modelinc[20]), 
			distribution.model = model$modeldesc$distribution, fixed.pars = as.list(pars))
	tmp =  xts(data[1:(N + fcreq)], index[1:(N + fcreq)])
	flt = .egarchfilter(data = tmp, spec = fspec, n.old = N)
	sigmafilter = flt@filter$sigma
	resfilter = flt@filter$residuals
	zfilter = flt@filter$z
	kappa = egarchKappa(ipars[idx["ghlambda",1],1], ipars[idx["shape",1],1], ipars[idx["skew",1],1], distribution)
	# forecast GARCH process
	seriesFor = sigmaFor = matrix(NA, ncol = n.roll+1, nrow = n.ahead)
	#seriesFor[1,] = fitted(flt)[(N+1):(N+n.roll+1)]
	#sigmaFor[1,]  = sigma(flt)[(N+1):(N+n.roll+1)]
	# n.ahead x n.roll +1 (n.ahead=1 generted by model)
	colnames(seriesFor) = colnames(sigmaFor) = as.character(index[N:(N+n.roll)])
	rownames(seriesFor) = rownames(sigmaFor) = paste("T+", 1:n.ahead, sep="")
	for(i in 1:(n.roll+1)){
		np = N + i - 1
		if(modelinc[1] > 0){
			mu = rep(ipars[idx["mu",1]:idx["mu",2], 1], N+i+n.ahead-1)
		} else{
			mu = rep(0, N+i+n.ahead-1)
		}
		omega = rep(ipars[idx["omega",1]:idx["omega",2], 1], N+i+n.ahead-1)
		# no look-ahead
		h = c(sigmafilter[1:(N+i-1)], rep(0, n.ahead))
		epsx = c(resfilter[1:(N+i-1)], rep(0, n.ahead))
		x = c(data[1:(N+i-1)], rep(0, n.ahead))
		z = c(zfilter[1:(N+i-1)], rep(0, n.ahead))
		# forecast of externals is provided outside the system
		mxfi = mxf[1:(N+i-1+n.ahead), , drop = FALSE]
		vxfi = vxf[1:(N+i-1+n.ahead), , drop = FALSE]
		ans = .negarchforecast(ipars, modelinc, idx, mu, omega, kappa, mxfi, vxfi, h, epsx, z, data = x, N = np, n.ahead)
		sigmaFor[,i] = ans$h
		seriesFor[,i] = ans$x
	}
	
	fcst = list()
	fcst$n.ahead = n.ahead
	fcst$N = N+ns
	fcst$n.start = ns
	fcst$n.roll = n.roll
	fcst$sigmaFor = sigmaFor
	fcst$seriesFor = seriesFor
	model$modeldata$sigma = flt@filter$sigma
	model$modeldata$residuals = flt@filter$residuals
	ans = new("uGARCHforecast",
			forecast = fcst,
			model = model)
	return(ans)
}	
.negarchforecast = function(ipars, modelinc, idx, mu, omega, kappa, mxfi, vxfi, h, epsx, z, data, N, n.ahead)
{
	if(modelinc[15]>0){
		omega = omega + vxfi%*%t(matrix(ipars[idx["vxreg",1]:idx["vxreg",2],1], ncol = modelinc[15]))
	}
	for(i in 1:n.ahead){
		if(modelinc[9]>0){
			h[N+i] = omega[N+i] + sum(ipars[idx["beta",1]:idx["beta",2],1]*log(h[N+i-(1:modelinc[9])]^2))
		} else{
			h[N+i] = omega[N+i]
		}
		if(modelinc[8]>0){
			for(j in 1:modelinc[8]){
				if (i-j > 0){
					s = 0 #unconditional
				} else{
					# first step ahead z is active else z=0
					s = ipars[idx["alpha",1]+j-1,1] * z[N+i-j] + ipars[idx["gamma",1]+j-1,1]*(abs(z[N+i-j]) - kappa)
				}
				h[N+i] = h[N+i] + s
			}
		}
		h[N+i] = sqrt(exp(h[N+i]))
	}
	# forecast arma process
	if(modelinc[4]>0){
		res = arfimaf(ipars, modelinc[1:21], idx, mu, mxfi, h, epsx, z, data, N, n.ahead)
	} else{
		res = armaf(ipars, modelinc[1:21], idx, mu, mxfi, h, epsx, z, data, N, n.ahead)
	}
	return(list(h = h[(N+1):(N+n.ahead)], x = res[(N+1):(N+n.ahead)]))
}

#---------------------------------------------------------------------------------
# 2nd dispatch method for forecast
.egarchforecast2 = function(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0, 
		external.forecasts = list(mregfor = NULL, vregfor = NULL), ...)
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
	model = spec@model
	ns = out.sample
	if( n.roll > ns ) stop("\nugarchforecast-->error: n.roll must not be greater than out.sample!")
	N = Nor - ns
	distribution = model$modeldesc$distribution
	
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
	xreg = .forcregressors(model, external.forecasts$mregfor, external.forecasts$vregfor, n.ahead, Nor, out.sample = ns, n.roll)
	mxf = xreg$mxf
	vxf = xreg$vxf
	
	# filter data (check external regressor data - must equal length of origData)
	fcreq = ifelse(ns >= (n.ahead+n.roll), n.ahead+n.roll, ns)
	fspec = ugarchspec(variance.model = list(model = "eGARCH", 
					garchOrder = c(modelinc[8], modelinc[9]), submodel = NULL, 
					external.regressors = vxf[1:(N + fcreq), , drop = FALSE]), 
			mean.model = list(armaOrder = c(modelinc[2], modelinc[3]),
					include.mean = modelinc[1], 
					archm = ifelse(modelinc[5]>0,TRUE,FALSE), archpow = modelinc[5], arfima = modelinc[4], 
					external.regressors = mxf[1:(N + fcreq), , drop = FALSE], archex = modelinc[20]), 
			distribution.model = model$modeldesc$distribution, fixed.pars = as.list(pars))
	tmp =  xts(data[1:(N + fcreq)], index[1:(N + fcreq)])
	flt = .egarchfilter(data = tmp, spec = fspec, n.old = N)
	sigmafilter = flt@filter$sigma
	resfilter = flt@filter$residuals
	zfilter = flt@filter$z
	kappa = egarchKappa(ipars[idx["ghlambda",1],1], ipars[idx["shape",1],1], ipars[idx["skew",1],1], distribution)
	# forecast GARCH process
	seriesFor = sigmaFor = matrix(NA, ncol = n.roll+1, nrow = n.ahead)
	#seriesFor[1,] = fitted(flt)[(N+1):(N+n.roll+1)]
	#sigmaFor[1,]  = sigma(flt)[(N+1):(N+n.roll+1)]
	# n.ahead x n.roll +1 (n.ahead=1 generted by model)
	colnames(seriesFor) = colnames(sigmaFor) = as.character(index[N:(N+n.roll)])
	rownames(seriesFor) = rownames(sigmaFor) = paste("T+", 1:n.ahead, sep="")
	for(i in 1:(n.roll+1)){
		np = N + i - 1
		if(modelinc[1] > 0){
			mu = rep(ipars[idx["mu",1]:idx["mu",2], 1], N+i+n.ahead-1)
		} else{
			mu = rep(0, N+i+n.ahead-1)
		}
		omega = rep(ipars[idx["omega",1]:idx["omega",2], 1], N+i+n.ahead-1)
		# no look-ahead
		h = c(sigmafilter[1:(N+i-1)], rep(0, n.ahead))
		epsx = c(resfilter[1:(N+i-1)], rep(0, n.ahead))
		x = c(data[1:(N+i-1)], rep(0, n.ahead))
		z = c(zfilter[1:(N+i-1)], rep(0, n.ahead))
		# forecast of externals is provided outside the system
		mxfi = mxf[1:(N+i-1+n.ahead), , drop = FALSE]
		vxfi = vxf[1:(N+i-1+n.ahead), , drop = FALSE]
		ans = .negarchforecast(ipars, modelinc, idx, mu, omega, kappa, mxfi, vxfi, h, epsx, z, data = x, N = np, n.ahead)
		sigmaFor[,i] = ans$h
		seriesFor[,i] = ans$x
	}
	fcst = list()
	fcst$n.ahead = n.ahead
	fcst$N = N+ns
	fcst$n.start = ns
	fcst$n.roll = n.roll
	fcst$sigmaFor = sigmaFor
	fcst$seriesFor = seriesFor
	model$modeldata$sigma = flt@filter$sigma
	model$modeldata$residuals = flt@filter$residuals
	ans = new("uGARCHforecast",
			forecast = fcst,
			model = model)
	return(ans)
}
#---------------------------------------------------------------------------------
# SECTION eGARCH simulate
#---------------------------------------------------------------------------------
.egarchsim = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, startMethod = 
				c("unconditional","sample"), presigma = NA, prereturns = NA, 
		preresiduals = NA, rseed = NA, custom.dist = list(name = NA, distfit = NA), 
		mexsimdata = NULL, vexsimdata = NULL, ...)
{
	if( (n.sim+n.start) < 1000 && m.sim > 100 ){
		ans = .egarchsim2(fit = fit, n.sim = n.sim, n.start = n.start, m.sim = m.sim, 
				startMethod = startMethod, presigma = presigma, prereturns = prereturns, 
				preresiduals = preresiduals, rseed = rseed, custom.dist = custom.dist, 
				mexsimdata = mexsimdata, vexsimdata = vexsimdata)
	} else{
		ans = .egarchsim1(fit = fit, n.sim = n.sim, n.start = n.start, m.sim = m.sim, 
				startMethod = startMethod, presigma = presigma, prereturns = prereturns, 
				preresiduals = preresiduals, rseed = rseed, custom.dist = custom.dist, 
				mexsimdata = mexsimdata, vexsimdata = vexsimdata)
	}
	return( ans )
}

.egarchsim1 = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, startMethod = 
				c("unconditional","sample"), presigma = NA, prereturns = NA, 
		preresiduals = NA, rseed = NA, custom.dist = list(name = NA, distfit = NA), 
		mexsimdata = NULL, vexsimdata = NULL)
{
	if(fit@model$modelinc[4]>0){
		if(n.start<fit@model$modelinc[3]){
			warning("\nugarchsim-->warning: n.start>=MA order for arfima model...automatically setting.")
			n.start = fit@model$modelinc[3]
		}
	}
	# some checks
	if(is.na(rseed[1])){
		sseed = as.integer(runif(1,0,as.integer(Sys.time())))
	} else{
		if(length(rseed) != m.sim) sseed = as.integer(rseed[1]) else sseed = rseed[1:m.sim]
	}
	# Enlarge Series:
	# need to allow for arfima case:
	n = n.sim + n.start
	startMethod = startMethod[1]
	data = fit@model$modeldata$data
	N = length(as.numeric(data))
	data = data[1:(N - fit@model$n.start)]
	N = length(as.numeric(data))
	m = fit@model$maxOrder
	resids = fit@fit$residuals
	sigma = fit@fit$sigma
	model = fit@model
	modelinc = model$modelinc
	idx = model$pidx
	ipars = fit@fit$ipars
	distribution = model$modeldesc$distribution
	
	# check if necessary the external regressor forecasts provided first
	xreg = .simregressors(model, mexsimdata, vexsimdata, N, n, m.sim, m)	
	mexsim = xreg$mexsimlist
	vexsim = xreg$vexsimlist
	kappa = egarchKappa(ipars[idx["ghlambda",1],1], ipars[idx["shape",1],1], ipars[idx["skew",1],1], distribution)
	if(N < n.start){
		startmethod[1] = "unconditional"
		warning("\nugarchsim-->warning: n.start greater than length of data...using unconditional start method...\n")
	}
	
	# Random Samples from the Distribution
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
	if(startMethod == "unconditional"){
		z = rbind(matrix(0, nrow = m, ncol = m.sim), z)
	} else{
		z = rbind(matrix(tail(fit@fit$z, m), nrow = m, ncol = m.sim), z) 
	}
	# create the presample information
	if(!is.na(presigma[1])){
		presigma = as.vector(presigma)
		if(length(presigma)<m) stop(paste("\nugarchsim-->error: presigma must be of length ", m, sep=""))
	}
	if(!is.na(prereturns[1])){
		prereturns = as.vector(prereturns)
		if(length(prereturns)<m) stop(paste("\nugarchsim-->error: prereturns must be of length ", m, sep=""))
	}
	if(!is.na(preresiduals[1])){
		preresiduals = as.vector(preresiduals)
		if(length(preresiduals)<m) stop(paste("\nugarchsim-->error: preresiduals must be of length ", m, sep=""))
		preres = matrix(preresiduals, nrow = m)
	}
	if(is.na(presigma[1])){
		if(startMethod[1] == "unconditional"){
			hEst = uncvariance(fit)^0.5
			presigma = as.numeric(rep(hEst, m))}
		else{
			presigma  = tail(sigma, m)
		}
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
	h = c(presigma^2, rep(0, n))
	x = c(prereturns, rep(0, n))
	constm = matrix(ipars[idx["mu",1]:idx["mu",2], 1], ncol = m.sim, nrow = n + m)
	
	# outpus matrices
	sigmaSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	seriesSim = matrix(0, ncol = m.sim, nrow = n.sim)
	residSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	z[is.na(z) | is.nan(z) | !is.finite(z)] = 0
	
	for(i in 1:m.sim){
		if(is.na(preresiduals[1])){
			if(startMethod[1] == "unconditional"){
				preres = as.numeric(z[1:m, i])*presigma
			} else{
				preres = tail(resids, m)
			}
		}
		res = c(preres, rep(0, n))
		ans1 = try(.C("egarchsimC", model = as.integer(modelinc[1:21]), 
						pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
						meanz = as.double(kappa), h = as.double(h), z = as.double(z[,i]), 
						res = as.double(res), vexdata = as.double(vexsim[[i]]), 
						T = as.integer(n+m), m = as.integer(m), 
						PACKAGE = "rugarch"), silent = TRUE)
		if(inherits(ans1, "try-error")) stop("\nugarchsim-->error: error in calling C function....\n")
		
		sigmaSim[,i] = ans1$h[(n.start + m + 1):(n+m)]^(1/2)
		residSim[,i] = ans1$res[(n.start + m + 1):(n+m)]
		if(modelinc[6]>0){
			mxreg = matrix( ipars[idx["mxreg",1]:idx["mxreg",2], 1], ncol = modelinc[6] )
			if(modelinc[20]==0){
				constm[,i] = constm[,i] + mxreg %*%t( matrix( mexsim[[i]], ncol = modelinc[6] ) )
			} else{
				if(modelinc[20] == modelinc[6]){
					constm[,i] = constm[,i] + mxreg %*%t( matrix( mexsim[[i]]*sqrt(ans1$h), ncol = modelinc[6] ) )
				} else{
					constm[,i] = constm[,i] + mxreg[,1:(modelinc[6]-modelinc[20]),drop=FALSE] %*%t( matrix( mexsim[[i]][,1:(modelinc[6]-modelinc[20]),drop=FALSE], ncol = modelinc[6]-modelinc[20] ) )
					constm[,i] = constm[,i] + mxreg[,(modelinc[6]-modelinc[20]+1):modelinc[6],drop=FALSE] %*%t( matrix( mexsim[[i]][,(modelinc[6]-modelinc[20]+1):modelinc[6],drop=FALSE]*sqrt(ans1$h), ncol = modelinc[20] ) )
				}
			}
		}
		
		if(modelinc[5]>0) constm[,i] = constm[,i] + ipars[idx["archm",1]:idx["archm",2], 1]*(sqrt(ans1$h)^modelinc[5])
		
		if(modelinc[4]>0){
			fres = c(ans1$res[(m+1):(n+m)], if(modelinc[3]>0) rep(0, modelinc[3]) else NULL)
			ans2 = .arfimaxsim(modelinc[1:21], ipars, idx, constm[1:n, i], fres, T = n)
			seriesSim[,i] = head(ans2$series, n.sim)
		} else{
			ans2 = .armaxsim(modelinc[1:21], ipars, idx, constm[,i],  x, ans1$res, T = n + m, m)
			seriesSim[,i] = ans2$x[(n.start + m + 1):(n+m)]
		}
	}
	sim = list(sigmaSim = sigmaSim, seriesSim = seriesSim, residSim = residSim)
	model$modeldata$sigma = sigma
	sol = new("uGARCHsim",
			simulation = sim,
			model = model,
			seed = as.integer(sseed))
	return(sol)
}

.egarchsim2 = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, startMethod = 
				c("unconditional","sample"), presigma = NA, prereturns = NA, 
		preresiduals = NA, rseed = NA, custom.dist = list(name = NA, distfit = NA), 
		mexsimdata = NULL, vexsimdata = NULL)
{
	if(fit@model$modelinc[4]>0){
		if(n.start<fit@model$modelinc[3]){
			warning("\nugarchsim-->warning: n.start>=MA order for arfima model...automatically setting.")
			n.start = fit@model$modelinc[3]
		}
	}
	if(is.na(rseed[1])){
		sseed = as.integer(runif(1,0,as.integer(Sys.time())))
	} else{
		if(length(rseed) != m.sim) sseed = as.integer(rseed[1]) else sseed = rseed[1:m.sim]
	}
	# Enlarge Series:
	# need to allow for arfima case:
	n = n.sim + n.start
	startMethod = startMethod[1]
	data = fit@model$modeldata$data
	N = fit@model$modeldata$T
	data = data[1:N]
	m = fit@model$maxOrder
	resids = fit@fit$residuals
	sigma = fit@fit$sigma
	model = fit@model
	modelinc = model$modelinc
	idx = model$pidx
	ipars = fit@fit$ipars
	distribution = model$modeldesc$distribution
	
	# check if necessary the external regressor forecasts provided first
	xreg = .simregressors(model, mexsimdata, vexsimdata, N, n, m.sim, m)	
	mexsim = xreg$mexsimlist
	vexsim = xreg$vexsimlist
	kappa = egarchKappa(ipars[idx["ghlambda",1],1], ipars[idx["shape",1],1], ipars[idx["skew",1],1], distribution)
	
	if(N < n.start){
		startmethod[1] = "unconditional"
		warning("\nugarchsim-->warning: n.start greater than length of data...using unconditional start method...\n")
	}
	
	# Random Samples from the Distribution
	if(length(sseed) == 1){
		zmatrix = data.frame(dist = model$modeldesc$distribution, lambda = ipars[idx["ghlambda",1], 1], 
				skew = ipars[idx["skew",1], 1], shape = ipars[idx["shape",1], 1],  n = n * m.sim, seed = sseed[1])
		z = .custzdist(custom.dist, zmatrix, m.sim, n)
	} else{
		zmatrix = data.frame(dist = rep(model$modeldesc$distribution, m.sim), lambda = rep(ipars[idx["ghlambda",1], 1], m.sim), 
				skew = rep(ipars[idx["skew",1], 1], m.sim), shape = rep(ipars[idx["shape",1], 1], m.sim), 
				n = rep(n, m.sim), seed = sseed)
		z = .custzdist(custom.dist, zmatrix, m.sim, n)
	}
	if(startMethod == "unconditional"){
		z = rbind(matrix(0, nrow = m, ncol = m.sim), z)
	} else{
		z = rbind(matrix(tail(fit@fit$z, m), nrow = m, ncol = m.sim), z) 
	}
	
	
	# create the presample information
	if(!is.na(presigma[1])){
		presigma = as.vector(presigma)
		if(length(presigma)<m) stop(paste("\nugarchsim-->error: presigma must be of length ", m, sep=""))
	}
	if(!is.na(prereturns[1])){
		prereturns = as.vector(prereturns)
		if(length(prereturns)<m) stop(paste("\nugarchsim-->error: prereturns must be of length ", m, sep=""))
	}
	if(!is.na(preresiduals[1])){
		preresiduals = as.vector(preresiduals)
		if(length(preresiduals)<m) stop(paste("\nugarchsim-->error: preresiduals must be of length ", m, sep=""))
		preres = matrix(preresiduals[1:m], nrow = m, ncol = m.sim)
	}
	
	if(is.na(presigma[1])){
		if(startMethod[1] == "unconditional"){
			hEst = uncvariance(fit)^0.5
			presigma = as.numeric(rep(hEst, m))
		} else{
			presigma  = tail(sigma, m)
		}
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
	h = matrix(c(presigma * presigma, rep(0, n)), nrow = n + m, ncol = m.sim)
	x = matrix(c(prereturns, rep(0, n)), nrow = n + m, ncol = m.sim)
	constm = matrix(ipars[idx["mu",1]:idx["mu",2],1], nrow = n + m, ncol = m.sim)
	
	# outpus matrices
	sigmaSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	seriesSim = matrix(0, ncol = m.sim, nrow = n.sim)
	residSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	z[is.na(z) | is.nan(z) | !is.finite(z)] = 0
	
	if(is.na(preresiduals[1])){
		if(startMethod[1] == "unconditional"){
			preres = matrix( z[1:m, 1:m.sim] * presigma, nrow = m, ncol = m.sim )
		} else{
			preres = matrix(tail(resids, m), nrow = m, ncol = m.sim)
		}
	}
	res =  rbind(preres, matrix(0, nrow = n, ncol = m.sim))
	if(modelinc[15]>0){
		vxreg = matrix( ipars[idx["vxreg",1]:idx["vxreg",2], 1], ncol = modelinc[15] )
		vxs = sapply(vexsim, FUN = function(x) vxreg%*%t(matrix(x, ncol = modelinc[15])))
	} else{
		vxs = matrix(0, nrow = m + n, ncol = m.sim)
	}
	
	ans = .Call("megarchsim", model = as.integer(modelinc[1:21]), 
			pars = as.numeric(ipars[,1]), idx = as.integer(idx[,1]-1), 
			meanz = as.numeric(kappa), h = h, z = z, res = res, vxs = vxs, 
			N = as.integer( c(m, n) ), PACKAGE = "rugarch")
	
	sigmaSim = matrix(sqrt( ans$h[(n.start + m + 1):(n+m), ] ), ncol = m.sim)
	residSim = matrix(ans$res[(n.start + m + 1):(n+m), ], ncol = m.sim)
	
	if(modelinc[6]>0){
		mxreg = matrix( ipars[idx["mxreg",1]:idx["mxreg",2], 1], ncol = modelinc[6] )
		if(modelinc[20]==0){
			mxs = sapply(mexsim, FUN = function(x) mxreg%*%t(matrix(x, ncol = modelinc[6])))
		} else{
			if(modelinc[20] == modelinc[6]){
				mxs = sapply(mexsim, FUN = function(x) mxreg%*%t(matrix(x*sqrt(ans$h), ncol = modelinc[6])))
			} else{
				mxs = sapply(mexsim, FUN = function(x) mxreg[,1:(modelinc[6]-modelinc[20]),drop=FALSE]%*%t(matrix(x[,1:(modelinc[6]-modelinc[20]),drop=FALSE], ncol = modelinc[6])))
				mxs = mxs + sapply(mexsim, FUN = function(x) mxreg[,(modelinc[6]-modelinc[20]+1):modelinc[6],drop=FALSE]%*%t(matrix(x[,(modelinc[6]-modelinc[20]+1):modelinc[6],drop=FALSE]*sqrt(ans$h), ncol = modelinc[20])))				
			}
		}
	} else{
		mxs = 0
	}
	if(modelinc[5]>0){
		imh = ipars[idx["archm",1],1]*(sqrt( ans$h )^modelinc[5])
	} else{
		imh = 0
	}
	
	constm = constm + mxs + imh
	if(modelinc[4]>0){
		#if(constant) constm[,i] = constm[,i]*(1-sum(ar))
		for(i in 1:m.sim){
			fres = c(ans$res[(m+1):(n+m), i], if(modelinc[3]>0) rep(0, modelinc[3]) else NULL)
			tmp = .arfimaxsim(modelinc[1:21], ipars, idx, constm[1:n, i], fres, T = n)
			seriesSim[,i] = head(tmp$series, n.sim)
		}
	} else{
		#if(constant) constm = constm * ( 1 - sum(ar) )
		tmp = .Call("marmaxsim", model = as.integer(modelinc[1:21]), 
				pars = as.numeric(ipars[,1]), idx = as.integer(idx[,1]-1), 
				mu = constm, x = x, res = ans$res, N = as.integer( c(m, n) ), 
				PACKAGE = "rugarch")
		seriesSim = matrix(tmp$x[(n.start + m + 1):(n+m), ], ncol = m.sim)
	}
	
	sim = list(sigmaSim = sigmaSim, seriesSim = seriesSim, residSim = residSim)
	model$modeldata$sigma = sigma
	sol = new("uGARCHsim",
			simulation = sim,
			model = model,
			seed = as.integer(sseed))
	return(sol)
}
#---------------------------------------------------------------------------------
# SECTION eGARCH path
#---------------------------------------------------------------------------------
.egarchpath = function(spec, n.sim = 1000, n.start = 0, m.sim = 1,
		presigma = NA, prereturns = NA, preresiduals = NA, rseed = NA, 
		custom.dist = list(name = NA, distfit = NA), mexsimdata = NULL, 
		vexsimdata = NULL, ...)
{
	if( (n.sim+n.start) < 1000 && m.sim > 100 ){
		ans = .egarchpath2(spec = spec, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
				presigma = presigma, prereturns = prereturns, preresiduals = preresiduals, 
				rseed = rseed, custom.dist = custom.dist, mexsimdata = mexsimdata, 
				vexsimdata = vexsimdata)
	} else{
		ans = .egarchpath1(spec = spec, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
				presigma = presigma, prereturns = prereturns, preresiduals = preresiduals, 
				rseed = rseed, custom.dist = custom.dist, mexsimdata = mexsimdata, 
				vexsimdata = vexsimdata)
	}
	return( ans )
}

.egarchpath1 = function(spec, n.sim = 1000, n.start = 0, m.sim = 1,
		presigma = NA, prereturns = NA, preresiduals = NA, rseed = NA, 
		custom.dist = list(name = NA, distfit = NA), mexsimdata = NULL, 
		vexsimdata = NULL)
{
	if(spec@model$modelinc[4]>0){
		if(n.start<spec@model$modelinc[3]){
			warning("\nugarchpath-->warning: n.start>=MA order for arfima model...automatically setting.")
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
		cat("\nugarchpath-->error: parameters names do not match specification\n")
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
	# Enlarge Series:
	n = n.sim + n.start
	m = model$maxOrder
	N = 0
	if(modelinc[6]>0) {
		mexdata = matrix(model$modeldata$mexdata, ncol = modelinc[6])
		N = dim(mexdata)[1]
	} else { mexdata = NULL }
	if(modelinc[15]>0) {
		vexdata = matrix(model$modeldata$vexdata, ncol = modelinc[15]) 
		N = dim(vexdata)[1]
	} else { vexdata = NULL }
	distribution = model$modeldesc$distribution	
	# check if necessary the external regressor forecasts provided first
	xreg = .simregressors(model, mexsimdata, vexsimdata, N, n, m.sim, m)	
	mexsim = xreg$mexsimlist
	vexsim = xreg$vexsimlist
	
	kappa = egarchKappa(ipars[idx["ghlambda",1],1], ipars[idx["shape",1],1], ipars[idx["skew",1],1], distribution)
	
	persist = sum(ipars[idx["beta",1]:idx["beta",2], 1])	
	if(persist >= 1) warning(paste("\nugarchpath->warning: persitence :", round(persist, 5), sep=""))
	
	# Random Samples from the Distribution
	if(length(sseed) == 1){
		zmatrix = data.frame(dist = distribution, lambda = ipars[idx["ghlambda",1], 1], 
				skew = ipars[idx["skew",1], 1], shape = ipars[idx["shape",1], 1], 
				n = n * m.sim, seed = sseed[1])
		z = .custzdist(custom.dist, zmatrix, m.sim, n)
	} else{
		zmatrix = data.frame(dist = rep(distribution, m.sim), lambda = rep(ipars[idx["ghlambda",1], 1], m.sim), 
				skew = rep(ipars[idx["skew",1], 1], m.sim), shape = rep(ipars[idx["shape",1], 1], m.sim), 
				n = rep(n, m.sim), seed = sseed)
		z = .custzdist(custom.dist, zmatrix, m.sim, n)
	}
	z = rbind(matrix(0, nrow = m, ncol = m.sim), z)
	
	# create the presample information
	if(!is.na(presigma[1])){
		presigma = as.vector(presigma)
		if(length(presigma)<m) stop(paste("\nugarchsim-->error: presigma must be of length ", m, sep=""))
	}
	
	if(!is.na(prereturns[1])){
		prereturns = as.vector(prereturns)
		if(length(prereturns)<m) stop(paste("\nugarchsim-->error: prereturns must be of length ", m, sep=""))
	}
	
	if(!is.na(preresiduals[1])){
		preresiduals = as.vector(preresiduals)
		if(length(preresiduals)<m) stop(paste("\nugarchsim-->error: preresiduals must be of length ", m, sep=""))
		preres = matrix(preresiduals, nrow = m)
	}
	
	if(is.na(presigma[1])){
		hEst = uncvariance(spec)^0.5
		presigma = as.numeric(rep(hEst, m))
	}
	if(is.na(prereturns[1])){
		prereturns = as.numeric(rep(uncmean(spec), times = m))
	}
	
	# input vectors/matrices
	h = c(presigma^2, rep(0, n))
	x = c(prereturns, rep(0, n))
	constm = matrix(ipars[idx["mu",1]:idx["mu",2],1], ncol = m.sim, nrow = n + m)
	
	# outpus matrices
	sigmaSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	seriesSim = matrix(0, ncol = m.sim, nrow = n.sim)
	residSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	
	for(i in 1:m.sim){
		if(is.na(preresiduals[1])){
			preres = as.numeric(z[1:m,i])*presigma
		}
		z[1:m, 1:m.sim] = preres[1:m]/presigma[1:m]
		z[is.na(z) | is.nan(z) | !is.finite(z)] = 0
		res = c(preres, rep(0, n))
		
		ans1 = try(.C("egarchsimC", model = as.integer(modelinc[1:21]), 
						pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
						meanz = as.double(kappa), h = as.double(h), 
						z = as.double(z[,i]), res = as.double(res),
						vexdata = as.double(vexsim[[i]]), T = as.integer(n+m), 
						m = as.integer(m), PACKAGE = "rugarch"), silent = TRUE)
		if(inherits(ans1, "try-error")) stop("\nugarchsim-->error: error in calling C function....\n")
		sigmaSim[,i] = ans1$h[(n.start + m + 1):(n+m)]^(1/2)
		residSim[,i] = ans1$res[(n.start + m + 1):(n+m)]
		if(modelinc[6]>0){
			mxreg = matrix( ipars[idx["mxreg",1]:idx["mxreg",2], 1], ncol = modelinc[6] )
			if(modelinc[20]==0){
				constm[,i] = constm[,i] + as.numeric( mxreg %*%t( matrix( mexsim[[i]], ncol = modelinc[6] ) ) )
			} else{
				if(modelinc[20] == modelinc[6]){
					constm[,i] = constm[,i] + as.numeric( mxreg %*%t( matrix( mexsim[[i]]*sqrt(ans1$h), ncol = modelinc[6] ) ) )
				} else{
					constm[,i] = constm[,i] + as.numeric( mxreg[,1:(modelinc[6]-modelinc[20]),drop=FALSE] %*%t( matrix( mexsim[[i]][,1:(modelinc[6]-modelinc[20]),drop=FALSE], ncol = modelinc[6]-modelinc[20] ) ) )
					constm[,i] = constm[,i] + as.numeric( mxreg[,(modelinc[6]-modelinc[20]+1):modelinc[6],drop=FALSE] %*%t( matrix( mexsim[[i]][,(modelinc[6]-modelinc[20]+1):modelinc[6],drop=FALSE]*sqrt(ans1$h), ncol = modelinc[20] ) ) )
				}
			}
		}
		if(modelinc[5]>0) constm[,i] = constm[,i] + ipars[idx["archm",1]:idx["archm",2], 1]*(sqrt(ans1$h)^modelinc[5])
		if(modelinc[4]>0){
			fres = c(ans1$res[(m+1):(n+m)], if(modelinc[3]>0) rep(0, modelinc[3]) else NULL)
			ans2 = .arfimaxsim(modelinc[1:21], ipars, idx, constm[1:n, i], fres, T = n)
			seriesSim[,i] = head(ans2$series, n.sim)
		} else{
			#if(constant) constm[,i] = constm[,i]*(1-sum(ar))
			ans2 = .armaxsim(modelinc[1:21], ipars, idx, constm[,i],  x, ans1$res, T = n + m, m)
			seriesSim[,i] = ans2$x[(n.start + m + 1):(n+m)]
		}
	}
	
	path = list(sigmaSim = sigmaSim, seriesSim = seriesSim, residSim = residSim)
	sol = new("uGARCHpath",
			path = path,
			model = model,
			seed = as.integer(sseed))
	return(sol)
}


.egarchpath2 = function(spec, n.sim = 1000, n.start = 0, m.sim = 1,
		presigma = NA, prereturns = NA, preresiduals = NA, rseed = NA, 
		custom.dist = list(name = NA, distfit = NA), mexsimdata = NULL, 
		vexsimdata = NULL)
{
	if(spec@model$modelinc[4]>0){
		if(n.start<spec@model$modelinc[3]){
			warning("\nugarchpath-->warning: n.start>=MA order for arfima model...automatically setting.")
			n.start = spec@model$modelinc[3]
		}
	}
	# some checks
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
		cat("\nugarchpath-->error: parameters names do not match specification\n")
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
	# Enlarge Series:
	n = n.sim + n.start
	m = spec@model$maxOrder
	N = 0
	if(modelinc[6]>0) {
		mexdata = matrix(model$modeldata$mexdata, ncol = modelinc[6])
		N = dim(mexdata)[1]
	} else { mexdata = NULL }
	if(modelinc[15]>0) {
		vexdata = matrix(model$modeldata$vexdata, ncol = modelinc[15]) 
		N = dim(vexdata)[1]
	} else { vexdata = NULL }
	distribution = model$modeldesc$distribution
	
	# check if necessary the external regressor forecasts provided first
	xreg = .simregressors(model, mexsimdata, vexsimdata, N, n, m.sim, m)	
	mexsim = xreg$mexsimlist
	vexsim = xreg$vexsimlist
	
	kappa = egarchKappa(ipars[idx["ghlambda",1],1], ipars[idx["shape",1],1], ipars[idx["skew",1],1], distribution)
	persist = sum(ipars[idx["beta",1]:idx["beta",2], 1])	
	if(persist >= 1) warning(paste("\nugarchpath->warning: persitence :", round(persist, 5), sep=""))
	
	# Random Samples from the Distribution
	if(length(sseed) == 1){
		zmatrix = data.frame(dist = distribution, lambda = ipars[idx["ghlambda",1], 1], 
				skew = ipars[idx["skew",1], 1], shape = ipars[idx["shape",1], 1], 
				n = n * m.sim, seed = sseed[1])
		z = .custzdist(custom.dist, zmatrix, m.sim, n)
	} else{
		zmatrix = data.frame(dist = rep(distribution, m.sim), lambda = rep(ipars[idx["ghlambda",1], 1], m.sim), 
				skew = rep(ipars[idx["skew",1], 1], m.sim), shape = rep(ipars[idx["shape",1], 1], m.sim), 
				n = rep(n, m.sim), seed = sseed)
		z = .custzdist(custom.dist, zmatrix, m.sim, n)
	}
	z = rbind(matrix(0, nrow = m, ncol = m.sim), z)
	
	# create the presample information
	if(!is.na(presigma[1])){
		presigma = as.vector(presigma)
		if(length(presigma)<m) stop(paste("\nugarchsim-->error: presigma must be of length ", m, sep=""))
	}
	
	if(!is.na(prereturns[1])){
		prereturns = as.vector(prereturns)
		if(length(prereturns)<m) stop(paste("\nugarchsim-->error: prereturns must be of length ", m, sep=""))
	}
	
	if(!is.na(preresiduals[1])){
		preresiduals = as.vector(preresiduals)
		if(length(preresiduals)<m) stop(paste("\nugarchsim-->error: preresiduals must be of length ", m, sep=""))
		preres = matrix(preresiduals[1:m], nrow = m, ncol = m.sim)
	}
	if(is.na(presigma[1])){
		hEst = uncvariance(spec)^0.5
		presigma = as.numeric(rep(hEst, m))
	}
	if(is.na(prereturns[1])){
		prereturns = as.numeric(rep(uncmean(spec), times = m))
	}
	# input vectors/matrices
	h = matrix(c(presigma * presigma, rep(0, n)), nrow = n + m, ncol = m.sim)
	x = matrix(c(prereturns, rep(0, n)), nrow = n + m, ncol = m.sim)
	constm = matrix(ipars[idx["mu",1]:idx["mu",2],1], nrow = n + m, ncol = m.sim)
	
	# outpus matrices
	sigmaSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	seriesSim = matrix(0, ncol = m.sim, nrow = n.sim)
	residSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	
	if(is.na(preresiduals[1])){
		preres = matrix( z[1:m, 1:m.sim] * presigma, nrow = m, ncol = m.sim )
	} else{
		preres = matrix( preresiduals, nrow = m, ncol = m.sim )
	}
	z[1:m, 1:m.sim] = preres[1:m, 1:m.sim]/presigma[1:m]
	z[is.na(z) | is.nan(z) | !is.finite(z)] = 0
	
	res =  rbind(preres, matrix(0, nrow = n, ncol = m.sim))
	# we'll do the external regressors first for speed.
	if(modelinc[15]>0){
		vxreg = matrix( ipars[idx["vxreg",1]:idx["vxreg",2], 1], ncol = modelinc[15] )
		vxs = sapply(vexsim, FUN = function(x) vxreg%*%t(matrix(x, ncol = modelinc[15])))
	} else{
		vxs = matrix(0, nrow = m + n, ncol = m.sim)
	}
	ans = .Call("megarchsim", model = as.integer(modelinc[1:21]), 
			pars = as.numeric(ipars[,1]), idx = as.integer(idx[,1]-1), 
			meanz = as.numeric(kappa), h = h, z = z, res = res, vxs = vxs, 
			N = as.integer( c(m, n) ), PACKAGE = "rugarch")
	
	sigmaSim = matrix(sqrt( ans$h[(n.start + m + 1):(n+m), ] ), ncol = m.sim)
	residSim = matrix(ans$res[(n.start + m + 1):(n+m), ], ncol = m.sim)
	
	if(modelinc[6]>0){
		mxreg = matrix( ipars[idx["mxreg",1]:idx["mxreg",2], 1], ncol = modelinc[6] )
		if(modelinc[20]==0){
			mxs = sapply(mexsim, FUN = function(x) mxreg%*%t(matrix(x, ncol = modelinc[6])))
		} else{
			if(modelinc[20] == modelinc[6]){
				mxs = sapply(mexsim, FUN = function(x) mxreg%*%t(matrix(x*sqrt( ans$h ), ncol = modelinc[6])))
			} else{
				mxs = sapply(mexsim, FUN = function(x) mxreg[,1:(modelinc[6]-modelinc[20]),drop=FALSE]%*%t(matrix(x[,1:(modelinc[6]-modelinc[20]),drop=FALSE], ncol = modelinc[6])))
				mxs = mxs + sapply(mexsim, FUN = function(x) mxreg[,(modelinc[6]-modelinc[20]+1):modelinc[6],drop=FALSE]%*%t(matrix(x[,(modelinc[6]-modelinc[20]+1):modelinc[6],drop=FALSE]*sqrt( ans$h ), ncol = modelinc[20])))				
			}
		}
	} else{
		mxs = 0
	}
	if(modelinc[5]>0){
		imh = ipars[idx["archm",1],1]*(sqrt( ans$h )^modelinc[5])
	} else{
		imh = 0
	}
	
	constm = constm + mxs + imh
	if(modelinc[4]>0){
		for(i in 1:m.sim){
			fres = c(ans$res[(m+1):(n+m), i], if(modelinc[3]>0) rep(0, modelinc[3]) else NULL)
			tmp = .arfimaxsim(modelinc[1:21], ipars, idx, constm[1:n, i], fres, T = n)
			seriesSim[,i] = head(tmp$series, n.sim)
		}
	} else{
		tmp = .Call("marmaxsim", model = as.integer(modelinc[1:21]), 
				pars = as.numeric(ipars[,1]), idx = as.integer(idx[,1]-1), 
				mu = constm, x = x, res = ans$res, N = as.integer( c(m, n) ), 
				PACKAGE = "rugarch")
		seriesSim = matrix(tmp$x[(n.start + m + 1):(n+m), ], ncol = m.sim)
	}
	path = list(sigmaSim = sigmaSim, seriesSim = seriesSim, residSim = residSim)
	sol = new("uGARCHpath",
			path = path,
			model = model,
			seed = as.integer(sseed))
	return(sol)
}