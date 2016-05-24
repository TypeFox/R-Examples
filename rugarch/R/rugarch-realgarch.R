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
# SECTION realGARCH fit
#---------------------------------------------------------------------------------
# [mu ar ma arfima im mxreg omega alpha beta vxreg eta11 eta21 delta lambda xi skew shape]

.realgarchfit = function(spec, data, out.sample = 0, solver = "solnp", solver.control = list(), 
		fit.control = list(stationarity = 1, fixed.se = 0, scale = 0, rec.init = 'all'), 
		numderiv.control = list(grad.eps=1e-4, grad.d=0.0001, grad.zero.tol=sqrt(.Machine$double.eps/7e-7),
				hess.eps=1e-4, hess.d=0.1, hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2),
		realizedVol = NULL)
{
	tic = Sys.time()
	if(is.null(solver.control$trace)) trace = 0 else trace = solver.control$trace
	
	if(is.null(fit.control$stationarity)) fit.control$stationarity = TRUE
	if(is.null(fit.control$fixed.se)) fit.control$fixed.se = FALSE
	if(is.null(fit.control$scale)){
		fit.control$scale = FALSE
	} else{
		if(fit.control$scale){
			warning("\nscaling not valid for realGARCH model.")
			fit.control$scale=0
		}
	}
	if(is.null(fit.control$rec.init)) fit.control$rec.init = 'all'
	mm = match(names(fit.control), c("stationarity", "fixed.se", "scale", "rec.init"))
	if(any(is.na(mm))){
		idx = which(is.na(mm))
		enx = NULL
		for(i in 1:length(idx)) enx = c(enx, names(fit.control)[idx[i]])
		warning(paste(c("unidentified option(s) in fit.control:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
	}
	# if we have external regressors in variance turn off scaling
	if(spec@model$modelinc[15] > 0) fit.control$scale = FALSE
	# if we have arch-in-mean turn off scaling
	if(spec@model$modelinc[5] > 0) fit.control$scale = FALSE
	# if there are fixed pars we do no allow scaling as there would be no way of mixing scaled
	# amd non scaled parameters	
	if(sum(spec@model$pars[,2]) > 0) fit.control$scale = FALSE

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
	
	if(is.null(realizedVol)){
		stop("\nugarchfit-->error: you must supply the realized volatility (realizedVol) for the realGARCH model\n")
	} else{
		if(!is(realizedVol, "xts")) stop("\nugarchfit-->error: realizedVol must be an xts object\n")
		realizedVolIndex = format(index(realizedVol), format="%Y-%m-%d")
		if(!all.equal(realizedVolIndex, format(xdata$index, format="%Y-%m-%d"))) stop("\nugarchfit-->error: realizedVol must match the time index of the data")
		if(ncol(realizedVol)>1) stop("\nugarchfit-->error: realizedVol should be a univariate series\n")
	}
	
	
	# arglist replaces the use of custom environments (resolves the problem of 
	# non-shared environment in parallel estimation, particularly windows)
	garchenv = new.env(hash = TRUE)
	arglist = list()
	###################
	arglist$garchenv <- garchenv
	arglist$pmode = 0
	model = spec@model
	modelinc = model$modelinc
	pidx = model$pidx
	arglist$realized = coredata(realizedVol)
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
	# store length of data for easy retrieval
	model$modeldata$T = T = length(as.numeric(data))
	dist = model$modeldesc$distribution
	# if(fit.control$scale) dscale = sd(as.numeric(data)) else dscale = 1
	dscale = 1
	# zdata = data/dscale
	zdata = data
	recinit = .checkrec(fit.control$rec.init, T)
	arglist$data = zdata
	arglist$recinit = recinit
	arglist$dscale = dscale
	arglist$model = model
	ipars = model$pars
	# Optimization Starting Parameters Vector & Bounds
	tmp = .garchstart(ipars, arglist)
	####################################################
	# Scale Omega in order to avoid problems: (25-01-2013)
	#arglist$omegascale = 1/1e8
	#tmp$ipars["omega",c(1,5,6)] = tmp$ipars["omega",c(1,5,6)]  *1e8
	####################################################
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
				return(ugarchfilter(data = xts(origdata, origindex), spec = spec, out.sample = out.sample, realizedVol = realizedVol))
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
	# assign solver constraints (solnp directly else exterior type penalty
	# for other solvers)
	if(fit.control$stationarity == 1 && modelinc[15] == 0){
		Ifn = .realgarchcon
		ILB = 0
		IUB = 0.999999
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
		if(modelinc[1] > 0) parscale["mu"] = abs(mean(as.numeric(zdata)))
		if(modelinc[7] > 0) parscale["omega"] = var(as.numeric(zdata))
		arglist$returnType = "llh"
		solution = .garchsolver(solver, pars = ipars[estidx, 1], fun = .realgarchLLH, 
				Ifn, ILB, IUB, gr = NULL, hessian = NULL, parscale = parscale, 
				control = solver.control, LB = ipars[estidx, 5], 
				UB = ipars[estidx, 6], ux = NULL, ci = NULL, mu = NULL, arglist)
		sol = solution$sol
		hess = solution$hess
		timer = Sys.time()-tic
		if(!is.null(sol$par)){
			ipars[estidx, 1] = sol$par
			if(modelinc[7]==0){
				# call it once more to get omega
				tmpx = .realgarchLLH(sol$par, arglist)
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
	# create a copy of ipars in case we need to change it below to calculate standard errors
	# which we will need to reset later (because for example, infocriteria uses estimated
	# parameters, not fixed.
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
		fit = .makefitmodel(garchmodel = "realGARCH", f = .realgarchLLH, T = T, m = m, 
				timer = timer, convergence = convergence, message = sol$message, 
				hess, arglist = arglist, numderiv.control = numderiv.control)
		model$modeldata$realizedVol = realizedVol
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
		########################################################################
		
		fit$z = fit$residuals/fit$sigma
	} else{
		fit$message = sol$message
		fit$convergence = 1
		model$modeldata$realizedVol = realizedVol
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
# SECTION sGARCH LLH
#---------------------------------------------------------------------------------
.realgarchLLH = function(pars, arglist)
{
	# prepare inputs
	data = arglist$data
	realized = arglist$realized
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
	rx = .arfimaxfilter(modelinc[1:21], pars = ipars[,1], idx = idx, mexdata = mexdata, h = hm, data = data, N = N)
	res = rx$res
	zrf = rx$zrf
	res[is.na(res) | !is.finite(res) | is.nan(res)] = 0
	mvar = ifelse(recinit$type==1, mean(res[1:recinit$n]*res[1:recinit$n]), backcastv(res, T, recinit$n))
	if(modelinc[15]>0) {
		mv = sum(apply(matrix(vexdata, ncol = modelinc[15]), 2, "mean")*ipars[idx["vxreg",1]:idx["vxreg",2],1])
	} else{
		mv = 0
	}
	hEst = mvar
	# variance targeting
	persist = .persistrealgarch1(ipars, idx, distribution)
	if(modelinc[7]==0){
		mvar2 = ifelse(!is.na(modelinc[22]), modelinc[22]/dscale, mvar)
		ipars[idx["omega",1],1] = log(mvar2) * max(1 - persist, 0.001) - mv
		assign("omega", ipars[idx["omega",1],1], garchenv)
	}
	# E[eres^2]=1
	if(is.na(hEst) | !is.finite(hEst) | is.nan(hEst)) hEst = var(data)
	# if we have external regressors in variance equation we cannot have
	# stationarity checks in likelihood. Stationarity conditions not valid for
	# solnp solver which implements them internally as constraints.
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
	# modelinc (1:21) since 22 is either NA or numeric and no longer needed (paased to ipars if used)
	ans = try( .C("realgarchfilterC", model = as.integer(modelinc[1:21]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					hEst = as.double(hEst), x = as.double(data), res = as.double(res), 
					mexdata = mexdata, vexdata = vexdata, 
					zrf = as.double(zrf), constm = double(T), condm = double(T), 
					m = as.integer(m), T = as.integer(T), h = double(T), 
					z = double(T), tau = double(T), r = as.double(realized)[1:T], 
					u = double(T), llh = double(1), LHT1P = double(T), 
					LHT = double(T), 
					PACKAGE = "rugarch"), silent = TRUE )
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
	epsx = res
	llh = ans$llh
	u = ans$u
	tau = ans$tau
	LHT1P = ans$LHT1P
	
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
			all = list(llh = llh, h = h, epsx = epsx, z = z, LHT = LHT, persistence = persist, u = u, tau = tau, 
					LHT1P = -LHT1P))
	return(ans)
}

#---------------------------------------------------------------------------------
# SECTION sGARCH filter
#---------------------------------------------------------------------------------
.realgarchfilter = function(spec, data, out.sample = 0, n.old = NULL, rec.init = 'all', realizedVol)
{
	# n.old is optional and indicates the length of the original dataseries (in
	# cases when this represents a dataseries augmented by newer data). The reason
	# for using this is so that the old and new datasets agree since the original
	# recursion uses the sum of the residuals to start the recursion and therefore
	# is influenced by new data. For a small augmentation the values converge after
	# x periods, but it is sometimes preferable to have this option so that there is
	# no forward looking information contaminating the study.
	if(missing(realizedVol)){
		stop("\nugarchfilter-->error: you must supply the realized volatility (realizedVol) for the realGARCH model\n")
	} else{
		if(!is(realizedVol, "xts")) stop("\nugarchfilter-->error: realizedVol must be an xts object\n")
		realizedVolIndex = format(index(realizedVol), format="%Y-%m-%d")
		if(!all.equal(realizedVolIndex, format(index(data), format="%Y-%m-%d"))) stop("\nugarchfit-->error: realizedVol must match the time index of the data")
		if(ncol(realizedVol)>1) stop("\nugarchfilter-->error: realizedVol should be a univariate series\n")
	}
	realized = coredata(realizedVol)
	# we are not going to extract the data just yet
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
	# NB Any changes made to the spec are not preserved once we apply set fixed
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
	rx = .arfimaxfilter(modelinc[1:21], ipars[,1], idx, mexdata = mexdata, h = 0, data = as.numeric(data), N = N)
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
	persist = .persistrealgarch1(ipars, idx, distribution)
	if(modelinc[7]==0){
		mvar2 = ifelse(!is.na(modelinc[22]), modelinc[22], mvar)
		ipars[idx["omega",1],1] = log(mvar2) * max(1 - persist, 0.001) - mv
	}
	if(modelinc[6]>0) mexdata = as.double(as.vector(mexdata)) else mexdata = double(1)
	if(modelinc[15]>0) vexdata = as.double(as.vector(vexdata)) else vexdata = double(1)
	
	ans = try( .C("realgarchfilterC", model = as.integer(modelinc[1:21]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					hEst = as.double(hEst), x = as.double(data), res = as.double(res), 
					mexdata = mexdata, vexdata = vexdata, 
					zrf = as.double(zrf), constm = double(T), condm = double(T), 
					m = as.integer(m), T = as.integer(T), h = double(T), 
					z = double(T), tau = double(T), r = as.double(realized)[1:T], 
					u = double(T), llh = double(1), LHT1P = double(T), 
					LHT = double(T), 
					PACKAGE = "rugarch"), silent = TRUE )
	
	filter = list()
	filter$residuals = res
	filter$LLH = -ans$llh
	filter$log.likelihoods = ans$LHT
	filter$partial.log.likelihoods = ans$LHT1P
	filter$u = ans$u
	filter$tau = ans$tau
	filter$persistence = persist
	filter$distribution = distribution
	filter$ipars = ipars
	model$modeldata$data = origdata
	model$modeldata$realizedVol = realizedVol
	model$modeldata$index = origindex
	model$modeldata$period = period
	model$n.start = out.sample
	filter$sigma = sqrt(ans$h)
	filter$z = filter$residuals/filter$sigma
	sol = new("uGARCHfilter",
			filter = filter,
			model = model)
	return(sol)
}

#---------------------------------------------------------------------------------
# SECTION realGARCH forecast
#---------------------------------------------------------------------------------
.realgarchforecast = function(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0, 
		external.forecasts = list(mregfor = NULL, vregfor = NULL), n.sim = 1000, returnDistribution = TRUE, 
		...)
{
	# will not allow to combine n.ahead and n.roll in this model
	fit    = fitORspec
	data   = fit@model$modeldata$data
	realized = coredata(fit@model$modeldata$realizedVol)[,1]
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
	xreg = .forcregressors(model, external.forecasts$mregfor, external.forecasts$vregfor, n.ahead, Nor, out.sample = ns, n.roll)
	mxf = xreg$mxf
	vxf = xreg$vxf
	
	# filter data (check external regressor data - must equal length of origData)
	fcreq = ifelse(ns >= (n.ahead+n.roll), n.ahead+n.roll, ns)
	fspec = ugarchspec(variance.model = list(model = "realGARCH", 
					garchOrder = c(modelinc[8], modelinc[9]), submodel = NULL, 
					external.regressors = vxf[1:(N + fcreq), , drop = FALSE]), 
			mean.model = list(armaOrder = c(modelinc[2], modelinc[3]),
					include.mean = modelinc[1], 
					archm = ifelse(modelinc[5]>0,TRUE,FALSE), archpow = modelinc[5], arfima = modelinc[4], 
					external.regressors = mxf[1:(N + fcreq), , drop = FALSE], archex = modelinc[20]), 
			distribution.model = model$modeldesc$distribution, fixed.pars = as.list(pars))
	
	tmp =  xts(data[1:(N + fcreq)], index[1:(N + fcreq)])	
	realtmp = xts(realized[1:(N + fcreq)], index[1:(N + fcreq)])	
	flt = .realgarchfilter(data = tmp, spec = fspec, n.old = N, realizedVol = realtmp)
	sigmafilter = flt@filter$sigma
	resfilter = flt@filter$residuals
	zfilter = flt@filter$z
	ufilter = flt@filter$u
	taufilter = flt@filter$tau
	seriesfilter = as.numeric(fitted(flt))
	# forecast GARCH process
	realizedFor = seriesFor = sigmaFor = matrix(NA, ncol = n.roll+1, nrow = n.ahead)
	colnames(realizedFor) = colnames(seriesFor) = colnames(sigmaFor) = as.character(index[N:(N+n.roll)])
	rownames(realizedFor) = rownames(sigmaFor) = paste("T+", 1:n.ahead, sep="")
	if(returnDistribution){
		realizedDF = sigmaDF = array(NA, dim = c(n.ahead, n.sim, n.roll+1), 
				dimnames=list(paste("T+", 1:n.ahead, sep=""), paste("Sim_", 1:n.sim, sep=""), as.character(index[N:(N+n.roll)])))
	}
	for(i in 1:(n.roll+1)){
		np = N + i - 1
		if(modelinc[1] > 0){
			mu = rep(ipars[idx["mu",1]:idx["mu",2], 1], N+i+n.ahead-1)
		} else{
			mu = rep(0, N+i+n.ahead-1)
		}
		omega = rep(ipars[idx["omega",1]:idx["omega",2], 1], N+i+n.ahead-1)
		h = c(sigmafilter[1:(N+i-1)], rep(0, n.ahead))
		r = c(as.numeric(realtmp[1:(N+i-1),1]), rep(0, n.ahead))
		epsx = c(resfilter[1:(N+i-1)], rep(0, n.ahead))
		x = c(data[1:(N+i-1)], rep(0, n.ahead))
		z = c(zfilter[1:(N+i-1)], rep(0, n.ahead))
		# forecast of externals is provided outside the system
		mxfi = mxf[1:(N+i-1+n.ahead), , drop = FALSE]
		vxfi = vxf[1:(N+i-1+n.ahead), , drop = FALSE]
		ans = .nrealgarchforecast(ipars, modelinc, idx, mu, mxfi, omega, vxfi, h, realized = r, epsx, z, data = x, 
				N = np, n.ahead, n.sim, returnD = returnDistribution)
		sigmaFor[,i] = ans$h
		seriesFor[,i] = ans$x
		realizedFor[,i] = ans$r
		if(returnDistribution){
			realizedDF[,,i] = ans$rforcdist
			sigmaDF[,,i] = ans$hforcdist
		}
	}
	fcst = list()
	fcst$n.ahead = n.ahead
	fcst$N = N+ns
	fcst$n.start = ns
	fcst$n.roll = n.roll
	fcst$sigmaFor = sigmaFor
	fcst$realizedFor = realizedFor
	fcst$seriesFor = seriesFor
	if(returnDistribution){
		fcst$sigmaDF = sigmaDF
		fcst$realizedDF = realizedDF
	} else{
		fcst$sigmaDF = NULL
		fcst$realizedDF = NULL
	}
	model$modeldata$sigma = flt@filter$sigma
	model$modeldata$residuals = flt@filter$residuals
	ans = new("uGARCHforecast",
			forecast = fcst,
			model = model)
	return(ans)
}

.taufn = function(z, tau1, tau2){ tau1*z + tau2*(z^2-1) }

.nrealgarchforecast = function(ipars, modelinc, idx, mu, mxfi, omega, vxfi, h, realized, epsx, z, 
		data, N, n.ahead, n.sim = 1000, returnD = FALSE)
{
	p = modelinc[9]
	q = modelinc[8]
	beta = ipars[idx["beta",1]:idx["beta",2],1]
	alpha = ipars[idx["alpha",1]:idx["alpha",2],1]
	xi = ipars[idx["xi",1]:idx["xi",2],1]
	delta = ipars[idx["delta",1]:idx["delta",2],1]
	lambda = ipars[idx["lambda",1]:idx["lambda",2],1]
	tau1 = ipars[idx["eta1",1]:idx["eta1",2],1]
	tau2 = ipars[idx["eta2",1]:idx["eta2",2],1]
	dist = c("norm", "snorm", "std", "sstd","ged", "sged", "nig", "ghyp", "jsu", "ghst")[modelinc[21]]
	if(modelinc[15]>0){
		omega = omega + vxfi%*%t(matrix(ipars[idx["vxreg",1]:idx["vxreg",2],1], ncol = modelinc[15]))
	}
	
	zs = lapply(1:n.ahead, function(i) rdist(dist, n = n.sim, mu = 0, sigma = 1, skew = ipars[idx["skew",1],1], shape = ipars[idx["shape",1],1],
						lambda = ipars[idx["ghlambda",1],1]))
	u = lapply(1:n.ahead, function(i) rnorm(n.sim, 0, lambda))
	Y = c(rev(tail(h[1:N]^2,p)), rev(tail(realized[1:N],q)))
	Y = log(Y)
	A = cbind(rbind(matrix(beta, nrow=1), if(p>1) cbind(diag(p-1), matrix(0, nrow=p-1,ncol=1)) else NULL, 
					matrix(delta*beta, nrow=1),
					if(q>1 && p>0) matrix(0, nrow=q-1, ncol=p) else NULL),
			rbind(matrix(alpha, nrow=1), if(p>1 && q>0) matrix(0, nrow=p-1,ncol=q) else NULL, 
					matrix(delta*alpha, nrow=1),
					if(q>1) cbind(diag(q-1), matrix(0, nrow=q-1,ncol=1)) else NULL))
	e = lapply(1:n.ahead, function(i) rbind(if(p>0) matrix(0, nrow=p, ncol=n.sim) else NULL, 
						matrix(.taufn(zs[[i]], tau1, tau2)+u[[i]], nrow=1, ncol=n.sim), 
						if(q>1) matrix(0, nrow=q-1, ncol=n.sim) else NULL))
	sigmaDFor = matrix(NA, nrow = n.ahead, ncol = n.sim)
	realizedDFor = matrix(NA, nrow = n.ahead, ncol = n.sim)
	xDFor = matrix(NA, ncol = n.ahead, nrow = n.sim)
	hFor = realizedFor = xFor = rep(0, n.ahead)
	n = ncol(e[[1]])
	D = matrix(0, ncol = n, nrow=nrow(e[[1]]))
	for(i in 1:n.ahead){
		# we can absorb vxfi into omega and make sure to adjust b at every iteration
		b = c(omega[N+i], if(p>1) rep(0, p-1) else NULL, xi+delta*omega[N+i], if(q>1) rep(0, q-1) else NULL)
		D = epsfn(A, e[[i]], b, i, D)
		tmp = matrix(A%^%i%*%Y, nrow = length(Y), ncol = ncol(D)) + D
		sigmaDFor[i,] = sqrt(exp(tmp[1,]))
		realizedDFor[i,] = exp(tmp[p+1,])
		h[N+i] = sqrt(exp(mean(tmp[1,])))
		realized[N+i] = exp(mean(tmp[p+1,]))
	}
	if(modelinc[4]>0){
		res = arfimaf(ipars, modelinc[1:21], idx, mu, mxfi, h, epsx, z, data, N, n.ahead)
	} else{
		res = armaf(ipars, modelinc[1:21], idx, mu, mxfi, h, epsx, z, data, N, n.ahead)
	}
	sol = list(h = h[(N+1):(N+n.ahead)], x = res[(N+1):(N+n.ahead)], r = realized[(N+1):(N+n.ahead)], 
			hforcdist = if(returnD) sigmaDFor else NULL,
			rforcdist = if(returnD) realizedDFor else NULL)
	return(sol)
}


epsfn = function(A, e, b, n.ahead, ed)
{
	n = ncol(e)
	for(i in 1:n){
		ed[,i] = ed[,i]+(A%^%(n.ahead-1) %*% (b+e[,i]))
	}
	return(ed)
}
#---------------------------------------------------------------------------------
# 2nd dispatch method for forecast
.realgarchforecast2 = function(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0, 
		external.forecasts = list(mregfor = NULL, vregfor = NULL), realizedVol = NULL, 
		n.sim = 1000, returnDistribution = TRUE, ...)
{
	# first we filter the data to get the results:
	if(is.null(realizedVol)){
		stop("\nugarchforecast-->error: you must supply the realized volatility (realizedVol) for the realGARCH model\n")
	} else{
		if(!is(realizedVol, "xts")) stop("\nugarchforecast-->error: realizedVol must be an xts object\n")
		realizedVolIndex = format(index(realizedVol), format="%Y-%m-%d")
		if(!all.equal(realizedVolIndex, format(index(data), format="%Y-%m-%d"))) stop("\nugarchfit-->error: realizedVol must match the time index of the data")
		if(ncol(realizedVol)>1) stop("\nugarchforecast-->error: realizedVol should be a univariate series\n")
	}
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
	realized = coredata(realizedVol)
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
	xreg = .forcregressors(model, external.forecasts$mregfor, external.forecasts$vregfor, n.ahead, Nor, out.sample = ns, n.roll)
	mxf = xreg$mxf
	vxf = xreg$vxf
	
	# filter data (check external regressor data - must equal length of origData)
	fcreq = ifelse(ns >= (n.ahead+n.roll), n.ahead+n.roll, ns)
	fspec = ugarchspec(variance.model = list(model = "realGARCH", 
					garchOrder = c(modelinc[8], modelinc[9]), submodel = NULL, 
					external.regressors = vxf[1:(N + fcreq), , drop = FALSE]), 
			mean.model = list(armaOrder = c(modelinc[2], modelinc[3]),
					include.mean = modelinc[1], 
					archm = ifelse(modelinc[5]>0,TRUE,FALSE), archpow = modelinc[5], arfima = modelinc[4], 
					external.regressors = mxf[1:(N + fcreq), , drop = FALSE], archex = modelinc[20]),
			distribution.model = model$modeldesc$distribution, fixed.pars = as.list(pars))
	tmp =  xts(data[1:(N + fcreq)], index[1:(N + fcreq)])	
	realtmp = xts(realized[1:(N + fcreq)], index[1:(N + fcreq)])	
	flt = .realgarchfilter(data = tmp, spec = fspec, n.old = N, realizedVol = realtmp)
	sigmafilter = flt@filter$sigma
	resfilter = flt@filter$residuals
	zfilter = flt@filter$z
	ufilter = flt@filter$u
	taufilter = flt@filter$tau
	seriesfilter = as.numeric(fitted(flt))
	# forecast GARCH process
	realizedFor = seriesFor = sigmaFor = matrix(NA, ncol = n.roll+1, nrow = n.ahead)
	colnames(realizedFor) = colnames(seriesFor) = colnames(sigmaFor) = as.character(index[N:(N+n.roll)])
	rownames(realizedFor) = rownames(sigmaFor) = paste("T+", 1:n.ahead, sep="")
	if(returnDistribution){
		realizedDF = sigmaDF = array(NA, dim = c(n.ahead, n.sim, n.roll+1), 
				dimnames=list(paste("T+", 1:n.ahead, sep=""), paste("Sim_", 1:n.sim, sep=""), as.character(index[N:(N+n.roll)])))
	}
	for(i in 1:(n.roll+1)){
		np = N + i - 1
		if(modelinc[1] > 0){
			mu = rep(ipars[idx["mu",1]:idx["mu",2], 1], N+i+n.ahead-1)
		} else{
			mu = rep(0, N+i+n.ahead-1)
		}
		omega = rep(ipars[idx["omega",1]:idx["omega",2], 1], N+i+n.ahead-1)
		h = c(sigmafilter[1:(N+i-1)], rep(0, n.ahead))
		r = c(as.numeric(realtmp[1:(N+i-1),1]), rep(0, n.ahead))
		epsx = c(resfilter[1:(N+i-1)], rep(0, n.ahead))
		x = c(data[1:(N+i-1)], rep(0, n.ahead))
		z = c(zfilter[1:(N+i-1)], rep(0, n.ahead))
		# forecast of externals is provided outside the system
		mxfi = mxf[1:(N+i-1+n.ahead), , drop = FALSE]
		vxfi = vxf[1:(N+i-1+n.ahead), , drop = FALSE]
		ans = .nrealgarchforecast(ipars, modelinc, idx, mu, mxfi, omega, vxfi, h, realized = r, epsx, z, data = x, 
				N = np, n.ahead, n.sim, returnD = returnDistribution)
		sigmaFor[,i] = ans$h
		seriesFor[,i] = ans$x
		realizedFor[,i] = ans$r
		if(returnDistribution){
			realizedDF[,,i] = ans$rforcdist
			sigmaDF[,,i] = ans$hforcdist
		}
	}
	fcst = list()
	fcst$n.ahead = n.ahead
	fcst$N = N+ns
	fcst$n.start = ns
	fcst$n.roll = n.roll
	fcst$sigmaFor = sigmaFor
	fcst$seriesFor = seriesFor
	fcst$realizedFor = realizedFor
	if(returnDistribution){
		fcst$sigmaDF = sigmaDF
		fcst$realizedDF = realizedDF
	} else{
		fcst$sigmaDF = NULL
		fcst$realizedDF = NULL
	}
	model$modeldata$sigma = flt@filter$sigma
	model$modeldata$realizedVol = realizedVol
	model$modeldata$residuals = flt@filter$residuals
	ans = new("uGARCHforecast",
			forecast = fcst,
			model = model)
	return(ans)
}
#---------------------------------------------------------------------------------
# SECTION realGARCH simulate
#---------------------------------------------------------------------------------
# DGP and conditional mean (\mu)
# y_{t+1} = \mu_{t+1} + \sigma_{t+1}z_{t+1}
# conditional variance
# \sigma^2_{t+1} = \omega+\beta_{t+1}\sigma^2_{t} + \alpha r_{t}
# measurement equation
# r_{t+1} = \xi + \delta sigma^2_{t+1} + \tau(z_{t+1}) + u_{t+1}
# 1. simulate z_{t+1}
# 2. Generate \sigma_{t+1}
# 3. simulate u_{t+1}
# 4. simulate r_{t+1} using \sigma_{t+1} and u_{t+1} and \tau(z_{t+1})
# 5. simulate y_{t+1} using \sigma_{t+1} and z_{t+1}
# Given a starting value for r_t we do not need any additional inputs for n.sim

.realgarchsim = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, startMethod = 
				c("unconditional","sample"), presigma = NA, prereturns = NA, 
		preresiduals = NA, rseed = NA, custom.dist = list(name = NA, distfit = NA), 
		mexsimdata = NULL, vexsimdata = NULL, prerealized = NA, ...)
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
	realized = coredata(fit@model$modeldata$realizedVol)[,1]
	N = length(as.numeric(data))
	data = data[1:(N - fit@model$n.start)]
	realized = realized[1:(N - fit@model$n.start)]
	N = length(as.numeric(data))
	m = fit@model$maxOrder
	resids = fit@fit$residuals
	sigma = fit@fit$sigma
	u = fit@fit$u
	tau = fit@fit$tau
	model = fit@model
	modelinc = model$modelinc
	idx = model$pidx
	ipars = fit@fit$ipars
	# check if necessary the external regressor forecasts provided first
	xreg = .simregressors(model, mexsimdata, vexsimdata, N, n, m.sim, m)	
	mexsim = xreg$mexsimlist
	vexsim = xreg$vexsimlist
	
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
	set.seed(sseed[1])
	usim = rbind(matrix(tail(u,m), nrow=m, ncol=m.sim), matrix(rnorm(m.sim*n, 0, ipars[idx["lambda",1]]), nrow=n, ncol=m.sim))
	
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
	if(!is.na(prerealized[1])){
		prerealized = as.vector(prerealized)
		if(length(prerealized)<m) stop(paste("\nugarchsim-->error: prerealized must be of length ", m, sep=""))
		prereal = matrix(prerealized, nrow = m)
	}
	if(is.na(presigma[1])){
		if(startMethod[1] == "unconditional"){
			hEst = uncvariance(fit)^(1/2)
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
	if(is.na(prerealized[1])){
		if(startMethod[1] == "unconditional"){
			xEst = exp(.uncrealgarchr(ipars, idx))
			prereal = as.numeric(rep(xEst, m))}
		else{
			prereal  = tail(realized, m)
		}
	}
	
	# input vectors/matrices
	h = c(presigma^2, rep(0, n))
	x = c(prereturns, rep(0, n))
	r = c(prereal, rep(0, n))
	constm = matrix(ipars[idx["mu",1]:idx["mu",2], 1], ncol = m.sim, nrow = n + m)
	
	# outpus matrices
	sigmaSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	realizedSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	tauSim =  matrix(0, ncol = m.sim, nrow = n.sim)
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
		ans1 = try(.C("realgarchsimC", model = as.integer(modelinc[1:21]), 
						pars = as.double(ipars[,1]), 
						idx = as.integer(idx[,1]-1), 
						res = as.double(res),
						vexdata = as.double(vexsim[[i]]),
						m = as.integer(m),
						T = as.integer(n+m),
						h = as.double(h), 
						z = as.double(z[,i]),
						tau = double(n+m),
						r = as.double(r),
						u = as.double(usim[,i]),
						PACKAGE = "rugarch"), silent = TRUE)
		
		if(inherits(ans1, "try-error")) stop("\nugarchsim-->error: error in calling C function....\n")
		
		sigmaSim[,i] = ans1$h[(n.start + m + 1):(n+m)]^(1/2)
		residSim[,i] = ans1$res[(n.start + m + 1):(n+m)]
		realizedSim[,i] = ans1$r[(n.start + m + 1):(n+m)]
		tauSim[,i] = ans1$tau[(n.start + m + 1):(n+m)]
		# ToDo: change to accomodate modelinc[20]
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
	sim = list(sigmaSim = sigmaSim, seriesSim = seriesSim, residSim = residSim, realizedSim = realizedSim, uSim = usim,
			tauSim = tauSim)
	model$modeldata$sigma = sigma
	sol = new("uGARCHsim",
			simulation = sim,
			model = model,
			seed = as.integer(sseed))
	return(sol)
}
#---------------------------------------------------------------------------------
# SECTION sGARCH path
#---------------------------------------------------------------------------------
.realgarchpath = function(spec, n.sim = 1000, n.start = 0, m.sim = 1,
		presigma = NA, prereturns = NA, preresiduals = NA, rseed = NA, 
		custom.dist = list(name = NA, distfit = NA), mexsimdata = NULL, 
		vexsimdata = NULL, prerealized = NA, ...)
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
	persist = persistence(spec)
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
	set.seed(sseed[1])
	usim = rbind(matrix(rep(0,m), nrow=m, ncol=m.sim), matrix(rnorm(m.sim*n, 0, ipars[idx["lambda",1]]), nrow=n, ncol=m.sim))
	
	# create the presample information
	if(!is.na(presigma[1])){
		presigma = as.vector(presigma)
		if(length(presigma)<m) stop(paste("\nugarchpath-->error: presigma must be of length ", m, sep=""))
	}
	
	if(!is.na(prereturns[1])){
		prereturns = as.vector(prereturns)
		if(length(prereturns)<m) stop(paste("\nugarchpath-->error: prereturns must be of length ", m, sep=""))
	}
	
	if(!is.na(preresiduals[1])){
		preresiduals = as.vector(preresiduals)
		if(length(preresiduals)<m) stop(paste("\nugarchpath-->error: preresiduals must be of length ", m, sep=""))
		preres = matrix(preresiduals, nrow = m)
	}
	if(!is.na(prerealized[1])){
		prerealized = as.vector(prerealized)
		if(length(prerealized)<m) stop(paste("\nugarchsim-->error: prerealized must be of length ", m, sep=""))
		prereal = matrix(prerealized, nrow = m)
	}
	if(is.na(presigma[1])){
		hEst = uncvariance(spec)^(1/2)
		presigma = as.numeric(rep(hEst, times = m))
	}
	if(is.na(prereturns[1])){
		prereturns = as.numeric(rep(uncmean(spec), times = m))
	}
	if(is.na(prerealized[1])){
		xEst = exp(.uncrealgarchr(ipars, idx))
		prereal = as.numeric(rep(xEst, m))
	}
	# input vectors/matrices
	h = c(presigma^2, rep(0, n))
	x = c(prereturns, rep(0, n))
	r = c(prereal, rep(0, n))
	constm = matrix(ipars[idx["mu",1]:idx["mu",2],1], ncol = m.sim, nrow = n + m)
	
	# outpus matrices
	sigmaSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	realizedSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	tauSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	seriesSim = matrix(0, ncol = m.sim, nrow = n.sim)
	residSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	
	for(i in 1:m.sim){
		if(is.na(preresiduals[1])){
			preres = as.numeric(z[1:m,i])*presigma
		}
		z[1:m, 1:m.sim] = preres[1:m]/presigma[1:m]
		z[is.na(z) | is.nan(z) | !is.finite(z)] = 0
		res = c(preres, rep(0, n))
		ans1 = try(.C("realgarchsimC", model = as.integer(modelinc[1:21]), 
						pars = as.double(ipars[,1]), 
						idx = as.integer(idx[,1]-1), 
						res = as.double(res),
						vexdata = as.double(vexsim[[i]]),
						m = as.integer(m),
						T = as.integer(n+m),
						h = as.double(h), 
						z = as.double(z[,i]),
						tau = double(n+m),
						r = as.double(r),
						u = as.double(usim[,i]),
						PACKAGE = "rugarch"), silent = TRUE)
		if(inherits(ans1, "try-error")) stop("\nugarchpath-->error: error in calling C function....\n")
		sigmaSim[,i] = ans1$h[(n.start + m + 1):(n+m)]^(1/2)
		realizedSim[,i] = ans1$r[(n.start + m + 1):(n+m)]
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
	
	path = list(sigmaSim = sigmaSim, seriesSim = seriesSim, residSim = residSim, realizedSim = realizedSim, uSim = usim,
			tauSim = tauSim)	
	sol = new("uGARCHpath",
			path = path,
			model = model,
			seed = as.integer(sseed))
	return(sol)
}
################################################################################
# special purpose rolling method for the mcsGARCH model
################################################################################

.rollfdensity.real = function(spec, data, n.ahead = 1, forecast.length = 500, 
		n.start = NULL, refit.every = 25, refit.window = c("recursive", "moving"), 
		window.size = NULL, solver = "hybrid", fit.control = list(), solver.control = list(),
		calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), cluster = NULL,
		keep.coef = TRUE, realizedVol = NULL, n.sim=1000,...)
{
	tic = Sys.time()
	if(is.null(realizedVol)){
		stop("\nugarchroll-->error: you must supply the realized volatility (realizedVol) for the realGARCH model\n")
	} else{
		if(!is(realizedVol, "xts")) stop("\nugarchroll-->error: realizedVol must be an xts object\n")
		realizedVolIndex = format(index(realizedVol), format="%Y-%m-%d")
		if(!all.equal(realizedVolIndex, format(index(data), format="%Y-%m-%d"))) stop("\nugarchroll-->error: realizedVol must match the time index of the data")
		if(ncol(realizedVol)>1) stop("\nugarchroll-->error: realizedVol should be a univariate series\n")
	}
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
	refit.window = refit.window[1]
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
	if( modelinc[15]>0 ){
		vexdata = spec@model$modeldata$vexdata
		vex = TRUE
	} else{
		vex = FALSE
		vexdata = NULL
	}
	if(n.ahead>1) stop("\nugarchroll:--> n.ahead>1 not supported...try again.")
	if(is.null(n.start)){
		if(is.null(forecast.length)) stop("\nugarchroll:--> forecast.length amd n.start are both NULL....try again.")
		n.start = T - forecast.length
	} else{
		forecast.length = T - n.start
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
	
	if( !is.null(cluster) ){
		clusterEvalQ(cl = cluster, library(rugarch))
		clusterExport(cluster, c("data", "index", "s","refit.every", "realizedVol","n.sim",
						"keep.coef", "shaped", "skewed", "ghyp", 
						"rollind", "spec", "out.sample", "mex", "vex",
						"solver", "solver.control", "fit.control"), envir = environment())
		if(mex) clusterExport(cluster, c("mexdata"), envir = environment())
		if(vex) clusterExport(cluster, c("vexdata"), envir = environment())
		tmp = parLapply(cl = cluster, 1:m, fun = function(i){
					if(mex) spec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
					if(vex) spec@model$modeldata$vexdata = vexdata[rollind[[i]],,drop=FALSE]
					fit = try(ugarchfit(spec, xts::xts(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i], 
									solver = solver, solver.control = solver.control, 
									fit.control = fit.control, realizedVol = realizedVol[rollind[[i]],1]), silent=TRUE)
					# 3 cases: General Error, Failure to Converge, Failure to invert Hessian (bad solution)
					if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
						ans = list(y = NA, cf = NA, converge = FALSE, lik = NA, plik = NA)
					} else{
						if(mex) fmex = tail(mexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fmex = NULL
						if(vex) fvex = tail(vexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fvex = NULL
						f = ugarchforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1, 
								external.forecasts = list(mregfor = fmex, vregfor = fvex), returnDistribution=FALSE,
								n.sim=n.sim)
						sig = as.numeric( sigma(f) )
						ret = as.numeric( fitted(f) )
						realized = as.numeric(f@forecast$realizedFor)
						if(shaped) shp = rep(coef(fit)["shape"], out.sample[i]) else shp = rep(0, out.sample[i])
						if(skewed) skw = rep(coef(fit)["skew"], out.sample[i]) else skw = rep(0, out.sample[i])
						if(ghyp) shpgig = rep(coef(fit)["ghlambda"], out.sample[i]) else shpgig = rep(0, out.sample[i])
						rlz = tail(data[rollind[[i]]], out.sample[i])
						realv = tail(as.numeric(realizedVol[rollind[[i]],1]), out.sample[i])
						# use xts for indexing the forecasts
						y = as.data.frame(cbind(ret, sig, skw, shp, shpgig, rlz, realv, realized))
						rownames(y) = tail(as.character(index[rollind[[i]]]), out.sample[i])
						colnames(y) = c("Mu", "Sigma", "Skew", "Shape", "Shape(GIG)", "Realized", "RVol", "RVolForecast")
						if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
						ans = list(y = y, cf = cf, converge = TRUE, lik = likelihood(fit), plik = sum(-fit@fit$partial.log.likelihoods))
					}
					return(ans)})
	} else{
		tmp = lapply(as.list(1:m), FUN = function(i){
					if(mex) spec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
					if(vex) spec@model$modeldata$vexdata = vexdata[rollind[[i]],,drop=FALSE]
					fit = try(ugarchfit(spec, xts(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i], 
									solver = solver, solver.control = solver.control, 
									fit.control = fit.control, realizedVol = realizedVol[rollind[[i]],1]), silent=TRUE)
					if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
						ans = list(y = NA, cf = NA, converge = FALSE, lik = NA, plik = NA)
					} else{
						if(mex) fmex = tail(mexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fmex = NULL
						if(vex) fvex = tail(vexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fvex = NULL
						f = ugarchforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1, 
								external.forecasts = list(mregfor = fmex, vregfor = fvex), returnDistribution=FALSE,
								n.sim = n.sim)
						sig = as.numeric( sigma(f) )
						ret = as.numeric( fitted(f) )
						realized = as.numeric(f@forecast$realizedFor)
						if(shaped) shp = rep(coef(fit)["shape"], out.sample[i]) else shp = rep(0, out.sample[i])
						if(skewed) skw = rep(coef(fit)["skew"], out.sample[i]) else skw = rep(0, out.sample[i])
						if(ghyp) shpgig = rep(coef(fit)["ghlambda"], out.sample[i]) else shpgig = rep(0, out.sample[i])
						rlz = tail(data[rollind[[i]]], out.sample[i])
						realv = tail(as.numeric(realizedVol[rollind[[i]],1]), out.sample[i])
						# use xts for indexing the forecasts
						y = as.data.frame(cbind(ret, sig, skw, shp, shpgig, rlz, realv, realized))
						rownames(y) = tail(as.character(index[rollind[[i]]]), out.sample[i])
						colnames(y) = c("Mu", "Sigma", "Skew", "Shape", "Shape(GIG)", "Realized", "RVol", "RVolForecast")
						if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
						ans = list(y = y, cf = cf, converge = TRUE, lik = likelihood(fit), plik = sum(-fit@fit$partial.log.likelihoods))
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
		model$realizedVol = realizedVol
		model$n.sim = n.sim
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
		model$rollind = rollind
		model$out.sample = out.sample
		forecast = tmp
		toc = Sys.time()-tic
		model$elapsed = toc
		ans = new("uGARCHroll",
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
		model$rollind = rollind
		model$out.sample = out.sample
		model$loglik = sapply(tmp, function(x) x$lik)
		model$partial.loglik = sapply(tmp, function(x) x$plik)
		forecast = list(VaR = VaR.matrix, density = forc)
	}
	toc = Sys.time()-tic
	model$elapsed = toc
	ans = new("uGARCHroll",
			model = model,
			forecast = forecast)
	return( ans )
}

.resumeroll.real = function(object, solver = "hybrid", fit.control = list(), 
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
		realizedVol = model$realizedVol
		n.sim = model$n.sim
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
		if( modelinc[15]>0 ){
			vexdata = spec@model$modeldata$vexdata
			vex = TRUE
		} else{
			vex = FALSE
			vexdata = NULL
		}
		n.ahead = model$n.ahead
		n.start = model$n.start
		forecast.length = model$forecast.length
		refit.every = model$refit.every
		refit.window = model$refit.window
		window.size = model$window.size
		if(n.ahead>1) stop("\nugarchroll:--> n.ahead>1 not supported...try again.")
		if(is.null(n.start)){
			if(is.null(forecast.length)) stop("\nugarchroll:--> forecast.length amd n.start are both NULL....try again.")
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
		if( !is.null(cluster) ){
			clusterEvalQ(cl = cluster, library(rugarch))
			clusterExport(cluster, c("data", "index","s","refit.every","realizedVol","n.sim",
							"keep.coef", "shaped", "skewed", "ghyp", 
							"rollind", "spec", "out.sample", "mex", "vex", 
							"noncidx", "solver", "solver.control", "fit.control"),
					envir = environment())
			if(mex) clusterExport(cluster, c("mexdata"), envir = environment())
			if(vex) clusterExport(cluster, c("vexdata"), envir = environment())
			tmp = parLapply(cl = cluster, as.list(noncidx), fun = function(i){
						if(mex) spec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
						if(vex) spec@model$modeldata$vexdata = vexdata[rollind[[i]],,drop=FALSE]
						fit = try(ugarchfit(spec, xts::xts(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i], 
										solver=solver, solver.control = solver.control, 
										fit.control = fit.control, realizedVol = realizedVol[rollind[[i]],1]), silent=TRUE)
						if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
							ans = list(y = NA, cf = NA, converge = FALSE, lik = NA, plik = NA)
						} else{
							if(mex) fmex = tail(mexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fmex = NULL
							if(vex) fvex = tail(vexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fvex = NULL
							f = ugarchforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1, 
									external.forecasts = list(mregfor = fmex, vregfor = fvex), returnDistribution=FALSE,
									n.sim=n.sim)
							sig = as.numeric( sigma(f) )
							ret = as.numeric( fitted(f) )
							realized = as.numeric(f@forecast$realizedFor)
							if(shaped) shp = rep(coef(fit)["shape"], out.sample[i]) else shp = rep(0, out.sample[i])
							if(skewed) skw = rep(coef(fit)["skew"], out.sample[i]) else skw = rep(0, out.sample[i])
							if(ghyp) shpgig = rep(coef(fit)["ghlambda"], out.sample[i]) else shpgig = rep(0, out.sample[i])
							rlz = tail(data[rollind[[i]]], out.sample[i])
							realv = tail(as.numeric(realizedVol[rollind[[i]],1]), out.sample[i])
							# use xts for indexing the forecasts
							y = as.data.frame(cbind(ret, sig, skw, shp, shpgig, rlz, realv, realized))
							rownames(y) = tail(as.character(index[rollind[[i]]]), out.sample[i])
							colnames(y) = c("Mu", "Sigma", "Skew", "Shape", "Shape(GIG)", "Realized", "RVol", "RVolForecast")
							if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
							ans = list(y = y, cf = cf, converge = TRUE, lik = likelihood(fit), plik = sum(-fit@fit$partial.log.likelihoods))
						}
						return(ans)})
		} else{
			tmp = lapply(as.list(noncidx), FUN = function(i){
						if(mex) spec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
						if(vex) spec@model$modeldata$vexdata = vexdata[rollind[[i]],,drop=FALSE]
						fit = try(ugarchfit(spec, xts(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i], 
										solver=solver, solver.control = solver.control, 
										fit.control = fit.control, realizedVol = realizedVol[rollind[[i]],1]), silent=TRUE)
						if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
							ans = list(y = NA, cf = NA, converge = FALSE, lik = NA, plik = NA)
						} else{
							if(mex) fmex = tail(mexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fmex = NULL
							if(vex) fvex = tail(vexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fvex = NULL
							f = ugarchforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1, 
									external.forecasts = list(mregfor = fmex, vregfor = fvex), returnDistribution=FALSE,
									n.sim = n.sim)
							sig = as.numeric( sigma(f) )
							ret = as.numeric( fitted(f) )
							realized = as.numeric(f@forecast$realizedFor)
							if(shaped) shp = rep(coef(fit)["shape"], out.sample[i]) else shp = rep(0, out.sample[i])
							if(skewed) skw = rep(coef(fit)["skew"], out.sample[i]) else skw = rep(0, out.sample[i])
							if(ghyp) shpgig = rep(coef(fit)["ghlambda"], out.sample[i]) else shpgig = rep(0, out.sample[i])
							rlz = tail(data[rollind[[i]]], out.sample[i])
							realv = tail(as.numeric(realizedVol[rollind[[i]],1]), out.sample[i])
							# use xts for indexing the forecasts
							y = as.data.frame(cbind(ret, sig, skw, shp, shpgig, rlz, realv, realized))
							rownames(y) = tail(as.character(index[rollind[[i]]]), out.sample[i])
							colnames(y) = c("Mu", "Sigma", "Skew", "Shape", "Shape(GIG)", "Realized", "RVol", "RVolForecast")
							if(keep.coef) cf = fit@fit$robust.matcoef else cf = NA
							ans = list(y = y, cf = cf, converge = TRUE, lik = likelihood(fit), plik = sum(-fit@fit$partial.log.likelihoods))
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
			model$realizedVol = realizedVol
			model$n.sim = n.sim
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
			model$rollind = rollind
			model$out.sample = out.sample
			forecast = forecast
			toc = Sys.time()-tic
			model$elapsed = toc
			ans = new("uGARCHroll",
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
			model$realizedVol = realizedVol
			model$n.sim = n.sim
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
			model$rollind = rollind
			model$out.sample = out.sample
			model$coef = cf
			model$loglik = sapply(tmp, function(x) x$lik)
			model$partial.loglik = sapply(tmp, function(x) x$plik)
			forecast = list(VaR = VaR.matrix, density = forc)
		}
		toc = Sys.time()-tic
		model$elapsed = toc
		ans = new("uGARCHroll",
				model = model,
				forecast = forecast)
	} else{
		# do nothing...all converged
		ans = object
	}
	return( ans )
}