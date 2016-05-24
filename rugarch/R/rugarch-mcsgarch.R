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
# t = daily index
# i = intraday bin (e.g 390 1 min indices)
# R[t,i] = ... + sqrt(RV[t] * S[t,i] * q[t,i]) * z[t,i]
#---------------------------------------------------------------------------------
# SECTION msGARCH fit
#---------------------------------------------------------------------------------
# [mu ar ma arfima im mxreg omega alpha beta gamma gamma11 gamma21 delta lambda vxreg skew shape dlamda aux aux aux aux]

.mcsgarchfit = function(spec, data, out.sample = 0, solver = "solnp", solver.control = list(), 
		fit.control = list(stationarity = 1, fixed.se = 0, scale = 0, rec.init = 'all'), 
		numderiv.control = list(grad.eps=1e-4, grad.d=0.0001, grad.zero.tol=sqrt(.Machine$double.eps/7e-7),
				hess.eps=1e-4, hess.d=0.1, hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2), DailyVar)
{
	tic = Sys.time()
	if(is.null(solver.control$trace)) trace = 0 else trace = solver.control$trace
	
	if(is.null(fit.control$stationarity)) fit.control$stationarity = TRUE
	if(is.null(fit.control$fixed.se)) fit.control$fixed.se = FALSE
	if(is.null(fit.control$scale)){
		fit.control$scale = FALSE
	} else{
		if(fit.control$scale) stop("\nscaling not valid for mcsGARCH model.")
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
	if(is.null(DailyVar)){
		stop("\nugarchfit-->error: you must supply the daily forecast variance (DailyVar) for the msGARCH model\n")
	} else{
		if(!is(DailyVar, "xts")) stop("\nugarchfit-->error: DailyVar must be an xts object\n")
		DailyVarIndex = format(index(DailyVar), format="%Y-%m-%d")
	}
	# we are not going to extract the data just yet
	UIndex = unique(format(index(data), format="%Y-%m-%d"))
	DIndex = format(index(data), format="%Y-%m-%d")
	RIndex = index(data)
	M = length(UIndex)
	matchD = all.equal(UIndex, DailyVarIndex)
	if(!matchD) stop("\nugarchfit-->error: DailyVar dates do not match the data dates (unique days).\n")
	Tb = lapply(1:M, function(i) RIndex[which(DIndex==UIndex[i])])
	DVar = lapply(1:M, function(i) rep(DailyVar[i], length(Tb[[i]])))
	# can't unlist a POSIXct object...need to manually concatentate (can't use 'c' with recursive option either)
	dTT = Tb[[1]]
	if(length(Tb)>1) for(i in 2:length(Tb)) dTT = c(dTT, Tb[[i]])
	DVar = xts(as.numeric(unlist(DVar)), dTT)
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
	# for the mcsGARCH model we need to work with xts
	data = xts(data, index)
	itime = .unique_intraday(data)
	DV = DVar[1:(n-n.start)]
	idx1 = .unique_time(data)
	idx2 = .stime(data)
	# unique intraday time is kept to be used for forecasting and simulation.

	# arglist replaces the use of custom environments (resolves the problem of 
	# non-shared environment in parallel estimation, particularly windows)
	garchenv = new.env(hash = TRUE)
	arglist = list()
	###################
	# needed for the msGARCH model
	arglist$idx1 = idx1
	arglist$idx2 = idx2
	arglist$DV = as.numeric(DV)
	###################
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
				return(ugarchfilter(data = xts(origdata, origindex), spec = spec, out.sample = out.sample, DailyVar = DailyVar))
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
		cb = .garchconbounds()
		Ifn = .sgarchcon
		ILB = cb$LB
		IUB = cb$UB
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
		solution = .garchsolver(solver, pars = ipars[estidx, 1], fun = .mcsgarchLLH, 
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
				tmpx = .mcsgarchLLH(sol$par, arglist)
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
		fit = .makefitmodel(garchmodel = "mcsGARCH", f = .mcsgarchLLH, T = T, m = m, 
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
		########################################################################
		# V (daily forecast Variance) and S (diurnal Variance) are aligned to the
		# original time index
		# model$idx1 = .unique_time(xts(origdata, origindex))
		# model$idx2 = .stime(xts(origdata, origindex))
		# itime == unique intraday intervals on which diurnal vol is based and for
		# use with forecasting and simulation
		model$dtime = itime
		idx1 = .unique_time(xts(origdata[1:T], origindex[1:T]))
		idx2 = .stime(xts(origdata[1:T], origindex[1:T]))
		model$dvalues = .diurnal_series(fit$residuals, as.numeric(DVar)[1:T], idx1)
		# DailyVar will be of length T+out.sample (since it does not depend on any endogenous
		# variables and we can safely store the full values)
		model$DailyVar = DVar
		# DiurnalVar will be of length T (because it depends on residuals)		
		model$DiurnalVar = xts(.diurnal_series_aligned(xts(fit$residuals, index), DVar[1:T], idx1, idx2), origindex[1:T])
		# adjust sigma (q = component volatility i.e. on deasonalized data)
		fit$q = fit$sigma
		fit$sigma = fit$q * sqrt(DVar[1:T]*model$DiurnalVar[1:T])
		fit$z = fit$residuals/fit$sigma
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
# SECTION sGARCH LLH
#---------------------------------------------------------------------------------
.mcsgarchLLH = function(pars, arglist)
{
	# prepare inputs
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
	rx = .arfimaxfilter(modelinc[1:21], pars = ipars[,1], idx = idx, mexdata = mexdata, h = hm, data = data, N = N)
	res = rx$res
	zrf = rx$zrf
	res[is.na(res) | !is.finite(res) | is.nan(res)] = 0
	# 1. Create the diurnal series (bins)
	# 2. Adjust res by denominator
	dseries = .diurnal_series_aligned(res, arglist$DV, arglist$idx1, arglist$idx2)
	eres = res/sqrt(arglist$DV*dseries)
	# sgarch persistence value
	kappa = 1
	persist = (sum(ipars[idx["alpha",1]:idx["alpha",2],1]) + sum(ipars[idx["beta",1]:idx["beta",2],1]))
	# recursion initialization
	mvar = ifelse(recinit$type==1, mean(eres[1:recinit$n]*eres[1:recinit$n]), backcastv(eres, T, recinit$n))
	if(modelinc[15]>0) {
		mv = sum(apply(matrix(vexdata, ncol = modelinc[15]), 2, "mean")*ipars[idx["vxreg",1]:idx["vxreg",2],1])
	} else{
		mv = 0
	}
	hEst = mvar
	# variance targeting
	if(modelinc[7]>0){
		ipars[idx["omega",1],1] = max(eps, ipars[idx["omega",1],1]) 
	} else{
		mvar2 = ifelse(!is.na(modelinc[22]), modelinc[22]/dscale, mvar)
		ipars[idx["omega",1],1] = mvar2 * (1 - persist) - mv
		assign("omega", ipars[idx["omega",1],1], garchenv)
	}
	# E[eres^2]=1
	if(is.na(hEst) | !is.finite(hEst) | is.nan(hEst)) hEst = (1 - persist)
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
	if(modelinc[15]>0) vexdata = as.double(as.vector(vexdata)) else vexdata = double(1)
	# modelinc (1:21) since 22 is either NA or numeric and no longer needed (paased to ipars if used)
	ans = try( .C("mcsgarchfilterC", model = as.integer(modelinc[1:21]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					hEst = as.double(hEst), res = as.double(eres), 
					e = as.double(eres*eres), s = as.double(dseries), 
					v = as.double(arglist$DV), vexdata = vexdata, m = as.integer(m), 
					T = as.integer(T), h = double(T), z = double(T), llh = double(1), 
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
					LHT = LHT, persistence = persist, dseries = dseries))
	return(ans)
}

#---------------------------------------------------------------------------------
# SECTION sGARCH filter
#---------------------------------------------------------------------------------
.mcsgarchfilter = function(spec, data, out.sample = 0, n.old = NULL, rec.init = 'all', DailyVar)
{
	# n.old is optional and indicates the length of the original dataseries (in
	# cases when this represents a dataseries augmented by newer data). The reason
	# for using this is so that the old and new datasets agree since the original
	# recursion uses the sum of the residuals to start the recursion and therefore
	# is influenced by new data. For a small augmentation the values converge after
	# x periods, but it is sometimes preferable to have this option so that there is
	# no forward looking information contaminating the study.
	if(missing(DailyVar)){
		stop("\nugarchfilter-->error: you must supply the daily forecast variance (DailyVar) for the msGARCH model\n")
	} else{
		if(!is(DailyVar, "xts")) stop("\nugarchfilter-->error: DailyVar must be an xts object\n")
		DailyVarIndex = format(index(DailyVar), format="%Y-%m-%d")
	}
	# we are not going to extract the data just yet
	UIndex = unique(format(index(data), format="%Y-%m-%d"))
	DIndex = format(index(data), format="%Y-%m-%d")
	RIndex = index(data)
	M = length(UIndex)
	matchD = all.equal(UIndex, DailyVarIndex)
	if(!matchD) stop("\nugarchfilter-->error: DailyVar dates do not match the data dates (unique days).\n")
	Tb = lapply(1:M, function(i) RIndex[which(DIndex==UIndex[i])])
	DVar = lapply(1:M, function(i) rep(DailyVar[i], length(Tb[[i]])))
	# can't unlist a POSIXct object...need to manually concatentate (can't use 'c' with recursive option either)
	dTT = Tb[[1]]
	if(length(Tb)>1) for(i in 2:length(Tb)) dTT = c(dTT, Tb[[i]])
	DVar = xts(as.numeric(unlist(DVar)), dTT)	
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
	# for the msGARCH model we need to work with xts
	data = xts(data, index)
	itime = .unique_intraday(data)
	DV = DVar[1:T]
	idx1 = .unique_time(data[1:Nx])
	idx2 = .stime(data)
	# unique intraday time is kept to be used for forecasting and simulation.
	
	
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
	kappa = 1
	persist = (sum(ipars[idx["alpha",1]:idx["alpha",2],1]) + sum(ipars[idx["beta",1]:idx["beta",2],1]))
	rx = .arfimaxfilter(modelinc[1:21], ipars[,1], idx, mexdata = mexdata, h = 0, data = as.numeric(data), N = N)
	res = rx$res
	zrf = rx$zrf
	
	# 1. Create the diurnal series (bins)
	# 2. Adjust res by denominator
	dseries = .diurnal_series_aligned(res, as.numeric(DV), idx1, idx2)
	eres = res/sqrt(as.numeric(DV)*dseries)
	if(!is.null(n.old)){
		rx2 = .arfimaxfilter(modelinc[1:21], ipars[,1], idx, mexdata = mexdata[1:Nx, , drop = FALSE], h = 0, data = origdata[1:Nx], N = c(m, Nx))
		res2 = rx2$res
		xdseries = .diurnal_series_aligned(res2, as.numeric(DV), idx1, idx2)
		xeres = res2/sqrt(as.numeric(DV[1:Nx])*xdseries[1:Nx])
		mvar = ifelse(recinit$type==1, mean(xeres[1:recinit$n]*xeres[1:recinit$n]), backcastv(xeres, Nx, recinit$n))
	} else{
		mvar = ifelse(recinit$type==1, mean(eres[1:recinit$n]*eres[1:recinit$n]), backcastv(eres, T, recinit$n))
	}
	
	hEst = mvar
	if(modelinc[15]>0) {
		mv = sum(apply(matrix(vexdata, ncol = modelinc[15]), 2, "mean")*ipars[idx["vxreg",1]:idx["vxreg",2],1])
	} else{
		mv = 0
	}
	if(modelinc[7]>0){
		ipars[idx["omega",1],1] = max(eps, ipars[idx["omega",1],1]) 
	} else{
		mvar2 = ifelse(!is.na(modelinc[22]), modelinc[22], mvar)
		ipars[idx["omega",1],1] = mvar2 * (1 - persist) - mv		
	}
	
	if(modelinc[6]>0) mexdata = as.double(as.vector(mexdata)) else mexdata = double(1)
	if(modelinc[15]>0) vexdata = as.double(as.vector(vexdata)) else vexdata = double(1)
	
	ans = try( .C("mcsgarchfilterC", model = as.integer(modelinc[1:21]), 
					pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
					hEst = as.double(hEst), res = as.double(eres), 
					e = as.double(eres*eres), s = as.double(dseries), 
					v = as.double(DV), vexdata = vexdata, m = as.integer(m), 
					T = as.integer(T), h = double(T), z = double(T), llh = double(1), 
					LHT = double(T), 
					PACKAGE = "rugarch"), silent = TRUE )
	
	filter = list()
	filter$residuals = res
	filter$LLH = -ans$llh
	filter$log.likelihoods = ans$LHT
	filter$persistence = persist
	filter$distribution = distribution
	filter$ipars = ipars
	model$modeldata$data = origdata
	model$modeldata$index = origindex
	model$modeldata$period = period
	model$n.start = out.sample
	# V (daily forecast Variance) and S (diurnal Variance) are aligned to the
	# original time index
	# set idx1 to only the n.old value.
	model$idx1 = .unique_time(xts(origdata[1:Nx], origindex[1:Nx]))
	model$idx2 = .stime(xts(origdata, origindex))
	# itime == unique intraday intervals on which diurnal vol is based and for
	# use with forecasting and simulation
	model$itime = itime
	model$DailyVar = DVar
	model$DiurnalVar = xts(.diurnal_series_aligned(res, as.numeric(DVar), model$idx1, model$idx2), origindex)
	filter$q = sqrt(ans$h)
	filter$sigma = filter$q * sqrt(DVar[1:T]*model$DiurnalVar[1:T])
	filter$z = filter$residuals/filter$sigma
	
	sol = new("uGARCHfilter",
			filter = filter,
			model = model)
	return(sol)
}

#---------------------------------------------------------------------------------
# SECTION mcsGARCH forecast
#---------------------------------------------------------------------------------
# Recipe for a GARCH forecast:
# 1. split forecast into arma and garch forecasts.
# 2. take into account garch-in-mean, arfima, exogenous regressors
# 3. check kappa and persistence
# 4. allow for rolling forecast based on filtering out.sample data from fit stage

# n.roll signifies whether to filter/roll the sigma else will use unconditional
# i.e. n.ahead=10 with n.roll=1 means that the first forecast is based on the
# previous value whereas all subsequent forecasts are based on the unconditional
# expectation formula
# with n.roll = n means that while n.ahead < n.roll we filter the data using the coef
# of the garch model to obtain filtered sigma estimates which are used to roll
# the forecast. This requires that during the fit process the out.sample option
# was used and that (n.roll) < out.sample (otherwise we revert to the unconditional
# expectation formula for the long run sigma.

## mcs GARCH has a lot of extra processing because of the time/date issues
# involved
# All pre-processed values must conform to the n.ahead by n.roll matrix, 
# but unlike other GARCH models, we now need to know what n.ahead is, in terms
# of time/date
.mcsgarchforecast = function(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0, 
		external.forecasts = list(mregfor = NULL, vregfor = NULL), DailyVar, ...)
{
	fit    = fitORspec
	data   = fit@model$modeldata$data
	Nor    = length(as.numeric(data))
	index  = fit@model$modeldata$index
	period = fit@model$modeldata$period
	# check DailyVar forecast provided with what is available from the fitted object
	# (if out.sample was used it may not be needed).
	DiurnalVar = fit@model$DiurnalVar
	inDailyVar = fit@model$DailyVar
	lastDate = format(tail(index(inDailyVar), 1), "%Y-%m-%d")
	dtime = fit@model$dtime
	dvalues = fit@model$dvalues
	
	# prepare the diurnal, daily vols
	if(fit@model$n.start>0){
		if( (n.ahead+n.roll)<=fit@model$n.start ){
			# we don't require external DailyVaR
			# completely in-the-sample
			DVar = fit@model$DailyVar
			DiurnalVar = NULL
		} else{
			# mixed in and out sample
			needT = (n.ahead+n.roll) - fit@model$n.start
			outD = ftseq(T0 = as.POSIXct(tail(index, 1)), 
					length.out = needT, by = period, 
					interval = fit@model$dtime, 
					exclude.weekends = TRUE)
			Dmatch = match(format(outD, "%H:%M:%S"), dtime)
			D2 = xts(dvalues[Dmatch], outD)
			D1 = xts(dvalues[match(format(index, "%H:%M:%S"), dtime)], index)
			DiurnalVar = c(D1, D2)
			DVar = .intraday2daily(fit@model$DailyVar)
			# Check to see whether we need DailyVar Forecast
			U = unique(format(index(DiurnalVar), "%Y-%m-%d"))
			Y = unique(c(format(index(DVar), "%Y-%m-%d"), if(!missing(DailyVar)) format(index(DailyVar),  "%Y-%m-%d") else NULL))
			DV = c(as.numeric(DVar), if(!missing(DailyVar)) as.numeric(DailyVar) else NULL)
			test = match(U, Y)
			if(any(is.na(test))){
				idx = which(is.na(test))
				stop(paste(c("DailyVar requires forecasts for: ", U[idx],"...resubmit."), sep="",collapse=" "))
			} else{
				# create the DailyVar
				M = length(Y)
				RIndex = index(DiurnalVar)
				UIndex = unique(format(index(DiurnalVar), format="%Y-%m-%d"))
				DIndex = format(index(DiurnalVar), format="%Y-%m-%d")
				Tb = lapply(1:M, function(i) RIndex[which(DIndex==UIndex[i])])
				DVar = lapply(1:M, function(i) rep(DV[i], length(Tb[[i]])))
				# can't unlist a POSIXct object...need to manually concatentate (can't use 'c' with recursive option either)
				dTT = Tb[[1]]
				if(length(Tb)>1) for(i in 2:length(Tb)) dTT = c(dTT, Tb[[i]])
				DVar = xts(as.numeric(unlist(DVar)), dTT)
			}
		}
	} else{
		# completely out of the sample
		outD = ftseq(T0 = as.POSIXct(tail(index, 1)), 
				length.out = n.ahead+n.roll, by = period, 
				interval = fit@model$dtime, 
				exclude.weekends = TRUE)
		Dmatch = match(format(outD, "%H:%M:%S"), dtime)
		D2 = xts(dvalues[Dmatch], outD)
		D1 = xts(dvalues[match(format(index, "%H:%M:%S"), dtime)], index)
		DiurnalVar = c(D1, D2)
		DVar = .intraday2daily(fit@model$DailyVar)
		# Check to see whether we need DailyVar Forecast
		U = unique(format(index(DiurnalVar), "%Y-%m-%d"))
		Y = unique(c(format(index(DVar), "%Y-%m-%d"), if(!missing(DailyVar)) format(index(DailyVar), "%Y-%m-%d") else NULL))
		DV = c(as.numeric(DVar), if(!missing(DailyVar)) as.numeric(DailyVar) else NULL)
		test = match(U, Y)
		if(any(is.na(test))){
			idx = which(is.na(test))
			stop(paste(c("DailyVar requires forecasts for: ", U[idx],"...resubmit."), sep="",collapse=" "))
		} else{
			# create the DailyVar
			M = length(Y)
			RIndex = index(DiurnalVar)
			UIndex = unique(format(index(DiurnalVar), format="%Y-%m-%d"))
			DIndex = format(index(DiurnalVar), format="%Y-%m-%d")
			Tb = lapply(1:M, function(i) RIndex[which(DIndex==UIndex[i])])
			DVar = lapply(1:M, function(i) rep(DV[i], length(Tb[[i]])))
			# can't unlist a POSIXct object...need to manually concatentate (can't use 'c' with recursive option either)
			dTT = Tb[[1]]
			if(length(Tb)>1) for(i in 2:length(Tb)) dTT = c(dTT, Tb[[i]])
			DVar = xts(as.numeric(unlist(DVar)), dTT)
		}
	}
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
	fspec = ugarchspec(variance.model = list(model = "mcsGARCH", 
					garchOrder = c(modelinc[8], modelinc[9]), submodel = NULL, 
					external.regressors = vxf[1:(N + fcreq), , drop = FALSE]), 
			mean.model = list(armaOrder = c(modelinc[2], modelinc[3]),
					include.mean = modelinc[1], 
					archm = ifelse(modelinc[5]>0,TRUE,FALSE), archpow = modelinc[5], arfima = modelinc[4], 
					external.regressors = mxf[1:(N + fcreq), , drop = FALSE], archex = modelinc[20]), 
			distribution.model = model$modeldesc$distribution, fixed.pars = as.list(pars))
	tmp =  xts(data[1:(N + fcreq)], index[1:(N + fcreq)])
	DailyV = .intraday2daily(DVar[1:(N + fcreq)])
	flt = .mcsgarchfilter(data = tmp, spec = fspec, n.old = N, DailyVar = DailyV)
	# For the case that the forecast is completey in the out.sample period, we
	# need to filtered DiurnalVar which extends to the whole period since the
	# one from the uGARCHfit object only extends to T (Totalobs-out.sample).
	if(is.null(DiurnalVar)) DiurnalVar = flt@model$DiurnalVar
	sigmafilter = as.numeric(sigma(flt))
	qfilter = flt@filter$q
	resfilter = as.numeric(residuals(flt))
	zfilter = as.numeric(flt@filter$z)
	eresfilter = resfilter/sqrt(as.numeric(DVar[1:(N + fcreq)])*as.numeric(DiurnalVar[1:(N + fcreq)]))
	# forecast GARCH process
	qFor = seriesFor = sigmaFor = matrix(NA, ncol = n.roll+1, nrow = n.ahead)
	colnames(qFor) = colnames(seriesFor) = colnames(sigmaFor) = as.character(index[N:(N+n.roll)])
	rownames(qFor) = rownames(seriesFor) = rownames(sigmaFor) = paste("T+", 1:n.ahead, sep="")
	for(i in 1:(n.roll+1)){
		np = N + i - 1
		if(modelinc[1] > 0){
			mu = rep(ipars[idx["mu",1]:idx["mu",2], 1], N+i+n.ahead-1)
		} else{
			mu = rep(0, N+i+n.ahead-1)
		}
		omega = rep(ipars[idx["omega",1]:idx["omega",2], 1], N+i+n.ahead-1)
		h = c(sigmafilter[1:(N+i-1)], rep(0, n.ahead))
		q = c(qfilter[1:(N+i-1)], rep(0, n.ahead))
		res = c(resfilter[1:(N+i-1)], rep(0, n.ahead))
		eres = c(eresfilter[1:(N+i-1)], rep(0, n.ahead))	
		x = c(data[1:(N+i-1)], rep(0, n.ahead))
		z = c(zfilter[1:(N+i-1)], rep(0, n.ahead))
		# forecast of externals is provided outside the system
		mxfi = mxf[1:(N+i-1+n.ahead), , drop = FALSE]
		vxfi = vxf[1:(N+i-1+n.ahead), , drop = FALSE]
		ans = .nmcsgarchforecast(ipars, modelinc, idx, mu, omega, mxfi, vxfi, q, h, eres, res, z, DVar, DiurnalVar, 
				data = x, N = np, n.ahead)
		sigmaFor[,i] = ans$h
		qFor[,i] = ans$q
		seriesFor[,i] = ans$x
	}
	
	fcst = list()
	fcst$n.ahead = n.ahead
	fcst$N = N+ns
	fcst$n.start = ns
	fcst$n.roll = n.roll
	fcst$sigmaFor = sigmaFor
	fcst$seriesFor = seriesFor
	fcst$qFor = qFor
	model$modeldata$sigma = flt@filter$sigma
	model$modeldata$residuals = flt@filter$residuals
	ans = new("uGARCHforecast",
			forecast = fcst,
			model = model)
	return(ans)
}

.nmcsgarchforecast = function(ipars, modelinc, idx, mu, omega, mxfi, vxfi, q, h, 
		eres, res, z, DVar, DiurnalVar, data, N, n.ahead)
{
	if(modelinc[15]>0){
		omega = omega + vxfi%*%t(matrix(ipars[idx["vxreg",1]:idx["vxreg",2],1], ncol = modelinc[15]))
	}
	for(i in 1:n.ahead){
		if(modelinc[9]>0){
			q[N+i] = omega[N+i] + sum(ipars[idx["beta",1]:idx["beta",2],1]*q[N+i-(1:modelinc[9])]^2)
		} else{
			q[N+i] = omega[N+i]
		}
		if(modelinc[8]>0){
			for (j in 1:modelinc[8]){
				if (i-j > 0){
					# Equation 8 of Engle/Sokalska
					s = q[N + i - j]^2
				} else{ 
					s = eres[N + i - j]^2
				}
				q[N+i] = q[N+i] + ipars[idx["alpha",1]+j-1,1] * s
			}
		}
		q[N+i] = sqrt(q[N+i])
		h[N+i] = q[N+i]*sqrt(as.numeric(DVar)[N+i]*as.numeric(DiurnalVar)[N+i])
	}
	
	if(modelinc[4]>0){
		xres = arfimaf(ipars, modelinc[1:21], idx, mu, mxfi, h, res, z, data, N, n.ahead)
	} else{
		xres = armaf(ipars, modelinc[1:21], idx, mu, mxfi, h, res, z, data, N, n.ahead)
	}
	return(list(h = h[(N+1):(N+n.ahead)], q = q[(N+1):(N+n.ahead)], x = xres[(N+1):(N+n.ahead)]))
}

#---------------------------------------------------------------------------------
# 2nd dispatch method for forecast
.mcsgarchforecast2 = function(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0, 
		external.forecasts = list(mregfor = NULL, vregfor = NULL), ...)
{
	stop("\nugarchforecast with specification object not available for mcsGARCH model")
}
#---------------------------------------------------------------------------------
# SECTION mcsGARCH simulate
#---------------------------------------------------------------------------------
# DailyVar is the simulated daily variance as an n by m xts object, where n by m
# contains the days within n.sim by m.sim.
.mcsgarchsim = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, startMethod = 
				c("unconditional","sample"), presigma = NA, prereturns = NA, 
		preresiduals = NA, rseed = NA, custom.dist = list(name = NA, distfit = NA), 
		mexsimdata = NULL, vexsimdata = NULL, DailyVar, ...)
{
	if(fit@model$modelinc[4]>0){
		if(n.start<fit@model$modelinc[3]){
			warning("\nugarchsim-->warning: n.start>=MA order for arfima model...automatically setting.")
			n.start = fit@model$modelinc[3]
		}
	}
	# 1. Create Diurnal Var (n.sim)
	T = fit@model$modeldata$T
	
	T0 = fit@model$modeldata$index[T]
	dtime = fit@model$dtime
	dvalues = fit@model$dvalues
	D = ftseq(T0, length.out = n.sim+n.start, by = fit@model$modeldata$period, interval = dtime)
	Dmatch = match(format(D, "%H:%M:%S"), dtime)
	DiurnalVar = xts(dvalues[Dmatch], D)
	# 2. Transform DailyVar into intraday (n.sim by m.sim)
	if(missing(DailyVar)) stop("\nDailyVar cannot be missing for the mcsGARCH model.")
	if(!is(DailyVar, "xts")) stop("\nDailyVar must be an xts object of daily variance simulations for the n.sim period")
	
	DailyVarOld = .intraday2daily(fit@model$DailyVar)
	Unique1 =  unique(format(index(DiurnalVar), "%Y-%m-%d"))
	Unique2 =  unique(format(index(DailyVarOld), "%Y-%m-%d"))
	# find the required dates
	ReqDates = setdiff(Unique1, Unique2)
	Unique3  = format(index(DailyVar), "%Y-%m-%d")
	if(any(is.na(match(ReqDates, Unique3)))){
		stop(paste(c("The required dates for the DailyVar are:", ReqDates), sep="", collapse=" "))
	}
	if(NCOL(DailyVar)!=m.sim){
		if(NCOL(DailyVar)==1){
			warning("\nDailyVar not equal to m.sim...will replicate to use the same for all independent simulations (m.sim)")
			DVar = matrix(NA, ncol = m.sim, nrow = n.sim+n.start)
			# need to align old Daily Var with new Daily Var
			DailyVarOld = .intraday2daily(fit@model$DailyVar)
			DailyVar = c(DailyVarOld, DailyVar)
			Dx = .daily2intraday(DiurnalVar, DailyVar)
			Dy = tail(coredata(Dx), n.sim+n.start)
			for(i in 1:m.sim) DVar[,i] = Dy
		} else{
			stop("\nNCOL(DailyVar) is greater than 1 and less than m.sim...resubmit something which makes more sense.")
		}
	} else{
		DVar = matrix(NA, ncol = m.sim, nrow = n.sim+n.start)
		DailyVarOld = .intraday2daily(fit@model$DailyVar)
		for(i in 1:m.sim){
			DailyVarx = c(DailyVarOld, DailyVar[,i])
			Dx = .daily2intraday(DiurnalVar, DailyVarx)
			DVar[,i] = tail(coredata(Dx), n.sim+n.start)
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
	# default to sample based method for the startMethod
	if(length(startMethod)==2) startMethod = startMethod[2]
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
	if(startMethod == "unconditional"){
		z = rbind(matrix(0, nrow = m, ncol = m.sim), z)
	} else{
		z = rbind(matrix(tail(fit@fit$z, m), nrow = m, ncol = m.sim), z) 
	}
	if(!is.na(preresiduals[1])){
		preresiduals = as.vector(preresiduals)
		if(length(preresiduals)<m) stop(paste("\nugarchsim-->error: preresiduals must be of length ", m, sep=""))
		preres = preresiduals
	} else{
		preres = tail(fit@fit$residuals, m)
	}
	# create the presample information
	if(!is.na(prereturns[1])){
		prereturns = as.vector(prereturns)
		if(length(prereturns)<m) stop(paste("\nugarchsim-->error: prereturns must be of length ", m, sep=""))
	} else{
		if(startMethod[1] == "unconditional"){
			prereturns = as.numeric(rep(uncmean(fit), m))
		}
		else{
			prereturns = tail(data, m)
		}
	}
	
	if(!is.null(list(...)$preq)){
		preq = tail(list(...)$preq^2, m)
	} else{
		preq = tail(as.numeric(fit@fit$q)^2, m)
	}
	eres = fit@fit$residuals/sqrt(as.numeric(fit@model$DiurnalVar[1:N])*as.numeric(fit@model$DailyVar[1:N]))
	eres = c(tail(eres, m), rep(0, n))
	# input vectors/matrices
	h = c(preq, rep(0, n))
	x = c(prereturns, rep(0, n))
	constm = matrix(ipars[idx["mu",1]:idx["mu",2], 1], ncol = m.sim, nrow = n + m)
	
	# outpus matrices
	qSim = sigmaSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	seriesSim = matrix(0, ncol = m.sim, nrow = n.sim)
	residSim =  matrix(0, ncol = m.sim, nrow = n.sim)
	z[is.na(z) | is.nan(z) | !is.finite(z)] = 0
	# eres

	for(i in 1:m.sim){
		ans1 = try(.C("mcsgarchsimC", model = as.integer(modelinc[1:21]), pars = as.double(ipars[,1]), idx = as.integer(idx[,1]-1), 
						h = as.double(h), z = as.double(z[,i]), eres = as.double(eres), e = as.double(eres*eres),
						vexdata = as.double(vexsim[[i]]), T = as.integer(n+m), m = as.integer(m), PACKAGE = "rugarch"), silent = TRUE)
		if(inherits(ans1, "try-error")) stop("\nugarchsim-->error: error in calling C function....\n")
		qSim[,i] = ans1$h[(n.start + m + 1):(n+m)]^(1/2)
		sigmaSim[,i] = qSim[,i]*sqrt(DVar[,i]*as.numeric(DiurnalVar))
		res = c(preres, ans1$eres[-c(1:m)]*sqrt(DVar[,i]*as.numeric(DiurnalVar)))
		residSim[,i] = res[(n.start + m + 1):(n+m)]
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
		# not allowed in model
		#if(modelinc[5]>0) constm[,i] = constm[,i] + ipars[idx["archm",1]:idx["archm",2], 1]*(sqrt(ans1$h)^modelinc[5])
		if(modelinc[4]>0){
			fres = c(res[(m+1):(n+m)], if(modelinc[3]>0) rep(0, modelinc[3]) else NULL)
			ans2 = .arfimaxsim(modelinc[1:21], ipars, idx, constm[1:n, i], fres, T = n)
			seriesSim[,i] = head(ans2$series, n.sim)
		} else{
			ans2 = .armaxsim(modelinc[1:21], ipars, idx, constm[,i],  x, res, T = n + m, m)
			seriesSim[,i] = ans2$x[(n.start + m + 1):(n+m)]
		}
	}
	sim = list(qSim = qSim, sigmaSim = sigmaSim, seriesSim = seriesSim, residSim = residSim)
	sim$DiurnalVar = DiurnalVar
	sim$DailyVar = DVar
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
.mcsgarchpath = function(spec, n.sim = 1000, n.start = 0, m.sim = 1,
		presigma = NA, prereturns = NA, preresiduals = NA, rseed = NA, 
		custom.dist = list(name = NA, distfit = NA), mexsimdata = NULL, 
		vexsimdata = NULL, ...)
{
	stop("\nugarchpath method not available for mcsGARCH model")
}
################################################################################
# special purpose rolling method for the mcsGARCH model
################################################################################

.rollfdensity.mcs = function(spec, data, n.ahead = 1, forecast.length = 500, 
		n.start = NULL, refit.every = 25, refit.window = c("recursive", "moving"), 
		window.size = NULL, solver = "hybrid", fit.control = list(), solver.control = list(),
		calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), cluster = NULL,
		keep.coef = TRUE, DailyVar, ...)
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
	###########################################################################
	# Check DailyVar
	if(is.null(DailyVar)){
		stop("\nugarchroll-->error: you must supply the daily forecast variance (DailyVar) for the msGARCH model\n")
	} else{
		if(!is(DailyVar, "xts")) stop("\nugarchroll-->error: DailyVar must be an xts object\n")
		DailyVarIndex = format(index(DailyVar), format="%Y-%m-%d")
	}
	# we are not going to extract the data just yet
	UIndex = unique(format(index(data), format="%Y-%m-%d"))
	DIndex = format(index(data), format="%Y-%m-%d")
	RIndex = index(data)
	M = length(UIndex)
	matchD = all.equal(UIndex, DailyVarIndex)
	if(!matchD) stop("\nugarchroll-->error: DailyVar dates do not match the data dates (unique days).\n")
	
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
	DailyV = lapply(rollind, FUN = function(x){
				dindex = unique(format(index[x], "%Y-%m-%d"))
				DailyVar[dindex]
			})
	# distribution
	distribution = spec@model$modeldesc$distribution
	if(any(distribution==c("snorm","sstd","sged","jsu","nig","ghyp","ghst"))) skewed = TRUE else skewed = FALSE
	if(any(distribution==c("std","sstd","ged","sged","jsu","nig","ghyp","ghst"))) shaped = TRUE else shaped = FALSE
	if(any(distribution==c("ghyp"))) ghyp = TRUE else ghyp = FALSE
	
	if( !is.null(cluster) ){
		clusterEvalQ(cl = cluster, library(rugarch))
		clusterExport(cluster, c("data", "index", "s","refit.every", 
						"keep.coef", "shaped", "skewed", "ghyp", 
						"rollind", "spec", "out.sample", "mex", "vex",
						"solver", "solver.control", "fit.control", "DailyV"), envir = environment())
		if(mex) clusterExport(cluster, c("mexdata"), envir = environment())
		if(vex) clusterExport(cluster, c("vexdata"), envir = environment())
		tmp = parLapply(cl = cluster, 1:m, fun = function(i){
					if(mex) spec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
					if(vex) spec@model$modeldata$vexdata = vexdata[rollind[[i]],,drop=FALSE]
					fit = try(ugarchfit(spec, xts::xts(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i], 
									solver = solver, solver.control = solver.control, 
									fit.control = fit.control, DailyVar = DailyV[[i]]), silent=TRUE)
					# 3 cases: General Error, Failure to Converge, Failure to invert Hessian (bad solution)
					if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
						ans = list(y = NA, cf = NA, q = NA, DiurnalVar = NA, DailyVar = NA, converge = FALSE)
					} else{
						if(mex) fmex = tail(mexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fmex = NULL
						if(vex) fvex = tail(vexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fvex = NULL
						# since the forecast is in the out.sample we don't need to supply DailyVar
						f = ugarchforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1, 
								external.forecasts = list(mregfor = fmex, vregfor = fvex))
						sig = as.numeric( sigma(f) )
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
						ans = list(y = y, cf = cf, q = f@forecast$qFor, DiurnalVar = f@model$DiurnalVar, 
								DailyVar = f@model$DailyVar, converge = TRUE)
					}
					return(ans)})
	} else{
		tmp = lapply(as.list(1:m), FUN = function(i){
					if(mex) spec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
					if(vex) spec@model$modeldata$vexdata = vexdata[rollind[[i]],,drop=FALSE]
					fit = try(ugarchfit(spec, xts(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i], 
									solver = solver, solver.control = solver.control, 
									fit.control = fit.control, DailyVar = DailyV[[i]]), silent=TRUE)
					if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
						ans = list(y = NA, cf = NA, q = NA, DiurnalVar = NA, DailyVar = NA, converge = FALSE)
					} else{
						if(mex) fmex = tail(mexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fmex = NULL
						if(vex) fvex = tail(vexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fvex = NULL
						f = ugarchforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1, 
								external.forecasts = list(mregfor = fmex, vregfor = fvex))
						sig = as.numeric( sigma(f) )
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
						ans = list(y = y, cf = cf, q = f@forecast$qFor, DiurnalVar = f@model$DiurnalVar, 
								DailyVar = f@model$DailyVar, converge = TRUE)
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
		model$rollind = rollind
		model$out.sample = out.sample
		model$DailyVar = DailyVar
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
		forecast = list(VaR = VaR.matrix, density = forc)
	}
	toc = Sys.time()-tic
	model$elapsed = toc
	ans = new("uGARCHroll",
			model = model,
			forecast = forecast)
	return( ans )
}

.resumeroll.mcs = function(object, solver = "hybrid", fit.control = list(), 
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
		DailyVar = object@model$DailyVar
		DailyV = lapply(rollind, FUN = function(x){
					dindex = unique(format(index[x], "%Y-%m-%d"))
					DailyVar[dindex]
				})
		# distribution
		distribution = spec@model$modeldesc$distribution
		if(any(distribution==c("snorm","sstd","sged","jsu","nig","ghyp","ghst"))) skewed = TRUE else skewed = FALSE
		if(any(distribution==c("std","sstd","ged","sged","jsu","nig","ghyp","ghst"))) shaped = TRUE else shaped = FALSE
		if(any(distribution==c("ghyp"))) ghyp = TRUE else ghyp = FALSE
		if( !is.null(cluster) ){
			clusterEvalQ(cl = cluster, library(rugarch))
			clusterExport(cluster, c("data", "index","s","refit.every",
							"keep.coef", "shaped", "skewed", "ghyp", 
							"rollind", "spec", "out.sample", "mex", "vex", 
							"noncidx", "solver", "solver.control", "fit.control", "DailyV"),
					envir = environment())
			if(mex) clusterExport(cluster, c("mexdata"), envir = environment())
			if(vex)  clusterExport(cluster, c("vexdata"), envir = environment())
			tmp = parLapply(cl = cluster, as.list(noncidx), fun = function(i){
						if(mex) spec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
						if(vex) spec@model$modeldata$vexdata = vexdata[rollind[[i]],,drop=FALSE]
						fit = try(ugarchfit(spec, xts::xts(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i], 
										solver=solver, solver.control = solver.control, 
										fit.control = fit.control, DailyVar = DailyV[[i]]), silent=TRUE)
						if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
							ans = list(y = NA, cf = NA, q = NA, DiurnalVar = NA, DailyVar = NA, converge = FALSE)
						} else{
							if(mex) fmex = tail(mexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fmex = NULL
							if(vex) fvex = tail(vexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fvex = NULL
							f = ugarchforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1, 
									external.forecasts = list(mregfor = fmex, vregfor = fvex))
							sig = as.numeric( sigma(f) )
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
							ans = list(y = y, cf = cf, q = f@forecast$qFor, DiurnalVar = f@model$DiurnalVar, 
									DailyVar = f@model$DailyVar, converge = TRUE)
						}
						return(ans)})
		} else{
			tmp = lapply(as.list(noncidx), FUN = function(i){
						if(mex) spec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
						if(vex) spec@model$modeldata$vexdata = vexdata[rollind[[i]],,drop=FALSE]
						fit = try(ugarchfit(spec, xts(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i], 
										solver=solver, solver.control = solver.control, 
										fit.control = fit.control), silent=TRUE)
						if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
							ans = list(y = NA, cf = NA, q = NA, DiurnalVar = NA, DailyVar = NA, converge = FALSE)
						} else{
							if(mex) fmex = tail(mexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fmex = NULL
							if(vex) fvex = tail(vexdata[rollind[[i]],,drop=FALSE], out.sample[i]) else fvex = NULL
							f = ugarchforecast(fit, n.ahead = 1, n.roll = out.sample[i]-1, 
									external.forecasts = list(mregfor = fmex, vregfor = fvex))
							sig = as.numeric( sigma(f) )
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
							ans = list(y = y, cf = cf, q = f@forecast$qFor, DiurnalVar = f@model$DiurnalVar, 
									DailyVar = f@model$DailyVar, converge = TRUE)
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
			model$rollind = rollind
			model$out.sample = out.sample
			model$DailyVar = DailyVar
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

################################################################################
# Special Purpose Functions for the intraday multiplicative component sGARCH model
################################################################################
.daily2intraday = function(x, v)
{
	Z = as.POSIXct(format(index(x), format="%Y-%m-%d"))
	Y = xts(rep(0, NROW(x)), index(x))
	T = index(v)
	for(i in 1:length(T)){
		idx = which(Z==T[i])
		Y[idx] = v[i]
	}
	return(Y)
}

.intraday2daily = function(v)
{
	idx = unique(format(index(v), "%Y-%m-%d"))
	idx1 = format(index(v),"%Y-%m-%d")
	idx2 = sapply(idx, function(x) min(which(idx1==x)))
	ans = xts(as.numeric(v[idx2]), as.POSIXct(idx))
	return(ans)
}
.unique_intraday = function(x)
{
	Z = format(index(x), format="%H:%M:%S")
	return(sort(unique(Z)))
}
.unique_time = function(x)
{
	Z = format(index(x), format="%H:%M:%S")
	Y = sort(unique(Z))
	idx = vector(mode="list", length=length(Y))
	for(i in 1:length(Y)){
		idx[[i]] = which(Z==Y[i])
	}
	return(idx)
}
.stime = function(x)
{
	Z = format(index(x), format="%H:%M:%S")
	Y = sort(unique(Z))
	U = as.POSIXct(format(index(x), format="%Y-%m-%d"))
	UX = unique(U)
	idx = vector(mode="list", length=length(UX))
	for(i in 1:length(UX)){
		tmp = which(U==UX[i])
		idx[[i]] = tmp[match(Y, Z[tmp])]
	}
	xidx = lapply(idx, function(i) ifelse(!is.na(i), 1, NA))
	return(xidx)
}
# idx2 = .stime(residuals)
# idx1 = .unique_time(residuals)
# v = .daily2intraday(residuals, dailyvar)
# x = residuals
.diurnal_series_aligned = function(x, v, idx1, idx2)
{
	s = sapply(idx1, function(i) mean((x[i]^2)/v[i]))
	sx = as.numeric(unlist(sapply(idx2, function(x) na.omit(s*x))))
	return(sx)
}

.diurnal_series = function(x, v, idx1)
{
	s = sapply(idx1, function(i) mean(x[i]^2/v[i]))
	return(s)
}