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

# this function implements a simulation based method for generating standard errors and the
# bootstrapped predictive density for the n.ahead forecasts of a garch model, 
# taking into account parameter uncertainty partly based on the paper by Pascual, Romo and
# Ruiz (2006) : Computational Statistics and Data Analysis 50
# "Bootstrap prediction for returns and volatilities in GARCH models"
# [PRR paper]

# the conditional bootstrap (Partial) does not consider parameter uncertainty
# the Full model does (but is more expensive since we need to simulate the parameter distribution)


.ugarchbootfit = function(fitORspec, data = NULL, method = c("Partial", "Full"),
		sampling = c("raw", "kernel", "spd"), 
		spd.options = list(upper = 0.9, lower = 0.1, type = "pwm", kernel = "normal"),
		n.ahead = 10, n.bootfit = 100, n.bootpred = 500, out.sample = 0, 
		rseed = NA, solver = "solnp", solver.control = list(), fit.control = list(), 
		external.forecasts =  list(mregfor = NULL, vregfor = NULL), 
		mexsimdata = NULL, vexsimdata = NULL, cluster = NULL, verbose = FALSE)
{
	method = tolower(method)
	ans = switch(method,
			partial = .ub1p1(fitORspec, data = data, sampling = sampling, spd.options = spd.options,
					n.ahead = n.ahead, n.bootfit = n.bootfit, n.bootpred = n.bootpred, rseed = rseed, 
					solver.control = solver.control, fit.control = fit.control, 
					external.forecasts =  external.forecasts, mexsimdata = mexsimdata, 
					vexsimdata = vexsimdata, cluster = cluster, verbose = verbose),
			full = .ub1f1(fitORspec, data = data, sampling = sampling, spd.options = spd.options, 
					n.ahead = n.ahead, n.bootfit = n.bootfit, n.bootpred = n.bootpred, rseed = rseed, 
					solver = solver, solver.control = solver.control, fit.control = fit.control, 
					external.forecasts =  external.forecasts, mexsimdata = mexsimdata, 
					vexsimdata = vexsimdata, cluster = cluster, verbose = verbose))
	return(ans)
}

.ugarchbootspec = function(fitORspec, data = NULL, method = c("Partial", "Full"), 
		sampling = c("raw", "kernel", "spd"), 
		spd.options = list(upper = 0.9, lower = 0.1, type = "pwm", kernel = "normal"),
		n.ahead = 10, n.bootfit = 100, n.bootpred = 500, out.sample = 0, rseed = NA, 
		solver = "solnp", solver.control = list(), fit.control = list(), 
		external.forecasts =  list(mregfor = NULL, vregfor = NULL), 
		mexsimdata = NULL, vexsimdata = NULL, cluster = NULL, verbose = FALSE)
{
	method = tolower(method)
	ans = switch(method,
			partial = .ub1p2(fitORspec, data = data, sampling = sampling, spd.options = spd.options, 
					n.ahead = n.ahead, n.bootfit = n.bootfit, n.bootpred = n.bootpred, 
					out.sample = out.sample, rseed = rseed, solver.control = solver.control, 
					fit.control = fit.control, external.forecasts =  external.forecasts, 
					mexsimdata = mexsimdata, vexsimdata = vexsimdata,
					cluster = cluster, verbose = verbose),
			full = .ub1f2(fitORspec, data = data, sampling = sampling, spd.options = spd.options, 
					n.ahead = n.ahead, n.bootfit = n.bootfit, n.bootpred = n.bootpred, 
					out.sample = out.sample, rseed = rseed, solver = solver, 
					solver.control = solver.control, fit.control = fit.control, 
					external.forecasts =  external.forecasts, mexsimdata = mexsimdata, 
					vexsimdata = vexsimdata, cluster = cluster, verbose = verbose))
	return(ans)
}

# method from a fit object
.ub1f1 = function(fitORspec, data = NULL, sampling = c("raw", "kernel", "spd"), 
		spd.options = list(upper = 0.9, lower = 0.1, type = "pwm", kernel = "normal"), 
		n.ahead = 10, n.bootfit = 100, n.bootpred = 500, rseed = NA, solver = "solnp", 
		solver.control = list(), fit.control = list(), 
		external.forecasts =  list(mregfor = NULL, vregfor = NULL), 
		mexsimdata = NULL, vexsimdata = NULL, cluster = NULL, verbose = FALSE)
{
	# inputs:
	# n.bootfit		: the number of simulations for which we will generate the parameter 
	#				uncertainty distribution (i.e. no of refits)
	# n.bootpred	: the number of simulations / n.ahead to use (for averaging to obtain forecast
	#				density)
	# out.sample	: data to keep for out of sample comparison purposes
	# n.ahead		: the horizon of the bootstrap density forecast
	sampling = tolower(sampling[1])
	if(is.na(match(sampling, c("raw", "kernel", "spd")))) stop("\nsampling method not recognized")
	fit = fitORspec
	model = fit@model
	vmodel = fit@model$modeldesc$vmodel
	m = model$maxOrder
	data = model$modeldata$data
	N = model$modeldata$T
	ns = model$n.start
	if(is.na(rseed[1])){
		sseed = as.integer(runif(n.bootpred + n.bootfit, 0, 65000))
	} else{
		if(length(rseed) < n.bootpred){
			stop("\nugarchboot-->error: seed length must equal n.bootpred + n.bootfit for full method\n")
		} else {
			sseed = rseed
		}
	}
	sseed1 = sseed[1:n.bootfit]
	sseed2 = sseed[(n.bootfit+1):(n.bootpred + n.bootfit)]
	# generate paths of equal length to data based on empirical re-sampling of z
	# p.2296 equation (5)
	fz = fit@fit$z
	if(sampling!="raw"){
		zfit = .bfit(fz, sampling, spd.options)
	} else{
		zfit = NULL
	}
	
	empz = matrix(0, ncol = n.bootfit, nrow = N)
	empz = apply(as.data.frame(1:n.bootfit), 1, FUN=function(i) {
				set.seed(sseed1[i]);
				.bsample(fz, N, sampling, zfit);})
	
	if(ns > 0) {
		spec = getspec(fit)
		setfixed(spec) <- as.list(coef(fit))
		realized.x = fit@model$modeldata$data[(N+1):(N+ns)]
		filtered.s = tail(ugarchfilter(data = .numeric2xts(fit@model$modeldata$data), spec = spec)@filter$sigma, fit@model$n.start)
	} else{
		realized.x = NULL
		filtered.s = NULL
	}
	
	# presigma uses the same starting values as the original fit
	# -> in paper they use alternatively the unconditional long run sigma (P.2296 paragraph 2
	# "...marginal variance..."
	paths = ugarchsim(fit, n.sim = N, m.sim = n.bootfit, 
			presigma = tail(fit@fit$sigma, m), 
			prereturns = tail(model$modeldata$data[1:N], m), 
			preresiduals = tail(residuals(fit), m), 
			startMethod = "sample", 
			custom.dist = list(name = "sample", distfit = as.matrix(empz)),
			rseed = sseed1, mexsimdata = mexsimdata, vexsimdata = vexsimdata)
	fitlist = vector(mode="list", length = n.bootfit)
	simseries = fitted(paths)
	spec = getspec(fit)
	# help the optimization with good starting parameters
	spec@model$start.pars = as.list(coef(fit))
	nx = NCOL(simseries)
	
	# get the distribution of the parameters (by fitting to the new path data)
	#-------------------------------------------------------------------------
	if( verbose ) cat("\nfitting stage...")
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("simseries", "spec", "solver", 
						"solver.control", "fit.control"), envir = environment())
		fitlist = parLapply(cluster, as.list(1:nx), fun = function(i){
					.safefit(spec = spec, data = .numeric2xts(simseries[,i]), 
							out.sample = 0, solver = solver, fit.control = fit.control, 
							solver.control = solver.control)
				})
	} else{
		fitlist = lapply(as.list(1:nx), FUN = function(i){
					.safefit(spec = spec, data = .numeric2xts(simseries[,i]), 
							out.sample = 0, solver = solver, fit.control = fit.control, 
							solver.control = solver.control)
				})
	}
	# check for any non convergence and remove
	if(any(sapply(fitlist, FUN=function(x) is.null(coef(x))))){
		exc = which(sapply(fitlist,FUN=function(x) is.null(coef(x))))
		fitlist = fitlist[-exc]
		n.bootfit = n.bootfit - length(exc)
		# in case something went very wrong:
		if(n.bootfit == 0) stop("\nugarchboot-->error: the fitting routine failed. No convergene at all!\n", call. = FALSE)
	}
	if( verbose ) cat("done!\n")
	# extract the coefficient distribution and generate spec list to feed to path function
	coefdist = matrix(NA, ncol = length(coef(fit)), nrow = n.bootfit)
	speclist = vector(mode = "list", length = n.bootfit)
	for(i in 1:n.bootfit){
		coefdist[i,]  = coef(fitlist[[i]])
		speclist[[i]] = spec
		speclist[[i]]@model$fixed.pars = as.list(coef(fitlist[[i]]))
	}
	colnames(coefdist) = names(coef(fit))
	#-------------------------------------------------------------------------
	# generate path based forecast values
	# for each path we generate n.bootpred vectors of resampled data of length n.ahead
	# Equation (6) in the PRR paper (again using z from original fit)
	#-------------------------------------------------------------------------
	empzlist = vector(mode = "list", length = n.bootfit)
	for(i in 1:n.bootfit){
		empz = matrix(0, ncol = n.bootpred, nrow = n.ahead)
		empz = apply(as.data.frame(1:n.bootpred), 1, FUN=function(i) {
					set.seed(sseed2[i]);
					.bsample(fz, n.ahead, sampling, zfit);})
		empzlist[[i]] = matrix(empz, ncol = n.bootpred)
	}
	# we start all forecasts from the last value of sigma based on the original series but
	# the pertrubed parameters as in equation (7) in PRR paper.
	if( !is.null(cluster) ){
		nx = length(speclist)
		clusterExport(cluster, c("speclist", "data", "m", "N"), envir = environment())
		xtmp = parLapply(cluster, as.list(1:n.bootfit), fun = function(i){
					ans = .sigmat(spec = speclist[[i]], origdata = .numeric2xts(data[1:N]), m)
					st = ans$st
					if(vmodel=="csGARCH"){
						qx = ans$qx
					} else{
						qx = rep(NA, m)
					}
					rx = ans$res
					return(list(st = st, rx = rx, qx = qx))
				})
		st = matrix(sapply(xtmp, FUN = function(x) x$st), ncol = n.bootfit)
		qx = matrix(sapply(xtmp, FUN = function(x) x$qx), ncol = n.bootfit)
		rx = matrix(sapply(xtmp, FUN = function(x) x$rx), ncol = n.bootfit)
	} else{
		xtmp = lapply(speclist, FUN=function(x){
					ans = .sigmat(spec = x, origdata = .numeric2xts(data[1:N]), m)
					st = ans$st
					if(vmodel=="csGARCH"){
						qx = ans$qx
					} else{
						qx = rep(NA, m)
					}
					rx = ans$res
					return(list(st = st, rx = rx, qx = qx))
		})
		st = matrix(sapply(xtmp, FUN = function(x) x$st), ncol = n.bootfit)
		qx = matrix(sapply(xtmp, FUN = function(x) x$qx), ncol = n.bootfit)
		rx = matrix(sapply(xtmp, FUN = function(x) x$rx), ncol = n.bootfit)
	}
	forcseries = NULL
	forcsigma  = NULL
	tmp = vector(mode = "list", length = n.bootfit)
	if( verbose ) cat("\nprediction stage...")
	xdat = tail( data[1:N], m )
	if( !is.null(cluster) ){
		clusterExport(cluster, c("fitlist", "n.ahead", "n.bootpred", "n.bootfit", 
						"st", "rx", "qx", "xdat", "sseed", "empzlist","mexsimdata", "vexsimdata"), 
				envir = environment())
		tmp = parLapply(cluster, as.list(1:n.bootfit), fun = function(i){
					try(.quicksimulate(fitlist[[i]], n.sim = n.ahead, 
							m.sim = n.bootpred, presigma = st[,i], prereturns = xdat, 
							preresiduals = rx[,i],
							n.start = 0, rseed = sseed[(n.bootfit+1):(n.bootfit + n.bootpred)], 
							custom.dist = list(name = "sample", distfit = as.matrix(empzlist[[i]])), 
							mexsimdata = mexsimdata, vexsimdata = vexsimdata, preq = qx[,i]), silent = TRUE)
				})
	} else{
		tmp = lapply(as.list(1:n.bootfit), FUN = function(i){
					try(.quicksimulate(fit = fitlist[[i]], n.sim = n.ahead, 
					m.sim = n.bootpred, presigma = st[,i], prereturns = xdat,  
					preresiduals = rx[,i],	
					n.start = 0, rseed = sseed[(n.bootfit+1):(n.bootfit + n.bootpred)], 
					custom.dist = list(name = "sample", distfit = as.matrix(empzlist[[i]])), 
					mexsimdata = mexsimdata, vexsimdata = vexsimdata, preq = qx[,i]), silent=TRUE)
		})
	}
	# we will have (n.bootfit x n.bootpred) x n.ahead matrix
	# eliminate errors:
	idxer = which(sapply(tmp, function(x) is(x, "try-error")))
	if(length(idxer)>0){
		tmp = tmp[-idxer]
		n.bootfit = length(tmp)
	}
	# we will have (n.bootfit x n.bootpred) x n.ahead matrix
	forcseries = lapply(tmp, FUN = function(x) x[1:n.ahead, ,drop = F])
	meanseriesfit = t(sapply(forcseries, FUN = function(x) apply(t(x), 2, "mean")))
	forcseries = matrix(unlist(forcseries), nrow = n.ahead, ncol = n.bootpred * n.bootfit, byrow = FALSE)
	forcseries = t(forcseries)
	forcsigma = lapply(tmp, FUN = function(x) x[(n.ahead+1):(2*n.ahead), ,drop = F])
	meansigmafit = t(sapply(forcsigma, FUN = function(x) apply(t(x), 2, "mean")))
	forcsigma = matrix(unlist(forcsigma), nrow = n.ahead, ncol = n.bootpred * n.bootfit, byrow = FALSE)
	forcsigma = t(forcsigma)
	# do some cleanup
	rm(tmp)
	rm(fitlist)
	rm(empzlist)
	gc(verbose = FALSE)	
	if( verbose ) cat("done!\n")
	#-------------------------------------------------------------------------

	# now we have the bootstrapped distribution of n.ahead forecast values
	# original forecast
	forc = ugarchforecast(fitORspec = fit, n.ahead = n.ahead, n.roll = 0, external.forecasts = external.forecasts)
	model$truecoef = coef(fit)
	model$modeldata$realized.x = realized.x
	model$modeldata$filtered.s = filtered.s
	model$n.ahead = n.ahead
	model$n.bootfit = n.bootfit
	model$n.bootpred = n.bootpred
	model$modeldata$meanfit.x = meanseriesfit
	model$modeldata$meanfit.s = meansigmafit
	model$seeds = sseed
	model$type = "full"
	model$indexT = model$modeldata$index[N]
	ans = new("uGARCHboot",
			fseries = forcseries,
			fsigma = forcsigma,
			bcoef = as.data.frame(coefdist),
			model = model,
			forc = forc)
	return(ans)
}


# method from a spec object
.ub1f2 = function(fitORspec, data = NULL, sampling = c("raw", "kernel", "spd"), 
		spd.options = list(upper = 0.9, lower = 0.1, type = "pwm", kernel = "normal"), 
		n.ahead = 10, n.bootfit = 100, n.bootpred = 500, out.sample = 0, rseed = NA, 
		solver = "solnp", solver.control = list(), fit.control = list(),  
		external.forecasts =  list(mregfor = NULL, vregfor = NULL), 
		mexsimdata = NULL, vexsimdata = NULL, cluster = NULL, verbose = FALSE)
{
	sampling = tolower(sampling[1])
	if(is.na(match(sampling, c("raw", "kernel", "spd")))) stop("\nsampling method not recognized")
	spec = fitORspec
	model = spec@model
	vmodel = model$modeldesc$vmodel
	if(is.null(data))
		stop("\nugarchboot-->error: data must be supplied if SPEC object used.", call. = FALSE)	
	flt = ugarchfilter(data = data, spec = spec, out.sample = out.sample)
	xdata = flt@model$modeldata$data
	m = spec@model$maxOrder
	N = flt@model$modeldata$T
	if (is.na(rseed[1])){
		sseed = as.integer(runif(n.bootpred + n.bootfit,0,65000))
	} else{
		if(length(rseed) < n.bootpred){
			stop("\nugarchboot-->error: seed length must equal n.bootpred + n.bootfit for full method\n")
		} else {
			sseed = rseed
		}
	}
	sseed1 = sseed[1:n.bootfit]
	sseed2 = sseed[(n.bootfit+1):(n.bootpred + n.bootfit)]
	
	if(out.sample > 0) {
		flt2 = ugarchfilter(data = data, spec = spec, out.sample = 0)
		realized.x = tail(xdata, out.sample)
		filtered.s = tail(flt2@filter$sigma, out.sample)
	} else{
		realized.x = NULL
		filtered.s = NULL
	}
	# generate paths of equal length to data based on empirical re-sampling of z
	# p.2296 equation (5)
	# use of filter when using spec/this will also check whether the fixed.pars are correctly
	# specified
	# ::length(data)-out.sample
	fz = flt@filter$z
	
	empz = matrix(0, ncol = n.bootfit, nrow = N)
	if(sampling!="raw"){
		zfit = .bfit(fz, sampling, spd.options)
	} else{
		zfit = NULL
	}
	
	empz = matrix(0, ncol = n.bootfit, nrow = N)
	empz = apply(as.data.frame(1:n.bootfit), 1, FUN=function(i) {
				set.seed(sseed1[i]);
				.bsample(fz, N, sampling, zfit);})	
	# presigma uses the same starting values as the original fit
	# -> in paper they use alternatively the unconditional long run sigma (P.2296 paragraph 2
	# "...marginal variance...")
	if(spec@model$modeldesc$vmodel=="csGARCH"){
		preq = tail(flt@filter$q, m)
	} else{
		preq = NULL
	}
	paths = ugarchpath(spec, n.sim = N, m.sim = n.bootfit, presigma = tail(flt@filter$sigma, m), 
			prereturns = tail(xdata[1:N], m), preresiduals = tail(residuals(flt), m), rseed = sseed1,
			n.start = 0, custom.dist = list(name = "sample", distfit = as.matrix(empz)), 
			mexsimdata = mexsimdata, vexsimdata = vexsimdata, preq = preq)
	fitlist = vector(mode="list", length = n.bootfit)
	simseries = fitted(paths)
	nx = NCOL(simseries)
	
	# help the optimization with good starting parameters
	spex = spec
	setfixed(spex) <- list(NA)
	setstart(spex) <- spec@model$fixed.pars
	# get the distribution of the parameters (by fitting to the new path data)
	#-------------------------------------------------------------------------
	if( verbose ) cat("\nfitting stage...")
	if( !is.null(cluster) ){
			clusterEvalQ(cluster, library(rugarch))
			nx = NCOL(simseries)
			clusterExport(cluster, c("simseries", "spex", "solver", "out.sample", 
							"solver.control", "fit.control"), envir = environment())
			fitlist = parLapply(cluster, as.list(1:nx), fun = function(i){
						.safefit(spec = spex, data = .numeric2xts(simseries[,i]), 
								out.sample = 0, solver = solver, fit.control = fit.control, 
								solver.control = solver.control)
					})
	} else{
		fitlist = lapply(as.list(1:nx), FUN = function(i){
					.safefit(spec = spex, data = .numeric2xts(simseries[,i]), 
							out.sample = 0, solver = solver, fit.control = fit.control, 
							solver.control = solver.control)
				})
	}
	# check for any non convergence and remove
	if(any(sapply(fitlist,FUN=function(x) is.null(coef(x))))){
		exc = which(sapply(fitlist,FUN=function(x) is.null(coef(x))))
		fitlist = fitlist[-exc]
		n.bootfit = n.bootfit - length(exc)
		# in case something went very wrong:
		if(n.bootfit == 0) stop("\nugarchboot-->error: the fitting routine failed. No convergence at all!\n", call. = FALSE)
	}
	if( verbose ) cat("done!\n")
	
	# extract the coefficient distribution and generate spec list to feed to path function
	coefdist = matrix(NA, ncol = length(coef(flt)), nrow = n.bootfit)
	speclist = vector(mode = "list", length = n.bootfit)
	for(i in 1:n.bootfit){
		coefdist[i,] = coef(fitlist[[i]])
		speclist[[i]] = spex
		speclist[[i]]@model$start.pars = NULL
		speclist[[i]]@model$fixed.pars = as.list(coef(fitlist[[i]]))
	}
	colnames(coefdist) = names(coef(flt))
	
	# generate path based forecast values
	# for each path we generate n.bootpred vectors of resampled data of length n.ahead
	# Equation (6) in the PRR paper (again using z from original fit)
	empzlist = vector(mode = "list", length = n.bootfit)
	for(i in 1:n.bootfit){
		empz = matrix(0, ncol = n.bootpred, nrow = n.ahead)
		empz = apply(as.data.frame(1:n.bootpred), 1, FUN=function(i) {
					set.seed(sseed2[i]);
					.bsample(fz, n.ahead, sampling, zfit);})
		empzlist[[i]] = matrix(empz, ncol = n.bootpred)
	}
	
	# we start all forecasts from the last value of sigma based on the original series but
	# the pertrubed parameters as in equation (7) in PRR paper.
	if( !is.null(cluster) ){
		nx = length(speclist)
		clusterExport(cluster, c("speclist", "xdata", "m", "N"), envir = environment())
		xtmp = parLapply(cluster, as.list(1:n.bootfit), fun = function(i){
					ans = .sigmat(spec = speclist[[i]], origdata = .numeric2xts(xdata[1:N]), m)
					st = ans$st
					if(vmodel=="csGARCH"){
						qx = ans$qx
					} else{
						qx = rep(NA, m)
					}
					rx = ans$res
					return(list(st = st, rx = rx, qx = qx))
				})
		
		st = matrix(sapply(xtmp, FUN = function(x) x$st), ncol = n.bootfit)
		qx = matrix(sapply(xtmp, FUN = function(x) x$qx), ncol = n.bootfit)
		rx = matrix(sapply(xtmp, FUN = function(x) x$rx), ncol = n.bootfit)
	} else{
		xtmp = lapply(speclist, FUN=function(x){
					ans = .sigmat(spec = x, origdata = .numeric2xts(xdata[1:N]), m)
					st = ans$st
					if(vmodel=="csGARCH"){
						qx = ans$qx
					} else{
						qx = rep(NA, m)
					}
					rx = ans$res
					return(list(st = st, rx = rx, qx = qx))
				})
		st = matrix(sapply(xtmp, FUN = function(x) x$st), ncol = n.bootfit)
		qx = matrix(sapply(xtmp, FUN = function(x) x$qx), ncol = n.bootfit)
		rx = matrix(sapply(xtmp, FUN = function(x) x$rx), ncol = n.bootfit)
	}
	forcseries = NULL
	forcsigma  = NULL
	forcseries = NULL
	forcsigma  = NULL
	tmp = vector(mode = "list", length = n.bootfit)
	if( verbose ) cat("\nprediction stage...")
	
	xdat =  tail(xdata[1:N], m)
	if( !is.null(cluster) ){
		clusterExport(cluster, c("fitlist", "n.ahead", "n.bootpred", "n.bootfit", 
						"st", "rx", "qx", "xdat", "sseed", "empzlist", "mexsimdata", "vexsimdata"), 
				envir = environment())
		tmp = parLapply(cluster, as.list(1:n.bootfit), fun = function(i){
					try(.quicksimulate(fitlist[[i]], n.sim = n.ahead, 
							m.sim = n.bootpred, presigma = st[,i], prereturns = xdat,
							preresiduals = rx[,i], n.start = 0, 
							rseed = sseed[(n.bootfit+1):(n.bootfit + n.bootpred)], 
							custom.dist = list(name = "sample", distfit = as.matrix(empzlist[[i]])), 
							mexsimdata = mexsimdata, vexsimdata = vexsimdata, preq = qx[,i]), silent = TRUE)
				})
	} else{
		tmp = lapply(as.list(1:n.bootfit), FUN = function(i){
					try(.quicksimulate(fitlist[[i]], n.sim = n.ahead, m.sim = n.bootpred, 
							presigma = st[,i], prereturns = xdat,  
							preresiduals = rx[,i], n.start = 0, 
							rseed = sseed[(n.bootfit+1):(n.bootfit + n.bootpred)], 
							custom.dist = list(name = "sample", distfit = as.matrix(empzlist[[i]])), 
							mexsimdata = mexsimdata, vexsimdata = vexsimdata, preq = qx[,i]), silent = TRUE)
				})
	}
	# we will have (n.bootfit x n.bootpred) x n.ahead matrix
	# eliminate errors:
	idxer = which(sapply(tmp, function(x) is(x, "try-error")))
	if(length(idxer)>0){
		tmp = tmp[-idxer]
		n.bootfit = length(tmp)
	}
	forcseries = lapply(tmp, FUN = function(x) x[1:n.ahead, ,drop = F])
	meanseriesfit = t(sapply(forcseries, FUN = function(x) apply(t(x), 2, "mean")))
	forcseries = matrix(unlist(forcseries), nrow = n.ahead, ncol = n.bootpred * n.bootfit, byrow = FALSE)
	forcseries = t(forcseries)
	forcsigma = lapply(tmp, FUN = function(x) x[(n.ahead+1):(2*n.ahead), ,drop = F])
	meansigmafit = t(sapply(forcsigma, FUN = function(x) apply(t(x), 2, "mean")))
	forcsigma = matrix(unlist(forcsigma), nrow = n.ahead, ncol = n.bootpred * n.bootfit, byrow = FALSE)
	forcsigma = t(forcsigma)
	# do some cleanup
	rm(empzlist)
	rm(tmp)
	rm(fitlist)
	gc(verbose = FALSE)
	if( verbose ) cat("done!\n")
	#-------------------------------------------------------------------------
	# now we have the bootstrapped distribution of n.ahead forecast values
	# original forecast
	forc = ugarchforecast(fitORspec = spec, data = head(data, N), n.ahead = n.ahead, n.roll = 0, external.forecasts = external.forecasts)	
	model$truecoef = coef(flt)
	model$modeldata$realized.x = realized.x
	model$modeldata$filtered.s = filtered.s
	model$n.ahead = n.ahead
	model$n.bootfit = n.bootfit
	model$n.bootpred = n.bootpred
	model$modeldata$meanfit.x = meanseriesfit
	model$modeldata$meanfit.s = meansigmafit
	model$seeds = sseed
	model$type = "full"
	model$indexT = flt@model$modeldata$index[N]	
	ans = new("uGARCHboot",
			fseries = forcseries,
			fsigma = forcsigma,
			bcoef = as.data.frame(coefdist),
			model = model,
			forc = forc)
	return(ans)
}

# Partial method (very fast but does not take into account the parameter uncertainty)
.ub1p1 = function(fitORspec, data = NULL, sampling = c("raw", "kernel", "spd"), 
		spd.options = list(upper = 0.9, lower = 0.1, type = "pwm", kernel = "normal"), 
		n.ahead = 10, n.bootfit = 100, n.bootpred = 500, rseed = NA, 
		solver.control = list(), fit.control = list(),
		external.forecasts =  list(mregfor = NULL, vregfor = NULL), 
		mexsimdata = NULL, vexsimdata = NULL,  cluster = NULL, verbose = FALSE)
{
	sampling = tolower(sampling[1])
	if(is.na(match(sampling, c("raw", "kernel", "spd")))) stop("\nsampling method not recognized")
	fit = fitORspec
	model = fit@model
	m = fit@model$maxOrder
	ns = fit@model$n.start
	data = xts(fit@model$modeldata$data, fit@model$modeldata$index)
	xdata = fit@model$modeldata$data
	N = length(xdata) - ns
	spec = getspec(fit)
	if(is.na(rseed[1])){
		sseed = as.integer(runif(n.bootpred,0,65000))
	} else{
		if(length(rseed) < n.bootpred){
			stop("\nugarchboot-->error: seed length must equal n.bootpred for partial method\n")
		} else {
			sseed = rseed
		}
	}
	# generate path based forecast values
	# for each path we generate n.bootpred vectors of resampled data of length n.ahead
	# Equation (6) in the PRR paper
	fz = fit@fit$z
	
	if(sampling!="raw"){
		zfit = .bfit(fz, sampling, spd.options)
	} else{
		zfit = NULL
	}
	
	empz = matrix(apply(as.data.frame(1:n.bootpred), 1, 
					FUN = function(i) {
						set.seed(sseed[i]);
						.bsample(fz, n.ahead, sampling, zfit);}), 
			ncol = n.bootpred)
	
	#empz[1:m,] = matrix(rep(tail(fit@fit$z, m), n.bootpred), nrow = m)
	
	# we start all forecasts from the last value of sigma based on the original series
	setfixed(spec) <- as.list(coef(fit))
	#st = .sigmat(spec, data[1:N], m)
	
	if(ns > 0){
		realized.x = data[(N+1):(N+ns)]
		filtered.s = tail(sigma(ugarchfilter(data = data, spec = spec)), ns)
	} else{
		realized.x = NULL
		filtered.s = NULL
	}
	
	forcseries = matrix(NA, ncol = n.ahead, nrow = n.bootpred)
	forcsigma  = matrix(NA, ncol = n.ahead, nrow = n.bootpred)
	sim = vector(mode = "list", length = 1)
	sim = ugarchsim(fit, n.sim = n.ahead, m.sim = n.bootpred, 
				presigma = tail(fit@fit$sigma, m), 
				prereturns = tail(data[1:N], m), 
				preresiduals = tail(residuals(fit), m),
				rseed = sseed, n.start = 0, startMethod = "sample", 
				custom.dist = list(name = "sample", distfit = empz),
				mexsimdata = mexsimdata, vexsimdata = vexsimdata)
	# we transpose to get n.boot x n.ahead
	forcseries = t(sim@simulation$seriesSim)
	forcsigma =  t(sim@simulation$sigmaSim)
	
	# now we have the bootstrapped distribution of n.ahead forecast values
	# original forecast
	forc = ugarchforecast(fitORspec = fit, n.ahead = n.ahead, n.roll = 0, external.forecasts = external.forecasts)
	coefdist = as.data.frame(coef(fit))
	model$truecoef = coef(fit)
	model$modeldata$realized.x = realized.x
	model$modeldata$filtered.s = filtered.s
	model$n.ahead = n.ahead
	model$n.bootfit = n.bootfit
	model$n.bootpred = n.bootpred
	model$type = "partial"
	model$indexT = model$modeldata$index[N]
	ans = new("uGARCHboot",
			fseries = (forcseries),
			fsigma = (forcsigma),
			bcoef = coefdist,
			model = model,
			forc = forc)
	return(ans)
}

# using spec method
.ub1p2 = function(fitORspec, data = NULL, sampling = c("raw", "kernel", "spd"), 
		spd.options = list(upper = 0.9, lower = 0.1, type = "pwm", kernel = "normal"), 
		n.ahead = 10, n.bootfit = 100, n.bootpred = 500, out.sample = 0, rseed = NA, 
		solver.control = list(), fit.control = list(), 
		external.forecasts =  list(mregfor = NULL, vregfor = NULL), 
		mexsimdata = NULL, vexsimdata = NULL, cluster = NULL, verbose = FALSE)
{
	sampling = tolower(sampling[1])
	if(is.na(match(sampling, c("raw", "kernel", "spd")))) stop("\nsampling method not recognized")
	spec = fitORspec
	if(is.null(data))
		stop("\nugarchboot-->error: data must be supplied if SPEC object used.", call. = FALSE)
	model = spec@model
	flt = ugarchfilter(data = data, spec = spec, out.sample = out.sample)
	tmp = .extractdata(data)
	xdata = tmp$data
	index = tmp$index
	period = tmp$period
	ns = out.sample
	N = length(xdata) - out.sample
	sigma = flt@filter$sigma
	m = spec@model$maxOrder
	
	if (is.na(rseed[1])){
		sseed = as.integer(runif(n.bootpred,0,65000))
	} else{
		if(length(rseed) < n.bootpred){
			stop("\nugarchboot-->error: seed length must equal n.bootpred for partial method\n", call. = FALSE)
		} else {
			sseed = rseed
		}
	}
	# generate path based forecast values
	fz = flt@filter$z
	if(sampling!="raw"){
		zfit = .bfit(fz, sampling, spd.options)
	} else{
		zfit = NULL
	}
	
	empz = matrix(apply(as.data.frame(1:n.bootpred), 1, 
					FUN = function(i) {
						set.seed(sseed[i]);
						.bsample(fz, n.ahead, sampling, zfit);}), 
			ncol = n.bootpred)
	
	if(ns > 0) {
		realized.x = xdata[(N+1):(N+ns)]
		filtered.s = tail(sigma(ugarchfilter(data = data, spec = spec)), ns)
	} else{
		realized.x = NULL
		filtered.s = NULL
	}
	
	# we start all forecasts from the last value of sigma based on the original series
	#st = .sigmat(spec, data, m)
	
	if(spec@model$modeldesc$vmodel=="csGARCH"){
		preq = tail(flt@filter$q, m)
	} else{
		preq = NULL
	}
	forcseries = matrix(NA, ncol = n.ahead, nrow = n.bootpred)
	forcsigma  = matrix(NA, ncol = n.ahead, nrow = n.bootpred)
	sim = ugarchpath(spec, n.sim = n.ahead, m.sim = n.bootpred, presigma = tail(sigma, m),
				prereturns = tail(xdata[1:N], m), preresiduals = tail(flt@filter$residuals, m),
				rseed = sseed, n.start = 0, custom.dist = list(name = "sample", 
						distfit = as.matrix(empz)), mexsimdata = mexsimdata, 
				vexsimdata = vexsimdata, preq = preq)
	# we transpose to get n.boot x n.ahead
	forcseries = t(sim@path$seriesSim)
	forcsigma =  t(sim@path$sigmaSim)
	
	# now we have the bootstrapped distribution of n.ahead forecast values
	# original forecast
	forc = ugarchforecast(fitORspec = spec, data = head(data, N), n.ahead = n.ahead, n.roll = 0, external.forecasts = external.forecasts) 
	
	model$truecoef = coef(flt)
	model$modeldata$realized.x = realized.x
	model$modeldata$filtered.s = filtered.s
	model$n.ahead = n.ahead
	model$n.bootfit = n.bootfit
	model$n.bootpred = n.bootpred
	coefdist = data.frame(NULL)
	model$type = "partial"
	model$indexT = flt@model$modeldata$index[N]
	ans = new("uGARCHboot",
			fseries = (forcseries),
			fsigma = (forcsigma),
			bcoef = coefdist,
			model = model,
			forc = forc)
	return(ans)
}



.sigmat = function(spec, origdata, m)
{
	flt = ugarchfilter(data = origdata, spec = spec)
	st = tail(flt@filter$sigma, m)
	res = tail(flt@filter$residuals, m)
	if(spec@model$modeldesc$vmodel=="csGARCH") qx = tail(flt@filter$q, m) else qx = NULL
	return(list(st = st, res = res, qx = qx))
}

.quicksimulate = function(fit, n.sim, m.sim, presigma = NA, prereturns = NA,  
		preresiduals = NA, n.start = 0, rseed = NA, 
		custom.dist = list(name = "sample", distfit = NULL), mexsimdata = NULL, 
		vexsimdata = NULL, preq = NULL)
{
	ans = ugarchsim(fit = fit, n.sim = n.sim, m.sim = m.sim, presigma = presigma, 
			prereturns = prereturns, preresiduals = preresiduals, 
			n.start = n.start, rseed = rseed, custom.dist = custom.dist, 
			mexsimdata = mexsimdata, vexsimdata = vexsimdata, preq = preq)
	ret = rbind(ans@simulation$seriesSim, ans@simulation$sigmaSim)
	return(ret)
}

.bsample = function(Z, n, method, fit)
{
	return(switch(tolower(method),
			"raw" = sample(Z, n, replace = TRUE),
			"kernel" = rkde(n, fit),
			"spd" = rspd(n, fit)))
}

.bfit = function(Z, method, spd.options){
	return(switch(tolower(method),
					"kernel" = .kfit(Z),
					"spd" = .sfit(Z, spd.options)))
}

.kfit = function(Z){
	H.pi = hpi(Z)
	fit = kde(Z, h=H.pi)
	return(fit)
}

.sfit = function(Z, spd.options){
	mm = match(names(spd.options), c("upper", "lower", "type", "kernel"))
	if(any(is.na(mm))){
		idx = which(is.na(mm))
		enx = NULL
		for(i in 1:length(idx)) enx = c(enx, names(spd.options)[idx[i]])
		warning(paste(c("unidentified option(s) in spd.options:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
	}
	if(is.null(spd.options$upper)) upper = 0.9 else upper = spd.options$upper
	if(is.null(spd.options$upper)) lower = 0.1 else lower = spd.options$lower
	if(is.null(spd.options$type)) type = "pwm" else type  = spd.options$type[1]
	if(is.null(spd.options$kernel)) kernel = "normal" else kernel  = spd.options$kernel[1]
	fit = spdfit(Z, upper = upper, lower = lower, type = type, kernelfit = kernel)
	return(fit)
}