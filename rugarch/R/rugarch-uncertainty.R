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

# inference of the parameters distribution via simulation

.ugarchdistribution = function(fitORspec, n.sim = 2000, n.start = 1, 
		m.sim = 100,  recursive = FALSE, recursive.length = 6000, recursive.window = 1000,
		presigma = NA, prereturns = NA, preresiduals = NA, rseed = NA,
		custom.dist = list(name = NA, distfit = NA), mexsimdata = NULL, vexsimdata = NULL, 
		fit.control = list(), solver = "solnp", solver.control = list(), cluster = NULL, ...)
{
	# recursive = FALSE, recursive.maxlength = 6000, recursive.window = 1000,
	# simulate
	# fit
	# extract distribution
	# look at the joint for persistence (apply persistence)
	# look at the joint for arma
	# create plot functions and report
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
	if(is(fitORspec, "uGARCHfit")){
		for(i in 1:nwindows){
			sim = ugarchsim(fitORspec, n.sim = n.sim + (i-1)*recursive.window, 
					n.start = n.start, m.sim = m.sim, presigma = presigma, prereturns = prereturns, 
					preresiduals = preresiduals, rseed = rseed, custom.dist = custom.dist, 
					mexsimdata = mexsimdata, vexsimdata = vexsimdata)
			swindow[[i]]$simseries = fitted(sim)
			swindow[[i]]$seed = sim@seed
		}
		fixpars = as.list(coef(fitORspec))
		truecoef = fitORspec@fit$robust.matcoef
		spec = getspec(fitORspec)
	}
	# simulate series paths
	if(is(fitORspec, "uGARCHspec")){

		for(i in 1:nwindows){
			sim  = ugarchpath(fitORspec, n.sim =  n.sim + (i-1)*recursive.window, 
					n.start = n.start, m.sim = m.sim, presigma = presigma, prereturns = prereturns, 
					preresiduals = preresiduals, rseed = rseed, custom.dist = custom.dist, 
					mexsimdata = mexsimdata, vexsimdata = vexsimdata)
			swindow[[i]]$simseries = fitted(sim)
			swindow[[i]]$seed = sim@seed
			
		}
		spec = fitORspec
		setfixed(spec)<-list(NA)
		fixpars = fitORspec@model$fixed.pars
		truecoef = as.matrix(cbind(unlist(fitORspec@model$fixed.pars),rep(0,length(fixpars)),
						rep(10, length(fixpars)),rep(0,length(fixpars))))
	}
	# fit to simulated series (using starting parameters)
	#setstart(spec) = fixpars

	fitlist = vector( mode = "list", length = m.sim )
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("spec", "swindow", "solver", "fit.control", "solver.control"), 
				envir = environment())
		for(i in 1:nwindows){
			clusterExport(cluster, c("i"), envir = environment())
			nx = NCOL(swindow[[i]]$simseries)
			rwindow[[i]]$fitlist = parLapply(cluster, as.list(1:nx), fun = function(j){
						.fitandextract(spec, .numeric2xts(swindow[[i]]$simseries[,j]), 
								out.sample = 0, solver = solver, fit.control = fit.control, 
								solver.control = solver.control)
					})
			}
	} else{
		for(i in 1:nwindows){
			rwindow[[i]]$fitlist = apply(swindow[[i]]$simseries, 2,
					FUN = function(x){
						.fitandextract(spec, .numeric2xts(x), out.sample = 0, solver = solver, 
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
		reslist[[j]]$persist = 	rep(NA, length = m.sim)
		reslist[[j]]$vlongrun = rep(NA, length = m.sim)
		reslist[[j]]$mlongrun = rep(NA, length = m.sim)
		reslist[[j]]$simmaxdata  = matrix(NA, ncol = 3, nrow = m.sim)
		reslist[[j]]$simmindata  = matrix(NA, ncol = 3, nrow = m.sim)
		reslist[[j]]$simmeandata  = matrix(NA, ncol = 3, nrow = m.sim)
		reslist[[j]]$simmomdata = matrix(NA, ncol = 2, nrow = m.sim)
		reslist[[j]]$convergence = 	rep(1, length = m.sim)
		reslist[[j]]$seeds = 	rep(1, length = m.sim)
		for(i in 1:m.sim){
			if(rwindow[[j]]$fitlist[[i]]$convergence!=0) next()
			reslist[[j]]$simcoef[i, ] = rwindow[[j]]$fitlist[[i]]$simcoef
			reslist[[j]]$simcoefse[i,] = rwindow[[j]]$fitlist[[i]]$simcoefse
			reslist[[j]]$likelist[i] = 	rwindow[[j]]$fitlist[[i]]$llh
			reslist[[j]]$persist[i] = 	rwindow[[j]]$fitlist[[i]]$persist
			reslist[[j]]$vlongrun[i] = 	rwindow[[j]]$fitlist[[i]]$vlongrun
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
			recursive = recursive, recursive.length = recursive.length, 
			recursive.window = recursive.window, nwindows = nwindows)
	ans = new("uGARCHdistribution",
			dist = reslist,
			truecoef = truecoef,
			model = fitORspec@model)
	
	return(ans)
}

.fitandextract = function(spec, x, out.sample = 0,  solver = "solnp", fit.control = list(), solver.control = list())
{
	dist = list()
	fit = .safefit(spec, x, out.sample = 0, solver = solver, fit.control = fit.control, solver.control = solver.control)
	if( is.null(fit) || fit@fit$convergence == 1 || !is( fit, "uGARCHfit" ) || any( is.na( coef( fit ) ) ) ){
		dist$convergence = 1
		return(dist)
	}
	dist = list()
	dist$simcoef = coef(fit)
	dist$simcoefse = fit@fit$robust.matcoef[, 2]
	dist$llh = likelihood(fit)
	dist$persist = 	persistence(fit)
	dist$vlongrun = uncvariance(fit)
	dist$mlongrun = uncmean(fit)
	tmp = cbind(fit@model$modeldata$data[1:fit@model$modeldata$T], as.numeric(fitted(fit)), as.numeric(sigma(fit)))
	dist$maxdata = unname(apply(tmp, 2, "max"))
	dist$mindata = unname(apply(tmp, 2, "min"))
	dist$meandata = unname(apply(tmp, 2, "mean"))
	dist$momdata = c(.kurtosis(tmp[,1]), .skewness(tmp[,1]))
	dist$convergence = fit@fit$convergence
	return(dist)
}


.rmse = function(est, act)
{
	exc = unique(which(is.na(est), arr.ind = TRUE)[,1])
	if(length(exc)>0) est = est[-exc, , drop = FALSE]
	n = dim(est)[1]
	diff = est - repmat(t(act), n, 1)
	ans = apply(diff, 2, FUN = function(x) sum((x)^2)/(n-1) )
	return(sqrt(ans))
}

.exprmseroc = function(window)
{
	sqrt(window[-1] / window[1])
}

.kurtosis = function(x)
{
	sum((x-mean(x))^4/var(x)^2)/length(x) - 3
}

.skewness = function(x)
{
	sum((x-mean(x))^3/sqrt(var(x))^3)/length(x)
}