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

# master branch
.rollfdensity = function(spec, data, n.ahead = 1, forecast.length = 500, 
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
		clusterExport(cluster, c("data", "index", "s","refit.every", 
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
									fit.control = fit.control), silent=TRUE)
					# 3 cases: General Error, Failure to Converge, Failure to invert Hessian (bad solution)
					if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
						ans = list(y = NA, cf = NA, converge = FALSE, loglik = NA)
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
						ans = list(y = y, cf = cf, converge = TRUE, loglik = likelihood(fit))
					}
					return(ans)})
	} else{
		tmp = lapply(as.list(1:m), FUN = function(i){
					if(mex) spec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
					if(vex) spec@model$modeldata$vexdata = vexdata[rollind[[i]],,drop=FALSE]
					fit = try(ugarchfit(spec, xts(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i], 
									solver = solver, solver.control = solver.control, 
									fit.control = fit.control), silent=TRUE)
					if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
						ans = list(y = NA, cf = NA, converge = FALSE, loglik = NA)
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
						ans = list(y = y, cf = cf, converge = TRUE, loglik = likelihood(fit))
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
		model$loglik = sapply(tmp, function(x) x$loglik)
		forecast = list(VaR = VaR.matrix, density = forc)
	}
	toc = Sys.time()-tic
	model$elapsed = toc
	ans = new("uGARCHroll",
			model = model,
			forecast = forecast)
	return( ans )
}

.resumeroll1 = function(object, solver = "hybrid", fit.control = list(), 
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
							"noncidx", "solver", "solver.control", "fit.control"),
					envir = environment())
			if(mex) clusterExport(cluster, c("mexdata"), envir = environment())
			if(vex) clusterExport(cluster, c("vexdata"), envir = environment())
			tmp = parLapply(cl = cluster, as.list(noncidx), fun = function(i){
						if(mex) spec@model$modeldata$mexdata = mexdata[rollind[[i]],,drop=FALSE]
						if(vex) spec@model$modeldata$vexdata = vexdata[rollind[[i]],,drop=FALSE]
						fit = try(ugarchfit(spec, xts::xts(data[rollind[[i]]], index[rollind[[i]]]), out.sample = out.sample[i], 
										solver=solver, solver.control = solver.control, 
										fit.control = fit.control), silent=TRUE)
						if(inherits(fit, 'try-error') || convergence(fit)!=0 || is.null(fit@fit$cvar)){
							ans = list(y = NA, cf = NA, converge = FALSE, loglik = NA)
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
							ans = list(y = y, cf = cf, converge = TRUE, loglik = likelihood(fit))
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
							ans = list(y = NA, cf = NA, converge = FALSE, loglik = NA)
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
							ans = list(y = y, cf = cf, converge = TRUE, loglik = likelihood(fit))
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
			model$loglik = sapply(tmp, function(x) x$loglik)
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


# forecast performance measures
.ugarchrollreport = function(object, type = "VaR", VaR.alpha = 0.01, conf.level = 0.95)
{
	if(!is.null(object@model$noncidx)) stop("\nObject contains non-converged estimation windows.")
	switch(type,
			VaR = .rollVaRreport1(object, VaR.alpha, conf.level),
			fpm = .rollfpmreport1(object))
	invisible(object)
}

.rollfpmreport1 = function(object)
{
	vmodel = object@model$spec@model$modeldesc$vmodel
	vsubmodel = object@model$spec@model$modeldesc$vsubmodel
	cat(paste("\nGARCH Roll Mean Forecast Performance Measures", sep = ""))
	cat(paste("\n---------------------------------------------", sep = ""))
	cat(paste("\nModel\t\t: ", vmodel, sep = ""))
	if(vmodel == "fGARCH"){
		cat(paste("\nSubModel : ", vsubmodel, sep = ""))
	}
	cat(paste("\nNo.Refits\t: ", object@model$n.refits, sep = ""))
	cat(paste("\nNo.Forecasts: ", NROW(object@forecast$density), sep = ""))
	cat("\n\n")
	tmp = fpm(object)
	print(signif(tmp, 4))
	cat("\n\n")
}

.rollVaRreport1 = function(object, VaR.alpha = 0.01, conf.level = 0.95)
{
	VaR.alpha = VaR.alpha[1]
	v.a = object@model$VaR.alpha
	if(!object@model$calculate.VaR) stop("\nplot-->error: VaR was not calculated for this object\n", call.=FALSE)
	if(!any(v.a==VaR.alpha[1])) stop("\nplot-->error: VaR.alpha chosen is invalid for the object\n", call.=FALSE)
		dvar = object@forecast$VaR
		m = NROW(dvar)
		idx = match(VaR.alpha, v.a)
		.VaRreport(if(is.null(object@model$datanames)) "" else object@model$datanames, object@model$spec@model$modeldesc$vmodel, 
				object@model$spec@model$modeldesc$distribution, 
				p = VaR.alpha, actual = as.numeric(dvar[,"realized"]), 
				VaR = dvar[, idx], 
				conf.level = conf.level)
	invisible(object)
}

.embed = function(data, k, by = 1, ascending = FALSE) 
{
	# n = no of time points, k = number of columns
	# by = increment. normally =1 but if =b calc every b-th point 
	# ascending If TRUE, points passed in ascending order else descending.
	# Note that embed(1:n,k) corresponds to embedX(n,k,by=1,rev=TRUE)
	# e.g. embedX(10,3)
	if(is.null(dim(data)[1])) n<-length(data) else n<-dim(data)[1]
	s <- seq(1,n-k+1,by)
	lens <- length(s)
	cols <- if (ascending) 1:k else k:1
	return(matrix(data[s + rep(cols,rep(lens,k))-1],lens))
}
