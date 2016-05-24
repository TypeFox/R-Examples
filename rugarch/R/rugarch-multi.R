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

.multispecall = function( speclist ){
	model = unlist( strsplit(class(speclist[[1]]), "spec") )
	if( model == "ARFIMA" ){
		ans = .multispecarfima( speclist )
	} else{
		ans = .multispecgarch( speclist )
	}
	return( ans )
}

.multispecgarch = function( speclist )
{
	# first create a spec which goes through validation process
	tp = 1
	if( !all(unlist(lapply(speclist, FUN = function(x) is(x, "uGARCHspec"))) ) ){
		stop("\nNot a valid list of univariate GARCH specs.")
	}
	# then check type
	n = length(speclist)
	for(i in 2:n){
		modelnames1 = rownames( speclist[[i]]@model$pars[speclist[[i]]@model$pars[,3]==1, ] )
		modelnames2 = rownames( speclist[[i-1]]@model$pars[speclist[[i-1]]@model$pars[,3]==1, ] )
		if(length(modelnames1) != length(modelnames2))
		{
			tp  = 0
			break()
		} else{
			if(!all(modelnames1 == modelnames2))
			{
				tp  = 0
				break()
			}
		}
	}
	if(tp) type = "equal" else type = "unequal"
	if(type=="unequal"){
		# mcsGARCH and realGARCH cannot be in unequal specification. Either all the same or none.
		mod = sapply(speclist, function(x) x@model$modeldesc$vmodel)
		if(any(mod=="mcsGARCH")) stop("\nmultispec-->error: cannot have unequal spec containing mcsGARCH model.\n")
		if(any(mod=="realGARCH")) stop("\nmultispec-->error: cannot have unequal spec containing realGARCH model.\n")
	}
	ans = new("uGARCHmultispec",
			spec = speclist,
			type = type)
	return(ans)
}

# a multifit function possible utilizing parallel execution returning a fitlist
# object
.multifitgarch = function(multispec, data, out.sample = 0, solver = "solnp", 
		solver.control = list(), fit.control = list(stationarity = 1, fixed.se = 0, scale = 0, 
				rec.init = "all"), cluster = NULL, ...)
{
	n = length(multispec@spec)
	if(is.null(data)) stop("\nmultifit GARCH-->error: multifit requires a data object", call. = FALSE)
	if(!is.xts(data) & !is.matrix(data) & !is.data.frame(data)) stop("\nmultifit GARCH-->error: multifit only supports xts, matrix or data.frame objects for the data", call. = FALSE)
	asset.names = colnames(data)
	if(dim(data)[2] != n)
		stop("\nmultifit GARCH-->error: speclist length not equal to data length", call. = FALSE)
	fitlist = vector(mode = "list", length = n)
	if(length(out.sample) == 1 | length(out.sample) < n) out.sample = rep(out.sample, n)
	if(multispec@type=="equal"){
		mod = multispec@spec[[1]]@model$modeldesc$vmodel
		if(mod=="realGARCH"){
			realVol = list(...)$realizedVol
			if(is.null(realVol)) stop("\nmultifit-->error: realGARCH model requires realizedVol xts matrix.")
			if(!is.xts(realVol)) stop("\nmultifit-->error: realizedVol must be an xts matrix.")
			if(ncol(realVol)!=n) stop("\nmultifit-->error: realizedVol must have the same number of columns as the data.")
			if(nrow(realVol)!=nrow(data)) stop("\nmultifit-->error: realizedVol must have the same number of rows as the data.")
		}
		if(mod=="mcsGARCH"){
			DailyVar = list(...)$DailyVar
			if(is.null(DailyVar)) stop("\nmultifit-->error: mcsGARCH model requires DailyVar xts matrix.")
			if(!is.xts(DailyVar)) stop("\nmultifit-->error: DailyVar must be an xts matrix.")
			if(ncol(DailyVar)!=n) stop("\nmultifit-->error: DailyVar must have the same number of columns as the data.")
		}
	} else{
		mod = "X"
	}
	##################
	# Parallel Execution Prelim Checks
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("multispec", "data", "out.sample", "solver", 
						"solver.control", "fit.control"), envir = environment())
		if(mod=="realGARCH"){
			clusterExport(cluster, "realVol", envir = environment())
			fitlist = parLapply(cluster, as.list(1:n), fun = function(i){
						ans<-try(ugarchfit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
								out.sample = out.sample[i], solver = solver, 
								solver.control = solver.control, fit.control = fit.control, 
								realizedVol = realVol[,i]), silent=TRUE)
						if(inherits(ans, 'try-error')){
							ans<-try(ugarchfit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
											out.sample = out.sample[i], solver = "gosolnp", 
											fit.control = fit.control, 
											realizedVol = realVol[,i]), silent=TRUE)
						}
						if(convergence(ans)==1){
							ans<-try(ugarchfit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
											out.sample = out.sample[i], solver = "gosolnp", 
											fit.control = fit.control, 
											realizedVol = realVol[,i]), silent=TRUE)
						}
						return(ans)
					})
		} else if(mod=="mcsGARCH"){
			clusterExport(cluster, "DailyVar", envir = environment())
			fitlist = parLapply(cluster, as.list(1:n), fun = function(i){
						ans<-try(ugarchfit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
										out.sample = out.sample[i], solver = solver, 
										solver.control = solver.control, fit.control = fit.control, 
										DailyVar = DailyVar[,i]), silent=TRUE)
						if(inherits(ans, 'try-error')){
							ans<-try(ugarchfit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
											out.sample = out.sample[i], solver = "gosolnp", 
											fit.control = fit.control, 
											DailyVar = DailyVar[,i]), silent=TRUE)
						}
						if(convergence(ans)==1){
							ans<-try(ugarchfit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
											out.sample = out.sample[i], solver = "gosolnp", 
											fit.control = fit.control, 
											DailyVar = DailyVar[,i]), silent=TRUE)
						}
						return(ans)
					})
		} else{		
			fitlist = parLapply(cluster, as.list(1:n), fun = function(i){
						ans<-try(ugarchfit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
										out.sample = out.sample[i], solver = solver, 
										solver.control = solver.control, fit.control = fit.control), silent=TRUE)
						if(inherits(ans, 'try-error')){
							ans<-try(ugarchfit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
											out.sample = out.sample[i], solver = "gosolnp", 
											fit.control = fit.control), silent=TRUE)
						}
						if(convergence(ans)==1){
							ans<-try(ugarchfit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
											out.sample = out.sample[i], solver = "gosolnp", 
											fit.control = fit.control), silent=TRUE)
						}
						return(ans)
					})
		}
	} else{
		if(mod=="realGARCH"){
			fitlist = lapply(as.list(1:n), FUN = function(i){
						ans<-try(ugarchfit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
										out.sample = out.sample[i], solver = solver, 
										solver.control = solver.control, fit.control = fit.control, 
										realizedVol = realVol[,i]), silent=TRUE)
						if(inherits(ans, 'try-error')){
							ans<-try(ugarchfit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
											out.sample = out.sample[i], solver = "gosolnp", 
											fit.control = fit.control, 
											realizedVol = realVol[,i]), silent=TRUE)
						}
						if(convergence(ans)==1){
							ans<-try(ugarchfit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
											out.sample = out.sample[i], solver = "gosolnp", 
											fit.control = fit.control, 
											realizedVol = realVol[,i]), silent=TRUE)
						}
						return(ans)
					})
		} else if(mod=="mcsGARCH"){
			fitlist = lapply(as.list(1:n), FUN = function(i){
						ans<-try(ugarchfit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
										out.sample = out.sample[i], solver = solver, 
										solver.control = solver.control, fit.control = fit.control, 
										DailyVar = DailyVar[,i]), silent=TRUE)
						if(inherits(ans, 'try-error')){
							ans<-try(ugarchfit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
											out.sample = out.sample[i], solver = "gosolnp", 
											fit.control = fit.control, 
											DailyVar = DailyVar[,i]), silent=TRUE)
						}
						if(convergence(ans)==1){
							ans<-try(ugarchfit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
											out.sample = out.sample[i], solver = "gosolnp", 
											fit.control = fit.control, 
											DailyVar = DailyVar[,i]), silent=TRUE)
						}
						return(ans)
					})
		} else{			
			fitlist = lapply(as.list(1:n), FUN = function(i){
						ans<-try(ugarchfit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
										out.sample = out.sample[i], solver = solver, 
										solver.control = solver.control, fit.control = fit.control), silent=TRUE)
						if(inherits(ans, 'try-error')){
							ans<-try(ugarchfit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
											out.sample = out.sample[i], solver = "gosolnp", 
											fit.control = fit.control), silent=TRUE)
						}
						if(convergence(ans)==1){
							ans<-try(ugarchfit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
											out.sample = out.sample[i], solver = "gosolnp", 
											fit.control = fit.control), silent=TRUE)
						}
						return(ans)
					})
		}
	}
	# converged: print
	desc = list()
	desc$type = multispec@type
	desc$asset.names = asset.names
	ans = new("uGARCHmultifit",
			fit  = fitlist,
			desc = desc)
	return(ans)
}

.multifiltergarch1 = function(multifitORspec, data = NULL, out.sample = 0, 
		n.old = NULL, rec.init = 'all', cluster = NULL, ...)
{
	fitlist = multifitORspec
	n = length(fitlist@fit)
	if(is.null(data)) 
		stop("\nmultifilter GARCH-->error: multifilter requires a data object", call. = FALSE)
	if(!is.xts(data) & !is.matrix(data) & !is.data.frame(data)) 
		stop("\nmultifilter GARCH-->error: multifilter only supports xts, matrix or data.frame objects for the data", call. = FALSE)
	if(dim(data)[2] != n)
		stop("\nmultifilter GARCH-->error: fitlist length not equal to data length", call. = FALSE)
	asset.names = colnames(data)
	filterlist = vector(mode = "list", length = n)
	if(length(out.sample) == 1 | length(out.sample) < n) out.sample = rep(out.sample, n)
	if(length(rec.init) == 1 | length(rec.init) < n) rec.init = rep(rec.init, n)
	
	mod = fitlist@fit[[1]]@model$modeldesc$vmodel
	if(mod=="realGARCH"){
		realVol = list(...)$realizedVol
		if(is.null(realVol)) stop("\nmultifilter-->error: realGARCH model requires realizedVol xts matrix.")
		if(!is.xts(realVol)) stop("\nmultifilter-->error: realizedVol must be an xts matrix.")
		if(ncol(realVol)!=n) stop("\nmultifilter-->error: realizedVol must have the same number of columns as the data.")
		if(nrow(realVol)!=nrow(data)) stop("\nmultifilter-->error: realizedVol must have the same number of rows as the data.")
	}
	if(mod=="mcsGARCH"){
		DailyVar = list(...)$DailyVar
		if(is.null(DailyVar)) stop("\nmultifilter-->error: mcsGARCH model requires DailyVar xts matrix.")
		if(!is.xts(DailyVar)) stop("\nmultifilter-->error: DailyVar must be an xts matrix.")
		if(ncol(DailyVar)!=n) stop("\nmultifilter-->error: DailyVar must have the same number of columns as the data.")
	}
	specx = vector(mode = "list", length = n)
	for(i in 1:n){
		specx[[i]] = getspec(fitlist@fit[[i]])
		specx[[i]]@model$fixed.pars = as.list(coef(fitlist@fit[[i]]))
	}
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("specx", "data", "out.sample", "n.old", "rec.init"), envir = environment())
		if(mod=="realGARCH"){
			clusterExport(cluster, "realVol", envir = environment())
			filterlist = parLapply(cluster, as.list(1:n), fun = function(i){
						ugarchfilter(data = data[, i, drop = FALSE], spec = specx[[i]], 
								out.sample =  out.sample[i], n.old = n.old, rec.init = rec.init[i],
								realizedVol = realVol[,i])
					})
		} else if(mod=="mcsGARCH"){
			filterlist = parLapply(cluster, as.list(1:n), fun = function(i){
						ugarchfilter(data = data[, i, drop = FALSE], spec = specx[[i]], 
								out.sample =  out.sample[i], n.old = n.old, rec.init = rec.init[i],
								DailyVar = DailyVar[,i])
					})
		} else{
			filterlist = parLapply(cluster, as.list(1:n), fun = function(i){
						ugarchfilter(data = data[, i, drop = FALSE], spec = specx[[i]], 
								out.sample =  out.sample[i], n.old = n.old, rec.init = rec.init[i])
					})
		}
	} else{
		if(mod=="realGARCH"){
			filterlist = lapply(as.list(1:n), FUN = function(i){
						ugarchfilter(data = data[, i, drop = FALSE], spec = specx[[i]], 
								out.sample =  out.sample[i], n.old = n.old, rec.init = rec.init[i],
								realizedVol = realVol[,i])
					})
			
		} else if(mod=="mcsGARCH"){
			filterlist = lapply(as.list(1:n), FUN = function(i){
						ugarchfilter(data = data[, i, drop = FALSE], spec = specx[[i]], 
								out.sample =  out.sample[i], n.old = n.old, rec.init = rec.init[i],
								DailyVar = DailyVar[,i])
					})
			
		} else{
			filterlist = lapply(as.list(1:n), FUN = function(i){
						ugarchfilter(data = data[, i, drop = FALSE], spec = specx[[i]], 
								out.sample =  out.sample[i], n.old = n.old, rec.init = rec.init[i])
					})
		}
	}
	desc = list()
	desc$type = "equal"
	desc$asset.names = asset.names
	ans = new("uGARCHmultifilter",
			filter = filterlist,
			desc = desc)
	return(ans)
}

.multifiltergarch2 = function(multifitORspec, data = NULL, out.sample = 0, 
		n.old = NULL, rec.init = 'all', cluster = NULL, ...)
{
	speclist = multifitORspec
	n = length(speclist@spec)
	if(is.null(data)) 
		stop("\nmultifilter GARCH-->error: multifilter requires a data object", call. = FALSE)
	if(!is.xts(data) & !is.matrix(data) & !is.data.frame(data)) 
		stop("\nmultifilter GARCH-->error: multifilter only supports xts, matrix or data.frame objects for the data", call. = FALSE)
	if(dim(data)[2] != n)
		stop("\nmultifilter GARCH-->error: multispec length not equal to data length", call. = FALSE)
	asset.names = colnames(data)
	filterlist = vector(mode = "list", length = n)
	if(length(out.sample) == 1 | length(out.sample) < n) out.sample = rep(out.sample, n)
	if(length(rec.init) == 1 | length(rec.init) < n) rec.init = rep(rec.init, n)
	if(speclist@type=="equal"){
		mod = speclist@spec[[1]]@model$modeldesc$vmodel
		if(mod=="realGARCH"){
			realVol = list(...)$realizedVol
			if(is.null(realVol)) stop("\nmultifilter-->error: realGARCH model requires realizedVol xts matrix.")
			if(!is.xts(realVol)) stop("\nmultifilter-->error: realizedVol must be an xts matrix.")
			if(ncol(realVol)!=n) stop("\nmultifilter-->error: realizedVol must have the same number of columns as the data.")
			if(nrow(realVol)!=nrow(data)) stop("\nmultifilter-->error: realizedVol must have the same number of rows as the data.")
		}
		if(mod=="mcsGARCH"){
			DailyVar = list(...)$DailyVar
			if(is.null(DailyVar)) stop("\nmultifilter-->error: mcsGARCH model requires DailyVar xts matrix.")
			if(!is.xts(DailyVar)) stop("\nmultifilter-->error: DailyVar must be an xts matrix.")
			if(ncol(DailyVar)!=n) stop("\nmultifilter-->error: DailyVar must have the same number of columns as the data.")
		}
	} else{
		mod = "X"
	}
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("speclist", "data", "out.sample", "n.old", "rec.init"), envir = environment())
		if(mod=="realGARCH"){
			clusterExport(cluster, "realVol", envir = environment())
			filterlist = parLapply(cluster, as.list(1:n), fun = function(i){
						ugarchfilter(data = data[, i, drop = FALSE], spec = speclist@spec[[i]], 
								out.sample =  out.sample[i], n.old = n.old, rec.init = rec.init[i], 
								realizedVol = realVol[,i])
					})
		} else if(mod=="mcsGARCH"){
			clusterExport(cluster, "DailyVar", envir = environment())
			filterlist = parLapply(cluster, as.list(1:n), fun = function(i){
						ugarchfilter(data = data[, i, drop = FALSE], spec = speclist@spec[[i]], 
								out.sample =  out.sample[i], n.old = n.old, rec.init = rec.init[i],
								DailyVar = DailyVar[,i])
					})
		} else{
			filterlist = parLapply(cluster, as.list(1:n), fun = function(i){
						ugarchfilter(data = data[, i, drop = FALSE], spec = speclist@spec[[i]], 
								out.sample =  out.sample[i], n.old = n.old, rec.init = rec.init[i])
					})
		}
	} else{
		if(mod=="realGARCH"){
			filterlist = lapply(as.list(1:n), FUN = function(i){
						ugarchfilter(data = data[, i, drop = FALSE], spec = speclist@spec[[i]], 
								out.sample =  out.sample[i], n.old = n.old, rec.init = rec.init[i],
								realizedVol = realVol[,i])
					})
		} else if(mod=="mcsGARCH"){
			filterlist = lapply(as.list(1:n), FUN = function(i){
						ugarchfilter(data = data[, i, drop = FALSE], spec = speclist@spec[[i]], 
								out.sample =  out.sample[i], n.old = n.old, rec.init = rec.init[i],
								DailyVar = DailyVar[,i])
					})
		} else{
			filterlist = lapply(as.list(1:n), FUN = function(i){
						ugarchfilter(data = data[, i, drop = FALSE], spec = speclist@spec[[i]], 
								out.sample =  out.sample[i], n.old = n.old, rec.init = rec.init[i])
					})	
		}
	}
	# converged: print
	desc = list()
	desc$type = speclist@type
	desc$asset.names = asset.names
	
	ans = new("uGARCHmultifilter",
			filter  = filterlist,
			desc = desc)
	return(ans)
}

.multiforecastgarch1 = function(multifitORspec, data = NULL, n.ahead = 1, 
		n.roll = 0, out.sample = 0, external.forecasts = list(mregfor = NULL, vregfor = NULL), 
		cluster = NULL, ...)
{
	# only need to account for mcsGARCH and only partially
	multifit = multifitORspec
	n = length(multifit@fit)
	asset.names = multifit@desc$asset.names
	forecastlist = vector(mode = "list", length = n)
	
	mod = multifit@fit[[1]]@model$modeldesc$vmodel
	if(mod=="mcsGARCH"){
		DailyVar = list(...)$DailyVar
		if(is.null(DailyVar)) includeD = FALSE else includeD = TRUE
	}
	
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("multifit", "n.ahead", "n.roll", "external.forecasts"), envir = environment())
		if(mod=="mcsGARCH"){
			if(includeD) clusterExport(cluster, c("includeD","DailyVar"), envir = environment())
			forecastlist = parLapply(cluster, as.list(1:n), fun = function(i){
						ugarchforecast(fitORspec = multifit@fit[[i]], data = NULL, n.ahead = n.ahead, 
								n.roll = n.roll, external.forecasts = external.forecasts,
								if(includeD) DailyVar = DailyVar[,i] else DailyVar = NULL)
					})
		} else{
			forecastlist = parLapply(cluster, as.list(1:n), fun = function(i){
						ugarchforecast(fitORspec = multifit@fit[[i]], data = NULL, n.ahead = n.ahead, 
								n.roll = n.roll, external.forecasts = external.forecasts)
					})		
		}
		
	} else{
		if(mod=="mcsGARCH"){
			forecastlist = lapply(as.list(1:n), FUN = function(i){
						ugarchforecast(fitORspec = multifit@fit[[i]], data = NULL, n.ahead = n.ahead, 
								n.roll = n.roll, external.forecasts = external.forecasts,
								if(includeD) DailyVar = DailyVar[,i] else DailyVar = NULL)
					})
		} else{
			forecastlist = lapply(as.list(1:n), FUN = function(i){
						ugarchforecast(fitORspec = multifit@fit[[i]], data = NULL, n.ahead = n.ahead, 
								n.roll = n.roll, external.forecasts = external.forecasts)
					})
		}
	}
	desc = list()
	desc$type = "equal"
	desc$asset.names = asset.names
	ans = new("uGARCHmultiforecast",
			forecast  = forecastlist,
			desc = desc)
	return(ans)
}

.multiforecastgarch2 = function(multifitORspec, data = NULL, n.ahead = 1, n.roll = 0, out.sample = 0,
		external.forecasts = list(mregfor = NULL, vregfor = NULL), cluster =NULL, ...)
{
	multispec = multifitORspec
	n = length(multispec@spec)
	if(is.null(data)) 
		stop("\nmultiforecast GARCH-->error: multiforecast with multiple spec requires a data object", call. = FALSE)
	if(!is.xts(data) & !is.matrix(data) & !is.data.frame(data)) 
		stop("\nmultiforecast GARCH-->error: multiforecast only supports xts, matrix or data.frame objects for the data", call. = FALSE)
	if(dim(data)[2] != n)
		stop("\nmultiforecast GARCH-->error: multispec length not equal to data length", call. = FALSE)
	asset.names = colnames(data)
	forecastlist = vector(mode = "list", length = n)
	if(is.null(out.sample)) out.sample = 0
	if(length(out.sample) == 1) out.sample = rep(out.sample, n)
	if(length(out.sample) !=n ) 
		stop("\nmultiforecast GARCH-->error: out.sample length not equal to data length", call. = FALSE)
	
	if(multispec@type=="equal"){
		mod = multispec@spec[[1]]@model$modeldesc$vmodel
		if(mod=="realGARCH"){
			realVol = list(...)$realizedVol
			if(is.null(realVol)) stop("\nmultiforecast GARCH-->error: realGARCH model requires realizedVol xts matrix.")
			if(!is.xts(realVol)) stop("\nmultiforecast GARCH-->error: realizedVol must be an xts matrix.")
			if(ncol(realVol)!=n) stop("\nmultiforecast GARCH-->error: realizedVol must have the same number of columns as the data.")
			if(nrow(realVol)!=nrow(data)) stop("\nmultiforecast GARCH-->error: realizedVol must have the same number of rows as the data.")
		}
		if(mod=="mcsGARCH"){
			stop("\nugarchforecast (and multiforecast) with specification object not available for mcsGARCH model")
		}
	} else{
		mod = "X"
	}
	
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("multispec", "data", "n.ahead", "n.roll", 
						"out.sample", "external.forecasts"), envir = environment())
		if(mod=="realGARCH"){
			clusterExport(cluster, "realVol", envir = environment())
			forecastlist = parLapply(cluster, as.list(1:n), fun = function(i){
						ugarchforecast(fitORspec = multispec@spec[[i]], 
								data = data[, i, drop = FALSE], n.ahead = n.ahead, n.roll = n.roll, 
								out.sample = out.sample[i], external.forecasts = external.forecasts,
								realizedVol = realVol[,i])
					})
		} else{
			forecastlist = parLapply(cluster, as.list(1:n), fun = function(i){
						ugarchforecast(fitORspec = multispec@spec[[i]], 
								data = data[, i, drop = FALSE], n.ahead = n.ahead, n.roll = n.roll, 
								out.sample = out.sample[i], external.forecasts = external.forecasts)
					})			
		}
		
	} else{
		if(mod=="realGARCH"){
			forecastlist = lapply(as.list(1:n), FUN = function(i){
						ugarchforecast(fitORspec = multispec@spec[[i]], 
								data = data[, i, drop = FALSE], n.ahead = n.ahead, n.roll = n.roll, 
								out.sample = out.sample[i], external.forecasts = external.forecasts, 
								realizedVol = realVol[,i])
					})
		} else{
			forecastlist = lapply(as.list(1:n), FUN = function(i){
						ugarchforecast(fitORspec = multispec@spec[[i]], 
								data = data[, i, drop = FALSE], n.ahead = n.ahead, n.roll = n.roll, 
								out.sample = out.sample[i], external.forecasts = external.forecasts, ...)
					})	
		}
		
	}
	desc = list()
	desc$type = multispec@type
	desc$asset.names = asset.names
	
	ans = new("uGARCHmultiforecast",
			forecast  = forecastlist,
			desc = desc)
	return(ans)
}