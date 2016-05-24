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

# a multifit function possible utilizing parallel execution returning a fitlist
# object
.multifitarfima = function(multispec, data, out.sample = 0, solver = "solnp", 
		solver.control = list(), fit.control = list(fixed.se = 0, scale = 0), 
		cluster = NULL, ...)
{
	n = length(multispec@spec)
	if(is.null(data)) 
		stop("\nmultifit ARFIMA-->error: multifit requires a data object", call. = FALSE)
	if(!is.matrix(data) & !is.data.frame(data)) 
		stop("\nmultifit ARFIMA-->error: multifit only supports matrix or data.frame objects for the data", call. = FALSE)
	if(is.matrix(data)) data = as.data.frame(data)
	asset.names = colnames(data)
	if(dim(data)[2] != n) 
		stop("\nmultifit ARFIMA-->error: speclist length not equal to data length", call. = FALSE)
	fitlist = vector(mode = "list", length = n)
	if(length(out.sample) == 1 | length(out.sample) < n) out.sample = rep(out.sample, n)
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("multispec", "data", "out.sample", "solver", 
						"solver.control", "fit.control"), envir = environment())
		fitlist = parLapply(cluster, as.list(1:n), fun = function(i){
					arfimafit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
							out.sample = out.sample[i], solver = solver, 
							solver.control = solver.control, fit.control = fit.control)
				})
	} else{
		fitlist = lapply(as.list(1:n), FUN = function(i){
					arfimafit(spec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
							out.sample = out.sample[i], solver = solver, 
							solver.control = solver.control, fit.control = fit.control)
				})
	}
	# converged: print
	desc = list()
	desc$type = multispec@type
	desc$asset.names = asset.names
	ans = new("ARFIMAmultifit",
			fit  = fitlist,
			desc = desc)
	return(ans)
}

.multifilterarfima1 = function(multifitORspec, data = NULL, out.sample = 0, 
		n.old = NULL, rec.init = "all", cluster = NULL, ...)
{
	fitlist = multifitORspec
	n = length(fitlist@fit)
	if(is.null(data)) 
		stop("\nmultifilter ARFIMA-->error: multifilter requires a data object", call. = FALSE)
	if(!is.matrix(data) & !is.data.frame(data)) 
		stop("\nmultifilter ARFIMA-->error: multifilter only supports matrix or data.frame objects for the data", call. = FALSE)
	if(dim(data)[2] != n)
		stop("\nmultifilter ARFIMA-->error: fitlist length not equal to data length", call. = FALSE)
	if(is.matrix(data)) data = as.data.frame(data)
	asset.names = colnames(data)
	filterlist = vector(mode = "list", length = n)
	if(length(out.sample) == 1 | length(out.sample) < n) out.sample = rep(out.sample, n)
	
	specx = vector(mode = "list", length = n)
	for(i in 1:n){
		specx[[i]] = getspec(fitlist@fit[[i]])
		specx[[i]]@model$fixed.pars = as.list(coef(fitlist@fit[[i]]))
	}
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("specx", "data", "out.sample", "n.old"), envir = environment())
		filterlist = parLapply(cluster, as.list(1:n), fun = function(i){
					arfimafilter(data = data[, i, drop = FALSE], spec = specx[[i]], 
							out.sample =  out.sample[i], n.old = n.old)
				})
	} else{
		filterlist = lapply(as.list(1:n), FUN = function(i){
					arfimafilter(data = data[, i, drop = FALSE], spec = specx[[i]], 
							out.sample =  out.sample[i], n.old = n.old)
				})
	}
	
	desc = list()
	desc$type = "equal"
	desc$asset.names = asset.names
	
	ans = new("ARFIMAmultifilter",
			filter = filterlist,
			desc = desc)
	return(ans)
}

.multifilterarfima2 = function(multifitORspec, data = NULL, out.sample = 0, 
		n.old = NULL, rec.init = "all", cluster = NULL, ...)
{
	speclist = multifitORspec
	n = length(speclist@spec)
	if(is.null(data)) 
		stop("\nmultifilter ARFIMA-->error: multifilter requires a data object", call. = FALSE)
	if(!is.matrix(data) & !is.data.frame(data)) 
		stop("\nmultifilter ARFIMA-->error: multifilter only supports matrix or data.frame objects for the data", call. = FALSE)
	if(dim(data)[2] != n)
		stop("\nmultifilter ARFIMA-->error: multispec length not equal to data length", call. = FALSE)
	if(is.matrix(data)) data = as.data.frame(data)
	asset.names = colnames(data)
	filterlist = vector(mode = "list", length = n)
	if(length(out.sample) == 1 | length(out.sample) < n) out.sample = rep(out.sample, n)
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("speclist", "data", "out.sample", "n.old"), envir = environment())
		filterlist = parLapply(cluster, as.list(1:n), fun = function(i){
					arfimafilter(spec = speclist@spec[[i]], data = data[, i, drop = FALSE], 
							out.sample = out.sample[i], n.old = n.old)
				})
	} else{
		filterlist = lapply(as.list(1:n), FUN = function(i){
					arfimafilter(spec = speclist@spec[[i]], data = data[, i, drop = FALSE], 
							out.sample = out.sample[i], n.old = n.old)
				})
	}
	
	# converged: print
	desc = list()
	desc$type = speclist@type
	desc$asset.names = asset.names
	
	ans = new("ARFIMAmultifilter",
			filter  = filterlist,
			desc = desc)
	return(ans)
}

.multiforecastarfima1 = function(multifitORspec, data = NULL, n.ahead = 1, n.roll = 0, out.sample = 0,
		external.forecasts = list(mregfor = NULL), cluster = NULL, ...)
{
	multifit = multifitORspec
	n = length(multifit@fit)
	asset.names = multifit@desc$asset.names
	forecastlist = vector(mode = "list", length = n)
	
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("multifit", "n.ahead", "n.roll", "external.forecasts"), envir = environment())
		forecastlist = parLapply(cluster, as.list(1:n), fun = function(i){
					arfimaforecast(fitORspec = multifit@fit[[i]], data = NULL, 
							n.ahead = n.ahead, n.roll = n.roll, external.forecasts = external.forecasts)
				})
	} else{
		forecastlist = lapply(as.list(1:n), FUN = function(i){
					arfimaforecast(fitORspec = multifit@fit[[i]], data = NULL,
							n.ahead = n.ahead, n.roll = n.roll, external.forecasts = external.forecasts)
				})
	}
	
	desc = list()
	desc$type = "equal"
	desc$asset.names = asset.names
	ans = new("ARFIMAmultiforecast",
			forecast  = forecastlist,
			desc = desc)
	return(ans)
}

.multiforecastarfima2 = function(multifitORspec, data = NULL, n.ahead = 1, n.roll = 0, out.sample = 0,
		external.forecasts = list(mregfor = NULL), cluster = NULL, ...)
{
	multispec = multifitORspec
	n = length(multispec@spec)
	if(is.null(data)) 
		stop("\nmultiforecast ARFIMA-->error: multiforecast with multiple spec requires a data object", call. = FALSE)
	if(!is.matrix(data) & !is.data.frame(data)) 
		stop("\nmultiforecast GARCH-->error: multiforecast only supports matrix or data.frame objects for the data", call. = FALSE)
	if(is.matrix(data)) data = as.data.frame(data)
	if(dim(data)[2] != n)
		stop("\nmultiforecast ARFIMA-->error: multispec length not equal to data length", call. = FALSE)
	asset.names = colnames(data)
	forecastlist = vector(mode = "list", length = n)
	if(is.null(out.sample)) out.sample = 0
	if(length(out.sample) == 1) out.sample = rep(out.sample, n)
	if(length(out.sample) !=n ) stop("\nmultiforecast ARFIMA-->error: out.sample length not equal to data length", call. = FALSE)
	
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, library(rugarch))
		clusterExport(cluster, c("multispec", "data", "n.ahead", "n.roll", 
						"out.sample", "external.forecasts"), envir = environment())
		forecastlist = parLapply(cluster, as.list(1:n), fun = function(i){
					arfimaforecast(fitORspec = multispec@spec[[i]], data = data[, i, drop = FALSE],
								n.ahead = n.ahead, n.roll = n.roll, out.sample = out.sample[i],
								external.forecasts = external.forecasts)
					})
	} else{
		forecastlist = lapply(as.list(1:n), FUN = function(i){
					arfimaforecast(fitORspec = multispec@spec[[i]], data = data[, i, drop = FALSE], 
							n.ahead = n.ahead, n.roll = n.roll, out.sample = out.sample[i], 
							external.forecasts = external.forecasts)
				})
	}
	
	desc = list()
	desc$type = multispec@type
	desc$asset.names = asset.names
	
	ans = new("ARFIMAmultiforecast",
			forecast  = forecastlist,
			desc = desc)
	return(ans)
}