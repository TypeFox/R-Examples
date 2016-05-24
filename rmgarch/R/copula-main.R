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

#################################################################################
# Fit, Filter, Forecast and Simulation
#-------------------------------------
.cgarchspec = function(uspec, VAR = FALSE, robust = FALSE, lag = 1, lag.max = NULL, 
				lag.criterion = c("AIC", "HQ", "SC", "FPE"), 
				external.regressors = NULL, robust.control = list("gamma" = 0.25, 
						"delta" = 0.01, "nc" = 10, "ns" = 500), dccOrder = c(1,1), 
				asymmetric = FALSE, 
				distribution.model = list(copula = c("mvnorm", "mvt"), 
						method = c("Kendall", "ML"), time.varying = FALSE, 
						transformation = c("parametric", "empirical", "spd")),
		start.pars = list(), fixed.pars = list())
{
	.eps = .Machine$double.eps
	VAR.opt = list()
	if(is.null(VAR)){
		VAR.opt$VAR = FALSE 
	} else{
		VAR.opt$VAR = as.logical(VAR)
	}
	if(is.null(robust)){
		VAR.opt$robust = FALSE 
	} else{
		VAR.opt$robust = as.logical(robust)
	}
	if(is.null(lag)){
		VAR.opt$lag = 1
	} else{
		VAR.opt$lag = as.integer(lag)
	}
	if(is.null(lag.max)){
		VAR.opt$lag.max = NULL 
	} else{
		VAR.opt$lag.max = as.integer(min(1, lag.max))
	}
	if(is.null(lag.criterion)){
		VAR.opt$lag.criterion = "AIC" 
	} else{
		VAR.opt$lag.criterion = lag.criterion[1]
	}
	if(is.null(external.regressors)){
		VAR.opt$external.regressors = NULL 
	} else{
		VAR.opt$external.regressors = external.regressors
	}
	rc = list("gamma" = 0.25, "delta" = 0.01, "nc" = 10, "ns" = 500)
	rcmatch = match(names(robust.control), c("gamma", "delta", "nc", "ns"))
	if(length(rcmatch[!is.na(rcmatch)]) > 0){
		rx = which(!is.na(rcmatch))
		rc[rcmatch[!is.na(rcmatch)]] = robust.control[rx]
	}
	VAR.opt$robust.control = rc
	
	modeldata = list()
	modeldesc = list()
	m = length(uspec@spec)
	modelinc = rep(0, 11)
	names(modelinc) = c("var", "mvmxreg", "C", "dcca", "dccb", "dccg", "mshape", 
			"mskew", "aux", "aux", "aux")
	
	if(is.null(distribution.model$copula)){
		distribution = "mvnorm" 
	} else{
		distribution = tolower(distribution.model$copula)
	}
	distribution = distribution[1]
	valid.distributions = c("mvnorm", "mvt")
	if(!any(distribution == valid.distributions)) 
		stop("\nInvalid Copula Distribution Choice\n", call. = FALSE)
	modeldesc$distribution = distribution
	if(distribution == "mvt") modelinc[7] = 1
	if(is.null(distribution.model$method)){
		method = "kendall"
	} else{
		method = tolower(distribution.model$method)
	}
	method = method[1]
	valid.methods = c("kendall", "ml")
	if(!any(method == valid.methods)){
		if(distribution == "mvt") 
			warning("\nInvalid Rho Method Estimation Choice\n", call. = FALSE)
		method = "kendall"
	}
	modeldesc$cor.method = toupper(method)
	
	if(is.null(distribution.model$time.varying)){
		timecopula = FALSE
	} else{
		timecopula = as.logical(distribution.model$time.varying)
	}
	timecopula = timecopula[1]
	modeldesc$timecopula = timecopula
	if(!timecopula && tolower(method[1]) == "ml") modelinc[3] = ((m*m - m)/2)
	
	
	if(is.null(distribution.model$transformation)){
		transformation = "parametric" 
	} else{
		transformation = tolower(distribution.model$transformation)
	}
	transformation = transformation[1]
	valid.transformations = c("parametric", "empirical", "spd")
	if(!any(transformation == valid.transformations)) 
		stop("\nInvalid Copula Transformation Choice\n", call. = FALSE)
	modeldesc$transformation = transformation
	
	if(timecopula){
		if(is.null(dccOrder)){
			modelinc[4:5] = 1
		} else{
			modelinc[4] = as.integer(dccOrder[1])
			modelinc[5] = as.integer(dccOrder[2])
		}
		if(is.null(asymmetric)){
			modelinc[6] = 0
		} else{
			if(as.logical(asymmetric)) modelinc[6] = modelinc[4]
		}
	} else{
		modelinc[4:6] = 0
	}
	
	if( VAR ){
		if(is.null(VAR.opt$lag)){
			modelinc[1] = 1 
		} else{
			modelinc[1] = as.integer( VAR.opt$lag )
		}
		if(!is.null(VAR.opt$external.regressors)){
			if(!is.matrix(VAR.opt$external.regressors)) 
				stop("\nexternal.regressors must be a matrix.")
			modelinc[2] = dim(VAR.opt$external.regressors)[2]
			modeldata$mexdata = VAR.opt$external.regressors
		} else{
			modeldata$mexdata = NULL
		}
	}
	varmodel = list()
	umodel = vector(mode ="list")
	if( modelinc[1]>0 ){
		varmodel$robust = VAR.opt$robust
		varmodel$lag.max = VAR.opt$lag.max
		varmodel$lag.criterion = VAR.opt$lag.criterion
		varmodel$robust.control = VAR.opt$robust.control
		umodel$modelinc = matrix(0, ncol = m, nrow = 21)
		rownames(umodel$modelinc) = names(uspec@spec[[1]]@model$modelinc[1:21])
		umodel$modeldesc = list()
		umodel$vt = sapply(uspec@spec, FUN = function(x) x@model$modelinc[22])
		umodel$modeldesc$vmodel = vector(mode = "character", length = m)
		umodel$modeldesc$vsubmodel = vector(mode = "character", length = m)
		umodel$start.pars = umodel$fixed.pars = vector(mode = "list", length = m)
		umodel$modeldesc$distribution = vector(mode = "character", length = m)
		umodel$modeldata = list()
		umodel$modeldata$vexdata = vector(mode = "list", length = m)
		for(i in 1:m){
			# zero the mean equation since we are using VAR
			umodel$modeldesc$vmodel[i] = uspec@spec[[i]]@model$modeldesc$vmodel
			umodel$modeldesc$vsubmodel[i] = ifelse(is.null(uspec@spec[[i]]@model$modeldesc$vsubmodel), 
					"GARCH", uspec@spec[[i]]@model$modeldesc$vsubmodel)
			umodel$modeldesc$distribution[i] = uspec@spec[[i]]@model$modeldesc$distribution
			umodel$modelinc[,i] = uspec@spec[[i]]@model$modelinc[1:21]
			umodel$modelinc[1:6,i] = 0
			umodel$modeldata$vexdata[[i]] = if(is.null(uspec@spec[[i]]@model$modeldata$vexdata)) NA else uspec@spec[[i]]@model$modeldata$vexdata
			umodel$start.pars[[i]] = if(is.null(uspec@spec[[i]]@model$start.pars)) NA else uspec@spec[[i]]@model$start.pars
			umodel$fixed.pars[[i]] = if(is.null(uspec@spec[[i]]@model$fixed.pars)) NA else uspec@spec[[i]]@model$fixed.pars
			umodel$modeldata$mexdata[[i]] = NA
		}
	} else{
		varmodel$lag.max = 1
		varmodel$lag.criterion = "HQ"
		umodel$modelinc = matrix(0, ncol = m, nrow = 21)
		rownames(umodel$modelinc) = names(uspec@spec[[1]]@model$modelinc[1:21])
		umodel$modeldesc = list()
		# variance targeting custom
		umodel$vt = sapply(uspec@spec, FUN = function(x) x@model$modelinc[22])
		umodel$modeldesc$vmodel = vector(mode = "character", length = m)
		umodel$modeldesc$vsubmodel = vector(mode = "character", length = m)
		umodel$start.pars = umodel$fixed.pars = vector(mode = "list", length = m)
		umodel$modeldesc$distribution = vector(mode = "character", length = m)
		umodel$modeldata = list()
		umodel$modeldata$mexdata = vector(mode = "list", length = m)
		umodel$modeldata$vexdata = vector(mode = "list", length = m)
		for(i in 1:m){
			umodel$modeldesc$vmodel[i] = uspec@spec[[i]]@model$modeldesc$vmodel
			umodel$modeldesc$vsubmodel[i] = ifelse(is.null(uspec@spec[[i]]@model$modeldesc$vsubmodel),"GARCH",uspec@spec[[i]]@model$modeldesc$vsubmodel)
			umodel$modeldesc$distribution[i] = uspec@spec[[i]]@model$modeldesc$distribution
			umodel$modelinc[,i] = uspec@spec[[i]]@model$modelinc[1:21]
			umodel$modeldata$mexdata[[i]] = if(is.null(uspec@spec[[i]]@model$modeldata$mexdata)) NA else uspec@spec[[i]]@model$modeldata$mexdata
			umodel$modeldata$vexdata[[i]] = if(is.null(uspec@spec[[i]]@model$modeldata$vexdata)) NA else uspec@spec[[i]]@model$modeldata$vexdata
			umodel$start.pars[[i]] = if(is.null(uspec@spec[[i]]@model$start.pars)) NA else uspec@spec[[i]]@model$start.pars
			umodel$fixed.pars[[i]] = if(is.null(uspec@spec[[i]]@model$fixed.pars)) NA else uspec@spec[[i]]@model$fixed.pars
		}
	}
	modelinc[11] = which(c("mvnorm", "mvt") == distribution)
	
	maxdccOrder = max(dccOrder)
	maxgarchOrder =  max( sapply(uspec@spec, FUN = function(x) max(max(x@model$modelinc[2:3]), max(x@model$modelinc[8:9]) ) ) )
	maxOrder = max( maxdccOrder, maxgarchOrder )
	if(modelinc[1]>0){
		maxgarchOrder = max(c(maxgarchOrder, modelinc[1]))
	}
	modeldesc$dccmodel = ifelse(modelinc[6]>0, "ADCC", "DCC")
	
	pars = matrix(0, ncol = 6, nrow = 6)
	colnames(pars) = c("Level", "Fixed", "Include", "Estimate", "LB", "UB")
	pidx = matrix(NA, nrow = 6, ncol = 2)
	colnames(pidx) = c("begin", "end")
	rownames(pidx) =  c("C", "dcca", "dccb", "dccg", "mshape", "mskew")
	
	pos = 1
	pos.matrix = matrix(0, ncol = 3, nrow = 6)
	colnames(pos.matrix) = c("start", "stop", "include")
	rownames(pos.matrix) = c("C", "dcca", "dccb", "dccg", "mshape", "mskew")
	
	for(i in 1:6){
		if( modelinc[2+i] > 0 ){
			pos.matrix[i,1:3] = c(pos, pos+modelinc[2+i]-1, 1)
			pos = max(pos.matrix[1:i,2]+1)
		}
	}
	
	mm = sum(modelinc[3:8])
	mm = mm - length( which(modelinc[c(3:8)]>0) )
	pars = matrix(0, ncol = 6, nrow = 6 + mm)
	colnames(pars) = c("Level", "Fixed", "Include", "Estimate", "LB", "UB")
	pidx = matrix(NA, nrow = 6, ncol = 2)
	colnames(pidx) = c("begin", "end")
	rownames(pidx) =  c("C", "dcca", "dccb", "dccg", "mshape", "mskew")
	fixed.names = names(fixed.pars)
	start.names = names(start.pars)
	fixed.pars = unlist(fixed.pars)
	start.pars = unlist(start.pars)
	
	pnames = NULL
	nx = 0
	pn = 1
	pidx[1,1] = 1
	if(pos.matrix[1,3] == 1){
		cvm = 1:((m*m - m)/2)
		pn = length( seq(pos.matrix[1,1], pos.matrix[1,2], by = 1) )
		for(i in 1:pn){
			nnx = paste("C", cvm[i], sep="")
			pars[(nx+i), 1] = 0
			if(any(substr(start.names, 1, nchar(nnx))==nnx)){
				nix = which(start.names == nnx)
				pars[(nx+i), 1] = start.pars[nix]
			}
			pars[(nx+i), 3] = 1
			pars[(nx+i), 5] = -0.99
			pars[(nx+i), 6] =  0.99
			if(any(substr(fixed.names, 1, nchar(nnx))==nnx)){
				nix = which(fixed.names == nnx)
				pars[(nx+i), 1] = fixed.pars[nix]
				pars[(nx+i), 2] = 1
			} else{
				pars[(nx+i), 4] = 1
			}
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "C")
	}
	pidx[1,2] = pn
	nx = pn
	pn = 1
	pidx[2,1] =  nx+1
	
	if(pos.matrix[2,3] == 1){
		pn = length( seq(pos.matrix[2,1], pos.matrix[2,2], by = 1) )
		for(i in 1:pn){
			nnx = paste("dcca", i, sep="")
			pars[(nx+i), 1] = 0.05/pn
			if(any(substr(start.names, 1, nchar(nnx))==nnx)){
				nix = which(start.names == nnx)
				pars[(nx+i), 1] = start.pars[nix]
			}
			pars[(nx+i), 3] = 1
			pars[(nx+i), 5] = .eps
			pars[(nx+i), 6] = 1
			if(any(substr(fixed.names, 1, nchar(nnx))==nnx)){
				nix = which(fixed.names == nnx)
				pars[(nx+i), 1] = fixed.pars[nix]
				pars[(nx+i), 2] = 1
			} else{
				pars[(nx+i), 4] = 1
			}
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "dcca")
	}
	pidx[2,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[3,1] = nx+1
	
	if(pos.matrix[3,3] == 1){
		pn = length( seq(pos.matrix[3,1], pos.matrix[3,2], by = 1) )
		for(i in 1:pn){
			nnx = paste("dccb", i, sep="")
			pars[(nx+i), 1] = 0.9/pn
			if(any(substr(start.names, 1, nchar(nnx))==nnx)){
				nix = which(start.names == nnx)
				pars[(nx+i), 1] = start.pars[nix]
			}			
			pars[(nx+i), 3] = 1
			pars[(nx+i), 5] = .eps
			pars[(nx+i), 6] = 1
			if(any(substr(fixed.names, 1, nchar(nnx))==nnx)){
				nix = which(fixed.names == nnx)
				pars[(nx+i), 1] = fixed.pars[nix]
				pars[(nx+i), 2] = 1
			} else{
				pars[(nx+i), 4] = 1
			}
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "dccb")
	}
	pidx[3,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[4,1] = nx+1
	
	if(pos.matrix[4,3] == 1){
		pn = length( seq(pos.matrix[4,1], pos.matrix[4,2], by = 1) )
		for(i in 1:pn){
			nnx = paste("dccg", i, sep="")
			pars[(nx+i), 1] = 0.05/pn
			if(any(substr(start.names, 1, nchar(nnx))==nnx)){
				nix = which(start.names == nnx)
				pars[(nx+i), 1] = start.pars[nix]
			}		
			pars[(nx+i), 3] = 1
			pars[(nx+i), 5] = .eps
			pars[(nx+i), 6] = 1
			if(any(substr(fixed.names, 1, nchar(nnx))==nnx)){
				nix = which(fixed.names == nnx)
				pars[(nx+i), 1] = fixed.pars[nix]
				pars[(nx+i), 2] = 1
			} else{
				pars[(nx+i), 4] = 1
			}
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "dccg")
	}
	pidx[4,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[5,1] = nx+1

	if(modelinc[7]<=1){
		if(pos.matrix[5,3]==1){
			pars[nx+pn, 3] = 1
			pars[nx+pn, 1] = 5
			pars[(nx+pn), 5] = 4
			pars[(nx+pn), 6] = 50
			if(any(substr(start.names, 1, 6) == "mshape")) pars[nx+pn, 1] = start.pars["mshape"]
			if(any(substr(fixed.names, 1, 6) == "mshape")){
				pars[nx+pn,2] = 1
				pars[nx+pn, 1] = fixed.pars["mshape"]
			} else{
				pars[nx+pn,4] = 1
			}
		}
		pnames = c(pnames, "mshape")
	} else{
		if(pos.matrix[5,3] == 1){
			pn = length( seq(pos.matrix[5,1], pos.matrix[5,2], by = 1) )
			for(i in 1:pn){
				nnx = paste("mshape", i, sep="")
				pars[(nx+i), 1] = 5
				if(any(substr(start.names, 1, nchar(nnx))==nnx)){
					nix = which(start.names == nnx)
					pars[(nx+i), 1] = start.pars[nix]
				}		
				pars[(nx+i), 3] = 1
				pars[(nx+i), 5] = 4
				pars[(nx+i), 6] = 50
				if(any(substr(fixed.names, 1, nchar(nnx))==nnx)){
					nix = which(fixed.names == nnx)
					pars[(nx+i), 1] = fixed.pars[nix]
					pars[(nx+i), 2] = 1
				} else{
					pars[(nx+i), 4] = 1
				}
				pnames = c(pnames, nnx)
			}
		} else{
			pnames = c(pnames, "mshape")
		}
	}
	pidx[5,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[6,1] = nx+1
	
	if(modelinc[8]<=1){
		if(pos.matrix[6,3]==1){
			pars[nx+pn, 3] = 1
			pars[nx+pn, 1] = 0.5
			pars[(nx+pn), 5] = -1
			pars[(nx+pn), 6] = 1
			if(any(substr(start.names, 1, 5) == "mskew")) pars[nx+pn, 1] = start.pars["mskew"]
			if(any(substr(fixed.names, 1, 5) == "mskew")){
				pars[nx+pn,2] = 1
				pars[nx+pn, 1] = fixed.pars["mskew"]
			} else{
				pars[nx+pn,4] = 1
			}
		}
		pnames = c(pnames, "mskew")
	} else{
		if(pos.matrix[6,3] == 1){
			pn = length( seq(pos.matrix[6,1], pos.matrix[6,2], by = 1) )
			for(i in 1:pn){
				nnx = paste("mskew", i, sep="")
				pars[(nx+i), 1] = 3
				pars[(nx+i), 5] = -1
				pars[(nx+i), 6] = 1
				if(any(substr(start.names, 1, nchar(nnx))==nnx)){
					nix = which(start.names == nnx)
					pars[(nx+i), 1] = start.pars[nix]
				}		
				if(any(substr(fixed.names, 1, nchar(nnx))==nnx)){
					nix = which(fixed.names == nnx)
					pars[(nx+i), 1] = fixed.pars[nix]
					pars[(nx+i), 2] = 1
				} else{
					pars[(nx+i), 4] = 1
				}
				pnames = c(pnames, nnx)
			}
		} else{
			pnames = c(pnames, "mskew")
		}
	}
	pidx[6,2] = nx+pn
	rownames(pars) = pnames
	
	model = list(modelinc = modelinc, modeldesc = modeldesc, modeldata = modeldata, 
			varmodel = varmodel, pars = pars, start.pars = start.pars, 
			fixed.pars = fixed.pars, maxgarchOrder = maxgarchOrder, 
			maxdccOrder = maxdccOrder, pos.matrix = pos.matrix, pidx = pidx)
	
	ans = new("cGARCHspec",
			model = model,
			umodel = umodel)
	return(ans)
}



.cgarchfit = function(spec, data, spd.control = list(lower = 0.1, upper = 0.9, 
				type = "pwm", kernel = "epanech"),  
		fit.control = list(eval.se = TRUE, stationarity = TRUE, scale = FALSE), 
		solver = "solnp", solver.control = list(), out.sample = 0, cluster = NULL, fit = NULL, 
		VAR.fit = NULL, realizedVol = NULL, ...)
{
	type = ifelse(spec@model$modeldesc$timecopula, "dynamic", "static")
	ans = switch(type, 
			dynamic = .cgarchfit.dynamic(spec = spec, data = data, spd.control = spd.control, 
					fit.control = fit.control, solver = solver, solver.control = solver.control, 
					out.sample = out.sample, cluster = cluster, fit = fit, VAR.fit = VAR.fit, 
					realizedVol = realizedVol, ...),
			static = .cgarchfit.static(spec = spec, data = data, spd.control = spd.control, 
					fit.control = fit.control, solver = solver, solver.control = solver.control, 
					out.sample = out.sample, cluster = cluster, fit = fit, VAR.fit = VAR.fit, 
					realizedVol = realizedVol, ...))
	return( ans )
}


.cgarchfit.dynamic = function(spec, data, spd.control = list(lower = 0.1, upper = 0.9, 
				type = "pwm", kernel = "epanech"), 
		fit.control = list(eval.se = TRUE, stationarity = TRUE, scale = FALSE), 
		solver = "solnp", solver.control = list(), out.sample = 0, cluster = NULL, 
		fit = NULL, VAR.fit = NULL, realizedVol = NULL, ...)
{
	tic = Sys.time()
	.eps = .Machine$double.eps
	model = spec@model
	umodel = spec@umodel
	ufit.control = list()
	if(is.null(fit.control$stationarity)){
		ufit.control$stationarity = TRUE
	} else {
		ufit.control$stationarity = fit.control$stationarity
		fit.control$stationarity = NULL
	}
	if(is.null(fit.control$scale)){
		ufit.control$scale = TRUE 
	} else{
		ufit.control$scale = fit.control$scale
		fit.control$scale = NULL
	}
	if(is.null(fit.control$eval.se)) fit.control$eval.se = TRUE
	
	if(length(solver)==2){
		garch.solver = 	solver[1]
		solver = solver[2]
	} else{
		garch.solver = solver[1]
	}
	solver = match.arg(tolower(solver)[1], c("solnp", "nlminb", "lbfgs","gosolnp"))
	#-----------------------------------------------------------------------------------
	# Data Extraction
	m = dim(data)[2]
	if( is.null( colnames(data) ) ) cnames = paste("Asset_", 1:m, sep = "") else cnames = colnames(data)
	colnames(umodel$modelinc) = cnames
	
	xdata = .extractmdata(data)
	if(!is.numeric(out.sample)) 
		stop("\ndccfit-->error: out.sample must be numeric\n")
	if(as.numeric(out.sample) < 0) 
		stop("\ndccfit-->error: out.sample must be positive\n")
	n.start = round(out.sample, 0)
	n = NROW(xdata$data)
	if( (n-n.start) < 100) 
		stop("\ndccfit-->error: function requires at least 100 data\n points to run\n")
	data  	= xdata$data
	index 	= xdata$index
	period  = xdata$period
	
	# save the data to the model spec
	model$modeldata$data = data
	model$modeldata$index = index
	model$modeldata$period = period
	T = model$modeldata$T = n - n.start
	model$modeldata$n.start = n.start
	model$modeldata$asset.names = cnames
	#-----------------------------------------------------------------------------------
	transformation = model$modeldesc$transformation
	#-----------------------------------------------------------------------------------
	# VAR model
	if( model$modelinc[1]>0 ){
		tmp = mvmean.varfit(model = model, data = data, VAR.fit = VAR.fit, T = T, 
				out.sample = out.sample, cluster = cluster)
		model = tmp$model
		zdata = tmp$zdata
		mu = tmp$mu
		varcoef = tmp$varcoef
		p = tmp$p
		N = tmp$N
	} else{
		zdata = data
		ex = NULL
	}
	T = dim(zdata)[1] - out.sample
	#-----------------------------------------------------------------------------------
	# Univariate GARCH fit
	# create a multispec list, check for variance targeting (vt)
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, umodel$vt)
	
	# Check if a pre-fitted uGARCHmultifit object provided, else estimate
	if( !is.null(fit) && is(fit, "uGARCHmultifit") ){
		# check VAR and fit:
		if(model$modelinc[1]>0){
			for(i in 1:m){
				if(sum(fit@fit[[i]]@model$modelinc[1:6])>0)
					stop("\nThe user supplied fit object has a non-null mean specification but VAR already chosen for mean filtration!!!")
			}
		}
		fitlist = fit
		if(spec@model$modelinc[1]>0) model$mu = mu else model$mu = fitted(fitlist)
		model$residuals = res = residuals(fitlist)
		model$sigma = sig = sigma(fitlist)
	} else{
		fitlist = multifit(multispec = mspec, data = xts(zdata, index), out.sample = n.start, 
				solver = garch.solver, solver.control = solver.control, 
				fit.control = ufit.control, cluster = cluster, realizedVol = realizedVol, ...)
		converge = sapply(fitlist@fit, FUN = function(x) x@fit$convergence)
		if( any( converge == 1 ) ){
			pr = which(converge != 1)
			cat("\nNon-Converged:\n")
			print(pr)
			cat("\ndccfit-->error: convergence problem in univariate fit...")
			cat("\n...returning uGARCHmultifit object instead...check and resubmit...")
			return( fitlist )
		}
		if(spec@model$modelinc[1]>0) model$mu = mu else model$mu = fitted(fitlist)
		model$residuals = res = residuals(fitlist)
		model$sigma = sig = sigma(fitlist)
	}
	
	stdresid = res/sig
	#-----------------------------------------------------------------------------------
	# Create the Full Model Pars and Indices
	modelinc =  model$modelinc
	# create full par matrix
	midx = .fullinc2(modelinc, umodel)
	mpars = midx*0
	# This is the estimated parameter index
	eidx = .estindfn(midx, mspec, model$pars)
	
	unipars = lapply(fitlist@fit, FUN = function(x) x@fit$ipars[x@fit$ipars[,3]==1,1])
	if(is.list(unipars)){
		for(i in 1:length(unipars)){
			uninames = names(unipars[[i]])
			mpars[uninames, i] = unipars[[i]]
		}
	} else{
		uninames = rownames(unipars)
		mpars[uninames, 1:NCOL(unipars)] = unipars
	}
	# add include pars from DCC spec (includes the fixed pars)
	mpars[which(midx[,m+1]==1), m+1] = as.numeric( model$pars[model$pars[,3]==1,1] )
	
	# DCC parameters
	ipars = model$pars
	LB 	= ipars[,5]
	UB 	= ipars[,6]
	estidx = as.logical( ipars[,4] )
	npars = sum(estidx)
		
	# Transformation of z to [0,1] domain
	ans = switch(transformation,
			parametric = .pparametric(fitlist, stdresid),
			empirical = .pempirical(stdresid),
			spd = .pspd(stdresid, spd.control))
	
	if(transformation == "spd"){
		ures = ans$ures
		sfit = ans$sfit
	} else{
		ures = ans
		sfit = NULL
	}
	# make small tail adjustments in order to avoid problem with the quantile 
	# functions in optimization
	if(any(ures > 0.99999)){
		xn = which(ures > 0.99999)
		ures[xn] = 0.99999
	}
	if(any(ures < .eps)){
		xn = which(ures < (1.5*.eps))
		ures[xn] = .eps
	}
	
	model$spd.control = spd.control
	mgarchenv = new.env(hash = TRUE)
	arglist = list()
	arglist$mgarchenv = mgarchenv
	arglist$verbose = FALSE
	arglist$spd.control = spd.control
	arglist$cluster = cluster
	arglist$eval.se = fit.control$eval.se
	arglist$solver = solver
	arglist$fit.control = fit.control
	arglist$cnames = cnames
	arglist$m = m
	arglist$T = T
	arglist$data = zdata
	arglist$index = index
	arglist$realizedVol = realizedVol
	arglist$model = model
	arglist$fitlist = fitlist
	arglist$umodel = umodel
	arglist$midx = midx
	arglist$eidx = eidx
	arglist$mpars = mpars
	arglist$ipars = ipars
	arglist$estidx = estidx
	arglist$dccN = npars
	arglist$stdresid = stdresid
	arglist$ures = ures
	
	if(any(ipars[,2]==1)){
		if(npars == 0){
			xspex = spec
			for(i in 1:m) xspex@umodel$fixed.pars[[i]] = as.list(fitlist@fit[[i]]@model$pars[fitlist@fit[[i]]@model$pars[,3]==1,1])
			return(cgarchfilter(spec = xspex, data = xts(data, index), out.sample = out.sample, 
							filter.control = list(), spd.control = spd.control, 
							cluster = cluster, VAR.fit = VAR.fit, realizedVol = realizedVol))
		} else{
			# with some parameters fixed we extract them (to be rejoined at end)
			# so that they do not enter the solver
			use.solver = 1
		}
	} else{
		use.solver = 1
	}
	# start counter
	assign("rmgarch_llh", 1, envir = mgarchenv)
	ILB = 0
	IUB = 1
	# Asymmetric Spec has different constraints
	if(model$modelinc[6]> 0) Ifn = copula.adcccon else Ifn = copula.dcccon
	if( solver == "solnp" | solver == "gosolnp" ) fit.control$stationarity = FALSE else fit.control$stationarity = TRUE
	arglist$fit.control = fit.control
	
	# get
	if( use.solver )
	{
		arglist$returnType = "llh"
		solution = switch(model$modeldesc$distribution,
				mvnorm = .copulasolver(solver, pars = ipars[estidx, 1], 
						fun = copula.tvnormalLLH1, Ifn, ILB, IUB, gr = NULL, 
						hessian = NULL, control = solver.control, 
						LB = ipars[estidx, 5], UB = ipars[estidx, 6], 
						arglist = arglist),
				mvt = .copulasolver(solver, pars = ipars[estidx, 1], 
						fun = copula.tvstudentLLH1, Ifn, ILB, IUB, gr = NULL, 
						hessian = NULL, control = solver.control, 
						LB = ipars[estidx, 5], UB = ipars[estidx, 6], 
						arglist = arglist))
		sol = solution$sol
		hess = solution$hess
		timer = Sys.time()-tic
		convergence = sol$convergence
		mpars[which(eidx[,(m+1)]==1, arr.ind = TRUE),m+1] = sol$pars
		ipars[estidx, 1] = sol$pars
		arglist$mpars = mpars
		arglist$ipars = ipars
	} else{
		hess = NULL
		timer = Sys.time()-tic
		convergence = 1
		sol = list()
		sol$message = "all parameters fixed"
	}
	# add some tail data of the univariate fits to be used with simulation 
	# (we assume that noone is going to add more than 50 lags in a univariate fit!
	mfit = list()
	if( convergence == 0 ){
		arglist$returnType = "ALL"
		mfit = switch(model$modeldesc$distribution,
				mvnorm = copula.tvnormalLLH2(mpars[which(eidx==1, arr.ind = TRUE)], 
						arglist = arglist),
				mvt = copula.tvstudentLLH2(mpars[which(eidx==1, arr.ind = TRUE)], 
						arglist = arglist))
		mfit$tailusigma = tail(sigma(fitlist), 50)
		mfit$tailuresids = tail(residuals(fitlist), 50)
		mfit$tailuret = tail(data, 50)
		
		nderiv = switch(model$modeldesc$distribution,
				mvnorm = .cgarchmakefitmodel1(f = copula.tvnormalLLH2, arglist = arglist, 
						timer = timer, message = sol$message, fname = "copula.tvnormalLLH2"),
				mvt = .cgarchmakefitmodel1(f = copula.tvstudentLLH2, arglist = arglist, 
						timer = timer, message = sol$message, fname = "copula.tvstudentLLH2"))
		
		mfit = c(mfit, nderiv)
		# make model list to return some usefule information which
		model$mpars = mpars
		model$ipars = ipars
		model$pars[,1] = ipars[,1]
		model$midx = midx
		model$eidx = eidx
		model$umodel = umodel
	} else{	
		mfit$convergence = 1
	}
	mfit$realizedVol = realizedVol
	mfit$timer = Sys.time() - tic
	# Add the spd fit as this is needed for simulation
	model$sfit = sfit
	rm(mgarchenv)
	ans = new("cGARCHfit",
			mfit = mfit,
			model = model)
	return(ans)
}

.cgarchfit.static = function(spec, data, spd.control = list(lower = 0.1, upper = 0.9, 
				type = "pwm", kernel = "epanech"),  
		fit.control = list(eval.se = TRUE, stationarity = TRUE, scale = FALSE), 
		solver = "solnp", solver.control = list(), out.sample = 0, cluster = NULL, 
		fit = NULL, VAR.fit = NULL, realizedVol = NULL, ...)
{
	tic = Sys.time()
	.eps = .Machine$double.eps
	model = spec@model
	umodel = spec@umodel
	ufit.control = list()
	if(is.null(fit.control$stationarity)){
		ufit.control$stationarity = TRUE 
	} else {
		ufit.control$stationarity = fit.control$stationarity
		fit.control$stationarity = NULL
	}
	if(is.null(fit.control$scale)){
		ufit.control$scale = TRUE 
	} else{
		ufit.control$scale = fit.control$scale
		fit.control$scale = NULL
	}
	if(is.null(fit.control$eval.se)) fit.control$eval.se = TRUE
	
	if(length(solver)==2){
		garch.solver = 	solver[1]
		solver = solver[2]
	} else{
		garch.solver = solver[1]
	}
	solver = match.arg(tolower(solver)[1], c("solnp", "nlminb", "lbfgs","gosolnp"))
	#-----------------------------------------------------------------------------------
	# Data Extraction
	m = dim(data)[2]
	if( is.null( colnames(data) ) ) cnames = paste("Asset_", 1:m, sep = "") else cnames = colnames(data)
	colnames(umodel$modelinc) = cnames
	
	xdata = .extractmdata(data)
	if(!is.numeric(out.sample)) 
		stop("\ndccfit-->error: out.sample must be numeric\n")
	if(as.numeric(out.sample) < 0) 
		stop("\ndccfit-->error: out.sample must be positive\n")
	n.start = round(out.sample, 0)
	n = dim(xdata$data)[1]
	if( (n-n.start) < 100) 
		stop("\ndccfit-->error: function requires at least 100 data\n points to run\n")
	data  = xdata$data
	index = xdata$index
	period = xdata$period
	
	# save the data to the model spec
	model$modeldata$data = data
	model$modeldata$index = index
	model$modeldata$period = period
	T = model$modeldata$T = n - n.start
	model$modeldata$n.start = n.start
	model$modeldata$asset.names = cnames
	
	#-----------------------------------------------------------------------------------
	transformation = model$modeldesc$transformation
	
	#-----------------------------------------------------------------------------------
	# VAR model
	if( model$modelinc[1]>0 ){
		tmp = mvmean.varfit(model = model, data = data, VAR.fit = VAR.fit, T = T, 
				out.sample = out.sample, cluster = cluster)
		model = tmp$model
		zdata = tmp$zdata
		mu = tmp$mu
		varcoef = tmp$varcoef
		p = tmp$p
		N = tmp$N
	} else{
		zdata = data
		ex = NULL
	}
	T = dim(zdata)[1] - out.sample
	#-----------------------------------------------------------------------------------
	# Univariate GARCH fit (only subtract the ARMA order)
	# create a multispec list
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, umodel$vt)
	
	if( !is.null(fit) && is(fit, "uGARCHmultifit") ){
		# check VAR and fit:
		if(model$modelinc[1]>0){
			for(i in 1:m){
				if(sum(fit@fit[[i]]@model$modelinc[1:6])>0)
					stop("\nThe user supplied fit object has a non-null mean specification but VAR already chosen for mean filtration!!!")
			}
		}
		fitlist = fit
		if(spec@model$modelinc[1]>0) model$mu = mu else model$mu = fitted(fitlist)
		model$residuals = res = residuals(fitlist)
		model$sigma = sig = sigma(fitlist)
	} else{
		fitlist = multifit(multispec = mspec, data = xts(zdata, index), out.sample = n.start, 
				solver = garch.solver, solver.control = solver.control, 
				fit.control = ufit.control, cluster = cluster, realizedVol = realizedVol, ...)
		converge = sapply(fitlist@fit, FUN = function(x) x@fit$convergence)
		if( any( converge == 1 ) ){
			pr = which(converge != 1)
			cat("\nNon-Converged:\n")
			print(pr)
			cat("\ndccfit-->error: convergence problem in univariate fit...")
			cat("\n...returning uGARCHmultifit object instead...check and resubmit...")
			return( fitlist )
		}
		if(spec@model$modelinc[1]>0) model$mu = mu else model$mu = fitted(fitlist)
		model$residuals = res = residuals(fitlist)
		model$sigma = sig = sigma(fitlist)
	}
	stdresid = res/sig
	#-----------------------------------------------------------------------------------
	# Create the Full Model Pars and Indices
	modelinc =  model$modelinc
	# create full par matrix
	midx = .fullinc2(modelinc, umodel)
	mpars = midx*0
	# This is the estimated parameter index
	eidx = .estindfn(midx, mspec, model$pars)
	unipars = lapply(fitlist@fit, FUN = function(x) x@fit$ipars[x@fit$ipars[,3]==1,1])
	if(is.list(unipars)){
		for(i in 1:length(unipars)){
			uninames = names(unipars[[i]])
			mpars[uninames, i] = unipars[[i]]
		}
	} else{
		uninames = rownames(unipars)
		mpars[uninames, 1:NCOL(unipars)] = unipars
	}
	# add include pars from DCC spec (includes the fixed pars)
	mpars[which(midx[,m+1]==1), m+1] = as.numeric( model$pars[model$pars[,3]==1,1] )
	
	ipars = model$pars
	LB 	= ipars[,5]
	UB 	= ipars[,6]
	estidx = as.logical( ipars[,4] )
	npars = sum(estidx)
	# Transformation of z to [0,1] domain
	ans = switch(transformation,
			parametric = .pparametric(fitlist, stdresid),
			empirical = .pempirical(stdresid),
			spd = .pspd(stdresid, spd.control))
	
	if(transformation == "spd"){
		ures = ans$ures
		sfit = ans$sfit
	} else{
		ures = ans
		sfit = NULL
	}
	
	# make small tail adjustments in order to avoid problem with the quantile 
	# functions in optimization
	if(any(ures > 0.99999)){
		xn = which(ures > 0.99999)
		ures[xn] = 0.99999
	}
	if(any(ures < .eps)){
		xn = which(ures < (1.5*.eps))
		ures[xn] = .eps
	}
	
	model$spd.control = spd.control
	mgarchenv = new.env(hash = TRUE)
	arglist = list()
	arglist$mgarchenv = mgarchenv
	arglist$verbose = FALSE
	arglist$spd.control = spd.control
	arglist$cluster = cluster
	arglist$eval.se = fit.control$eval.se
	arglist$solver = solver
	arglist$fit.control = fit.control
	arglist$cnames = cnames
	arglist$m = m
	arglist$T = T 	
	arglist$data = zdata
	arglist$index = index
	arglist$realizedVol = realizedVol
	arglist$model = model
	arglist$fitlist = fitlist
	arglist$umodel = umodel
	arglist$midx = midx
	arglist$eidx = eidx
	arglist$mpars = mpars
	arglist$ipars = ipars
	arglist$estidx = estidx
	arglist$stdresid = stdresid
	arglist$ures = ures
	arglist$npars = npars
	
	if(any(ipars[,2]==1)){
		if(npars == 0){
			use.solver = 0
		} else{
			# with some parameters fixed we extract them (to be rejoined at end)
			# so that they do not enter the solver
			use.solver = 1
		}
	} else{
		use.solver = 1
	}
	# start counter
	# get
	if( use.solver )
	{
		arglist$returnType = "llh"
		solution = try(
				switch(model$modeldesc$distribution,
				mvnorm = .copulasolver(solver, pars = ipars[estidx, 1], 
						fun = copula.normalLLH1, Ifn = NULL, ILB = NULL, IUB = NULL, 
						gr = NULL, hessian = NULL, control = solver.control, 
						LB = ipars[estidx, 5], UB = ipars[estidx, 6], 
						arglist = arglist),
				mvt = .copulasolver(solver, pars = ipars[estidx, 1], 
						fun = copula.studentLLH1, Ifn = NULL, ILB = NULL, IUB = NULL, 
						gr = NULL, hessian = NULL, control = solver.control, 
						LB = ipars[estidx, 5], UB = ipars[estidx, 6], 
						arglist = arglist) ), silent = TRUE)
		if(inherits(solution, "try-error") | solution$sol$convergence == 1){
			solver = "lbfgs"
			solver.control = list(trace = 1)
			solution = switch(model$modeldesc$distribution,
					mvnorm = .copulasolver(solver, pars = ipars[estidx, 1], 
							fun = copula.normalLLH1, Ifn = NULL, ILB = NULL, 
							IUB = NULL, gr = NULL, hessian = NULL, 
							control = solver.control, LB = ipars[estidx, 5], 
							UB = ipars[estidx, 6], arglist = arglist),
					mvt = .copulasolver(solver, pars = ipars[estidx, 1], 
							fun = copula.studentLLH1, Ifn = NULL, ILB = NULL, 
							IUB = NULL, gr = NULL, hessian = NULL, 
							control = solver.control, LB = ipars[estidx, 5], 
							UB = ipars[estidx, 6], arglist = arglist))
		}
		sol = solution$sol
		hess = solution$hess
		timer = Sys.time()-tic
		convergence = sol$convergence
		mpars[which(eidx[,(m+1)]==1, arr.ind = TRUE),m+1] = sol$pars
		ipars[estidx, 1] = sol$pars
		arglist$mpars = mpars
		arglist$ipars = ipars
	} else{
		hess = NULL
		timer = Sys.time()-tic
		convergence = 1
		sol = list()
		sol$message = "all parameters fixed"
		arglist$returnType = "ALL"
		mfit = switch(model$modeldesc$distribution,
				mvnorm = copula.normalLLH2(mpars[which(eidx==1, arr.ind = TRUE)], arglist = arglist),
				mvt = copula.studentLLH2(mpars[which(eidx==1, arr.ind = TRUE)], arglist = arglist))
		nderiv = switch(model$modeldesc$distribution,
				mvnorm = .cgarchmakefitmodel2(f = copula.normalLLH2, arglist = arglist, 
						timer = timer, message = sol$message, fname = "copula.normalLLH2"),
				mvt = .cgarchmakefitmodel2(f = copula.studentLLH2, arglist = arglist, 
						timer = timer, message = sol$message, fname = "copula.studentLLH2"))
		mfit = c(mfit, nderiv)
		
		mfit$tailusigma = tail(sigma(fitlist), 25)
		mfit$tailuresids = tail(residuals(fitlist), 25)
		mfit$tailuret = tail(data, 25)
		# make model list to return some usefule information which
		model$mpars = mpars
		model$ipars = ipars
		model$pars[,1] = ipars[,1]
		model$midx = midx
		model$eidx = eidx
		model$umodel = umodel
	}
	# add some tail data of the univariate fits to be used with simulation 
	# (we assume that noone is going to add more than 50 lags in a univariate fit!
	if( convergence == 0 ){
		arglist$returnType = "ALL"
		mfit = switch(model$modeldesc$distribution,
				mvnorm = copula.normalLLH2(mpars[which(eidx==1, arr.ind = TRUE)], arglist = arglist),
				mvt = copula.studentLLH2(mpars[which(eidx==1, arr.ind = TRUE)], arglist = arglist))
		mfit$tailusigma = tail(sigma(fitlist), 50)
		mfit$tailuresids = tail(residuals(fitlist), 50)
		mfit$tailuret = tail(data, 50)

		nderiv = switch(model$modeldesc$distribution,
				mvnorm = .cgarchmakefitmodel2(f = copula.normalLLH2, arglist = arglist, 
						timer = timer, message = sol$message, fname = "copula.normalLLH2"),
				mvt = .cgarchmakefitmodel2(f = copula.studentLLH2, arglist = arglist, 
						timer = timer, message = sol$message, fname = "copula.studentLLH2"))
		
		mfit = c(mfit, nderiv)
		# make model list to return some usefule information which
		model$mpars = mpars
		model$ipars = ipars
		model$pars[,1] = ipars[,1]
		model$midx = midx
		model$eidx = eidx
		model$umodel = umodel
		mfit$timer = Sys.time() - tic
	} else{
		mfit$convergence = 1
		mfit$timer = Sys.time() - tic
	}
	model$sfit = sfit
	mfit$realizedVol = realizedVol
	ans = new("cGARCHfit",
			mfit = mfit,
			model = model)
	return(ans)
}


.cgarchfilter = function(spec, data, out.sample = 0, filter.control = list(n.old = NULL), 
		spd.control = list(lower = 0.1, upper = 0.9, type = "pwm", kernel = "epanech"), 
		cluster = NULL, varcoef = NULL, realizedVol = NULL, ...)
{
	type = ifelse(spec@model$modeldesc$timecopula, "dynamic", "static")
	ans = switch(type, 
			dynamic = .cgarchfilter.dynamic(spec = spec, data = data, out.sample = out.sample,
					filter.control = filter.control, spd.control = spd.control, 
					cluster = cluster, varcoef = varcoef, realizedVol = realizedVol, ...),
			static = .cgarchfilter.static(spec = spec, data = data, out.sample = out.sample,
					filter.control = filter.control, spd.control = spd.control, 
					cluster = cluster, varcoef = varcoef, realizedVol = realizedVol, ...))
	return( ans )
}

.cgarchfilter.static = function(spec, data, out.sample = 0, filter.control = list(n.old = NULL), 
		spd.control = list(lower = 0.1, upper = 0.9, type = "pwm", kernel = "epanech"), 
		cluster = NULL, varcoef = NULL, realizedVol = NULL, ...)
{
	tic = Sys.time()
	model = spec@model
	umodel = spec@umodel
	n.old = filter.control$n.old
	#-----------------------------------------------------------------------------------
	# Data Extraction
	m = dim(data)[2]
	if( is.null( colnames(data) ) ) cnames = paste("Asset_", 1:m, sep = "") else cnames = colnames(data)
	colnames(umodel$modelinc) = cnames
	
	xdata = .extractmdata(data)
	if(!is.numeric(out.sample)) 
		stop("\ndccfilter-->error: out.sample must be numeric\n")
	if(as.numeric(out.sample) < 0) 
		stop("\ndccfilter-->error: out.sample must be positive\n")
	n.start = round(out.sample, 0)
	n = dim(xdata$data)[1]
	if( (n-n.start) < 100) 
		stop("\ndccfilter-->error: function requires at least 100 data\n points to run\n")
	data  = xdata$data
	index = xdata$index
	period = xdata$period
	# save the data to the model spec
	model$modeldata$data = data
	model$modeldata$index = index
	model$modeldata$period = period
	T = model$modeldata$T = n - n.start
	model$modeldata$n.start = n.start
	model$modeldata$asset.names = cnames
	#-----------------------------------------------------------------------------------
	transformation = model$modeldesc$transformation
	if(is.null(filter.control$n.old)) n.old = T
	#-----------------------------------------------------------------------------------
	# VAR model
	if( spec@model$modelinc[1]>0 ){
		tmp = mvmean.varfilter(model = model, data = data, varcoef = varcoef, 
				T = T, out.sample = out.sample)
		model = tmp$model
		zdata = tmp$zdata
		mu = tmp$mu
		p = tmp$p
		N = tmp$N
	} else{
		zdata = data
		ex = NULL
	}
	T = dim(zdata)[1] - out.sample
	#-----------------------------------------------------------------------------------
	if(is.null(filter.control$n.old)) n.old = T
	#-----------------------------------------------------------------------------------
	# Univariate GARCH filter
	if(model$modelinc[1]>0){
		for(i in 1:m){
			if(sum(umodel$modelinc[1:6,i])>0)
				stop("\nThe user supplied univariate spec object has a non-null mean specification but VAR already chosen for mean filtration!!!")
		}
	}
	
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, umodel$vt)
	
	filterlist = multifilter(multifitORspec = mspec, data = xts(zdata, index), out.sample = out.sample, 
			cluster = cluster, n.old = n.old, realizedVol = realizedVol, ...)
	
	if(spec@model$modelinc[1]>0) model$mu = mu else model$mu = fitted(filterlist)
	model$residuals = res = residuals(filterlist)
	model$sigma = sig = sigma(filterlist)
	stdresid = res/sig
	N = dim(stdresid)[1]
	
	modelinc =  model$modelinc
	# create full par matrix
	midx = .fullinc2(modelinc, umodel)
	mpars = midx*0
	eidx = midx
	unipars = sapply(filterlist@filter, FUN = function(x) x@filter$ipars[x@filter$ipars[,3]==1,1])
	if(is.list(unipars)){
		for(i in 1:length(unipars)){
			uninames = names(unipars[[i]])
			mpars[uninames, i] = unipars[[i]]
		}
	} else{
		uninames = rownames(unipars)
		mpars[uninames, 1:NCOL(unipars)] = unipars
	}
	# add include pars from DCC spec (includes the fixed pars)
	mpars[which(midx[,m+1]==1, arr.ind = TRUE), m+1] = as.numeric( model$pars[model$pars[,3]==1,1] )
	
	# DCC parameters
	ipars = model$pars
	estidx = as.logical( ipars[,3] )
	npars = sum(estidx)
		
	model$spd.control = spd.control
	
	arglist = list()
	arglist$verbose = FALSE
	arglist$spd.control = spd.control
	arglist$cluster = cluster
	arglist$cnames = cnames
	arglist$m = m
	arglist$T = T
	arglist$data = zdata
	arglist$index = index
	arglist$realizedVol = realizedVol
	arglist$model = model
	arglist$filterlist = filterlist
	arglist$umodel = umodel
	arglist$midx = midx
	arglist$eidx = eidx
	arglist$mpars = mpars
	arglist$ipars = ipars
	arglist$estidx = estidx
	arglist$stdresid = stdresid
	arglist$npars = npars
	arglist$n.old = n.old
	arglist$filter.control = filter.control
	
	mfilter = list()
	arglist$returnType = "ALL"
	mfilter = switch(model$modeldesc$distribution,
				mvnorm = copula.normalLLH3(arglist = arglist),
				mvt = copula.studentLLH3(arglist = arglist))
	mfilter$tailusigma = tail(sigma(filterlist), 50)
	mfilter$tailuresids = tail(residuals(filterlist), 50)
	mfilter$tailuret = tail(data, 50)
	
	Rt = mfilter$R
	Ht = array( 0, dim = c(m, m, N) )
	stdresid = matrix(0, nrow = N, ncol = m)
	
	if( !is.null(cluster) ){
			clusterExport(cluster, c("sig", "Rt", "res"), envir = environment())
			tmp = parLapply(cluster, as.list(1:N), fun = function(i){
						tmph = diag( sig[i, ] ) %*% Rt %*% diag( sig[i, ] )
						zz = eigen( tmph )
						sqrtzz = ( zz$vectors %*% diag( sqrt( zz$values ) ) %*% solve( zz$vectors ) )
						tmpz = as.numeric( res[i, ] %*% solve( sqrtzz ) )
						return( list( H = tmph, Z = tmpz ) )
					})
		for(i in 1:N){
			Ht[,,i] = tmp[[i]]$H
			stdresid[i,] = tmp[[i]]$Z
		}
	} else{
		tmp = lapply(as.list(1:N), FUN = function(i){
					tmph = diag( sig[i, ] ) %*% Rt %*% diag( sig[i, ] )
					zz = eigen( tmph )
					sqrtzz = ( zz$vectors %*% diag( sqrt( zz$values ) ) %*% solve( zz$vectors ) )
					tmpz = as.numeric( res[i, ] %*% solve( sqrtzz ) )
					return( list( H = tmph, Z = tmpz ) )
				})
		for(i in 1:N){
			Ht[,,i] = tmp[[i]]$H
			stdresid[i,] = tmp[[i]]$Z
		}
	}
	mfilter$H = Ht
	allnames = NULL
	for(i in 1:m){
		allnames = c(allnames, paste("[",cnames[i],"].", rownames(midx[midx[,i]==1,i, drop = FALSE]), sep = ""))
	}
	garchnames = allnames
	dccnames = rownames(midx[midx[,m+1]==1,m+1, drop = FALSE])
	if(!is.null(dccnames)){
		allnames = c(garchnames, paste("[Joint]", rownames(midx[midx[,m+1]==1,m+1, drop = FALSE]), sep = ""))
	} else{
		allnames = garchnames
	}
	mfilter$coef = mpars[which(midx==1, arr.ind = TRUE)]
	names(mfilter$coef) = allnames
	mfilter$garchnames = garchnames
	mfilter$dccnames = dccnames
	mfilter$stdresid = stdresid
	# make model list to return some usefule information which
	model$mpars = mpars
	model$ipars = ipars
	model$pars[,1] = ipars[,1]
	model$midx = midx
	model$eidx = eidx
	model$umodel = umodel
	mfilter$realizedVol = realizedVol
	mfilter$timer = Sys.time() - tic
	#model$sfit = sfit
	
	ans = new("cGARCHfilter",
			mfilter = mfilter,
			model = model)
	return(ans)
}

# ToDo : change the likelihood to take into account the n.old argument for the cov calculation.
.cgarchfilter.dynamic = function(spec, data, out.sample = 0, filter.control = list(n.old = NULL), 
		spd.control = list(lower = 0.1, upper = 0.9, type = "pwm", kernel = "epanech"), 
		cluster = NULL, varcoef = NULL, realizedVol = NULL, ...)
{
	tic = Sys.time()
	model = spec@model
	umodel = spec@umodel
	n.old = filter.control$n.old
	#-----------------------------------------------------------------------------------
	# Data Extraction
	m = dim(data)[2]
	if( is.null( colnames(data) ) ) cnames = paste("Asset_", 1:m, sep = "") else cnames = colnames(data)
	colnames(umodel$modelinc) = cnames
	
	xdata = .extractmdata(data)
	if(!is.numeric(out.sample)) 
		stop("\ndccfilter-->error: out.sample must be numeric\n")
	if(as.numeric(out.sample) < 0) 
		stop("\ndccfilter-->error: out.sample must be positive\n")
	n.start = round(out.sample, 0)
	n = dim(xdata$data)[1]
	if( (n-n.start) < 100) 
		stop("\ndccfilter-->error: function requires at least 100 data\n points to run\n")
	data  = xdata$data
	index = xdata$index
	period = xdata$period
		
	# save the data to the model spec
	model$modeldata$data = data
	model$modeldata$index = index
	model$modeldata$period = period
	T = model$modeldata$T = n - n.start
	model$modeldata$n.start = n.start
	model$modeldata$asset.names = cnames
	#-----------------------------------------------------------------------------------
	transformation = model$modeldesc$transformation
	#-----------------------------------------------------------------------------------
	# VAR model
	if( spec@model$modelinc[1]>0 ){
		tmp = mvmean.varfilter(model = model, data = data, varcoef = varcoef, 
				T = T, out.sample = out.sample)
		model = tmp$model
		zdata = tmp$zdata
		mu = tmp$mu
		p = tmp$p
		N = tmp$N
	} else{
		zdata = data
		ex = NULL
	}
	T = dim(zdata)[1] - out.sample
	#-----------------------------------------------------------------------------------
	if(is.null(filter.control$n.old)) n.old = T
	#-----------------------------------------------------------------------------------
	# Univariate GARCH filter
	if(model$modelinc[1]>0){
		for(i in 1:m){
			if(sum(umodel$modelinc[1:6,i])>0)
				stop("\nThe user supplied univariate spec object has a non-null mean specification but VAR already chosen for mean filtration!!!")
		}
	}	
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, umodel$vt)
	
	filterlist = multifilter(multifitORspec = mspec, data = xts(zdata, index), out.sample = out.sample, 
			cluster = cluster, n.old = n.old, realizedVol = realizedVol, ...)
	
	if(spec@model$modelinc[1]>0) model$mu = mu else model$mu = fitted(filterlist)
	model$residuals = res = residuals(filterlist)
	model$sigma = sig = sigma(filterlist)
	stdresid = res/sig
	modelinc =  model$modelinc
	# create full par matrix
	midx = .fullinc2(modelinc, umodel)
	mpars = midx*0
	eidx = midx
	unipars = sapply(filterlist@filter, FUN = function(x) x@filter$ipars[x@filter$ipars[,3]==1,1])
	if(is.list(unipars)){
		for(i in 1:length(unipars)){
			uninames = names(unipars[[i]])
			mpars[uninames, i] = unipars[[i]]
		}
	} else{
		uninames = rownames(unipars)
		mpars[uninames, 1:NCOL(unipars)] = unipars
	}
	# add include pars from DCC spec (includes the fixed pars)
	mpars[which(midx[,m+1]==1, arr.ind = TRUE), m+1] = as.numeric( model$pars[model$pars[,3]==1,1] )
	
	# DCC parameters
	ipars = model$pars
	estidx = as.logical( ipars[,3] )
	npars = sum(estidx)
	
	model$spd.control = spd.control
	
	arglist = list()
	arglist$verbose = FALSE
	arglist$spd.control = spd.control
	arglist$cluster = cluster
	arglist$cnames = cnames
	arglist$m = m
	arglist$T = T
	arglist$zdata = zdata
	arglist$index = index
	arglist$realizedVol = realizedVol
	arglist$model = model
	arglist$filterlist = filterlist
	arglist$umodel = umodel
	arglist$midx = midx
	arglist$eidx = eidx
	arglist$mpars = mpars
	arglist$ipars = ipars
	arglist$estidx = estidx
	arglist$stdresid = stdresid
	arglist$npars = npars
	arglist$filter.control = filter.control
	arglist$n.old = n.old
	
	mfilter = list()
	arglist$returnType = "ALL"
	mfilter = switch(model$modeldesc$distribution,
			mvnorm = copula.tvnormalLLH3(arglist = arglist),
			mvt = copula.tvstudentLLH3(arglist = arglist))
	mfilter$tailusigma = tail(sigma(filterlist), 50)
	mfilter$tailuresids = tail(residuals(filterlist), 50)
	mfilter$tailuret = tail(data, 50)
	
	N = dim(stdresid)[1]
	Rt = mfilter$R
	Ht = array( 0, dim = c(m, m, N) )
	stdresid = matrix(0, nrow = N, ncol = m)
	
	if( !is.null(cluster) ){
			clusterExport(cluster, c("sig", "Rt", "res"), envir = environment())
			tmp = parLapply(cluster, as.list(1:N), fun = function(i){
						tmph = diag( sig[i, ] ) %*% Rt[[i]] %*% diag( sig[i, ] )
						zz = eigen( tmph )
						sqrtzz = ( zz$vectors %*% diag( sqrt( zz$values ) ) %*% solve( zz$vectors ) )
						tmpz = as.numeric( res[i, ] %*% solve( sqrtzz ) )
						return( list( H = tmph, Z = tmpz ) )
					})
		for(i in 1:N){
			Ht[,,i] = tmp[[i]]$H
			stdresid[i,] = tmp[[i]]$Z
		}
	} else{
		tmp = lapply(as.list(1:N), FUN = function(i){
					tmph = diag( sig[i, ] ) %*% Rt[[i]] %*% diag( sig[i, ] )
					zz = eigen( tmph )
					sqrtzz = ( zz$vectors %*% diag( sqrt( zz$values ) ) %*% solve( zz$vectors ) )
					tmpz = as.numeric( res[i, ] %*% solve( sqrtzz ) )
					return( list( H = tmph, Z = tmpz ) )
				})
		for(i in 1:N){
			Ht[,,i] = tmp[[i]]$H
			stdresid[i,] = tmp[[i]]$Z
		}
	}
	mfilter$H = Ht
	allnames = NULL
	for(i in 1:m){
		allnames = c(allnames, paste("[",cnames[i],"].", rownames(midx[midx[,i]==1,i, drop = FALSE]), sep = ""))
	}
	garchnames = allnames
	dccnames = rownames(midx[midx[,m+1]==1,m+1, drop = FALSE])
	allnames = c(garchnames, paste("[Joint]",dccnames, sep = ""))
	
	mfilter$coef = mpars[which(midx==1, arr.ind = TRUE)]
	names(mfilter$coef) = allnames
	mfilter$garchnames = garchnames
	mfilter$dccnames = dccnames
	mfilter$stdresid = stdresid
	# make model list to return some useful information which
	model$mpars = mpars
	model$ipars = ipars
	model$pars[,1] = ipars[,1]
	model$midx = midx
	model$eidx = eidx
	model$umodel = umodel
	mfilter$realizedVol = realizedVol
	mfilter$timer = Sys.time() - tic
	#model$sfit = sfit
	
	ans = new("cGARCHfilter",
			mfilter = mfilter,
			model = model)
	return(ans)
}

# allow returning only the scenario (return) matrix
.cgarchsim = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, 
		startMethod = c("unconditional", "sample"), presigma = NULL, 
		preresiduals = NULL, prereturns = NULL, preR = NULL, preQ = NULL, 
		preZ = NULL, rseed = NULL, mexsimdata = NULL, vexsimdata = NULL, 
		cluster = NULL, only.density = FALSE, prerealized = NULL, ...)
{
	timecopula = fit@model$modeldesc$timecopula
	if(timecopula){
		ans = .cgarchsim2(fit = fit, n.sim = n.sim, n.start = n.start, m.sim = m.sim, 
				startMethod = startMethod[1], presigma = presigma, preresiduals = preresiduals, 
				prereturns = prereturns, preR = preR, preQ = preQ, preZ = preZ, 
				rseed = rseed, mexsimdata = mexsimdata, vexsimdata = vexsimdata, 
				cluster = cluster, only.density = only.density, prerealized = prerealized)
	} else{
		ans = .cgarchsim1(fit = fit, n.sim = n.sim, n.start = n.start, m.sim = m.sim, 
				startMethod = startMethod[1], presigma = presigma, preresiduals = preresiduals, 
				prereturns = prereturns, preR = preR, rseed = rseed, 
				mexsimdata = mexsimdata, vexsimdata = vexsimdata, cluster = cluster, 
				only.density = only.density, prerealized = prerealized)
	}
	return(ans)
}

.cgarchsim1 = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, 
		startMethod = c("unconditional", "sample"), preresiduals = NULL, 
		presigma = NULL, prereturns = NULL, preR = NULL, rseed = NULL, 
		mexsimdata = NULL, vexsimdata = NULL, cluster = NULL, only.density = FALSE, 
		prerealized = NULL, ...)
{
	# first generate the copula random uniform numbers (static copula)
	# ures --> transform --> zres
	# pass zres to GARCH simulation
	# pass if needed to VAR model
	if( is.null(rseed) ){
		rseed = as.integer(runif(m.sim, 1, Sys.time()))
	} else {
		if(length(rseed) == 1) rseed = c(rseed[1], rseed[1]+seq_len(m.sim))
		rseed = as.integer( rseed[1:m.sim] )
	}
	
	model = fit@model
	umodel = model$umodel
	T = model$modeldata$T
	Data = model$modeldata$data[1:T, ]
	m = dim(Data)[2]
	
	nsim = n.sim + n.start
	mg = model$maxgarchOrder
	
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, 
			umodel$modeldesc$vsubmodel, umodel$modeldata$mexdata, 
			umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, NULL)
	for(i in 1:m){
		# variance targeting case
		if(model$midx["omega",i]==0){
			setfixed(mspec@spec[[i]]) = as.list(c("omega"=model$mpars["omega",i], model$mpars[which(model$midx[,i]==1), i]))
		} else{
			setfixed(mspec@spec[[i]]) = as.list(model$mpars[which(model$midx[,i]==1), i])    
		}
	}
	#setfixed(specx)<-as.list(model$mpars[model$midx[,m+1]==1,m+1])
	
	uncv = sapply(mspec@spec, FUN = function(x) uncvariance(x))
	
	if( !is.null(presigma) ){
		if( !is.matrix(presigma) ) 
			stop("\ncgarchsim-->error: presigma must be a matrix.")
		if( dim(presigma)[2] != m ) 
			stop("\ncgarchsim-->error: wrong column dimension for presigma.")
		if( dim(presigma)[1] != mg ) 
			stop(paste("\ncgarchsim-->error: wrong row dimension for presigma (need ", mg, " rows.", sep = ""))		
	} else{
		if(startMethod == "sample"){
			presigma = fit@mfit$tailusigma
		} else{
			presigma = matrix(uncv, ncol = m, nrow = 25)
		}
	}
	
	if(is.null(preR)){
		Rbar = rcor(fit)
	} else{
		if(dim(preR)[1] != m) 
			stop("\ncgarchsim-->error: wrong dimension for preR\n")
		if(dim(preR)[2] != m) 
			stop("\ncgarchsim-->error: wrong dimension for preR\n")
		if(any(diag(preR) != 1)) 
			stop("\ncgarchsim-->error: preR diagonals must be 1.\n")
		Rbar = preR
	}
	if( !is.null(prereturns) ){
		if( !is.matrix(prereturns) ) 
			stop("\ncgarchsim-->error: prereturns must be a matrix.")
		if( dim(prereturns)[2] != m ) 
			stop("\ncgarchsim-->error: wrong column dimension for prereturns.")
		if( dim(prereturns)[1] != mg ) 
			stop(paste("\ncgarchsim-->error: wrong row dimension for prereturns (need ", mg, " rows.", sep = ""))
	} else{
		prereturns = tail(model$modeldata$data[1:T, ], 25)
	}
	
	if(fit@model$umodel$modeldesc$vmodel[1]=="realGARCH"){
		if( !is.null(prerealized) ){
			if( !is.matrix(prerealized) ) 
				stop("\ncgarchsim-->error: prerealized must be a matrix.")
			if( dim(prerealized)[2] != m ) 
				stop("\ncgarchsim-->error: wrong column dimension for prerealized.")
			if( dim(prerealized)[1] != mg ) 
				stop(paste("\ncgarchsim-->error: wrong row dimension for prerealized (need ", mg, " rows.", sep = ""))
		} else{
			# might want to include the option for unconditional
		    # value of realized (unexported in rugarch)
			prerealized = tail(fit@mfit$realizedVol[1:T,], 25)
		}
	} else{
		prerealized = matrix(NA, ncol = m, nrow = 25)
	}
	
	ures = .sample.copula(model, Qbar = NULL, preQ = NULL, Rbar = Rbar, 
			Nbar = NULL, preZ = NULL, n.sim = n.sim, n.start = n.start, 
			m.sim = m.sim, rseed = rseed, cluster = cluster)
	
	# ures = [0,1] copula random numbers
	# now transform back into margins for use in garch
	transformation = model$modeldesc$transformation
	distu = umodel$modeldesc$distribution
	
	zres = array(NA, dim = c(nsim, m, m.sim))
	# parallel:
	if( !is.null(cluster) ){
		ssfit = model$sfit
		minc = umodel$modelinc
		mpars = model$mpars
		sxres  = fit@mfit$stdresid
		clusterEvalQ(cluster, loadNamespace('rmgarch'))
		clusterEvalQ(cluster, loadNamespace('rugarch'))
		if(transformation == "spd"){
			clusterExport(cluster, c("ures", "mpars", "minc", "m", 
							"sxres", "ssfit", "transformation"), envir = environment())
		} else{
			clusterExport(cluster, c("ures", "mpars", "minc", "m", "sxres", 
							"transformation"), envir = environment())
		}
		clusterExport(cluster, c(".qparametric", ".qempirical",".qspd"), envir = environment())
		mtmp = parLapply(cluster, as.list(1:m.sim), fun = function(i){
					switch(transformation,
							parametric = .qparametric(matrix(ures[,,i], ncol = m), pars = mpars, modelinc = minc),
							empirical = .qempirical(matrix(ures[,,i], ncol = m), sxres),
							spd = .qspd(matrix(ures[,,i], ncol = m), sfit = ssfit))
				})
		for(i in 1:m.sim) zres[,,i] = mtmp[[i]]
	} else{
		mtmp = lapply(as.list(1:m.sim), FUN = function(i){
					switch(transformation,
							parametric = .qparametric(matrix(ures[,,i], ncol = m), 
									pars = model$mpars, modelinc = umodel$modelinc),
							empirical = .qempirical(matrix(ures[,,i], ncol = m), fit@mfit$stdresid),
							spd = .qspd(matrix(ures[,,i], ncol = m), sfit = model$sfit))
				})
		for(i in 1:m.sim) zres[,,i] = mtmp[[i]]
	}
	
	# Now we have the innovations which we feed back into the garch routine to simulate	
	simRes = simX = simR = simQ = simH = simSeries = vector(mode = "list", length = m.sim)
	
	if( !is.null(cluster) ){
		tailres = fit@mfit$tailuresids
		clusterEvalQ(cluster, loadNamespace('rugarch'))
		clusterExport(cluster, c("mspec", "n.sim", "n.start", "m.sim", 
						"startMethod", "zres", "presigma", "tailres", 
						"preresiduals", "prereturns", "model", 
						"mexsimdata", "vexsimdata","prerealized"), envir = environment())	
		simlist = parLapply(cluster, as.list(1:m), fun = function(i){
						maxx = mspec@spec[[i]]@model$maxOrder;
						htmp = rugarch::ugarchpath(mspec@spec[[i]], n.sim = n.sim + n.start, n.start = 0, m.sim = m.sim,
								custom.dist = list(name = "sample", distfit = matrix(zres[,i,1:m.sim], ncol = m.sim)),
								presigma = tail(presigma[,i], maxx), 
								preresiduals = if( is.null(preresiduals) ) tail(tailres[,i], maxx) else tail(preresiduals[,i], maxx), 
								prereturns = if(model$modelinc[1]==0) tail(prereturns[,i], maxx) else NA,
								mexsimdata = if(model$modelinc[1]==0) mexsimdata[[i]] else NULL, vexsimdata = vexsimdata[[i]],
								prerealized = tail(prerealized[,i], maxx))
						h = matrix(tail(htmp@path$sigmaSim^2, n.sim), nrow = n.sim)
						x = matrix(htmp@path$seriesSim,  nrow = n.sim + n.start)
						return(list(h = h, x = x))
					})
		H = array(NA, dim = c(n.sim, m, m.sim))
		tmpH = array(NA, dim = c(m, m, n.sim))
		for(i in 1:n.sim) H[i,,] = t(sapply(simlist, FUN = function(x) as.numeric(x$h[i,])))
		for(i in 1:m.sim){
			for(j in 1:n.sim){
				tmpH[ , , j] = sqrt(diag( H[j, , i]) ) %*% Rbar %*% sqrt(diag( H[j, , i] ) )
			}
			simH[[i]] = tmpH
		}
		if(model$modelinc[1]>0){
			simxX = array(NA, dim = c(n.sim+n.start, m, m.sim))
			for(i in 1:m.sim) simxX[,,i] = sapply(simlist, FUN = function(x) as.numeric(x$x[,i]))
			simX = vector(mode = "list", length = m.sim)
			for(i in 1:m.sim) simX[[i]] = matrix(simxX[,,i], nrow = n.sim+n.start)
		} else{
			simxX = array(NA, dim = c(n.sim+n.start, m, m.sim))
			for(i in 1:m.sim) simxX[,,i] = sapply(simlist, FUN = function(x) as.numeric(x$x[,i]))
			simX = vector(mode = "list", length = m.sim)
			for(i in 1:m.sim) simX[[i]] = matrix(tail(matrix(simxX[,,i], ncol = m), n.sim), nrow = n.sim)
		}
	} else{
		simlist = vector(mode="list", length=m)
		for(i in 1:m){
			maxx = mspec@spec[[i]]@model$maxOrder
			htmp = ugarchpath(mspec@spec[[i]], n.sim = n.sim + n.start, n.start = 0, m.sim = m.sim,
					custom.dist = list(name = "sample", distfit = matrix(zres[,i,1:m.sim], ncol = m.sim)),
					presigma = tail(presigma[,i], maxx), 
					preresiduals = if( is.null(preresiduals) ) tail(fit@mfit$tailuresids[,i], maxx) else tail(preresiduals[,i], maxx), 
					prereturns = if(model$modelinc[1]==0) tail(prereturns[,i], maxx) else NA,
					mexsimdata = if(model$modelinc[1]==0) mexsimdata[[i]] else NULL, vexsimdata = vexsimdata[[i]],
					prerealized = tail(prerealized[,i], maxx))
			h = matrix(tail(htmp@path$sigmaSim^2, n.sim), nrow = n.sim)
			x = matrix(htmp@path$seriesSim,  nrow = n.sim + n.start)
			simlist[[i]] = list(h = h, x = x)
		}
		H = array(NA, dim = c(n.sim, m, m.sim))
		tmpH = array(NA, dim = c(m, m, n.sim))
		for(i in 1:n.sim) H[i,,] = t(sapply(simlist, FUN = function(x) as.numeric(x$h[i,])))
		for(i in 1:m.sim){
			for(j in 1:n.sim){
				tmpH[ , , j] = sqrt(diag( H[j, , i]) ) %*% Rbar %*% sqrt(diag( H[j, , i] ) )
			}
			simH[[i]] = tmpH
		}
		if(model$modelinc[1]>0){
			simxX = array(NA, dim = c(n.sim+n.start, m, m.sim))
			for(i in 1:m.sim) simxX[,,i] = sapply(simlist, FUN = function(x) as.numeric(x$x[,i]))
			simX = vector(mode = "list", length = m.sim)
			for(i in 1:m.sim) simX[[i]] = matrix(simxX[,,i], nrow = n.sim+n.start)
		} else{
			simxX = array(NA, dim = c(n.sim+n.start, m, m.sim))
			for(i in 1:m.sim) simxX[,,i] = sapply(simlist, FUN = function(x) as.numeric(x$x[,i]))
			simX = vector(mode = "list", length = m.sim)
			for(i in 1:m.sim) simX[[i]] = matrix(tail(matrix(simxX[,,i], ncol = m), n.sim), nrow = n.sim)
		}
	}
	
	if( model$modelinc[1]>0 ){
		simRes = simX
		simX = mvmean.varsim(model = model, Data = Data, res = simX, 
				mexsimdata = mexsimdata, prereturns = prereturns, m.sim = m.sim, 
				n.sim = n.sim, n.start = n.start, startMethod = startMethod, 
				cluster = cluster)
	} else{
		# reshape
		for(j in 1:m.sim) simX[[j]] = tail(simX[[j]], n.sim)
	}
	
	msim = list()
	if(only.density){
		msim$simX = simX
		msim$rseed = rseed
		model$n.sim = n.sim
		model$m.sim = m.sim
		model$n.start = n.start
		model$startMethod = startMethod[1]
	} else{
		msim$simH = simH
		msim$simR = simR
		msim$simQ = simQ
		msim$simX = simX
		msim$simRes = simRes
		msim$simZ = zres
		msim$rseed = rseed
		model$n.sim = n.sim
		model$m.sim = m.sim
		model$n.start = n.start
		model$startMethod = startMethod[1]
	}
	model$only.density = only.density
	
	# need to run gc and rm since 'memory accumulates here'
	if(exists("simX")) rm(simX)
	if(exists("simRes")) rm(simRes)
	if(exists("simH")) rm(simH)
	if(exists("H")) rm(H)
	if(exists("simlist")) rm(simlist)
	if(exists("tmpH")) rm(tmpH)
	if(exists("mtmp")) rm(mtmp)
	if(exists("simxX")) rm(simxX)
	if(exists("zres")) rm(zres)
	if(exists("ures")) rm(ures)
	
	gc(verbose = FALSE)
	
	ans = new("cGARCHsim",
			msim = msim,
			model = model)
	return( ans )

}

# time varying copula simulation
.cgarchsim2 = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, 
		startMethod = c("unconditional", "sample"), presigma = NULL, preresiduals = NULL, 
		prereturns = NULL, preR = NULL, preQ = NULL, preZ = NULL, rseed = NULL, 
		mexsimdata = NULL, vexsimdata = NULL, cluster = NULL, only.density = FALSE, 
		prerealized = NULL, ...)
{
	# first generate the copula random uniform numbers (static copula)
	# ures --> transform --> zres
	# pass zres to GARCH simulation
	# pass if needed to VAR model
	if( is.null(rseed) ){
		rseed = as.integer(runif(m.sim, 1, Sys.time()))
	} else {
		if(length(rseed) == 1) rseed = c(rseed[1], rseed[1]+seq_len(m.sim))
		rseed = as.integer( rseed[1:m.sim] )
	}
	model = fit@model
	umodel = model$umodel
	T = model$modeldata$T
	Data = model$modeldata$data[1:T, ]
	m = dim(Data)[2]
	
	nsim = n.sim + n.start
	mg = fit@model$maxgarchOrder
	
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, NULL)
	
	for(i in 1:m){
		# variance targeting case
		if(model$midx["omega",i]==0){
			setfixed(mspec@spec[[i]]) = as.list(c("omega"=model$mpars["omega",i], model$mpars[which(model$midx[,i]==1), i]))
		} else{
			setfixed(mspec@spec[[i]]) = as.list(model$mpars[which(model$midx[,i]==1), i])    
		}
	}		
	uncv = sapply(mspec@spec, FUN = function(x) uncvariance(x))
	
	if( !is.null(presigma) ){
		if( !is.matrix(presigma) ) 
			stop("\ncgarchsim-->error: presigma must be a matrix.")
		if( dim(presigma)[2] != m ) 
			stop("\ncgarchsim-->error: wrong column dimension for presigma.")
		if( dim(presigma)[1] != mg ) 
			stop(paste("\ncgarchsim-->error: wrong row dimension for presigma (need ", mg, " rows.", sep = ""))		
	} else{
		if(startMethod == "sample"){
			presigma = fit@mfit$tailusigma
		} else{
			presigma = matrix(uncv, ncol = m, nrow = 25)
		}
	}
	
	if(is.null(preR)){
		Rbar = last(rcor(fit), 1)[,,1]
	} else{
		if(dim(preR)[1] != m) 
			stop("\ncgarchsim-->error: wrong dimension for preR\n")
		if(dim(preR)[2] != m) 
			stop("\ncgarchsim-->error: wrong dimension for preR\n")
		if(any(diag(preR) != 1)) 
			stop("\ncgarchsim-->error: preR diagonals must be 1.\n")
		Rbar = preR
	}
	
	Qbar = fit@mfit$Qbar
	
	if(model$modelinc[6]>0){
		Nbar = fit@mfit$Nbar
	} else{
		Nbar = matrix(0, m, m)
	}
	
	
	if( !is.null(prereturns) ){
		if( !is.matrix(prereturns) ) 
			stop("\ncgarchsim-->error: prereturns must be a matrix.")
		if( dim(prereturns)[2] != m ) 
			stop("\ncgarchsim-->error: wrong column dimension for prereturns.")
		if( dim(prereturns)[1] != mg ) 
			stop(paste("\ncgarchsim-->error: wrong row dimension for prereturns (need ", mg, " rows.", sep = ""))
	} else{
		prereturns = tail(model$modeldata$data[1:T, ], 25)
	}
	
	if(fit@model$umodel$modeldesc$vmodel[1]=="realGARCH"){
		if( !is.null(prerealized) ){
			if( !is.matrix(prerealized) ) 
				stop("\ncgarchsim-->error: prerealized must be a matrix.")
			if( dim(prerealized)[2] != m ) 
				stop("\ncgarchsim-->error: wrong column dimension for prerealized.")
			if( dim(prerealized)[1] != mg ) 
				stop(paste("\ncgarchsim-->error: wrong row dimension for prerealized (need ", mg, " rows.", sep = ""))
		} else{
			# might want to include the option for unconditional
			# value of realized (unexported in rugarch)
			prerealized = tail(fit@mfit$realizedVol[1:T,], 25)
		}
	} else{
		prerealized = matrix(NA, ncol = m, nrow = 25)
	}
	
	# This step is KEY if we are to replicate the 1-ahead filter

	if(startMethod == "sample"){
		if(is.null(preZ)){
			preZ = matrix(tail(fit@mfit$Z, fit@model$maxdccOrder), ncol = m)
		} else{
			preZ = matrix(preZ[1,], ncol = m, nrow = fit@model$maxdccOrder)
		}
		if(is.null(preQ)){
			preQ = fit@mfit$Qt[[length(fit@mfit$Qt)]]
		} else{
			preQ = preQ
		}
	} else{
		if(is.null(preZ)){
			preZ = matrix(0, ncol = m, nrow = fit@model$maxdccOrder)
		} else{
			preZ = matrix(preZ[1,], ncol = m, nrow = fit@model$maxdccOrder)
		}
		if(is.null(preQ)){
			preQ = Rbar
		} else{
			preQ = preQ
		}		
	}
	ures = .sample.copula(model, Qbar = Qbar, preQ = preQ, Rbar = Rbar, Nbar = Nbar, 
			preZ = preZ, n.sim = n.sim, n.start = n.start, m.sim = m.sim, 
			rseed = rseed, cluster = cluster)
	# ures = [0,1] copula random numbers
	# now transform back into margins for use in garch	
	transformation = model$modeldesc$transformation
	distu = umodel$modeldesc$distribution
	
	zres = array(NA, dim = c(nsim, m, m.sim))
	# parallel:
	if( !is.null(cluster) ){
			ssfit = model$sfit
			minc = umodel$modelinc
			mpars = model$mpars
			sxres  = fit@mfit$stdresid
			clusterEvalQ(cluster, require('rmgarch'))
			if(transformation == "spd"){
				clusterExport(cluster, c("ures", "mpars", "minc", 
								"m", "sxres", "ssfit", "transformation"), 
						envir = environment())
			} else{
				clusterExport(cluster, c("ures", "mpars", "minc", 
								"m", "sxres", "ssfit", "transformation"), 
						envir = environment())
			}
			clusterExport(cluster, c(".qparametric", ".qempirical",".qspd"), envir = environment())
			mtmp = parLapply(cluster, as.list(1:m.sim), fun = function(i){
						switch(transformation,
								parametric = .qparametric(matrix(ures$Usim[,,i], ncol = m), pars = mpars, modelinc = minc),
								empirical = .qempirical(matrix(ures$Usim[,,i], ncol = m), sxres),
								spd = .qspd(matrix(ures$Usim[,,i], ncol = m), sfit = ssfit))
					})
			for(i in 1:m.sim) zres[,,i] = mtmp[[i]]
	} else{
		mtmp = lapply(as.list(1:m.sim), FUN = function(i) switch(transformation,
							parametric = .qparametric(matrix(ures$Usim[,,i], ncol = m), pars = model$mpars, modelinc = umodel$modelinc),
							empirical = .qempirical(matrix(ures$Usim[,,i], ncol = m), fit@mfit$stdresid),
							spd = .qspd(matrix(ures$Usim[,,i], ncol = m), sfit = model$sfit)))
		for(i in 1:m.sim) zres[,,i] = mtmp[[i]]
	}
	
	
	# Now we have the innovations which we feed back into the garch routine to simulate	
	simRes = simX = simR = simQ = simH = simSeries = vector(mode = "list", length = m.sim)
	#if( is.null(fit@mfit$vrmodel) ){
	#	if(is.null(prereturns)) prereturns = fit@mfit$tailuret
	#}
	
	if( !is.null(cluster) ){
			tailres = fit@mfit$tailuresids
			clusterEvalQ(cluster, loadNamespace('rugarch'))
			clusterExport(cluster, c("mspec", "n.sim", "n.start", "m.sim", 
							"startMethod", "zres", "presigma", "tailres",
							"preresiduals", "prereturns", "model", "mexsimdata", 
							"vexsimdata", "rseed","prerealized"), envir = environment())
			simR = ures$simR
			simlist = parLapply(cluster, as.list(1:m), fun = function(i){
						maxx = mspec@spec[[i]]@model$maxOrder;
						htmp = rugarch::ugarchpath(mspec@spec[[i]], n.sim = n.sim + n.start, n.start = 0, m.sim = m.sim,
								custom.dist = list(name = "sample", distfit = matrix(zres[,i,1:m.sim], ncol = m.sim)),
								presigma = tail(presigma[,i], maxx), 
								preresiduals = if( is.null(preresiduals) ) tail(tailres[,i], maxx) else tail(preresiduals[,i], maxx), 
								prereturns = if(model$modelinc[1]==0) tail(prereturns[,i], maxx) else NA,
								mexsimdata = if(model$modelinc[1]==0) mexsimdata[[i]] else NULL, vexsimdata = vexsimdata[[i]],
								prerealized = tail(prerealized[,i], maxx))
						h = matrix(tail(htmp@path$sigmaSim^2, n.sim), nrow = n.sim)
						x = matrix(htmp@path$seriesSim,  nrow = n.sim + n.start)
						return(list(h = h, x = x))
					})
			if(!only.density){	
				H = array(NA, dim = c(n.sim, m, m.sim))
				tmpH = array(NA, dim = c(m, m, n.sim))
				for(i in 1:n.sim) H[i,,] = t(sapply(simlist, FUN = function(x) as.numeric(x$h[i,])))
				for(i in 1:m.sim){
					for(j in 1:n.sim){
						tmpH[ , , j] = sqrt(diag( H[j, , i]) ) %*% simR[[i]][,,j] %*% sqrt(diag( H[j, , i] ) )
					}
					simH[[i]] = tmpH
				}
			}
			if(model$modelinc[1]>0){
				simxX = array(NA, dim = c(n.sim+n.start, m, m.sim))
				for(i in 1:m.sim) simxX[,,i] = sapply(simlist, FUN = function(x) as.numeric(x$x[,i]))
				simX = vector(mode = "list", length = m.sim)
				for(i in 1:m.sim) simX[[i]] = matrix(simxX[,,i], nrow = n.sim+n.start)
			} else{
				simxX = array(NA, dim = c(n.sim+n.start, m, m.sim))
				for(i in 1:m.sim) simxX[,,i] = sapply(simlist, FUN = function(x) as.numeric(x$x[,i]))
				simX = vector(mode = "list", length = m.sim)
				for(i in 1:m.sim) simX[[i]] = matrix(tail(matrix(simxX[,,i], ncol = m), n.sim), nrow = n.sim)
			}
	} else{
		simR = ures$simR
		#simQ = lapply(mtmp, FUN = function(x) if(is.matrix(x$Q)) array(x$Q, dim = c(m, m, n.sim)) else last(x$Q, n.sim))
		simlist = vector(mode="list", length=m)
		for(i in 1:m){		
			maxx = mspec@spec[[i]]@model$maxOrder
			htmp = ugarchpath(mspec@spec[[i]], n.sim = n.sim + n.start, n.start = 0, m.sim = m.sim,
					custom.dist = list(name = "sample", distfit = matrix(zres[,i,1:m.sim], ncol = m.sim)),
					presigma = tail(presigma[,i], maxx), 
					preresiduals = if( is.null(preresiduals) ) tail(fit@mfit$tailuresids[,i], maxx) else tail(preresiduals[,i], maxx), 
					prereturns = if(model$modelinc[1]==0) tail(prereturns[,i], maxx) else NA,
					mexsimdata = if(model$modelinc[1]==0) mexsimdata[[i]] else NULL, vexsimdata = vexsimdata[[i]],
					prerealized = tail(prerealized[,i], maxx))
			h = matrix(tail(htmp@path$sigmaSim^2, n.sim), nrow = n.sim)
			x = matrix(htmp@path$seriesSim,  nrow = n.sim + n.start)
			simlist[[i]] = list(h = h, x = x)
		}
		if(!only.density){
			H = array(NA, dim = c(n.sim, m, m.sim))
			tmpH = array(NA, dim = c(m, m, n.sim))
			for(i in 1:n.sim) H[i,,] = t(sapply(simlist, FUN = function(x) as.numeric(x$h[i,])))
			for(i in 1:m.sim){
				for(j in 1:n.sim){
					tmpH[ , , j] = sqrt(diag( H[j, , i]) ) %*% simR[[i]][,,j] %*% sqrt(diag( H[j, , i] ) )
				}
				simH[[i]] = tmpH
			}
		}
		if(model$modelinc[1]>0){
			simxX = array(NA, dim = c(n.sim+n.start, m, m.sim))
			for(i in 1:m.sim) simxX[,,i] = sapply(simlist, FUN = function(x) as.numeric(x$x[,i]))
			simX = vector(mode = "list", length = m.sim)
			for(i in 1:m.sim) simX[[i]] = matrix(simxX[,,i], nrow = n.sim+n.start)
		} else{
			simxX = array(NA, dim = c(n.sim+n.start, m, m.sim))
			for(i in 1:m.sim) simxX[,,i] = sapply(simlist, FUN = function(x) as.numeric(x$x[,i]))
			simX = vector(mode = "list", length = m.sim)
			for(i in 1:m.sim) simX[[i]] = matrix(tail(matrix(simxX[,,i], ncol = m), n.sim), nrow = n.sim)
		}
	}
	
	if( model$modelinc[1]>0 ){
		simRes = simX
		simX = mvmean.varsim(model = model, Data = Data, res = simX, 
				mexsimdata = mexsimdata, prereturns = prereturns, m.sim = m.sim, 
				n.sim = n.sim, n.start = n.start, startMethod = startMethod, 
				cluster = cluster)
	} else{
		# reshape
		for(j in 1:m.sim) simX[[j]] = tail(simX[[j]], n.sim)
	}
	msim = list()
	
	if(only.density){
		msim$simX = simX
		msim$rseed = rseed
		model$n.sim = n.sim
		model$m.sim = m.sim
		model$n.start = n.start
		model$startMethod = startMethod[1]
	} else{
		msim$simH = simH
		msim$simR = simR
		msim$simX = simX
		msim$simRes = simRes
		msim$simZ = zres
		msim$rseed = rseed
		model$n.sim = n.sim
		model$m.sim = m.sim
		model$n.start = n.start
		model$startMethod = startMethod[1]
	}
	model$only.density = only.density
	# need to run gc and rm since 'memory accumulates here'
	if(exists("simX")) rm(simX)
	if(exists("simRes")) rm(simRes)
	if(exists("simH")) rm(simH)
	if(exists("H")) rm(H)
	if(exists("simlist")) rm(simlist)
	if(exists("tmpH")) rm(tmpH)
	if(exists("mtmp")) rm(mtmp)
	if(exists("simxX")) rm(simxX)
	if(exists("zres")) rm(zres)
	if(exists("ures")) rm(ures)
	if(exists("simR")) rm(simR)
	
	gc(verbose = FALSE)
	
	ans = new("cGARCHsim",
			msim = msim,
			model = model)
	return( ans )
	
}

.fullinc2 = function(modelinc, umodel){
	m = dim(umodel$modelinc)[2]
	# the par matrix will have the max parameter vis Asset (if an asset does not have
	# the parameter it will be zero...for matrix estimation.
	# the index places the solver estimates into their place in the matrix
	###################################################################################
	#           Asset1 .... Asset2 .... AssetN ...   ........     Joint
	# mu
	# ar
	# ma
	# arfima
	# archm
	# mxreg
	# omega
	# alpha
	# beta
	# gamma
	# delta
	# lambda
	# eta1
	# eta2
	# vxreg
	# skew
	# shape
	# ghlambda
	# xi
	# C
	# dcca
	# dccb
	# dccg
	# mshape
	# mskew
	
	vecmax = rep(0, 19)
	names(vecmax) = rownames(umodel$modelinc[1:19,])
	vecmax = apply(umodel$modelinc, 1, FUN = function(x) max(x) )
	maxOrder = apply(umodel$modelinc, 2, FUN = function(x) max(c(x[2], x[3], x[8], x[9])))
	sumv = 19 + sum(pmax(1, vecmax[c(2,3,6,8,9,10,11,12,13,15,16)])) - 11
	tmpmat = matrix(0, ncol = m+1, nrow = sumv)
	nx = 0
	pnames = NULL
	# mu
	if(vecmax[1]>0){
		tmpmat[1, 1:m] = umodel$modelinc[1,]
	}
	nx = nx + max(1, vecmax[1])
	pnames = c(pnames, "mu")
	# ar
	if(vecmax[2]>0){
		for(i in 1:vecmax[2]){
			tmpmat[nx+i, 1:m] = as.integer( umodel$modelinc[2,] >= i)
			pnames = c(pnames, paste("ar", i, sep = ""))
		}
	} else{
		pnames = c(pnames, "ar")
	}
	nx = nx + max(1, vecmax[2])
	# ma
	if(vecmax[3]>0){
		for(i in 1:vecmax[3]){
			tmpmat[nx+i, 1:m] = as.integer( umodel$modelinc[3,] >= i)
			pnames = c(pnames, paste("ma", i, sep = ""))
		}
	} else{
		pnames = c(pnames, "ma")
	}
	nx = nx + max(1, vecmax[3])
	# arfima
	if(vecmax[4]>0){
		tmpmat[nx+1, 1:m] = umodel$modelinc[4, ]
	}
	nx = nx + max(1, vecmax[4])
	pnames = c(pnames, "arfima")
	# archm no allowed
	nx = nx + max(1, vecmax[5])
	pnames = c(pnames, "archm")
	# mxreg
	if(vecmax[6]>0){
		for(i in 1:vecmax[6]){
			tmpmat[nx+i, 1:m] = as.integer(umodel$modelinc[6, ] >= i)
			pnames = c(pnames, paste("mxreg", i, sep = ""))
		}
	} else{
		pnames = c(pnames, "mxreg")
	}
	nx = nx + max(1, vecmax[6])
	# omega
	if(vecmax[7]>0){
		tmpmat[nx+1, 1:m] = umodel$modelinc[7, ]
	}
	nx = nx + max(1, vecmax[7])
	pnames = c(pnames, "omega")
	# alpha
	if(vecmax[8]>0){
		for(i in 1:vecmax[8]){
			tmpmat[nx+i, 1:m] = as.integer( umodel$modelinc[8, ] >= i)
			pnames = c(pnames, paste("alpha", i, sep = ""))
		}
	} else{
		pnames = c(pnames, "alpha")
	}
	nx = nx + max(1, vecmax[8])
	# beta
	if(vecmax[9]>0){
		for(i in 1:vecmax[9]){
			tmpmat[nx+i, 1:m] = as.integer( umodel$modelinc[9, ] >= i)
			pnames = c(pnames, paste("beta", i, sep = ""))
		}
	} else{
		pnames = c(pnames, "beta")
	}
	nx = nx + max(1, vecmax[9])
	# gamma
	if(vecmax[10]>0){
		for(i in 1:vecmax[10]){
			tmpmat[nx+i, 1:m] = as.integer( umodel$modelinc[10, ] >= i)
			pnames = c(pnames, paste("gamma", i, sep = ""))
		}
	} else{
		pnames = c(pnames, "gamma")
	}
	nx = nx + max(1, vecmax[10])
	# eta1
	if(vecmax[11]>0){
		for(i in 1:vecmax[11]){
			tmpmat[nx+i, 1:m] = as.integer( umodel$modelinc[11, ] >= i)
			pnames = c(pnames, paste("eta1", i, sep = ""))
		}
	} else{
		pnames = c(pnames, "eta1")
	}
	nx = nx + max(1, vecmax[11])
	# eta2
	if(vecmax[12]>0){
		for(i in 1:vecmax[12]){
			tmpmat[nx+i, 1:m] = as.integer( umodel$modelinc[12, ] >= i)
			pnames = c(pnames, paste("eta2", i, sep = ""))
		}
	} else{
		pnames = c(pnames, "eta2")
	}
	nx = nx + max(1, vecmax[12])
	# delta
	if(vecmax[13]>0){
		tmpmat[nx+1, 1:m] = umodel$modelinc[13, ]
	}
	nx = nx + max(1, vecmax[13])
	pnames = c(pnames, "delta")
	# lambda
	if(vecmax[14]>0){
		tmpmat[nx+1, 1:m] = umodel$modelinc[14, ]
	}
	nx = nx + max(1, vecmax[14])
	pnames = c(pnames, "lambda")
	# vxreg
	if(vecmax[15]>0){
		for(i in 1:vecmax[15]){
			tmpmat[nx+i, 1:m] = as.integer( umodel$modelinc[15, ] >= i)
			pnames = c(pnames, paste("vxreg", i, sep = ""))
		}
	} else{
		pnames = c(pnames, "vxreg")
	}
	nx = nx + max(1, vecmax[15])
	# skew
	if(vecmax[16]>0){
		tmpmat[nx+1, 1:m] = umodel$modelinc[16, ]
	}
	nx = nx + max(1, vecmax[16])
	pnames = c(pnames, "skew")
	# shape
	if(vecmax[17]>0){
		tmpmat[nx+1, 1:m] = umodel$modelinc[17, ]
	}
	nx = nx + max(1, vecmax[17])
	pnames = c(pnames, "shape")
	# skew
	if(vecmax[18]>0){
		tmpmat[nx+1, 1:m] = umodel$modelinc[18, ]
	}
	nx = nx + max(1, vecmax[18])
	pnames = c(pnames, "ghlambda")
	
	if(vecmax[19]>0){
		tmpmat[nx+1, 1:m] = umodel$modelinc[19, ]
	}
	nx = nx + max(1, vecmax[19])
	pnames = c(pnames, "xi")
	
	sumdcc = 6 + sum(pmax(1, modelinc[c(3,4,5,6,7,8)])) - 6
	tmpmat = rbind(tmpmat, matrix(0, ncol = m+1, nrow = sumdcc))
	if(modelinc[3]>0){
		for(i in 1:modelinc[3]){
			tmpmat[nx+i, m+1] = 1
			pnames = c(pnames, paste("C", i, sep = ""))
		}
	} else{
		pnames = c(pnames, "C")
	}
	nx = nx + max(1, modelinc[3])
	
	if(modelinc[4]>0){
		for(i in 1:modelinc[4]){
			tmpmat[nx+i, m+1] = 1
			pnames = c(pnames, paste("dcca", i, sep = ""))
		}
	} else{
		pnames = c(pnames, "dcca")
	}
	
	nx = nx + max(1, modelinc[4])
	
	if(modelinc[5]>0){
		for(i in 1:modelinc[5]){
			tmpmat[nx+i, m+1] = 1
			pnames = c(pnames, paste("dccb", i, sep = ""))
		}
	} else{
		pnames = c(pnames, "dccb")
	}
	nx = nx + max(1, modelinc[5])
	
	if(modelinc[6]>0){
		for(i in 1:modelinc[6]){
			tmpmat[nx+i, m+1] = 1
			pnames = c(pnames, paste("dccg", i, sep = ""))
		}
	} else{
		pnames = c(pnames, "dccg")
	}
	nx = nx + max(1, modelinc[6])
	
	# we allow for possibility of vector valued shape and skew parameters for future expansion
	if(modelinc[7]>0){
		for(i in 1:modelinc[7]){
			tmpmat[nx+i, m+1] = 1
			if(modelinc[7]>1) pnames = c(pnames, paste("mshape", i, sep = "")) else pnames = c(pnames, "mshape")
		}
	} else{
		pnames = c(pnames, "mshape")
	}
	nx = nx + max(1, modelinc[7])
	
	if(modelinc[8]>0){
		for(i in 1:modelinc[8]){
			tmpmat[nx+i, m+1] = 1
			if(modelinc[8]>1) pnames = c(pnames, paste("mskew", i, sep = "")) else pnames = c(pnames, "mskew")
		}
	} else{
		pnames = c(pnames, "mskew")
	}
	# NOW we have and indicator matrix
	colnames(tmpmat) = c(paste("Asset", 1:m, sep = ""), "Joint")
	rownames(tmpmat) = pnames
	return(tmpmat)
}