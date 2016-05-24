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
.dccspec = function(uspec, VAR = FALSE, robust = FALSE, lag = 1, lag.max = NULL, 
				lag.criterion = c("AIC", "HQ", "SC", "FPE"), external.regressors = NULL, 
				robust.control = list("gamma" = 0.25, "delta" = 0.01, "nc" = 10, "ns" = 500), 
		dccOrder = c(1,1), asymmetric = FALSE, distribution = c("mvnorm", "mvt", "mvlaplace"), 
		start.pars = list(), fixed.pars = list())
{	
	VAR.opt = list()
	if(is.null(VAR)) VAR.opt$VAR = FALSE else VAR.opt$VAR = as.logical(VAR)
	if(is.null(robust)) VAR.opt$robust = FALSE else VAR.opt$robust = as.logical(robust)
	if(is.null(lag)) VAR.opt$lag = 1 else VAR.opt$lag = as.integer(lag)
	if(is.null(lag.max)) VAR.opt$lag.max = NULL else VAR.opt$lag.max = as.integer(min(1, lag.max))
	if(is.null(lag.criterion)) VAR.opt$lag.criterion = "AIC" else VAR.opt$lag.criterion = lag.criterion[1]
	if(is.null(external.regressors)) VAR.opt$external.regressors = NULL else VAR.opt$external.regressors = external.regressors
	
	rc = list("gamma" = 0.25, "delta" = 0.01, "nc" = 10, "ns" = 500)
	rcmatch = match(names(robust.control), c("gamma", "delta", "nc", "ns"))
	if(length(rcmatch[!is.na(rcmatch)]) > 0){
		rx = which(!is.na(rcmatch))
		rc[rcmatch[!is.na(rcmatch)]] = robust.control[rx]
	}
	VAR.opt$robust.control = rc
	.eps = .Machine$double.eps
	modeldata = list()
	modeldesc = list()
	m = length(uspec@spec)
	
	if(is.null(distribution)) distribution = "mvnorm"
	distribution = distribution[1]
	valid.distributions = c("mvnorm", "mvt", "mvlaplace")
	if(!any(distribution == valid.distributions)) stop("\nInvalid Distribution Choice\n", call. = FALSE)
	
	modelinc = rep(0, 10)
	names(modelinc) = c("var", "mvmxreg", "dcca", "dccb", "dccg", "mshape", "mskew", "aux", "aux", "aux")
	# if we ever introduce a distribution with vector shape all we need to do is increase modelinc[6] to
	# nl i.e. no of assets...same for skew
	if(distribution == "mvt") modelinc[6] = 1
	if(is.null(dccOrder)){
		modelinc[3] = 1
		modelinc[4] = 1
	} else{
		modelinc[3] = as.integer( dccOrder[1] )
		modelinc[4] = as.integer( dccOrder[2] )
	}
	# Asymemtry parameter follows ARCH parameter in lags
	if( asymmetric ) modelinc[5] = modelinc[3]
	if( VAR ){
		if(is.null(VAR.opt$lag)) modelinc[1] = 1 else modelinc[1] = as.integer( VAR.opt$lag )
		if(!is.null(VAR.opt$external.regressors)){
			if(!is.matrix(VAR.opt$external.regressors)) stop("\nexternal.regressors must be a matrix.")
			modelinc[2] = dim(VAR.opt$external.regressors)[2]
			modeldata$mexdata = VAR.opt$external.regressors
		} else{
			modeldata$mexdata = NULL
		}
	}
	modelinc[10] = which(c("mvnorm", "mvt", "mvlaplace") == distribution)
	maxdccOrder = max(dccOrder)
	#modeldata = modeldata
	modeldesc$distribution = distribution
	modeldesc$dccmodel = ifelse(asymmetric, "ADCC", "DCC")
	if( !is(uspec, "uGARCHmultispec") ) stop("\ndccspec-->error: uspec must be a uGARCHmultispec object")
	
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
			umodel$modeldesc$vsubmodel[i] = ifelse(is.null(uspec@spec[[i]]@model$modeldesc$vsubmodel),"GARCH",uspec@spec[[i]]@model$modeldesc$vsubmodel)
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
	
	maxgarchOrder = max( sapply(uspec@spec, FUN = function(x) x@model$maxOrder) )
	
	# Note: We run varxfit with the postpad option using "constant" so that
	# it resembles the ARMA-GARCH type functionality
	if(modelinc[1]>0){
		maxgarchOrder = max(c(maxgarchOrder, modelinc[1]))
	}
	
	################################################################
	# pars list
	pars = matrix(0, ncol = 6, nrow = 5)
	colnames(pars) = c("Level", "Fixed", "Include", "Estimate", "LB", "UB")
	pidx = matrix(NA, nrow = 5, ncol = 2)
	colnames(pidx) = c("begin", "end")
	rownames(pidx) =  c("dcca", "dccb", "dccg", "mshape", "mskew")
	
	pos = 1
	pos.matrix = matrix(0, ncol = 3, nrow = 5)
	colnames(pos.matrix) = c("start", "stop", "include")
	rownames(pos.matrix) = c("dcca", "dccb", "dccg", "mshape", "mskew")
	
	for(i in 1:5){
		if( modelinc[2+i] > 0 ){
			pos.matrix[i,1:3] = c(pos, pos+modelinc[2+i]-1, 1)
			pos = max(pos.matrix[1:i,2]+1)
		}
	}
	
	mm = sum(modelinc[3:7])
	mm = mm - length( which(modelinc[c(3:7)]>0) )
	pars = matrix(0, ncol = 6, nrow = 5 + mm)
	colnames(pars) = c("Level", "Fixed", "Include", "Estimate", "LB", "UB")
	pidx = matrix(NA, nrow = 5, ncol = 2)
	colnames(pidx) = c("begin", "end")
	rownames(pidx) =  c("dcca", "dccb", "dccg", "mshape", "mskew")
	fixed.names = names(fixed.pars)
	start.names = names(start.pars)
	fixed.pars = unlist(fixed.pars)
	start.pars = unlist(start.pars)
	
	pn = 1
	pnames = NULL
	nx = 0
	if(pos.matrix[1,3] == 1){
		pn = length( seq(pos.matrix[1,1], pos.matrix[1,2], by = 1) )
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
	pidx[1,1] = 1
	pidx[1,2] = pn
	nx = pn
	pn = 1
	pidx[2,1] =  nx+1
	if(pos.matrix[2,3] == 1){
		pn = length( seq(pos.matrix[2,1], pos.matrix[2,2], by = 1) )
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
	pidx[2,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[3,1] = nx+1
	
	if(pos.matrix[3,3] == 1){
		pn = length( seq(pos.matrix[3,1], pos.matrix[3,2], by = 1) )
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
	pidx[3,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[4,1] = nx+1
	if(modelinc[6]<=1){
		if(pos.matrix[4,3]==1){
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
		if(pos.matrix[4,3] == 1){
			pn = length( seq(pos.matrix[4,1], pos.matrix[4,2], by = 1) )
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
	pidx[4,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[5,1] = nx+1
	
	if(modelinc[7]<=1){
		if(pos.matrix[5,3]==1){
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
		if(pos.matrix[5,3] == 1){
			pn = length( seq(pos.matrix[5,1], pos.matrix[5,2], by = 1) )
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
	pidx[5,2] = nx+pn
	rownames(pars) = pnames
	#zf = match(fixed.names, rownames(pars))
	#if( length(zf)>0 ) pars[zf, 1] = unlist(fixed.pars)
	modeldesc$type = "2-step"
	model = list(modelinc = modelinc, modeldesc = modeldesc, modeldata = modeldata, varmodel = varmodel, 
			pars = pars, start.pars = start.pars, fixed.pars = fixed.pars, maxgarchOrder = maxgarchOrder, 
			maxdccOrder = maxdccOrder, pos.matrix = pos.matrix, pidx = pidx)
	model$DCC = ifelse(asymmetric, "aDCC", "DCC")
	
	ans = new("DCCspec",
			model = model,
			umodel = umodel)
	return(ans)
}

.dccfit = function(spec, data, out.sample = 0, solver = "solnp", 
		solver.control = list(), fit.control = list(eval.se = TRUE, 
				stationarity = TRUE, scale = FALSE), 
		cluster = NULL, fit = NULL, VAR.fit = NULL, verbose = FALSE, 
		realizedVol = NULL, ...)
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
	# allow first and second stage solvers to be defined
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
	data   = xdata$data
	index  = xdata$index
	period = xdata$period
	# save the data to the model spec
	model$modeldata$data = data
	model$modeldata$index = index
	model$modeldata$period = period
	T = model$modeldata$T = n - n.start
	model$modeldata$n.start = n.start
	model$modeldata$asset.names = cnames
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
		if(umodel$modeldesc$vmodel[1]=="realGARCH") plik = sapply(fitlist@fit, function(x) sum(-x@fit$partial.log.likelihoods)) else plik =  sapply(fitlist@fit, function(x) sum(-x@fit$log.likelihoods))
	} else{	
		fitlist = multifit(multispec = mspec, data = xts(zdata, index), out.sample = n.start, 
				solver = garch.solver, solver.control = solver.control, 
				fit.control = ufit.control, cluster = cluster, realizedVol = realizedVol)
		#assign(".fitlist", fitlist, envir = .GlobalEnv)
		#p = max(sapply(fitlist@fit, FUN = function(x) max(x@model$modelinc[2:3])))
		
		converge = sapply(fitlist@fit, FUN = function(x) x@fit$convergence)
		if( any( converge == 1 ) ){
			pr = which(converge != 1)
			cat("\nNon-Converged:\n")
			print(pr)
			cat("\ndccfit-->error: convergence problem in univariate fit...")
			cat("\n...returning uGARCHmultifit object instead...check and resubmit...")
			return( fitlist )
		}
		if(umodel$modeldesc$vmodel[1]=="realGARCH") plik = sapply(fitlist@fit, function(x) sum(-x@fit$partial.log.likelihoods)) else plik = sapply(fitlist@fit, function(x) sum(-x@fit$log.likelihoods))
		if(spec@model$modelinc[1]>0) model$mu = mu else model$mu = fitted(fitlist)
		model$residuals = res = residuals(fitlist)
		model$sigma = sig = sigma(fitlist)
	}
	
	stdresid = res/sig
	# In the VAR model we postpadded the first lag values with a constant.
	# In an ARMA model, the first lag values are constant.
	#-----------------------------------------------------------------------------------
	# Create the Full Model Pars and Indices
	modelinc =  model$modelinc
	# create full par matrix
	midx = .fullinc(modelinc, umodel)
	midx["omega",1:m]=1
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
	#-----------------------------------------------------------------------------------
	Qbar = cov(stdresid)
	# Take care of the Asymmetry Matrices
	if(modelinc[5]>0){
		Ibar = .asymI(stdresid)
		astdresid = Ibar*stdresid
		Nbar = cov(astdresid)
	} else{
		Ibar = .asymI(stdresid)
		# zero what is not needed
		astdresid = Ibar*stdresid*0
		Nbar = matrix(0, m, m)
	}
	H = sig^2
	#-----------------------------------------------------------------------------------
	# create a temporary environment to store values (deleted at end of function)
	mgarchenv = new.env(hash = TRUE)
	arglist = list()
	arglist$mgarchenv = mgarchenv
	arglist$verbose = verbose
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
	arglist$astdresid = astdresid
	arglist$Ibar = Ibar
	arglist$Qbar = Qbar
	arglist$Nbar = Nbar
	arglist$H = H
	#
	#-----------------------------------------------------------------------------------	
	
	# Check for fixed parameters in DCC
	if(any(ipars[,2]==1)){
		if(npars == 0){
			if(fit.control$eval.se==0) {
				warning("\ndccfit-->warning: all parameters fixed...returning dccfilter object instead\n")
				xspex = spec
				for(i in 1:m) xspex@umodel$fixed.pars[[i]] = as.list(fitlist@fit[[i]]@model$pars[fitlist@fit[[i]]@model$pars[,3]==1,1])
				return(dccfilter(spec = xspex, data = xts(data, index), out.sample = out.sample, 
								cluster = cluster, VAR.fit = VAR.fit, , realizedVol = realizedVol, ...))
			} else{
				# if all parameters are fixed but we require standard errors, we
				# skip the solver
				use.solver = 0
				ipars[ipars[,2]==1, 4] = 1
				ipars[ipars[,2]==1, 2] = 0
				arglist$pars = ipars
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
	assign("rmgarch_llh", 1, envir = mgarchenv)
	# get
	ILB = 0
	IUB = 1
	# Asymmetric Spec has different constraints
	if(model$modelinc[5]> 0) Ifn = .adcccon else Ifn = .dcccon
	if( solver == "solnp" | solver == "gosolnp") fit.control$stationarity = FALSE else fit.control$stationarity = TRUE
	arglist$fit.control = fit.control
	
	if( use.solver )
	{
		arglist$returnType = "llh"
		solution = switch(model$modeldesc$distribution,
				mvnorm = .dccsolver(solver, pars = ipars[estidx, 1], fun = normal.dccLLH1, Ifn, ILB, 
						IUB, gr = NULL, hessian = NULL, control = solver.control, 
						LB = ipars[estidx, 5], UB = ipars[estidx, 6], arglist = arglist),
				mvlaplace = .dccsolver(solver, pars = ipars[estidx, 1], fun = laplace.dccLLH1, Ifn, ILB, 
						IUB, gr = NULL, hessian = NULL, control = solver.control, 
						LB = ipars[estidx, 5], UB = ipars[estidx, 6], arglist = arglist),
				mvt = .dccsolver(solver, pars = ipars[estidx, 1], fun = student.dccLLH1, Ifn, ILB, 
						IUB, gr = NULL, hessian = NULL, control = solver.control, 
						LB = ipars[estidx, 5], UB = ipars[estidx, 6], arglist = arglist))
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
		convergence = 0
		sol = list()
		sol$message = "all parameters fixed"
	}
	
	fit = list()
	# check convergence else write message/return
	if( convergence == 0 ){
		fit = switch(model$modeldesc$distribution,
				mvnorm = .dccmakefitmodel(garchmodel = "dccnorm", f = normal.dccLLH2, 
						arglist = arglist, timer = 0, message = sol$message, fname = "normal.dccLLH2"),
				mvlaplace = .dccmakefitmodel(garchmodel = "dcclaplace", f = laplace.dccLLH2, 
						arglist = arglist, timer = 0, message = sol$message, fname = "laplace.dccLLH2"),		
				mvt =.dccmakefitmodel(garchmodel = "dccstudent", f = student.dccLLH2, 
						arglist = arglist, timer = 0, message = sol$message, fname = "student.dccLLH2"))
		fit$timer = Sys.time() - tic
	} else{
		fit$message = sol$message
		fit$convergence = 1
	}
	fit$Nbar = Nbar
	fit$Qbar = Qbar
	fit$realizedVol = realizedVol
	fit$plik = plik
	# make model list to return some usefule information which
	model$mpars = mpars
	model$ipars = ipars
	model$pars[,1] = ipars[,1]
	model$midx = midx
	model$eidx = eidx
	model$umodel = umodel
	ans = new("DCCfit",
			mfit = fit,
			model = model)
	return(ans)
}

.dccfilter = function(spec, data, out.sample = 0, filter.control = list(n.old = NULL), 
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
	data   = xdata$data
	index  = xdata$index
	period = xdata$period
	# save the data to the model spec
	model$modeldata$data = data
	model$modeldata$index = index
	model$modeldata$period = period
	T = model$modeldata$T = n - n.start
	model$modeldata$n.start = n.start
	model$modeldata$asset.names = cnames
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
			umodel$fixed.pars, NULL)
	
	filterlist = multifilter(multifitORspec = mspec, data = xts(zdata, index[1:nrow(zdata)]), 
			out.sample = out.sample, cluster = cluster, n.old = n.old, , realizedVol = realizedVol, ...)
	
	if(spec@model$modelinc[1]>0) model$mu = mu else model$mu = fitted(filterlist)
	model$residuals = res = residuals(filterlist)
	model$sigma = sig = sigma(filterlist)
	stdresid = res/sig
	
	if(is.null(filter.control$n.old)) dcc.old = dim(stdresid)[1] else dcc.old = n.old
	#-----------------------------------------------------------------------------------	
	# Create the Full Model Pars and Indices
	
	modelinc =  model$modelinc
	# create full par matrix
	midx = .fullinc(modelinc, umodel)
	mpars = midx*0
	eidx = midx
	# This is the estimated parameter index
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
	#-----------------------------------------------------------------------------------
	

	Qbar = cov(stdresid[1:dcc.old, ])
	# Take care of the Asymmetry Matrices
	if(modelinc[5]>0){
		Ibar = .asymI(stdresid)
		astdresid = Ibar*stdresid
		Nbar = cov(astdresid[1:dcc.old, ])
		
	} else{
		Ibar = .asymI(stdresid)
		# zero what is not needed
		astdresid = Ibar*stdresid*0
		Nbar = matrix(0, m, m)
	}
	H = sig^2
	
	#-----------------------------------------------------------------------------------
	mgarchenv = new.env(hash = TRUE)
	arglist = list()
	arglist$mgarchenv = mgarchenv
	arglist$verbose = FALSE
	arglist$cluster = cluster
	arglist$filter.control = filter.control
	arglist$cnames = cnames
	arglist$m = m
	arglist$T = T
	arglist$n.old = n.old
	arglist$dcc.old = dcc.old
	arglist$data = zdata
	arglist$index = index
	arglist$model = model
	arglist$filterlist = filterlist
	arglist$realizedVol = realizedVol
	arglist$umodel = umodel
	arglist$midx = midx
	arglist$eidx = eidx
	arglist$mpars = mpars
	arglist$ipars = ipars
	arglist$estidx = estidx
	arglist$dccN = npars
	arglist$stdresid = stdresid
	arglist$astdresid = astdresid
	arglist$Ibar = Ibar
	arglist$Qbar = Qbar
	arglist$Nbar = Nbar
	arglist$H = H
	assign("rmgarch_llh", 1, envir = mgarchenv)
	#-----------------------------------------------------------------------------------	
	
	filt = switch(model$modeldesc$distribution,
			mvnorm = .dccmakefiltermodel(garchmodel = "dccnorm", f = normalfilter.dccLLH2, 
					arglist = arglist, timer = 0, message = 0, fname = "normalfilter.dccLLH2"),
			mvlaplace = .dccmakefiltermodel(garchmodel = "dcclaplace", f = laplacefilter.dccLLH2, 
					arglist = arglist, timer = 0, message = 0, fname = "laplacefilter.dccLLH2"),
			mvt = .dccmakefiltermodel(garchmodel = "dccstudent", f = studentfilter.dccLLH2, 
					arglist = arglist, timer = 0, message = 0, fname = "studentfilter.dccLLH2"))
	
	model$mpars = mpars
	model$ipars = ipars
	model$pars[,1] = ipars[,1]
	model$midx = midx
	model$eidx = eidx
	model$umodel = umodel
	
	filt$Nbar = Nbar
	filt$Qbar = Qbar
	filt$realizedVol = realizedVol
	
	ans = new("DCCfilter",
			mfilter = filt,
			model = model)
	return(ans)
}

.dccforecast = function(fit, n.ahead = 1, n.roll = 0, external.forecasts = list(mregfor = NULL, vregfor = NULL), 
		cluster = NULL, ...)
{
	# checks for out.sample wrt n.roll
	model = fit@model
	modelinc = model$modelinc
	ns = model$modeldata$n.start
	if( n.roll > ns ) stop("n.roll must not be greater than out.sample!")
	if(n.roll>1 && n.ahead>1) stop("\ngogarchforecast-->error: n.ahead must be equal to 1 when using n.roll\n")
	if( fit@model$modelinc[5] > 0 && n.ahead > 1) stop("\ngogarchforecast-->error: asymmetric DCC specification only support n.ahead = 1 currently.\n")
	# checks for external forecasts
	tf = n.ahead + n.roll
	if( !is.null( external.forecasts$mregfor ) ){
		mregfor = external.forecasts$mregfor
		if( !is.matrix(mregfor) ) stop("\nmregfor must be a matrix.")
		if( dim(mregfor)[1] < tf ) stop("\nmregfor must have at least n.ahead + n.roll observations to be used")
		mregfor = mregfor[1:tf, , drop = FALSE]
	} else{
		mregfor = NULL
	}
	if( !is.null( external.forecasts$vregfor ) ){
		if( !is.matrix(vregfor) ) stop("\nvregfor must be a matrix.")
		if( dim(vregfor)[1] < tf ) stop("\nvregfor must have at least n.ahead + n.roll observations to be used")
		vregfor = vregfor[1:tf, , drop = FALSE]
	}
	if( modelinc[1]>0 ){
		if( modelinc[2] > 0 ){
			if( is.null(external.forecasts$mregfor ) ){
				warning("\nExternal Regressor Forecasts Matrix NULL...setting to zero...\n")
				mregfor = matrix(0, ncol = modelinc[2], nrow = (n.roll + n.ahead) )
			} else{
				if( dim(mregfor)[2] != modelinc[2] ) stop("\ndccforecast-->error: wrong number of external regressors!...", call. = FALSE)
				if( dim(mregfor)[1] < (n.roll + n.ahead) ) stop("\ndccforecast-->error: external regressor matrix has less points than requested forecast length (1+n.roll) x n.ahead!...", call. = FALSE)
			}
		} else{
			mregfor = NULL
		}
		if(n.roll>1 && n.ahead>1) stop("\ndccforecast-->error: n.ahead must be equal to 1 when using n.roll\n")
		
		if( n.ahead == 1 && (n.roll > ns) ) stop("\ndccforecast-->error: n.roll greater than out.sample!", call. = FALSE)
		#VARf = .varpredict(fit, n.ahead, n.roll, mregfor)$Mu
		mu = varxforecast(X = fit@model$modeldata$data, Bcoef = model$varcoef, p = modelinc[1], 
				out.sample = ns, n.ahead = n.ahead, n.roll = n.roll, mregfor = mregfor)
	} else{
		mu = NULL
	}
	exf = external.forecasts
	if(  modelinc[1] > 0 ){
		exf$mregfor = NULL
	}
	ans = .dccforecastm(fit, n.ahead = n.ahead, n.roll = n.roll, external.forecasts = exf, cluster = cluster, realizedVol = fit@mfit$realizedVol, ...)
	
	if(modelinc[1]==0) mu = ans$mu
	
	model$n.roll = n.roll
	model$n.ahead = n.ahead
	# keep H
	model$H = rcov(fit)
	mforecast = list( H = ans$H, R = ans$R, Q = ans$Q, Rbar = ans$Rbar, mu = mu )
	
	ans = new("DCCforecast",
			mforecast = mforecast,
			model = model)
	
	return( ans )
}


# n-ahead forecast based on the approximation method described in:
# ENGLE, R. and SHEPPARD (2001), K. Theoretical and empirical properties of
# dynamic conditional correlation multivariate garch. NBER Working Papers, No. 8554

.dccforecastm = function(fit, n.ahead = 1, n.roll = 10, external.forecasts = list(mregfor = NULL, vregfor = NULL), 
		cluster = NULL, realizedVol = NULL, ...)
{
	model = fit@model
	modelinc = model$modelinc
	umodel = model$umodel
	m = dim(umodel$modelinc)[2]
	ns = fit@model$modeldata$n.start
	# first check whether the data needs to be filtered by VAR:
	Data = fit@model$modeldata$data
	if(modelinc[1]>0){
		zdata = varxfilter(Data, p = model$modelinc[1], Bcoef = model$varcoef, 
				exogen = fit@model$modeldata$mexdata, postpad = c("constant"))$xresiduals
	} else{
		zdata = Data
	}
	fpars = lapply(1:m, FUN = function(i) fit@model$mpars[fit@model$midx[,i]==1,i])
	
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			fpars, NULL)
	filterlist = multifilter(multifitORspec = mspec, data = xts(zdata, fit@model$modeldata$index[1:nrow(zdata)]), out.sample = 0, 
			n.old = fit@model$modeldata$T, cluster = cluster, realizedVol = realizedVol)
	
	# all.equal(head(sigma(filterlist)), head(fit@model$sigma) )
	n.roll = n.roll + 1
	m = length(mspec@spec)
	out.sample = fit@model$modeldata$n.start
	mo = max(fit@model$maxdccOrder)
	# we augmented n.roll before, now subtract
	# forclist = multiforecast(multifitORspec = fitlist, n.ahead = n.ahead, n.roll = n.roll - 1)
	forclist = multiforecast(multifitORspec = mspec, data = xts(zdata, fit@model$modeldata$index[1:nrow(zdata)]), n.ahead = n.ahead, 
			out.sample = ns, n.roll = n.roll - 1, external.forecasts = external.forecasts, 
			cluster = cluster, realizedVol = realizedVol, ...)
	if(modelinc[1] == 0){
		mu = array(NA, dim=c(n.ahead, m, n.roll))
		f = lapply(forclist@forecast, function(x) fitted(x))
		for(i in 1:n.roll) mu[,,i] = matrix( sapply( f, function(x) x[,i] ), ncol = m)
	} else{
		mu = NULL
	}
	sig = sigma(filterlist)
	resid = residuals(filterlist)
	stdresid = resid/sig
	
	if(modelinc[5]>0){
		Ibar = .asymI(stdresid)
		astdresid = Ibar*stdresid		
	} else{
		Ibar = .asymI(stdresid)
		# zero what is not needed
		astdresid = Ibar*stdresid*0
	}
	T = dim(fit@mfit$H)[3]
	## call .nsdccforecast with rolling window
	Rbar = Rtfor = Htfor = Qtfor = vector(mode = "list", length = n.roll)
	Qstart = last( rcor(fit, type = "Q"), mo )
	Rstart = last( rcor(fit, type = "R"), mo )
	Hstart = last( rcov(fit), mo)
	
	f = lapply(forclist@forecast, function(x) sigma(x))
	for(i in 1:n.roll){
		xQbar = cov(stdresid[1:(T + i - 1), ])
		if(modelinc[5]>0) xNbar = cov(astdresid[1:(T + i - 1), ]) else xNbar = matrix(0, m, m)
		xstdresids = stdresid[(T - mo + i ):(T  +  i - 1), , drop = FALSE]
		xastdresids = astdresid[(T - mo + i ):(T  +  i - 1), , drop = FALSE]
		xfsig = matrix( sapply( f, function(x) x[,i] ), ncol = m)
		ans = .dccforecastn(model, Qstart, Rstart, Hstart, xQbar, xNbar, xstdresids, xastdresids, xfsig, n.ahead, mo)
		Rtfor[[i]] = ans$Rtfor
		Qtfor[[i]] = ans$Qtfor
		Htfor[[i]] = ans$Htfor
		Rbar[[i]] = ans$Rbar
		Qstart = last( rugarch:::.abind(Qstart, ans$Qtfor[, , 1]), mo )
		Rstart = last( rugarch:::.abind(Rstart, ans$Rtfor[, , 1]), mo )
		Hstart = last( rugarch:::.abind(Hstart, ans$Htfor[, , 1]), mo )
	}
	
	forc = list( H = Htfor, R = Rtfor, Q = Qtfor, Rbar = Rbar, mu = mu )
	
	return(forc)
}

.dccforecastn = function(model, Qstart, Rstart, Hstart, Qbar, Nbar, stdresids, 
		astdresids, fsig, n.ahead, mo)
{
	m = dim(Qbar)[1]
	modelinc = model$modelinc
	Qtfor = Rtfor = Htfor = array(NA, dim = c(m, m, n.ahead + mo))
	Qtfor[ , , 1:mo]  =  Qstart[, , 1:mo]
	Rtfor[ , , 1:mo]  =  Rstart[, , 1:mo]
	Htfor[ , , 1:mo] = Hstart[, , 1:mo]
	pars = model$ipars[,1]
	idx = model$pidx
	dccsum = sum(pars[idx["dcca",1]:idx["dcca",2]]) + sum(pars[idx["dccb",1]:idx["dccb",2]])
	Qt_1 = (1 - dccsum) * Qbar - sum(pars[idx["dccg",1]:idx["dccg",2]])*Nbar
	for(i in 1:n.ahead){
		Qtfor[, , mo + i] = Qt_1
		if( i == 1 ){
			if(modelinc[3]>0){
				for(j in 1:modelinc[3]){
					Qtfor[ , , mo + 1] = Qtfor[ , , mo + 1] + pars[idx["dcca",1]+j-1] * (stdresids[(mo + 1 - j), ] %*% t(stdresids[(mo + 1 - j), ]))
				}
			}
			if(modelinc[5]>0){
				for(j in 1:modelinc[5]){
					Qtfor[ , , mo + 1] = Qtfor[ , , mo + 1] + pars[idx["dccg",1]+j-1] * (astdresids[(mo + 1 - j), ] %*% t(astdresids[(mo + 1 - j), ]))
				}
			}
			if(modelinc[4]>0){
				for(j in 1:modelinc[4]){
					Qtfor[ , , mo + 1] = Qtfor[ , , mo + 1] + pars[idx["dccb",1]+j-1] * Qtfor[ , , mo + 1 - j]
				}
			}
			Qtmp = diag( 1/sqrt( diag(Qtfor[ , , mo + 1]) ) , m, m)
			Rtfor[ , , mo + 1] =  Qtmp %*% Qtfor[ , , mo + 1] %*% t(Qtmp)
			Dtmp = diag(fsig[1, ], m, m)
			Htfor[ , , mo + 1] = Dtmp %*% Rtfor[ , , mo + 1] %*% Dtmp
			# now the unconditional calculations
			Qt_1star = diag( 1/sqrt( diag(Qtfor[, , mo + 1]) ) , m, m)
			ER_1 = Qt_1star %*% Qtfor[, , mo + 1] %*% t(Qt_1star)
			Qbarstar = diag( 1/sqrt( diag(Qbar) ) , m, m)
			Rbar = Qbarstar %*% Qbar %*% t(Qbarstar)
		} else{
				Rtfor[, , mo + i] = (1 - dccsum^(i - 1) ) * Rbar + dccsum^(i - 1) * ER_1
				Dtmp = diag(fsig[i, ], m, m)
				Htfor[, , mo + i] = Dtmp %*% Rtfor[, , mo + i] %*% Dtmp
				Qtfor[, , mo + i] = Qtfor[, , mo + 1]
		}
	}
	return( list( Rtfor = Rtfor[, , -(1:mo), drop = FALSE], Htfor = Htfor[, , -(1:mo), drop = FALSE], 
					Qtfor = Qtfor[, , -(1:mo), drop = FALSE], Rbar = Rbar, ER_1 = ER_1) )
}

.dccsim.fit = function(fitORspec, n.sim = 1000, n.start = 0, m.sim = 1, 
		startMethod = c("unconditional", "sample"), presigma = NULL, 
		preresiduals = NULL, prereturns = NULL, preQ = NULL, preZ = NULL, 
		Qbar = NULL, Nbar = NULL, rseed = NULL, mexsimdata = NULL, 
		vexsimdata = NULL, cluster = NULL, VAR.fit = NULL, prerealized = NULL, ...)
{
	fit = fitORspec
	T = fit@model$modeldata$T
	Data = fit@model$modeldata$data[1:T,]
	m = dim(Data)[2]
	mo = fit@model$maxdccOrder
	mg = fit@model$maxgarchOrder
	
	startMethod = startMethod[1]
	if( is.null(rseed) ){
		rseed = as.integer(runif(1, 1, Sys.time())) 
	} else {
		if(length(rseed) == 1) rseed = as.integer(rseed[1]) else rseed = as.integer( rseed[1:m.sim] )
	}
	
	model = fit@model
	umodel = model$umodel
	modelinc = model$modelinc
	fpars = lapply(1:m, FUN = function(i) fit@model$mpars[fit@model$midx[,i]==1,i])
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			fpars, NULL)
	
	# Technically, we should pass an array of size (m, m, mo) for the
	# sample but how many use a dccOrder > 1? (and it would mean modifying
	# the C++ code to deal with arrays/lists for this input).
	if(startMethod == "sample"){
		if(is.null(preZ)){
			preZ = matrix(tail(residuals(fit)/sigma(fit), mo), ncol = m)
		} else{
			preZ = matrix(tail(preZ, 1), ncol = m, nrow = mo, byrow = TRUE)
		}
		if(is.null(preQ)){
			preQ = fit@mfit$Q[[length(fit@mfit$Q)]]
		} else{
			dcc.symcheck(preQ, m, d = NULL)
		}
		Rbar = preQ/(sqrt(diag(preQ)) %*% t(sqrt(diag(preQ))) )
	} else{
		if(is.null(preZ)){
			preZ = matrix(0, ncol = m, nrow = mo)
		} else{
			preZ = matrix(tail(preZ, 1), ncol = m, nrow = mo, byrow = TRUE)
		}
		Rbar = cor(Data)
		if(is.null(preQ)){
			preQ = Rbar
		} else{
			dcc.symcheck(preQ, m, d = NULL)
			Rbar = preQ/(sqrt(diag(preQ)) %*% t(sqrt(diag(preQ))) )
		}
	}
	if(is.null(Qbar)){
		Qbar = fit@mfit$Qbar
	} else{
		dcc.symcheck(Qbar, m, d = NULL)
	}
	if(model$modelinc[5]>0){
		if(is.null(Nbar)){
			Nbar = fit@mfit$Nbar
		} else{
			dcc.symcheck(Nbar, m, d = NULL)
		}
	} else{
		Nbar = matrix(0, m, m)
	}
	
	uncv = sapply(mspec@spec, FUN = function(x) uncvariance(x))
	
	if( !is.null(presigma) ){
		if( !is.matrix(presigma) ) 
			stop("\ndccsim-->error: presigma must be a matrix.")
		if( dim(presigma)[2] != m ) 
			stop("\ndccsim-->error: wrong column dimension for presigma.")
		if( dim(presigma)[1] != mg ) 
			stop(paste("\ndccsim-->error: wrong row dimension for presigma (need ", mg, " rows.", sep = ""))		
	} else{
		if(startMethod == "sample"){
			mx = max(sapply(mspec@spec, FUN = function(x) x@model$maxOrder))
			presigma = matrix(NA, ncol = m, nrow = mx)
			tmp = last(fit@mfit$H, mx)
			for(i in 1:mx) presigma[i,] = sqrt(diag(tmp[,,i]))
		}
	}
	
	if( !is.null(preresiduals) ){
		if( !is.matrix(preresiduals) ) 
			stop("\ndccsim-->error: preresiduals must be a matrix.")
		if( dim(preresiduals)[2] != m ) 
			stop("\ndccsim-->error: wrong column dimension for preresiduals.")
		if( dim(preresiduals)[1] != mg ) 
			stop(paste("\ndccsim-->error: wrong row dimension for preresiduals (need ", mg, " rows.", sep = ""))
	} else{
		if(startMethod == "sample"){
			mx = max(sapply(mspec@spec, FUN = function(x) x@model$maxOrder))
			preresiduals = matrix(NA, ncol = m, nrow = mx)
			tmp = tail(fit@model$residuals, mx)
			for(i in 1:mx) preresiduals[i,] = tmp[i,]
		}
	}
	
	if( !is.null(prereturns) ){
		if( !is.matrix(prereturns) ) 
			stop("\ndccsim-->error: prereturns must be a matrix.")
		if( dim(prereturns)[2] != m ) 
			stop("\ndccsim-->error: wrong column dimension for prereturns.")
		if( dim(prereturns)[1] != mg ) 
			stop(paste("\ndccsim-->error: wrong row dimension for prereturns (need ", mg, " rows.", sep = ""))
	} else{
		if(startMethod == "sample"){
			mx = max(sapply(mspec@spec, FUN = function(x) x@model$maxOrder))
			prereturns = matrix(NA, ncol = m, nrow = mx)
			tmp = tail(Data, mx)
			for(i in 1:mx) prereturns[i,] = tmp[i,]
		}
	}
	
	if(fit@model$umodel$modeldesc$vmodel[1]=="realGARCH"){
		if( !is.null(prerealized) ){
			if( !is.matrix(prerealized) ) 
				stop("\ndccsim-->error: prerealized must be a matrix.")
			if( dim(prerealized)[2] != m ) 
				stop("\ndccsim-->error: wrong column dimension for prerealized.")
			if( dim(prerealized)[1] != mg ) 
				stop(paste("\ndccsim-->error: wrong row dimension for prerealized (need ", mg, " rows.", sep = ""))
		} else{
			if(startMethod == "sample"){
				mx = max(sapply(mspec@spec, FUN = function(x) x@model$maxOrder))
				prerealized = matrix(NA, ncol = m, nrow = mx)
				tmp = tail(fit@mfit$realizedVol[1:T,], mx)
				for(i in 1:mx) prerealized[i,] = tmp[i,]
			}
		}
	} else{
		mx = max(sapply(mspec@spec, FUN = function(x) x@model$maxOrder))
		prerealized = matrix(NA, ncol = m, nrow = mx)
	}
	
	# switch distributions
	if(fit@model$modeldesc$distribution == "mvnorm"){
		if(length(rseed) == 1){
			set.seed( rseed )
			tmp = matrix(rnorm(m * (n.sim + n.start) * m.sim, 0, 1), ncol = m, nrow = n.sim+n.start)
			z = array(NA,  dim = c(n.sim + n.start + mo, m, m.sim))
			for(i in 1:m.sim) z[,,i] = rbind(preZ, tmp)
		} else{
			z = array(NA, dim = c(n.sim + n.start + mo, m, m.sim))
			for(i in 1:m.sim){
				set.seed( rseed[i] )
				z[,,i] = rbind(preZ, matrix(rnorm(m * (n.sim + n.start), 0, 1), nrow = n.sim + n.start, ncol = m))
			}
		}
	} else if(fit@model$modeldesc$distribution == "mvlaplace"){
		if(length(rseed) == 1){
			set.seed( rseed )
			tmp = matrix(rugarch:::rged(m * (n.sim + n.start) * m.sim, 0, 1, shape = 1), ncol = m, nrow = n.sim+n.start)
			z = array(NA,  dim = c(n.sim + n.start + mo, m, m.sim))
			for(i in 1:m.sim) z[,,i] = rbind(preZ, tmp)
		} else{
			z = array(NA, dim = c(n.sim + n.start + mo, m, m.sim))
			for(i in 1:m.sim){
				set.seed( rseed[i] )
				z[,,i] = rbind(preZ, matrix(rugarch:::rged(m * (n.sim + n.start), 0, 1, shape = 1), nrow = n.sim + n.start, ncol = m))
			}
		}
	} else{
		if(length(rseed) == 1){
			set.seed( rseed )
			tmp = matrix(rugarch:::rstd(m * (n.sim + n.start) * m.sim, 0, 1, shape = rshape(fit)), ncol = m, nrow = n.sim+n.start)
			z = array(NA,  dim = c(n.sim + n.start + mo, m, m.sim))
			for(i in 1:m.sim) z[,,i] = rbind(preZ, tmp)
		} else{
			z = array(NA, dim = c(n.sim + n.start + mo, m, m.sim))
			for(i in 1:m.sim){
				set.seed( rseed[i] )
				z[,,i] = rbind(preZ, matrix(rugarch:::rstd(m * (n.sim + n.start), 0, 1, shape = rshape(fit)), nrow = n.sim + n.start, ncol = m))
			}
		}
	}
	# ok now to expand rseed (but still allow for replication given initial seed)
	if(length(rseed) == 1){
		rseed = c(rseed, (1:m.sim)*(rseed+1)) 
	}
	
	
	simRes = simX = simR = simQ = simH = simSeries = vector(mode = "list", length = m.sim)
	
	if( !is.null(cluster) ){
		simH = vector(mode = "list", length = m.sim)
		simX = vector(mode = "list", length = m.sim)
		clusterEvalQ(cluster, require(rmgarch))
		clusterExport(cluster, c("model", "z", "preQ", "Rbar", 
						"Qbar", "Nbar", "mo", "n.sim", "n.start", "m", 
						"rseed",".dccsimf"), envir = environment())
		mtmp = parLapply(cluster, as.list(1:m.sim), fun = function(j){
					.dccsimf(model, Z = z[,,j], Qbar = Qbar, 
							preQ = preQ, Nbar = Nbar, Rbar = Rbar, mo = mo, 
							n.sim, n.start, m, rseed[j])
				})
		# need to pass m.sim and matrix of simulated z's to ugarchsim (speedup)
		simR = lapply(mtmp, FUN = function(x) if(is.matrix(x$R)) array(x$R, dim = c(m, m, n.sim)) else last(x$R, n.sim))
		simQ = lapply(mtmp, FUN = function(x) if(is.matrix(x$Q)) array(x$Q, dim = c(m, m, n.sim)) else last(x$Q, n.sim))
		#simZ = lapply(mtmp, FUN = function(x) x$Z)
		simZ = vector(mode = "list", length = m)
		for(i in 1:m) simZ[[i]] = sapply(mtmp, FUN = function(x) x$Z[,i])			
		clusterExport(cluster, c("fit", "n.sim", "n.start", "m.sim", 
						"startMethod", "simZ", "presigma", "preresiduals", 
						"prereturns", "mexsimdata", "vexsimdata", "prerealized"), 
				envir = environment())
		xtmp = parLapply(cluster, as.list(1:m), fun = function(j){
					maxx = mspec@spec[[j]]@model$maxOrder;
					htmp = ugarchpath(mspec@spec[[j]], n.sim = n.sim + n.start, n.start = 0, m.sim = m.sim,
							custom.dist = list(name = "sample", distfit = matrix(simZ[[j]][-(1:mo), ], ncol = m.sim)),
							presigma = if( is.null(presigma) ) NA else tail(presigma[,j], maxx), 
							preresiduals = if( is.null(preresiduals) ) NA else tail(preresiduals[,j], maxx), 
							prereturns = if( is.null(prereturns) | model$modelinc[1]>0 ) NA else tail(prereturns[,j], maxx),
							mexsimdata = if( model$modelinc[1]==0 ) mexsimdata[[j]] else NULL, 
							vexsimdata = vexsimdata[[j]], prerealized = tail(prerealized[,j], maxx));
					h = matrix(tail(htmp@path$sigmaSim^2, n.sim), nrow = n.sim);
					x = matrix(htmp@path$seriesSim,  nrow = n.sim + n.start);
					return(list(h = h, x = x))
				})
	} else{
		simH = vector(mode = "list", length = m.sim)
		simX = vector(mode = "list", length = m.sim)
		mtmp = lapply(as.list(1:m.sim), FUN = function(j){
					.dccsimf(model, Z = z[,,j], Qbar = Qbar, preQ = preQ, 
							Nbar = Nbar, Rbar = Rbar, mo = mo, n.sim, n.start, 
							m, rseed[j])
				})
		# need to pass m.sim and matrix of simulated z's to ugarchsim (speedup)
		simR = lapply(mtmp, FUN = function(x) if(is.matrix(x$R)) array(x$R, dim = c(m, m, n.sim)) else last(x$R, n.sim))
		simQ = lapply(mtmp, FUN = function(x) if(is.matrix(x$Q)) array(x$Q, dim = c(m, m, n.sim)) else last(x$Q, n.sim))
		#simZ = lapply(mtmp, FUN = function(x) x$Z)
		simZ = vector(mode = "list", length = m)
		for(i in 1:m) simZ[[i]] = sapply(mtmp, FUN = function(x) x$Z[,i])
		xtmp = lapply(as.list(1:m), FUN = function(j){
					maxx = mspec@spec[[j]]@model$maxOrder;
					htmp = ugarchpath(mspec@spec[[j]], n.sim = n.sim + n.start, n.start = 0, m.sim = m.sim,
							custom.dist = list(name = "sample", distfit = matrix(simZ[[j]][-(1:mo), ], ncol = m.sim)),
							presigma = if( is.null(presigma) ) NA else tail(presigma[,j], maxx), 
							preresiduals = if( is.null(preresiduals) ) NA else tail(preresiduals[,j], maxx), 
							prereturns = if( is.null(prereturns) | model$modelinc[1]>0 ) NA else tail(prereturns[,j], maxx),
							mexsimdata = if( model$modelinc[1]==0 ) mexsimdata[[j]] else NULL, 
							vexsimdata = vexsimdata[[j]], prerealized = tail(prerealized[,j], maxx));
					h = matrix(tail(htmp@path$sigmaSim^2, n.sim), nrow = n.sim);
					x = matrix(htmp@path$seriesSim,  nrow = n.sim + n.start);
					return(list(h = h, x = x))})
	}
	H = array(NA, dim = c(n.sim, m, m.sim))
	tmpH = array(NA, dim = c(m, m, n.sim))
	for(i in 1:n.sim) H[i,,] = t(sapply(xtmp, FUN = function(x) as.numeric(x$h[i,])))
	
	for(i in 1:m.sim){
		for(j in 1:n.sim){
			tmpH[ , , j] = diag(sqrt( H[j, , i]) ) %*% simR[[i]][ , , j] %*% diag(sqrt( H[j, , i] ) )
		}
		simH[[i]] = tmpH
	}
	if(model$modelinc[1]>0){
		simxX = array(NA, dim = c(n.sim+n.start, m, m.sim))
		for(i in 1:m.sim) simxX[,,i] = sapply(xtmp, FUN = function(x) as.numeric(x$x[,i]))
		simX = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim) simX[[i]] = matrix(simxX[,,i], nrow = n.sim+n.start)
	} else{
		simxX = array(NA, dim = c(n.sim+n.start, m, m.sim))
		for(i in 1:m.sim) simxX[,,i] = sapply(xtmp, FUN = function(x) as.numeric(x$x[,i]))
		simX = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim) simX[[i]] = matrix(tail(matrix(simxX[,,i], ncol = m), n.sim), nrow = n.sim)
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
	msim$simH = simH
	msim$simR = simR
	msim$simQ = simQ
	msim$simX = simX
	msim$simZ = simZ
	msim$simRes = simRes
	#msim$Z = stdresid
	msim$rseed = rseed
	model$n.sim = n.sim
	model$m.sim = m.sim
	model$n.start = n.start
	model$startMethod = startMethod[1]
	ans = new("DCCsim",
			msim = msim,
			model = model)
	
	return( ans )
}

.dccsim.spec = function(fitORspec, n.sim = 1000, n.start = 0, m.sim = 1, 
		startMethod = c("unconditional", "sample"), presigma = NULL, 
		preresiduals = NULL, prereturns = NULL, preQ = NULL, preZ = NULL, 
		Qbar = NULL, Nbar = NULL, rseed = NULL, mexsimdata = NULL, 
		vexsimdata = NULL, cluster = NULL, VAR.fit = NULL, 
		prerealized = NULL, ...)
{
	spec = fitORspec
	startMethod = startMethod[1]
	# check that VAR.fit is included if used
	if( spec@model$modelinc[1]>0 ){
		if( is.null(VAR.fit) ) stop("\ndccsim-->error: VAR.fit must not be NULL for VAR method when calling dccsim using spec!", call. = FALSE)
	}
	model = spec@model
	model$umodel = spec@umodel
	m = dim(spec@umodel$modelinc)[2]
	mo = spec@model$maxdccOrder
	mg = spec@model$maxgarchOrder
	model$modeldata$asset.names = paste("Asset", 1:m, sep = "")
	
	if( is.null(rseed) ){
		rseed = as.integer(runif(1, 1, Sys.time())) 
	} else {
		if(length(rseed) == 1) rseed = as.integer(rseed[1]) else rseed = as.integer( rseed[1:m.sim] )
	}
	if(is.null(preZ)){
		preZ = matrix(0, ncol = m, nrow = mo)
	} else{
		preZ = matrix(tail(preZ, 1), ncol = m, nrow = mo, byrow = TRUE)
	}
	if(is.null(preQ)){
		stop("\ndccsim-->error: preQ cannot be NULL when method uses spec!")		
	} else{
		dcc.symcheck(preQ, m, d = NULL)
		Rbar = preQ/(sqrt(diag(preQ)) %*% t(sqrt(diag(preQ))) )
	}
	if(is.null(Qbar)){
		stop("\ndccsim-->error: Qbar cannot be NULL when method uses spec!")		
	} else{
		dcc.symcheck(Qbar, m, d = NULL)
	}
	if(model$modelinc[5]>0){
		if(is.null(Nbar)){
			stop("\ndccsim-->error: Nbar cannot be NULL for aDCC when method uses spec!")
		} else{
			dcc.symcheck(Nbar, m, d = NULL)
		}
	} else{
		Nbar = matrix(0, m, m)
	}
	
	model = spec@model
	umodel = spec@umodel
	modelinc = model$modelinc
	
	midx = .fullinc(modelinc, umodel)
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, NULL)
	
	# uncvariance has uGARCHfit AND uGARCHspec dispatch methods
	uncv = sapply(mspec@spec, FUN = function(x) uncvariance(x))
	
	if( !is.null(presigma) ){
		if( !is.matrix(presigma) ) 
			stop("\ndccsim-->error: presigma must be a matrix.")
		if( dim(presigma)[2] != m ) 
			stop("\ndccsim-->error: wrong column dimension for presigma.")
		if( dim(presigma)[1] != mg ) 
			stop(paste("\ndccsim-->error: wrong row dimension for presigma (need ", mg, " rows.", sep = ""))		
	}
	
	if( !is.null(preresiduals) ){
		if( !is.matrix(preresiduals) ) 
			stop("\ndccsim-->error: preresiduals must be a matrix.")
		if( dim(preresiduals)[2] != m ) 
			stop("\ndccsim-->error: wrong column dimension for preresiduals.")
		if( dim(preresiduals)[1] != mg ) 
			stop(paste("\ndccsim-->error: wrong row dimension for preresiduals (need ", mg, " rows.", sep = ""))
	}
	
	if( !is.null(prereturns) ){
		if( !is.matrix(prereturns) ) 
			stop("\ndccsim-->error: prereturns must be a matrix.")
		if( dim(prereturns)[2] != m ) 
			stop("\ndccsim-->error: wrong column dimension for prereturns.")
		if( dim(prereturns)[1] != mg ) 
			stop(paste("\ndccsim-->error: wrong row dimension for prereturns (need ", mg, " rows.", sep = ""))
	}
	
	if(spec@umodel$modeldesc$vmodel[1]=="realGARCH"){
		if( !is.null(prerealized) ){
			if( !is.matrix(prerealized) ) 
				stop("\ndccsim-->error: prerealized must be a matrix.")
			if( dim(prerealized)[2] != m ) 
				stop("\ndccsim-->error: wrong column dimension for prerealized.")
			if( dim(prerealized)[1] != mg ) 
				stop(paste("\ndccsim-->error: wrong row dimension for prerealized (need ", mg, " rows.", sep = ""))
		}
	} else{
		prerealized = matrix(NA, ncol = m, nrow = mg)
	}
	
	
	# we use the GED distribution  with nu = 1 which corresponds to the Laplace case.
	if(model$modeldesc$distribution == "mvnorm"){
		if(length(rseed) == 1){
			set.seed( rseed )
			tmp = matrix(rnorm(m * (n.sim + n.start) * m.sim, 0, 1), ncol = m, nrow = n.sim+n.start)
			z = array(NA,  dim = c(n.sim + n.start + mo, m, m.sim))
			for(i in 1:m.sim) z[,,i] = rbind(preZ, tmp)
		} else{
			z = array(NA, dim = c(n.sim + n.start + mo, m, m.sim))
			for(i in 1:m.sim){
				set.seed( rseed[i] )
				z[,,i] = rbind(preZ, matrix(rnorm(m * (n.sim + n.start), 0, 1), nrow = n.sim + n.start, ncol = m))
			}
		}
	} else if(model$modeldesc$distribution == "mvlaplace"){
		if(length(rseed) == 1){
			set.seed( rseed )
			tmp = matrix(rugarch:::rged(m * (n.sim + n.start) * m.sim, 0, 1, shape = 1), , ncol = m, nrow = n.sim+n.start)
			z = array(NA,  dim = c(n.sim + n.start + mo, m, m.sim))
			for(i in 1:m.sim) z[,,i] = rbind(preZ, tmp)
		} else{
			z = array(NA, dim = c(n.sim + n.start + mo, m, m.sim))
			for(i in 1:m.sim){
				set.seed( rseed[i] )
				z[,,i] = rbind(preZ, matrix(rugarch:::rged(m * (n.sim + n.start), 0, 1, shape = 1), nrow = n.sim + n.start, ncol = m))
			}
		}
	} else{
		if(length(rseed) == 1){
			set.seed( rseed )
			tmp = matrix(rugarch:::rstd(m * (n.sim + n.start) * m.sim, 0, 1, shape = model$pars["mshape",1]), ncol = m, nrow = n.sim+n.start)
			z = array(NA,  dim = c(n.sim + n.start + mo, m, m.sim))
			for(i in 1:m.sim) z[,,i] = rbind(preZ, tmp)
		} else{
			z = array(NA, dim = c(n.sim + n.start + mo, m, m.sim))
			for(i in 1:m.sim){
				set.seed( rseed[i] )
				z[,,i] = rbind(preZ, matrix(rugarch:::rstd(m * (n.sim + n.start), 0, 1, shape = model$pars["mshape",1]), nrow = n.sim + n.start, ncol = m))
			}
		}
	}
	
	# ok now to expand rseed
	if(length(rseed) == 1){
		rseed = c(rseed, as.integer(runif(m.sim, 1, Sys.time()))) 
	}

	simRes = simX = simR = simQ = simH = simSeries = vector(mode = "list", length = m.sim)
	
	if( !is.null(cluster) ){
		simH = vector(mode = "list", length = m.sim)
		simX = vector(mode = "list", length = m.sim)
		clusterEvalQ(cluster, require(rmgarch))
		clusterExport(cluster, c("model", "z", "preQ", "Rbar", 
						"Qbar", "Nbar", "mo", "n.sim", "n.start", "m", 
						"rseed",".dccsimf"), envir = environment())
		mtmp = parLapply(cluster, as.list(1:m.sim), fun = function(j){
					.dccsimf(model, Z = z[,,j], Qbar = Qbar, 
							preQ = preQ, Nbar = Nbar, Rbar = Rbar, mo = mo, 
							n.sim, n.start, m, rseed[j])
				})
		# need to pass m.sim and matrix of simulated z's to ugarchsim (speedup)
		simR = lapply(mtmp, FUN = function(x) if(is.matrix(x$R)) array(x$R, dim = c(m, m, n.sim)) else last(x$R, n.sim))
		simQ = lapply(mtmp, FUN = function(x) if(is.matrix(x$Q)) array(x$Q, dim = c(m, m, n.sim)) else last(x$Q, n.sim))
		#simZ = lapply(mtmp, FUN = function(x) x$Z)
		
		simZ = vector(mode = "list", length = m)
		for(i in 1:m) simZ[[i]] = sapply(mtmp, FUN = function(x) x$Z[,i])
		clusterExport(cluster, c("mspec", "n.sim", "n.start", "m.sim", 
						"startMethod", "simZ", "presigma", "preresiduals", 
						"prereturns", "mexsimdata", "vexsimdata", "prerealized"), 
				envir = environment())
		xtmp = parLapply(cluster, as.list(1:m), fun = function(j){
					maxx = mspec@spec[[j]]@model$maxOrder;
					htmp = ugarchpath(mspec@spec[[j]], n.sim = n.sim + n.start, n.start = 0, m.sim = m.sim,
							custom.dist = list(name = "sample", distfit = matrix(simZ[[j]][-(1:mo), ], ncol = m.sim)),
							presigma = if( is.null(presigma) ) NA else tail(presigma[,j], maxx), 
							preresiduals = if( is.null(preresiduals) ) NA else tail(preresiduals[,j], maxx), 
							prereturns = if( is.null(prereturns) | model$modelinc[1]>0 ) NA else tail(prereturns[,j], maxx),
							mexsimdata = if( model$modelinc[1]==0 ) mexsimdata[[j]] else NULL, 
							vexsimdata = vexsimdata[[j]], prerealized = tail(prerealized[,j], maxx))
					h = matrix(tail(htmp@path$sigmaSim^2, n.sim), nrow = n.sim)
					x = matrix(htmp@path$seriesSim,  nrow = n.sim + n.start)
					xres = matrix(htmp@path$residSim,  nrow = n.sim + n.start)
					return(list(h = h, x = x, xres = xres))
				})
	} else{
		simH = vector(mode = "list", length = m.sim)
		simX = vector(mode = "list", length = m.sim)
		mtmp = lapply(as.list(1:m.sim), FUN = function(j) 
					.dccsimf(model, Z = z[,,j], Qbar = Qbar, preQ = preQ, Nbar = Nbar, 
							Rbar = Rbar, mo = mo, n.sim, n.start, m, rseed[j]))		
		# need to pass m.sim and matrix of simulated z's to ugarchsim (speedup)
		simR = lapply(mtmp, FUN = function(x) if(is.matrix(x$R)) array(x$R, dim = c(m, m, n.sim)) else last(x$R, n.sim))
		simQ = lapply(mtmp, FUN = function(x) if(is.matrix(x$Q)) array(x$Q, dim = c(m, m, n.sim)) else last(x$Q, n.sim))
		#simZ = lapply(mtmp, FUN = function(x) x$Z)
		simZ = vector(mode = "list", length = m)
		for(i in 1:m) simZ[[i]] = sapply(mtmp, FUN = function(x) x$Z[,i])			
		xtmp = lapply(as.list(1:m), FUN = function(j){
					maxx = mspec@spec[[j]]@model$maxOrder;
					htmp = ugarchpath(mspec@spec[[j]], n.sim = n.sim + n.start, n.start = 0, m.sim = m.sim,
							custom.dist = list(name = "sample", distfit = matrix(simZ[[j]][-(1:mo), ], ncol = m.sim)),
							presigma = if( is.null(presigma) ) NA else tail(presigma[,j], maxx), 
							preresiduals = if( is.null(preresiduals) ) NA else tail(preresiduals[,j], maxx), 
							prereturns = if( is.null(prereturns) | model$modelinc[1]>0 ) NA else tail(prereturns[,j], maxx),
							mexsimdata = if( model$modelinc[1]==0 ) mexsimdata[[j]] else NULL, vexsimdata = vexsimdata[[j]], 
							prerealized = tail(prerealized[,j], maxx));
					h = matrix(tail(htmp@path$sigmaSim^2, n.sim), nrow = n.sim);
					x = matrix(htmp@path$seriesSim,  nrow = n.sim + n.start);
					xres = matrix(htmp@path$residSim,  nrow = n.sim + n.start);
					return(list(h = h, x = x, xres = xres))
				})
	}
	H = array(NA, dim = c(n.sim, m, m.sim))
	tmpH = array(NA, dim = c(m, m, n.sim))
	for(i in 1:n.sim) H[i,,] = t(sapply(xtmp, FUN = function(x) as.numeric(x$h[i,])))	
	for(i in 1:m.sim){
		for(j in 1:n.sim){
			tmpH[ , , j] = diag(sqrt( H[j, , i]) ) %*% simR[[i]][ , , j] %*% diag(sqrt( H[j, , i] ) )
		}
		simH[[i]] = tmpH
	}
	if(model$modelinc[1]>0){
		simxX = array(NA, dim = c(n.sim+n.start, m, m.sim))
		for(i in 1:m.sim) simxX[,,i] = sapply(xtmp, FUN = function(x) as.numeric(x$x[,i]))
		simX = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim) simX[[i]] = matrix(simxX[,,i], nrow = n.sim+n.start)
	} else{
		simxX = array(NA, dim = c(n.sim+n.start, m, m.sim))
		for(i in 1:m.sim) simxX[,,i] = sapply(xtmp, FUN = function(x) as.numeric(x$x[,i]))
		simX = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim) simX[[i]] = matrix(tail(matrix(simxX[,,i], ncol = m), n.sim), nrow = n.sim)
	}
	simxRes = array(NA, dim = c(n.sim, m, m.sim))
	for(i in 1:n.sim) simxRes[i,,] = t(sapply(xtmp, FUN = function(x) as.numeric(x$xres[i,])))
	simRes = vector(mode = "list", length = m.sim)
	for(i in 1:m.sim) simRes[[i]] = matrix(simxRes[,,i], nrow = n.sim)
	
	if( model$modelinc[1]>0 ){
		model$varcoef = VAR.fit$Bcoef
		Data = VAR.fit$xfitted
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
	msim$simH = simH
	msim$simR = simR
	msim$simQ = simQ
	msim$simX = simX
	msim$simZ = simZ
	msim$simRes = simRes
	#msim$Z = stdresid
	msim$rseed = rseed
	msim$model$Data = NULL
	model$n.sim = n.sim
	model$m.sim = m.sim
	model$n.start = n.start
	model$startMethod = "unconditional"
	ans = new("DCCsim",
			msim = msim,
			model = model)	
	return( ans )
	
}

.rolldcc = function(spec, data, n.ahead = 1, forecast.length = 50, refit.every = 25, 
		n.start = NULL, refit.window = c("recursive", "moving"), window.size = NULL, 
		solver = "solnp", solver.control = list(), 
		fit.control = list(eval.se = TRUE, stationarity = TRUE, scale = FALSE), 
		cluster = NULL, save.fit = FALSE, save.wdir = NULL, realizedVol = NULL, ...)
{
	if(spec@model$DCC=="FDCC") stop("\nFDCC model rolling estimation not yet implemented.")
	model = spec@model
	verbose = FALSE
	model$umodel = spec@umodel
	if(is.null(solver.control$trace)) solver.control$trace = 0
	if(is.null(fit.control$stationarity)) fit.control$stationarity = TRUE
	if(is.null(fit.control$eval.se)) fit.control$eval.se = FALSE
	if(is.null(fit.control$scale)) fit.control$scale = FALSE
	
	mm = match(names(fit.control), c("stationarity", "eval.se", "scale"))
	if(any(is.na(mm))){
		idx = which(is.na(mm))
		enx = NULL
		for(i in 1:length(idx)) enx = c(enx, names(fit.control)[idx[i]])
		warning(paste(c("unidentified option(s) in fit.control:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
	}
	asset.names = colnames(data)
	xdata = .extractmdata(data)
	data = xts(xdata$data, xdata$index)
	index = xdata$index
	period = xdata$period
	if(is.null(fit.control$stationarity)) fit.control$stationarity = 1
	if(is.null(fit.control$fixed.se)) fit.control$fixed.se = 0
	T = dim(data)[1]
	
	if(n.ahead>1) 
		stop("\ndccroll:--> n.ahead>1 not supported...try again.")
	if(is.null(n.start)){
		if(is.null(forecast.length)) 
			stop("\ndccroll:--> forecast.length amd n.start are both NULL....try again.")
		n.start = T - forecast.length
	} else{
		forecast.length = T - n.start
	}
	if(T<=n.start) 
		stop("\ndccroll:--> start cannot be greater than length of data")
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
			if(window.size<100) stop("\ndccroll:--> window size must be greater than 100.")
			rollind = lapply(1:m, FUN = function(i) max(1, (s[i]-window.size-out.sample[i])):s[i])
		} else{
			rollind = lapply(1:m, FUN = function(i) (1+(i-1)*refit.every):s[i])
		}
	}
	WD <- getwd()
	if(is.null(save.wdir)){
		if (!is.null(WD)) setwd(WD)
	} else{
		ND = save.wdir
		if (!is.null(ND)) setwd(ND)
	}
	cf = lik = forc = vector(mode = "list", length = m)
	plik = vector(mode = "list", length = m)
	
	mspec = .makemultispec(model$umodel$modelinc, model$umodel$modeldesc$vmodel, model$umodel$modeldesc$vsubmodel, 
			model$umodel$modeldata$mexdata, model$umodel$modeldata$vexdata, model$umodel$start.pars, 
			model$umodel$fixed.pars, model$umodel$vt)
	# changed estimation strategy to catch non-converged or problematic solutions
	for(i in 1:m){
		if(!is.null(realizedVol)){
			mfit = multifit(mspec, data[rollind[[i]],], out.sample = out.sample[i], 
					solver = solver[1], fit.control = fit.control, cluster = cluster, 
					realizedVol = realizedVol[rollind[[i]],], solver.control = solver.control)
			k=1
			while(k==1){
				conv = sapply(mfit@fit, function(x) convergence(x))
				if(any(conv==1)){
					idx = which(conv==1)
					for(j in idx){ mfit@fit[[j]] = ugarchfit(mspec@spec[[j]], data[rollind[[i]],j], 
								out.sample = out.sample[i], solver = "gosolnp", fit.control = fit.control,
								realizedVol = realizedVol[rollind[[i]],j])
					}
				} else{
					k=0
				}
			}
			# for some reason realGARCH might converge to a problematic value
			k=1
			while(k==1){
				tmp = sapply(mfit@fit, function(x){
							L = try(likelihood(x), silent=TRUE)
							if(inherits(L, 'try-error') | !is.numeric(L)) L = 1e10
							L})
				conv=diff(log(abs(tmp)))
				if(any(conv>1)){
					idx = which(conv>1)+1
					for(j in idx){ 
						mfit@fit[[j]] = ugarchfit(mspec@spec[[j]], data[rollind[[i]],j], 
								out.sample = out.sample[i], 
								solver = "gosolnp", fit.control = fit.control,
								realizedVol = realizedVol[rollind[[i]],j])
					}
				} else{
					k=0
				}
			}
			mcfit = dccfit(spec, data[rollind[[i]],], out.sample = out.sample[i], 
					solver = solver, fit.control = fit.control, solver.control=solver.control,
					cluster = NULL, realizedVol = realizedVol[rollind[[i]],], 
					fit = mfit)
			plik[[i]] = mcfit@mfit$plik
		} else{
			mfit = multifit(mspec, data[rollind[[i]],], out.sample = out.sample[i], 
					solver = solver[1], fit.control = fit.control, solver.control = solver.control,
					cluster = cluster)
			k=1
			while(k==1){
				conv = sapply(mfit@fit, function(x) convergence(x))
				if(any(conv==1)){
					idx = which(conv==1)
					for(j in idx){ mfit@fit[[j]] = ugarchfit(mspec@spec[[j]], data[rollind[[i]],j], 
								out.sample = out.sample[i], solver = "gosolnp", fit.control = fit.control)
					}
				} else{
					k=0
				}
			}
			mcfit = dccfit(spec, data[rollind[[i]],], out.sample = out.sample[i], 
					solver = solver, fit.control = fit.control, solver.control = solver.control, 
					cluster = cluster, fit = mfit)
			plik[[i]] = mcfit@mfit$plik
		}
		cf[[i]]   = mcfit@model$mpars
		lik[[i]]  = likelihood(mcfit)
		forc[[i]] = dccforecast(mcfit, n.ahead = 1, n.roll = out.sample[i]-1, cluster = cluster)
		if(save.fit){
			eval(parse(text = paste("dccroll_",i,"=mcfit",sep = "")))
			eval(parse(text = paste("save(dccroll_",i,",file='dccroll_",i,".rda')",sep = "")))
		}
	}
	model$n.start = n.start
	model$n.refits = m
	model$refit.every = refit.every
	model$refit.window = refit.window
	model$window.size = window.size
	model$forecast.length = forecast.length
	model$n.start = n.start
	model$rollind = rollind
	model$out.sample = out.sample
	model$modeldata$asset.names = asset.names
	model$rollcoef = cf
	model$rolllik = lik
	model$index = index
	model$period = period
	model$data = xdata$data
	model$plik = plik
	ans = new("DCCroll",
			mforecast = forc,
			model = model)
	# return to the current directory
	setwd(WD)
	return(ans)
}

.dccsimf = function(model, Z, Qbar, preQ, Nbar, Rbar, mo, n.sim, n.start, m, rseed)
{
	modelinc = model$modelinc
	ipars = model$pars
	idx = model$pidx
	n = n.sim + n.start + mo
	set.seed(rseed[1]+1)
	# All distributions here make use of this since they are mixtures of normals
	NZ = matrix(rnorm(m * (n.sim+n.start+mo)), nrow = n.sim + n.start + mo, ncol = m)
	sumdcca = sum(ipars[idx["dcca",1]:idx["dcca",2],1])
	sumdccb = sum(ipars[idx["dccb",1]:idx["dccb",2],1])
	sumdcc = sumdcca + sumdccb
	sumdccg = sum(ipars[idx["dccg",1]:idx["dccg",2],1])
	res = switch(model$modeldesc$distribution,
			mvnorm = .Call( "dccsimmvn", model = as.integer(modelinc), pars = as.numeric(ipars[,1]), 
					idx = as.integer(idx[,1]-1), Qbar = as.matrix(Qbar),  preQ = as.matrix(preQ), 
					Rbar = as.matrix(Rbar), Nbar = as.matrix(Nbar), Z = as.matrix(Z), 
					NZ = as.matrix(NZ), epars = c(sumdcc, sumdccg, mo), PACKAGE = "rmgarch"),
			mvlaplace = .Call( "dccsimmvl", model = as.integer(modelinc), pars = as.numeric(ipars[,1]), 
					idx = as.integer(idx[,1]-1), Qbar = as.matrix(Qbar),  preQ = as.matrix(preQ), 
					Rbar = as.matrix(Rbar), Nbar = as.matrix(Nbar), Z = as.matrix(Z), 
					NZ = as.matrix(NZ), epars = c(sumdcc, sumdccg, mo), PACKAGE = "rmgarch"),
			mvt = .Call( "dccsimmvt", model = as.integer(modelinc), pars = as.numeric(ipars[,1]), 
					idx = as.integer(idx[,1]-1), Qbar = as.matrix(Qbar),  preQ = as.matrix(preQ), 
					Rbar = as.matrix(Rbar), Nbar = as.matrix(Nbar), Z = as.matrix(Z), 
					NZ = as.matrix(NZ), epars = c(sumdcc, sumdccg, mo), PACKAGE = "rmgarch"))
	Q = array(NA, dim = c(m, m, n.sim + n.start + mo))
	R = array(NA, dim = c(m, m, n.sim + n.start + mo))
	for(i in 1:(n.sim + n.start + mo)){
		R[,,i] = res[[2]][[i]]
		Q[,,i] = res[[1]][[i]]
	}
	ans = list( Q = Q, R = R, Z = res[[3]])
	return( ans )
}

.asymI = function(x){
	ans = (-sign(x)+1)/2
	# deal with cases which are zero
	ans[ans==0.5] = 0
	ans
}

.makemultispec = function(modelinc, vmodel, vsubmodel, mexdata, vexdata, spars, fpars, vt){
	m = dim(modelinc)[2]
	mspec = vector(mode = "list", length = m)
	dist = c("norm", "snorm", "std", "sstd","ged", "sged", "nig", "ghyp", "jsu", "ghst")
	for(i in 1:m){
		if(is.null(vt)){
			vtarget = FALSE
		} else{
			if(!is.na(vt[i])) vtarget = vt[i] else vtarget = ifelse(modelinc[7,i]==0, TRUE, FALSE)
		}
		mspec[[i]] = ugarchspec(variance.model = list(model = vmodel[i], garchOrder = modelinc[8:9,i], 
						submodel = vsubmodel[i], external.regressors = if(is.na(vexdata[[i]][1])) NULL else vexdata[[i]],
						variance.targeting = vtarget), 
				mean.model = list(armaOrder = modelinc[2:3,i], include.mean = as.logical(modelinc[1,i]), 
						archm = ifelse(modelinc[5,i]>0, TRUE, FALSE), archpow = modelinc[5,i], 
				arfima = modelinc[4,i], external.regressors = if(is.na(mexdata[[i]][1])) NULL else mexdata[[i]]), 
				distribution.model = dist[modelinc[21,i]], start.pars = if(is.na(spars[[i]][1])) NULL else spars[[i]], 
				fixed.pars = if(is.na(fpars[[i]][1])) NULL else fpars[[i]])
	}
	ans = multispec( mspec )
	return(ans)
}


# The challenge is in accomodating different lags for the garch and arma of 
# differing specs.
# We constrain the univariate distribution to a single family.

.fullinc = function(modelinc, umodel){
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
	
	sumdcc = 5 + sum(pmax(1, modelinc[c(3,4,5,6,7)])) - 5
	tmpmat = rbind(tmpmat, matrix(0, ncol = m+1, nrow = sumdcc))
	if(modelinc[3]>0){
		for(i in 1:modelinc[3]){
			tmpmat[nx+i, m+1] = 1
			pnames = c(pnames, paste("dcca", i, sep = ""))
		}
	} else{
		pnames = c(pnames, "dcca")
	}
	nx = nx + max(1, modelinc[3])
	
	if(modelinc[4]>0){
		for(i in 1:modelinc[4]){
			tmpmat[nx+i, m+1] = 1
			pnames = c(pnames, paste("dccb", i, sep = ""))
		}
	} else{
		pnames = c(pnames, "dccb")
	}
	
	nx = nx + max(1, modelinc[4])
	
	if(modelinc[5]>0){
		for(i in 1:modelinc[5]){
			tmpmat[nx+i, m+1] = 1
			pnames = c(pnames, paste("dccg", i, sep = ""))
		}
	} else{
		pnames = c(pnames, "dccg")
	}
	nx = nx + max(1, modelinc[5])
	
	# we allow for possibility of vector valued shape and skew parameters for future expansion
	if(modelinc[6]>0){
		for(i in 1:modelinc[6]){
			tmpmat[nx+i, m+1] = 1
			if(modelinc[6]>1) pnames = c(pnames, paste("mshape", i, sep = "")) else pnames = c(pnames, "mshape")
		}
	} else{
		pnames = c(pnames, "mshape")
	}
	nx = nx + max(1, modelinc[6])
	
	if(modelinc[7]>0){
		for(i in 1:modelinc[7]){
			tmpmat[nx+i, m+1] = 1
			if(modelinc[7]>1) pnames = c(pnames, paste("mskew", i, sep = "")) else pnames = c(pnames, "mskew")
		}
	} else{
		pnames = c(pnames, "mskew")
	}
	# NOW we have and indicator matrix
	colnames(tmpmat) = c(paste("Asset", 1:m, sep = ""), "Joint")
	rownames(tmpmat) = pnames
	return(tmpmat)
}

.estindfn = function(midx, mspec, dccpars){
	m = dim(midx)[2]-1
	eidx = midx*0
	rnx = rownames(midx)
	for(i in 1:m){
		um = mspec@spec[[i]]@model$pars
		zi = match(rownames(um), rnx)
		zi = zi[!is.na(zi)]
		zx = match(rnx, rownames(um))
		zx = zx[!is.na(zx)]
		eidx[zi,i] = um[zx,4]
	}
	zi = match(rownames(dccpars), rnx)
	eidx[zi,m+1] = dccpars[,4]
	return(eidx)
}