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
.fdccspec = function(uspec, VAR = FALSE, robust = FALSE, lag = 1, lag.max = NULL, 
		lag.criterion = c("AIC", "HQ", "SC", "FPE"), external.regressors = NULL, 
		robust.control = list("gamma" = 0.25, "delta" = 0.01, "nc" = 10, "ns" = 500), 
		fdccOrder=c(1,1), fdccindex = c(1,1,1,2,2,2,3,3,3), 
		distribution = c("mvnorm", "mvt", "mvlaplace"), 
		start.pars = list(), fixed.pars = list())
{
	# VAR checks
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
	
	## distribution checks
	if(is.null(distribution)) distribution = "mvnorm"
	distribution = distribution[1]
	valid.distributions = c("mvnorm", "mvt", "mvlaplace")
	if(!any(distribution == valid.distributions)) stop("\nInvalid Distribution Choice\n", call. = FALSE)
	
	## fdcc checks
	fdccOrder = as.integer(abs(fdccOrder))
	if(max(fdccOrder)>1) stop("\nFDCC model only allows (1,1) DCC dynamics and not higher.")
	fdccindex = as.integer(abs(fdccindex))
	if(any(fdccindex)==0) stop("\nfdccindex cannot contain zeros...index starts at 1!")
	lenF = length(unique(fdccindex))
	if(length(fdccindex)!=length(uspec@spec)) stop("\nfdccindex length not equal to length of uspec.")
	# allow for case of just 1 index
	if(!all(fdccindex==1) && any(diff(unique(fdccindex)))!=1) stop("\nfdccindex must use contiguous indexing, starting from one and incrementing by 1.")
	
	modelinc = rep(0, 10)
	names(modelinc) = c("var", "mvmxreg", "fdccC", "fdccA", "fdccB", "mshape", "mskew", "aux", "aux", "aux")
	if(distribution == "mvt") modelinc[6] = 1
	modelinc[3] = lenF
	modelinc[4] = fdccOrder[1]
	modelinc[5] = fdccOrder[2]
	# indicator matrix of groups
	sidx = matrix(0, ncol = lenF, nrow = length(fdccindex))
	for(i in 1:lenF) sidx[fdccindex==i,i] = 1
	# sidx %*% a
	# ai_j...i=sector, j = order
	# bi_j...
	# ci
	
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
	modeldesc$distribution = distribution
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
	pos = 1
	pos.matrix = matrix(0, ncol = 3, nrow = 4)
	colnames(pos.matrix) = c("start", "stop", "include")
	rownames(pos.matrix) = c("fdcca", "fdccb", "mshape", "mskew")
	
	# c("var", "mvmxreg", "fdccC", "fdccA", "fdccB", "mshape", "mskew", "aux", "aux", "aux")
	if(modelinc[4]>0) pos.matrix[1,1:3] = c(pos, pos+(modelinc[4]*modelinc[3])-1, 1)
	pos = pos + (modelinc[4]*modelinc[3])
	if(modelinc[5]>0) pos.matrix[2,1:3] = c(pos, pos+(modelinc[5]*modelinc[3])-1, 1)
	pos = pos + (modelinc[5]*modelinc[3])
	if(modelinc[6]>0) pos.matrix[3,1:3] = c(pos, pos, 1)
	pos = pos + modelinc[6]
	if(modelinc[7]>0) pos.matrix[4,1:3] = c(pos, pos, 1)
	
	# minimum parameter matrix = modelinc[3]+dcc1+dccb+mshape+mskew = modelinc[3]+4
	mm = max(1, (modelinc[4]*modelinc[3]))+max(1, (modelinc[5]*modelinc[3])) + max(1, modelinc[6]) + max(1, modelinc[6])
	pars = matrix(0, ncol = 6, nrow = mm)
	colnames(pars) = c("Level", "Fixed", "Include", "Estimate", "LB", "UB")
	pidx = matrix(NA, nrow = 4, ncol = 2)
	colnames(pidx) = c("begin", "end")
	rownames(pidx) =  c("fdcca", "fdccb", "mshape", "mskew")
	fixed.names = names(fixed.pars)
	start.names = names(start.pars)
	fixed.pars = unlist(fixed.pars)
	start.pars = unlist(start.pars)
	
	pn = 1
	pnames = NULL
	nx = 0
	pidx[1,1] = 1
	if(pos.matrix[1,3] == 1){
		pn = length( seq(pos.matrix[1,1], pos.matrix[1,2], by = 1) )
		xid = 0
		for(i in 1:modelinc[4]){
			for(j in 1:modelinc[3]){
				nnx = paste("fdcca[g=", j,",P=", i, "]", sep="")
				pars[(nx+xid+j), 1] = 0.01
				sp = na.omit(match(nnx, start.names))
				if(length(sp)==1) pars[(nx+xid+j), 1] = start.pars[sp]			
				pars[(nx+xid+j), 3] = 1
				pars[(nx+xid+j), 5] = 0
				pars[(nx+xid+j), 6] =  0.999
				sp = na.omit(match(nnx, fixed.names))
				if(length(sp)==1){
					pars[(nx+xid+j), 1] = fixed.pars[sp]
					pars[(nx+xid+j), 2] = 1
				} else{
					pars[(nx+xid+j), 4] = 1
				}
				pnames = c(pnames, nnx)
			}
			xid = xid + modelinc[3]
		}
	} else{
		pnames = c(pnames, "fdcca")
	}
	pidx[1,2] = nx + pn
	nx = nx + pn
	pn = 1
	pidx[2,1] = nx+1
	if(pos.matrix[2,3] == 1){
		pn = length( seq(pos.matrix[2,1], pos.matrix[2,2], by = 1) )
		xid = 0
		for(i in 1:modelinc[5]){
			for(j in 1:modelinc[3]){
				nnx = paste("fdccb[g=", j,",Q=", i, "]", sep="")
				pars[(nx+xid+j), 1] = 0.9
				sp = na.omit(match(nnx, start.names))
				if(length(sp)==1) pars[(nx+xid+j), 1] = start.pars[sp]			
				pars[(nx+xid+j), 3] = 1
				pars[(nx+xid+j), 5] = 0
				pars[(nx+xid+j), 6] = 1
				sp = na.omit(match(nnx, fixed.names))
				if(length(sp)==1){
					pars[(nx+xid+j), 1] = fixed.pars[sp]
					pars[(nx+xid+j), 2] = 1
				} else{
					pars[(nx+xid+j), 4] = 1
				}
				pnames = c(pnames, nnx)
			}
			xid = xid+modelinc[3]
		}
	} else{
		pnames = c(pnames, "fdccb")
	}
	pidx[2,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[3,1] = nx+1
	if(modelinc[6]<=1){
		if(pos.matrix[3,3]==1){
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
		if(pos.matrix[3,3] == 1){
			pn = length( seq(pos.matrix[3,1], pos.matrix[3,2], by = 1) )
			for(i in 1:pn){
				nnx = paste("mshape", i, sep="")
				pars[(nx+i), 1] = 5
				if(any(substr(start.names, 1, nchar(nnx))==nnx)){
					nix = which(substr(start.names, 1, nchar(nnx)) == nnx)
					pars[(nx+i), 1] = start.pars[nix]
				}		
				pars[(nx+i), 3] = 1
				pars[(nx+i), 5] = 4
				pars[(nx+i), 6] = 50
				if(any(substr(fixed.names, 1, nchar(nnx))==nnx)){
					nix = which(substr(fixed.names, 1, nchar(nnx)) == nnx)
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
	pidx[3,2] = nx+pn
	
	nx = nx + pn
	pn = 1
	pidx[4,1] = nx+1
	
	if(modelinc[7]<=1){
		if(pos.matrix[4,3]==1){
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
		if(pos.matrix[4,3] == 1){
			pn = length( seq(pos.matrix[4,1], pos.matrix[4,2], by = 1) )
			for(i in 1:pn){
				nnx = paste("mskew", i, sep="")
				pars[(nx+i), 1] = 3
				pars[(nx+i), 5] = -1
				pars[(nx+i), 6] = 1
				if(any(substr(start.names, 1, nchar(nnx))==nnx)){
					nix = which(substr(start.names, 1, nchar(nnx)) == nnx)
					pars[(nx+i), 1] = start.pars[nix]
				}		
				if(any(substr(fixed.names, 1, nchar(nnx))==nnx)){
					nix = which(substr(fixed.names, 1, nchar(nnx)) == nnx)
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
	pidx[4,2] = nx+pn
	rownames(pars) = pnames
	modeldesc$type = "2-step"
	model = list(modelinc = modelinc, modeldesc = modeldesc, modeldata = modeldata, varmodel = varmodel, 
			pars = pars, start.pars = start.pars, fixed.pars = fixed.pars, maxgarchOrder = maxgarchOrder, 
			sidx = sidx, fdccindex = fdccindex, pos.matrix = pos.matrix, pidx = pidx)
	model$DCC = "FDCC"	
	ans = new("DCCspec",
			model = model,
			umodel = umodel)
	return(ans)
}


.fdccfit = function(spec, data, out.sample = 0, solver = "solnp", 
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
	} else{		
		fitlist = multifit(multispec = mspec, data = xts(zdata, index), out.sample = n.start, 
				solver = garch.solver, solver.control = solver.control, 
				fit.control = ufit.control, cluster = cluster, realizedVol = realizedVol)
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
	Qbar = (t(stdresid) %*% stdresid)/nrow(stdresid)
	# In the VAR model we postpadded the first lag values with a constant.
	# In an ARMA model, the first lag values are also constant.
	#-----------------------------------------------------------------------------------
	# Create the Full Model Pars and Indices
	modelinc =  model$modelinc
	# create full par matrix
	tmp = .fullinc_fdcc(modelinc, umodel)
	midx = tmp$midx
	fdccAidx = tmp$fdccAidx
	fdccBidx = tmp$fdccBidx
	
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
	arglist$H = H
	arglist$Qbar = Qbar
	# FDCC specific
	arglist$fdccAidx = fdccAidx
	arglist$fdccBidx = fdccBidx
	arglist$sidx = model$sidx
	#-----------------------------------------------------------------------------------	
	# Check for fixed parameters in DCC
	if(any(ipars[,2]==1)){
		if(npars == 0){
			if(fit.control$eval.se==0) {
				warning("\nfdccfit-->warning: all parameters fixed...returning dccfilter object instead\n")
				xspex = spec
				for(i in 1:m) xspex@umodel$fixed.pars[[i]] = as.list(fitlist@fit[[i]]@model$pars[fitlist@fit[[i]]@model$pars[,3]==1,1])
				return(.fdccfilter(spec = xspex, data = xts(data, index), out.sample = out.sample, cluster = cluster, VAR.fit = VAR.fit,
								realizedVol = realizedVol))
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
	# calculate constraint length:
	if(modelinc[3]==1){
		ILB = 0
		IUB = 1
	} else{
		if(modelinc[4]>0) xx = matrix(ipars[arglist$model$pidx["fdcca", 1]:arglist$model$pidx["fdcca", 2], 1], ncol = modelinc[4]) else xx = matrix(0, nrow = modelinc[3], ncol = modelinc[5])
		if(modelinc[5]>0) yy = matrix(ipars[arglist$model$pidx["fdccb", 1]:arglist$model$pidx["fdccb", 2], 1], ncol = modelinc[5]) else yy = matrix(0, nrow = modelinc[3], ncol = modelinc[4])
		# each pairwise combination
		nc = NROW(expand.grid(xx, xx, yy, yy))
		ILB = rep(0, nc)
		IUB = rep(1, nc)
	}
	# Asymmetric Spec has different constraints
	Ifn = .fdcccon
	if( solver == "solnp" | solver == "gosolnp") fit.control$stationarity = FALSE else fit.control$stationarity = TRUE
	arglist$fit.control = fit.control
	
	if( use.solver )
	{
		arglist$returnType = "llh"
		solution = switch(model$modeldesc$distribution,
				mvnorm = .dccsolver(solver, pars = ipars[estidx, 1], fun = normal.fdccLLH1, Ifn, ILB, 
						IUB, gr = NULL, hessian = NULL, control = solver.control, 
						LB = ipars[estidx, 5], UB = ipars[estidx, 6], arglist = arglist))
				#mvlaplace = .dccsolver(solver, pars = ipars[estidx, 1], fun = laplace.fdccLLH1, Ifn, ILB, 
				#		IUB, gr = NULL, hessian = NULL, control = solver.control, 
				#		LB = ipars[estidx, 5], UB = ipars[estidx, 6], arglist = arglist),
				#mvt = .dccsolver(solver, pars = ipars[estidx, 1], fun = student.fdccLLH1, Ifn, ILB, 
				#		IUB, gr = NULL, hessian = NULL, control = solver.control, 
				#		LB = ipars[estidx, 5], UB = ipars[estidx, 6], arglist = arglist))
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
				mvnorm = .dccmakefitmodel(garchmodel = "fdccnorm", f = normal.fdccLLH2, 
						arglist = arglist, timer = 0, message = sol$message, fname = "normal.fdccLLH2"))
				#mvlaplace = .dccmakefitmodel(garchmodel = "fdcclaplace", f = laplace.fdccLLH2, 
				#		arglist = arglist, timer = 0, message = sol$message, fname = "laplace.fdccLLH2"),		
				#mvt =.dccmakefitmodel(garchmodel = "fdccstudent", f = student.fdccLLH2, 
				#		arglist = arglist, timer = 0, message = sol$message, fname = "student.fdccLLH2"))
		fit$timer = Sys.time() - tic
	} else{
		fit$message = sol$message
		fit$convergence = 1
	}
	# make model list to return some usefule information which
	model$mpars = mpars
	model$ipars = ipars
	model$pars[,1] = ipars[,1]
	model$midx = midx
	model$eidx = eidx
	model$umodel = umodel
	fit$Qbar = Qbar
	fit$realizedVol = realizedVol
	ans = new("DCCfit",
			mfit = fit,
			model = model)
	return(ans)
}

getfdccpars = function(pars, model){
	modelinc = model$modelinc
	pidx = model$pidx
	sidx = model$sidx
	if(modelinc[4]>0){
		A = matrix(pars[pidx["fdcca",1]:pidx["fdcca",2],1], ncol = modelinc[4], byrow=FALSE)
		A = apply(A, 2, function(x) sidx %*% x)
	} else{
		A = matrix(0)
	}
	if(modelinc[5]>0){
		B = matrix(pars[pidx["fdccb",1]:pidx["fdccb",2],1], ncol = modelinc[5], byrow=FALSE)
		B = apply(B, 2, function(x) sidx %*% x)
	} else{
		B = matrix(0)
	}
	# Because only the DCC (1,1) model is allowed, this is simplified:
	BB = (B %*% t(B))
	AA = (A %*% t(A))
	C = matrix(1, NROW(BB), NCOL(BB)) - BB - AA
	return(list(A = A, B = B, C = C))
}

.fdccfilter = function(spec, data, out.sample = 0, filter.control = list(n.old = NULL), 
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
			out.sample = out.sample, cluster = cluster, n.old = n.old,
			realizedVol = realizedVol)
	
	if(spec@model$modelinc[1]>0) model$mu = mu else model$mu = fitted(filterlist)
	model$residuals = res = residuals(filterlist)
	model$sigma = sig = sigma(filterlist)
	stdresid = res/sig
	
	if(is.null(filter.control$n.old)) dcc.old = dim(stdresid)[1] else dcc.old = n.old
	#-----------------------------------------------------------------------------------	
	# Create the Full Model Pars and Indices
	
	modelinc =  model$modelinc
	# create full par matrix
	tmp = .fullinc_fdcc(modelinc, umodel)
	midx = tmp$midx
	fdccAidx = tmp$fdccAidx
	fdccBidx = tmp$fdccBidx
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
	arglist$realizedVol = realizedVol
	arglist$model = model
	arglist$filterlist = filterlist
	arglist$umodel = umodel
	arglist$midx = midx
	arglist$eidx = eidx
	arglist$mpars = mpars
	arglist$ipars = ipars
	arglist$estidx = estidx
	arglist$dccN = npars
	arglist$stdresid = stdresid
	arglist$H = H
	arglist$Qbar = Qbar
	# FDCC specific
	arglist$fdccAidx = fdccAidx
	arglist$fdccBidx = fdccBidx
	arglist$sidx = model$sidx
	assign("rmgarch_llh", 1, envir = mgarchenv)
	#-----------------------------------------------------------------------------------	
	
	filt = switch(model$modeldesc$distribution,
			mvnorm = .dccmakefiltermodel(garchmodel = "fdccnorm", f = normalfilter.fdccLLH2, 
					arglist = arglist, timer = 0, message = 0, fname = "normalfilter.fdccLLH2"))
	
	model$mpars = mpars
	model$ipars = ipars
	model$pars[,1] = ipars[,1]
	model$midx = midx
	model$eidx = eidx
	model$umodel = umodel	
	filt$Qbar = Qbar
	filt$realizedVol = realizedVol
	filt$timer = Sys.time() - tic
	ans = new("DCCfilter",
			mfilter = filt,
			model = model)
	return(ans)
}

.fdccforecast = function(fit, n.ahead = 1, n.roll = 0, 
		external.forecasts = list(mregfor = NULL, vregfor = NULL), cluster = NULL, ...)
{
	# checks for out.sample wrt n.roll
	model = fit@model
	modelinc = model$modelinc
	ns = model$modeldata$n.start
	if( n.roll > ns ) stop("n.roll must not be greater than out.sample!")
	if(n.roll>1 && n.ahead>1) stop("\ngogarchforecast-->error: n.ahead must be equal to 1 when using n.roll\n")
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
		#if(n.ahead!=1)    stop("\ndccforecast-->error: only n.ahead==1 supported\n")
		if( n.roll > ns ) stop("\ndccforecast-->error: n.roll greater than out.sample!", call. = FALSE)
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
	ans = .fdccforecastm(fit, n.ahead = n.ahead, n.roll = n.roll, external.forecasts = exf, cluster = cluster, ...)
	
	if(modelinc[1]==0) mu = ans$mu
	
	model$n.roll = n.roll
	model$n.ahead = n.ahead
	# keep last 50 points of R and H for plotting
	model$R = last( rcor(fit), min(length(fit@mfit$R), 50) )
	model$H = last( rcov(fit), min(length(fit@mfit$R), 50) )
	mforecast = list( H = ans$H, R = ans$R, Q = ans$Q, mu = mu )
	
	ans = new("DCCforecast",
			mforecast = mforecast,
			model = model)
	
	return( ans )
}


# n-ahead forecast based on the approximation method described in:
# ENGLE, R. and SHEPPARD (2001), K. Theoretical and empirical properties of
# dynamic conditional correlation multivariate garch. NBER Working Papers, No. 8554

.fdccforecastm = function(fit, n.ahead = 1, n.roll = 10, 
		external.forecasts = list(mregfor = NULL, vregfor = NULL), 
		cluster = NULL, ...)
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
			n.old = fit@model$modeldata$T, cluster = cluster, realizedVol = fit@mfit$realizedVol)
	
	# all.equal(head(sigma(filterlist)), head(fit@model$sigma) )
	n.roll = n.roll + 1
	m = length(mspec@spec)
	out.sample = fit@model$modeldata$n.start
	mo = max(modelinc[4:5])
	# we augmented n.roll before, now subtract
	# forclist = multiforecast(multifitORspec = fitlist, n.ahead = n.ahead, n.roll = n.roll - 1)
	forclist = multiforecast(multifitORspec = mspec, data = xts(zdata, fit@model$modeldata$index[1:nrow(zdata)]), n.ahead = n.ahead, 
			out.sample = ns, n.roll = n.roll - 1, external.forecasts = external.forecasts, 
			cluster = cluster, realizedVol = fit@mfit$realizedVol, ...)
	# ToDo:
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
	T = dim(fit@mfit$H)[3]
	## call .nsdccforecast with rolling window
	Rbar = Rtfor = Htfor = Qtfor = vector(mode = "list", length = n.roll)
	Qstart = last( rcor(fit, type = "Q"), mo )
	Rstart = last( rcor(fit, type = "R"), mo )
	Hstart = last( rcov(fit), mo)
	
	f = lapply(forclist@forecast, function(x) sigma(x))
	for(i in 1:n.roll){
		xQbar = cov(stdresid[1:(T + i - 1), ])
		xstdresids = stdresid[(T - mo + i ):(T  +  i - 1), , drop = FALSE]
		xfsig = matrix( sapply( f, function(x) x[,i] ), ncol = m)
		ans = .fdccforecastn(model, Qstart, Rstart, Hstart, xQbar, xstdresids, xfsig, n.ahead, mo)
		Rtfor[[i]] = ans$Rtfor
		Qtfor[[i]] = ans$Qtfor
		Htfor[[i]] = ans$Htfor
		Qstart = last( rugarch:::.abind(Qstart, ans$Qtfor[, , 1]), mo )
		Rstart = last( rugarch:::.abind(Rstart, ans$Rtfor[, , 1]), mo )
		Hstart = last( rugarch:::.abind(Hstart, ans$Htfor[, , 1]), mo )
	}
	
	forc = list( H = Htfor, R = Rtfor, Q = Qtfor, mu = mu )
	
	return(forc)
}

.fdccforecastn = function(model, Qstart, Rstart, Hstart, Qbar, stdresids, fsig, n.ahead, mo)
{
	m = dim(Qbar)[1]
	modelinc = model$modelinc
	Qtfor = Rtfor = Htfor = array(NA, dim = c(m, m, n.ahead + mo))
	Qtfor[ , , 1:mo]  =  Qstart[, , 1:mo]
	Rtfor[ , , 1:mo]  =  Rstart[, , 1:mo]
	Htfor[ , , 1:mo] = Hstart[, , 1:mo]
	tmp = getfdccpars(model$ipars, model)
	A = tmp$A
	B = tmp$B
	BB = (B %*% t(B))
	AA = (A %*% t(A))
	C = matrix(1, NROW(BB), NCOL(BB)) - BB - AA
	Q = C  * Qbar
	Qstar = diag( 1/sqrt( diag(Q) ) , m, m)
	Rbar = Qstar %*% Q %*% t(Qstar)
	for(i in 1:n.ahead){
		if(i==1){
			Qtfor[, , mo + i] = Q
			if(modelinc[4]>0){
				for(j in 1:modelinc[4]){
					Qtfor[ , , mo + i] = Qtfor[ , , mo + i] + (A[,j]%*%t(A[,j])) * (stdresids[(mo + i - j), ] %*% t(stdresids[(mo + i - j), ]))
				}
			}
			if(modelinc[5]>0){
				for(j in 1:modelinc[5]){
					Qtfor[ , , mo + i] = Qtfor[ , , mo + i] + (B[,j]%*%t(B[,j])) * Qtfor[ , , mo + i - j]
				}
			}
			Qtmp = diag( 1/sqrt( diag(Qtfor[ , , mo + i]) ) , m, m)
			Rtfor[ , , mo + i] =  Qtmp %*% Qtfor[ , , mo + i] %*% t(Qtmp)
			Dtmp = diag(fsig[i, ], m, m)
			Htfor[ , , mo + i] = Dtmp %*% Rtfor[ , , mo + i] %*% Dtmp
		} else{
			Rtfor[ , , mo + i] = Rbar
			if(modelinc[4]>0){
				for(j in 1:modelinc[4]){
					Rtfor[ , , mo + i] = Rtfor[ , , mo + i] + (A[,j]%*%t(A[,j])) * (Rtfor[ , , mo + i - j]-Rbar)
				}
			}
			if(modelinc[5]>0){
				for(j in 1:modelinc[5]){
					Rtfor[ , , mo + i] = Rtfor[ , , mo + i] + (B[,j]%*%t(B[,j])) * (Rtfor[ , , mo + i - j]-Rbar)
				}
			}
			Dtmp = diag(fsig[i, ], m, m)
			Htfor[ , , mo + i] = Dtmp %*% Rtfor[ , , mo + i] %*% Dtmp
		}
	}
	return( list( Rtfor = Rtfor[, , -(1:mo), drop = FALSE], Htfor = Htfor[, , -(1:mo), drop = FALSE], 
					Qtfor = Qtfor[, , -(1:mo), drop = FALSE] ) )
}

.fdccforecastn2 = function(model, Qstart, Rstart, Hstart, Qbar, stdresids, fsig, n.ahead, mo)
{
	m = dim(Qbar)[1]
	modelinc = model$modelinc
	Qtfor = Rtfor = Htfor = array(NA, dim = c(m, m, n.ahead + mo))
	Qtfor[ , , 1:mo]  =  Qstart[, , 1:mo]
	Rtfor[ , , 1:mo]  =  Rstart[, , 1:mo]
	Htfor[ , , 1:mo] = Hstart[, , 1:mo]
	tmp = getfdccpars(model$ipars, model)
	A = tmp$A
	B = tmp$B
	BB = (B %*% t(B))
	AA = (A %*% t(A))
	C = matrix(1, NROW(BB), NCOL(BB)) - BB - AA
	Q = C * Qbar
	for(i in 1:n.ahead){
		if(i==1){
			Qtfor[, , mo + i] = Q
			if(modelinc[4]>0){
				for(j in 1:modelinc[4]){
					Qtfor[ , , mo + i] = Qtfor[ , , mo + i] + (A[,j]%*%t(A[,j])) * (stdresids[(mo + i - j), ] %*% t(stdresids[(mo + i - j), ]))
				}
			}
			if(modelinc[5]>0){
				for(j in 1:modelinc[5]){
					Qtfor[ , , mo + i] = Qtfor[ , , mo + i] + (B[,j]%*%t(B[,j])) * Qtfor[ , , mo + i - j]
				}
			}
			Qtmp = diag( 1/sqrt( diag(Qtfor[ , , mo + i]) ) , m, m)
			Rtfor[ , , mo + i] =  Qtmp %*% Qtfor[ , , mo + i] %*% t(Qtmp)
			Dtmp = diag(fsig[i, ], m, m)
			Htfor[ , , mo + i] = Dtmp %*% Rtfor[ , , mo + i] %*% Dtmp
		} else{
			Qtfor[, , mo + i] = Q
			if(modelinc[4]>0){
				for(j in 1:modelinc[4]){
					Qtfor[ , , mo + i] = Qtfor[ , , mo + i] + (A[,j]%*%t(A[,j])) * (Qtfor[ , , mo + i - j])
				}
			}
			if(modelinc[5]>0){
				for(j in 1:modelinc[5]){
					Qtfor[ , , mo + i] = Qtfor[ , , mo + i] + (B[,j]%*%t(B[,j])) * (Qtfor[ , , mo + i - j])
				}
			}
			Qtmp = diag( 1/sqrt( diag(Qtfor[ , , mo + i]) ) , m, m)
			Rtfor[ , , mo + i] =  Qtmp %*% Qtfor[ , , mo + i] %*% t(Qtmp)
			Dtmp = diag(fsig[i, ], m, m)
			Htfor[ , , mo + i] = Dtmp %*% Rtfor[ , , mo + i] %*% Dtmp
		}
	}
	return( list( Rtfor = Rtfor[, , -(1:mo), drop = FALSE], Htfor = Htfor[, , -(1:mo), drop = FALSE], 
					Qtfor = Qtfor[, , -(1:mo), drop = FALSE] ) )
}

.fdccsim.fit = function(fitORspec, n.sim = 1000, n.start = 0, m.sim = 1, 
		startMethod = c("unconditional", "sample"), presigma = NULL, 
		preresiduals = NULL, prereturns = NULL, preQ = NULL, preZ = NULL, 
		Qbar = NULL, rseed = NULL, mexsimdata = NULL, 
		vexsimdata = NULL, cluster = NULL, VAR.fit = NULL, prerealized = NULL, ...)
{
	fit = fitORspec
	T = fit@model$modeldata$T
	Data = fit@model$modeldata$data[1:T,]
	m = dim(Data)[2]
	mo = max(fit@model$modelinc[4:5])
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
	tmp = getfdccpars(model$ipars, model)
	A = tmp$A
	B = tmp$B
	BB = (B %*% t(B))
	AA = (A %*% t(A))
	C = matrix(1, NROW(BB), NCOL(BB)) - BB - AA
	# model, Z, A, B, C, Qbar, preQ, Rbar, mo, n.sim, n.start, m, rseed
	simRes = simX = simR = simQ = simH = simSeries = vector(mode = "list", length = m.sim)
	
	if( !is.null(cluster) ){
		simH = vector(mode = "list", length = m.sim)
		simX = vector(mode = "list", length = m.sim)
		clusterEvalQ(cluster, require(rmgarch))
		clusterExport(cluster, c("model", "z", "A", "B", "C", 
						"preQ", "Rbar", "Qbar", "mo", "n.sim", "n.start", "m", 
						"rseed",".fdccsimf"), envir = environment())
		mtmp = parLapply(cluster, as.list(1:m.sim), fun = function(j){
					.fdccsimf(model, Z = z[,,j], A = A, B = B, C = C, 
							Qbar = Qbar, preQ = preQ, Rbar = Rbar, mo = mo, 
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
						"prereturns", "mexsimdata", "vexsimdata","prerealized"), 
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
					.fdccsimf(model, Z = z[,,j], A = A, B = B, C = C, 
							Qbar = Qbar, preQ = preQ, Rbar = Rbar, mo = mo, 
							n.sim, n.start, m, rseed[j])
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

.fdccsim.spec = function(fitORspec, n.sim = 1000, n.start = 0, m.sim = 1, 
		startMethod = c("unconditional", "sample"), presigma = NULL, 
		preresiduals = NULL, prereturns = NULL, preQ = NULL, preZ = NULL, 
		Qbar = NULL, rseed = NULL, mexsimdata = NULL, 
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
	mo = max(spec@model$modelinc[4:5])
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
			tmp = rugarch:::rged(m * (n.sim + n.start) * m.sim, 0, 1, shape = 1)
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
			tmp = rugarch:::rstd(m * (n.sim + n.start) * m.sim, 0, 1, shape = model$mpars["mshape",1])
			z = array(NA,  dim = c(n.sim + n.start + mo, m, m.sim))
			for(i in 1:m.sim) z[,,i] = rbind(preZ, tmp)
		} else{
			z = array(NA, dim = c(n.sim + n.start + mo, m, m.sim))
			for(i in 1:m.sim){
				set.seed( rseed[i] )
				z[,,i] = rbind(preZ, matrix(rugarch:::rstd(m * (n.sim + n.start), 0, 1, shape = model$mpars["mshape",1]), nrow = n.sim + n.start, ncol = m))
			}
		}
	}
	
	# ok now to expand rseed
	if(length(rseed) == 1){
		rseed = c(rseed, as.integer(runif(m.sim, 1, Sys.time()))) 
	}
	
	tmp = getfdccpars(model$ipars, model)
	A = tmp$A
	B = tmp$B
	BB = (B %*% t(B))
	AA = (A %*% t(A))
	C = matrix(1, NROW(BB), NCOL(BB)) - BB - AA
	
	simRes = simX = simR = simQ = simH = simSeries = vector(mode = "list", length = m.sim)
	
	if( !is.null(cluster) ){
		simH = vector(mode = "list", length = m.sim)
		simX = vector(mode = "list", length = m.sim)
		clusterEvalQ(cluster, require(rmgarch))
		clusterExport(cluster, c("model", "z", "A", "B", "C", "preQ", 
						"Rbar", "Qbar", "mo", "n.sim", "n.start", "m", 
						"rseed",".fdccsimf"), envir = environment())
		mtmp = parLapply(cluster, as.list(1:m.sim), fun = function(j){
					.fdccsimf(model, Z = z[,,j], A = A, B = B, C = C, 
							Qbar = Qbar, preQ = preQ, Rbar = Rbar, mo = mo, 
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
						"prereturns", "mexsimdata", "vexsimdata","prerealized"), 
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
		mtmp = lapply(as.list(1:m.sim), FUN = function(j){
					.fdccsimf(model, Z = z[,,j], A = A, B = B, C = C, 
							Qbar = Qbar, preQ = preQ, Rbar = Rbar, mo = mo, 
							n.sim, n.start, m, rseed[j])
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


.fdccsimf = function(model, Z, A, B, C, Qbar, preQ, Rbar, mo, n.sim, n.start, m, rseed)
{
	modelinc = model$modelinc
	ipars = model$pars
	idx = model$pidx
	n = n.sim + n.start + mo
	set.seed(rseed[1]+1)
	# All distributions here make use of this since they are mixtures of normals
	NZ = matrix(rnorm(m * (n.sim+n.start+mo)), nrow = n.sim + n.start + mo, ncol = m)
	
	# SEXP model, SEXP A, SEXP B, SEXP C, SEXP Qbar, SEXP preQ, SEXP Rbar, SEXP Z, SEXP NZ
	res = switch(model$modeldesc$distribution,
			mvnorm = .Call("fdccsimmvn", 
					model = as.integer(modelinc), 
					A = A, 
					B = B, 
					C = C,
					Qbar = as.matrix(Qbar),  
					preQ = as.matrix(preQ), 
					Rbar = as.matrix(Rbar),  
					Z = as.matrix(Z), 
					NZ = as.matrix(NZ), 
					PACKAGE = "rmgarch")
	)
	Q = array(NA, dim = c(m, m, n.sim + n.start + mo))
	R = array(NA, dim = c(m, m, n.sim + n.start + mo))
	for(i in 1:(n.sim + n.start + mo)){
		R[,,i] = res[[2]][[i]]
		Q[,,i] = res[[1]][[i]]
	}
	ans = list( Q = Q, R = R, Z = res[[3]])
	return( ans )
}

.fullinc_fdcc = function(modelinc, umodel){
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
	# fdccC
	# fdccA
	# fdccB
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
	
	sumdcc = max(1, (modelinc[4]*modelinc[3]))+max(1, (modelinc[5]*modelinc[3])) + 
			max(1, modelinc[6]) + max(1, modelinc[6])
	
	tmpmat = rbind(tmpmat, matrix(0, ncol = m+1, nrow = sumdcc))
	# intercept
	if(modelinc[4]>0){
		xid = 0
		for(i in 1:modelinc[4]){
			for(j in 1:modelinc[3]){
				tmpmat[(nx+xid+j), m+1] = 1
				pnames = c(pnames, paste("fdcca[g=", j,",P=", i, "]", sep=""))
			}
			xid = xid + modelinc[3]
		}
	} else{
		pnames = c(pnames, "fdcca")
	}
	fdccAidx = (nx+1):max( (nx+1), (nx+modelinc[3]*modelinc[4]) )
	
	nx = nx + max(1, modelinc[4]*modelinc[3])
	
	if(modelinc[5]>0){
		xid = 0
		for(i in 1:modelinc[5]){
			for(j in 1:modelinc[3]){
				tmpmat[(nx+xid+j), m+1] = 1
				pnames = c(pnames, paste("fdccb[g=", j,",Q=", i, "]", sep=""))
			}
			xid = xid + modelinc[3]
		}
	} else{
		pnames = c(pnames, "fdccb")
	}
	fdccBidx = (nx+1):max( (nx+1), (nx+modelinc[3]*modelinc[5]) )
	
	nx = nx + max(1, modelinc[5]*modelinc[3])
	
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
	return(list(midx = tmpmat, fdccAidx = fdccAidx, fdccBidx = fdccBidx))
}


.fdcccon = function(pars, arglist){
	ipars = arglist$ipars
	estidx = arglist$estidx
	idx = arglist$model$pidx
	ipars[estidx, 1] = pars
	modelinc = arglist$model$modelinc
	if(modelinc[4]>0) fdccA = matrix(ipars[idx["fdcca", 1]:idx["fdcca", 2], 1], ncol = modelinc[4]) else fdccA = matrix(0, nrow = modelinc[3], ncol = modelinc[5])
	if(modelinc[5]>0) fdccB = matrix(ipars[idx["fdccb", 1]:idx["fdccb", 2], 1], ncol = modelinc[5]) else fdccB = matrix(0, nrow = modelinc[3], ncol = modelinc[4])
	if(modelinc[3]==1){
		ans = fdccA + fdccB
	} else{
		# each pairwise combination (only allow DCC(1,1))
		ans = apply(expand.grid(fdccA, fdccA, fdccB, fdccB), 1, function(x) x[1]*x[2] + x[3]*x[4])
	}
	return( as.numeric(ans) )
}