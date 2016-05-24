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

#----------------------------------------------------------------------------------
# univariate spec method
#----------------------------------------------------------------------------------
ugarchspec = function(variance.model = list(model = "sGARCH", garchOrder = c(1,1),
				submodel = NULL, external.regressors = NULL, variance.targeting = FALSE),
		mean.model = list(armaOrder = c(1,1), include.mean = TRUE, archm = FALSE,
				archpow = 1, arfima = FALSE, external.regressors = NULL, archex = FALSE),
		distribution.model = "norm", start.pars = list(), fixed.pars = list(), ...)
{
	UseMethod("ugarchspec")
}

# [mu ar ma arfima im mxreg omega alpha beta gamma gamma11 gamma21 delta lambda vxreg skew shape dlamda aux aux aux aux]

.expand.model = function(model){
	modelnames = NULL
	for(i in 1:21){
		if(model[i]>0){
			if(any(c(2,3,6,8,9,10,11,12,15) == i)){
				modelnames = c(modelnames, paste(names(model)[i], 1:model[i], sep = ""))
			} else{
				modelnames = c(modelnames, names(model)[i])
			}
		}
	}
	return( modelnames )
}

# Changelog:
# 06-12-2011 Added "archex" option in mean.model so that external regressor might
# be multiplied by conditional variance. If used must be integer and represents the
# number of series from the end of the supplied external regressors.

# modelinc[20] = archex
# modelinc[21] = distribution no
# modelinc[22] = custom variance targeting no

.ugarchspec = function(variance.model = list(model = "sGARCH", garchOrder = c(1,1),
				submodel = NULL, external.regressors = NULL, variance.targeting = FALSE),
		mean.model = list(armaOrder = c(1,1), include.mean = TRUE, archm = FALSE,
				archpow = 1, arfima = FALSE, external.regressors = NULL, archex = FALSE),
		distribution.model = "norm", start.pars = list(), fixed.pars = list())
{
	# some checks and preparation to be passed on to specific models by switch
	modelinc = rep(0, 22)
	# set the custom variance target to NA
	modelinc[22] = NA
	# [19] from aux to xi to accomodate realGARCH model (25-12-2014)
	names(modelinc) = c("mu", "ar", "ma", "arfima", "archm", "mxreg", "omega", "alpha",
			"beta", "gamma", "eta1", "eta2", "delta", "lambda", "vxreg", "skew", "shape",
			"ghlambda", "xi", "aux", "aux", "aux")
	modeldesc = list()
	modeldata = list()

	# check the option parameters specified and stop on error
	mm = match(names(mean.model), c("armaOrder", "include.mean", "archm", "archpow", "arfima", "external.regressors", "archex"))
	if(any(is.na(mm))){
		idx = which(is.na(mm))
		enx = NULL
		for(i in 1:length(idx)) enx = c(enx, names(mean.model)[idx[i]])
		warning(paste(c("unidentified option(s) in mean.model:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
	}
	vm = match(names(variance.model), c("model", "garchOrder", "submodel", "external.regressors", "variance.targeting"))
	if(any(is.na(vm))){
		idx = which(is.na(vm))
		enx = NULL
		for(i in 1:length(idx)) enx = c(enx, names(variance.model)[idx[i]])
		warning(paste(c("unidentified option(s) in variance.model:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
	}

	# distribution model
	if(is.null(distribution.model)) modeldesc$distribution = "norm"
	valid.distribution = c("norm", "snorm", "std", "sstd","ged", "sged", "nig", "ghyp", "jsu", "ghst")
	distribution = distribution.model

	if(!is.character(distribution[1]))
		stop("\nugarchspec-->error: cond.distribution argument must be a character")

	if(!any(distribution==valid.distribution))
		stop("\nugarchspec-->error: the cond.distribution does not appear to be a valid choice.")

	if(length(distribution)!=1) distribution = distribution[1]

	modeldesc$distribution = distribution
	modeldesc$distno = which(distribution == valid.distribution)
	di = .DistributionBounds(distribution)
	modelinc[16] = di$include.skew
	modelinc[17] = di$include.shape
	modelinc[18] = di$include.ghlambda
	# the last aux value is the distribution number
	modelinc[21] = modeldesc$distno
	# variance model:
	vmodel = list()

	valid.model = c("sGARCH", "eGARCH", "gjrGARCH", "tGARCH", "fGARCH", "iGARCH", "apARCH", "csGARCH", "mcsGARCH","realGARCH")
	if(is.null(variance.model$model)){
		modeldesc$vmodel = "sGARCH"
	} else{
		modeldesc$vmodel = variance.model$model[1]

		if(!is.character(modeldesc$vmodel))
			stop("\nugarchspec-->error: garch model argument must be a character.\n", call. = FALSE)

		if(!any(modeldesc$vmodel == valid.model))
			stop("\nugarchspec-->error: the garch model does not appear to be a valid choice.\n", call. = FALSE)

		if(modeldesc$vmodel == "fGARCH"){
			modeldesc$vsubmodel = variance.model$submodel
			valid.submodel = c("GARCH","TGARCH","AVGARCH","NGARCH","NAGARCH","APARCH","ALLGARCH","GJRGARCH")
			if(is.null(modeldesc$vsubmodel))
				stop("\nugarchspec-->error: NULL not allowed for the submodel when model is of type fGARCH.\n", call. = FALSE)
			if(!any(modeldesc$vsubmodel == valid.submodel))
				stop("\nugarchspec-->error: the fGARCH submodel does not appear to be a valid choice. See documentation for valid choices.\n", call. = FALSE)
		}
	}

	# depending on model include the additional parameters
	if(is.null(variance.model$garchOrder)){
		modelinc[8] = 1
		modelinc[9] = 1
	} else{
		modelinc[8] = variance.model$garchOrder[1]
		modelinc[9] = variance.model$garchOrder[2]
	}

	if( modeldesc$vmodel == "gjrGARCH" ) modelinc[10] = modelinc[8]
	if( modeldesc$vmodel == "eGARCH" ) modelinc[10] = modelinc[8]
	if( modeldesc$vmodel == "fiGARCH" ){
		# fractional parameter (FIGARCH): \delta
		# Hyperbolic Alpha Parameters (HYGARCH): \lambda
		modelinc[13] = 1
		# modelinc[14] = 1
	}

	if( modeldesc$vmodel == "apARCH" ){
		modelinc[10] = modelinc[8]
		modelinc[13] = 1
	}
	if( modeldesc$vmodel == "fGARCH" ){
		if(modeldesc$vsubmodel == "AVGARCH"){
			modelinc[12] = modelinc[11] = modelinc[8]
		}
		if(modeldesc$vsubmodel == "GJRGARCH") modelinc[11] = modelinc[8]
		if(modeldesc$vsubmodel == "TGARCH") modelinc[11] = modelinc[8]
		if(modeldesc$vsubmodel == "NGARCH") modelinc[14] = 1
		if(modeldesc$vsubmodel == "NAGARCH") modelinc[12] = modelinc[8]
		if(modeldesc$vsubmodel == "APARCH"){
			modelinc[14] = 1
			modelinc[11] = modelinc[8]
		}
		if(modeldesc$vsubmodel == "ALLGARCH"){
			modelinc[12] = modelinc[11] = modelinc[8]
			modelinc[14] = 1
		}
	}
	if( modeldesc$vmodel == "csGARCH" ){
		modelinc[12] = modelinc[11] = 1
		variance.model$variance.targeting = vmodel$variance.targeting = FALSE
	}
	if( modeldesc$vmodel == "iGARCH" && modelinc[9] == 0 ) 	stop("\nugarchspec-->error: the iGARCH model requires the GARCH beta parameter.\n", call. = FALSE)

	if( modeldesc$vmodel == "realGARCH"){
		# xi
		modelinc[19] = 1
		# \xi  + \delta \sigma _t^2 + {\eta _1}{z_t} + {\eta _2}\left( {z_t^2 - 1} \right) + {u_t}
		# lambda is now the variance of u_t
		modelinc[11:14] = 1
	}

	modeldata$vexdata = variance.model$external.regressors
	if( !is.null(variance.model$external.regressors) ) modelinc[15] = dim( variance.model$external.regressors )[2]

	if(is.null(variance.model$variance.targeting)){
		modelinc[7] = 1
	} else{
		if(is.logical(variance.model$variance.targeting)){
			modelinc[7] = as.integer( 1 - variance.model$variance.targeting )
		} else{
			modelinc[7] = 0
			modelinc[22] = as.numeric(variance.model$variance.targeting)
		}
	}

	# mean model:
	if(is.null(mean.model$armaOrder)){
		modelinc[2] = modelinc[3] = 1
	} else{
		modelinc[2] = mean.model$armaOrder[1]
		modelinc[3] = mean.model$armaOrder[2]
	}
	if(is.null(mean.model$include.mean)) modelinc[1] = 1 else modelinc[1] = as.integer( mean.model$include.mean )

	if(is.null(mean.model$archm) || !mean.model$archm){
		modelinc[5] = 0
	} else{
		if(modeldesc$vmodel=="mcsGARCH") stop("\narchm not supported by mcsGARCH model.")
		if(is.null(mean.model$archpow)) mean.model$archpow = 1
		modelinc[5] = as.integer( mean.model$archpow )
	}


	if(is.null(mean.model$arfima)) modelinc[4] = 0 else modelinc[4] = as.integer( mean.model$arfima )

	modeldata$mexdata = mean.model$external.regressors
	if( !is.null(mean.model$external.regressors) ) modelinc[6] = dim( mean.model$external.regressors )[2]

	if(is.null(mean.model$archex) || !mean.model$archex){
		modelinc[20] = 0
	} else{
		if(modeldesc$vmodel=="mcsGARCH") stop("\narchex not supported by mcsGARCH model.")
		modelinc[20] = as.integer( mean.model$archex )
		if(modelinc[6] == 0) stop("\narchex cannot be used without external.regressors!!\n", call. = FALSE)
		if(modelinc[6] < modelinc[20]) stop("\narchex cannot be greater than number of external.regressors!!\n", call. = FALSE)
	}

	maxOrder = max(modelinc[c(2,3,8,9)])
	modelnames = .expand.model(modelinc)


	fmodel = NULL
	if(modeldesc$vmodel == "fGARCH"){
		valid.submodels = c("GARCH","AVGARCH","NGARCH","NAGARCH","TGARCH","GJRGARCH","APARCH","ALLGARCH")
		submodel = modeldesc$vsubmodel
		if(!any(submodel == valid.submodels)) stop("not a valid fmodel name for fGARCH specification. See documentation.")
		fspec = .fgarchModel(submodel)
		fmodel$fpars = fspec$parameters
		# lambda, delta, fb, fc, fk
		fmodel$finclude = fspec$indicator
		fmodel$fgarchtype = fspec$garchtype
		fmodel$fbounds = .fmodelBounds(submodel)
	}


	pos = 1
	pos.matrix = matrix(0, ncol = 3, nrow = 21)
	colnames(pos.matrix) = c("start", "stop", "include")
	rownames(pos.matrix) = c("mu", "ar", "ma", "arfima", "archm", "mxreg", "omega", "alpha",
			"beta", "gamma", "eta1", "eta2", "delta", "lambda", "vxreg", "skew", "shape",
			"ghlambda", "xi", "aux", "aux")

	# check if there are starting or fixed
	# check that it is included in the optimization
	# check that it is needed in the model


	if(modeldesc$vmodel == "fGARCH"){
		for(i in 1:21){
			if( modelinc[i] > 0 ){
				if(i == 11 && fmodel$finclude[3]==1){
					pos.matrix[i,1:3] = c(pos, pos+modelinc[i]-1, 1)
					pos = max(pos.matrix[1:i,2]+1)
				} else if(i == 12 && fmodel$finclude[4]==1){
					pos.matrix[i,1:3] = c(pos, pos+modelinc[i]-1, 1)
					pos = max(pos.matrix[1:i,2]+1)
				} else if(i == 13 && fmodel$finclude[2]==1){
					pos.matrix[i,1:3] = c(pos, pos+modelinc[i]-1, 1)
					pos = max(pos.matrix[1:i,2]+1)
				} else if(i == 14 && fmodel$finclude[1]==1){
					pos.matrix[i,1:3] = c(pos, pos+modelinc[i]-1, 1)
					pos = max(pos.matrix[1:i,2]+1)
				} else{
					pos.matrix[i,1:3] = c(pos, pos+modelinc[i]-1, 1)
					pos = max(pos.matrix[1:i,2]+1)
				}
			}
		}
	} else{
		for(i in 1:21){
			if( modelinc[i] > 0 ){
				pos.matrix[i,1:3] = c(pos, pos+modelinc[i]-1, 1)
				pos = max(pos.matrix[1:i,2]+1)
			}
		}
	}
	nn = length(modelnames)
	modelmatrix = matrix(0, ncol = 3, nrow = nn)
	rownames(modelmatrix) = modelnames
	colnames(modelmatrix) = c("opt", "fixed", "start")
	fixed.names = names(fixed.pars)
	fp = charmatch(fixed.names, modelnames)

	if(!is.null(fixed.names) && any(!is.na(fp))){
		fixed = fp[!is.na(fp)]
		modelmatrix[fixed,2] = 1
		fz = charmatch(modelnames, fixed.names)
		fz = fz[!is.na(fz)]
		fixed.pars = fixed.pars[fz]
		names(fixed.pars) = fixed.names[fz]
	} else{
		fixed.pars = NULL
	}
	modelmatrix[,1] = 1 - modelmatrix[,2]
	start.names = names(start.pars)
	sp = charmatch(start.names, modelnames)
	if(!is.null(start.names) && any(!is.na(sp))){
		start = sp[!is.na(sp)]
		modelmatrix[start,3] = 1
		sz = charmatch(modelnames, start.names)
		sz = sz[!is.na(sz)]
		start.pars = start.pars[sz]
	} else{
		start.pars = NULL
	}


	##################################################################
	# Parameter Matrix
	mm = sum(modelinc[c(2,3,6,8,9,10,11,12,15)])
	mm = mm - length( which(modelinc[c(2,3,6,8,9,10,11,12,15)]>0) )
	pars = matrix(0, ncol = 6, nrow = 19 + mm)
	colnames(pars) = c("Level", "Fixed", "Include", "Estimate", "LB", "UB")
	pidx = matrix(NA, nrow = 19, ncol = 2)
	colnames(pidx) = c("begin", "end")
	rownames(pidx) =  c("mu", "ar", "ma", "arfima", "archm", "mxreg", "omega", "alpha",
			"beta", "gamma", "eta1", "eta2", "delta", "lambda", "vxreg", "skew", "shape",
			"ghlambda","xi")
	fixed.names = names(fixed.pars)
	pnames = NULL
	nx = 0
	if(pos.matrix[1,3]==1){
		pars[1, 3] = 1
		pars[1, 1] = 0
		if(any(substr(fixed.names, 1, 2)=="mu")) pars[1,2] = 1 else pars[1,4] = 1

	}
	pidx[1,1] = 1
	pidx[1,2] = 1
	pnames = c(pnames, "mu")
	nx = 1
	pn = 1
	pidx[2,1] = 2
	if(pos.matrix[2,3] == 1){
		pn = length( seq(pos.matrix[2,1], pos.matrix[2,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("ar", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "ar")
	}
	pidx[2,2] = 1+pn

	nx = nx + pn
	pn = 1
	pidx[3,1] = nx+1
	if(pos.matrix[3,3] == 1){
		pn = length( seq(pos.matrix[3,1], pos.matrix[3,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("ma", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "ma")
	}
	pidx[3,2] = nx+pn

	nx = nx + pn
	pn = 1
	pidx[4,1] = nx+1
	if(pos.matrix[4,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(substr(fixed.names, 1, 6)=="arfima")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pnames = c(pnames, "arfima")
	pidx[4,2] = nx+pn

	nx = nx + pn
	pn = 1
	pidx[5,1] = nx+1
	if(pos.matrix[5,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(substr(fixed.names, 1, 5)=="archm")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pnames = c(pnames, "archm")
	pidx[5,2] = nx+pn

	nx = nx + pn
	pn = 1
	pidx[6,1] = nx+1
	if(pos.matrix[6,3]==1){
		pn = length( seq(pos.matrix[6,1], pos.matrix[6,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("mxreg", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "mxreg")
	}
	pidx[6,2] = nx+pn

	nx = nx + pn
	pn = 1
	pidx[7,1] = nx+1
	if(pos.matrix[7,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(substr(fixed.names, 1, 5)=="omega")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pnames = c(pnames, "omega")
	pidx[7,2] = nx+pn

	nx = nx + pn
	pn = 1
	pidx[8,1] = nx+1
	if(pos.matrix[8,3]==1){
		pn = length( seq(pos.matrix[8,1], pos.matrix[8,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("alpha", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "alpha")
	}
	pidx[8,2] = nx+pn

	nx = nx + pn
	pn = 1
	pidx[9,1] = nx+1
	if(pos.matrix[9,3]==1){
		pn = length( seq(pos.matrix[9,1], pos.matrix[9,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("beta", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
		#-------------------------------------------
		# special consideration for the iGARCH model
		#-------------------------------------------
		if(modeldesc$vmodel == "iGARCH"){
			# last beta not estimated
			pars[nx+pn, 4] = 0
			nnx = paste("beta", pn, sep="")
			# do not allow the last beta to be fixed
			if(any(substr(fixed.names, 1, nchar(nnx))==nnx)) pars[(nx+pn), 2] = 0
		}
	} else{
		pnames = c(pnames, "beta")
	}
	pidx[9,2] = nx+pn

	nx = nx + pn
	pn = 1
	pidx[10,1] = nx+1

	if(pos.matrix[10,3]==1){
		pn = length( seq(pos.matrix[10,1], pos.matrix[10,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("gamma", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "gamma")
	}
	pidx[10,2] = nx+pn

	nx = nx + pn
	pn = 1
	pidx[11,1] = nx+1

	if(pos.matrix[11,3]==1){
		pn = length( seq(pos.matrix[11,1], pos.matrix[11,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("eta1", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "eta1")
	}
	pidx[11,2] = nx+pn

	nx = nx + pn
	pn = 1
	pidx[12,1] = nx+1

	if(pos.matrix[12,3]==1){
		pn = length( seq(pos.matrix[12,1], pos.matrix[12,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("eta2", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "eta2")
	}
	pidx[12,2] = nx+pn

	nx = nx + pn
	pn = 1
	pidx[13,1] = nx+1

	if(pos.matrix[13,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(substr(fixed.names, 1, 5)=="delta")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	} else{
		#-------------------------------------------
		# special consideration for the fGARCH model
		#-------------------------------------------
		if(modeldesc$vmodel == "fGARCH")
		{
			pars[nx+pn, 3] = 1
			pars[nx+pn, 1] = fmodel$fpars$delta
			pars[nx+pn, 4] = pars[nx+pn, 2] = 0
		}
	}
	pidx[13,2] = nx+pn

	pnames = c(pnames, "delta")

	nx = nx + pn
	pn = 1
	pidx[14,1] = nx+1

	if(pos.matrix[14,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(substr(fixed.names, 1, 6)=="lambda")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	} else{
		#-------------------------------------------
		# special consideration for the fGARCH model
		#-------------------------------------------
		if(modeldesc$vmodel == "fGARCH")
		{
			pars[nx+pn, 3] = 1
			pars[nx+pn, 1] = fmodel$fpars$lambda
			pars[nx+pn, 4] = pars[nx+pn, 2] = 0
		}
	}
	pidx[14,2] = nx+pn

	pnames = c(pnames, "lambda")

	nx = nx + pn
	pn = 1
	pidx[15,1] = nx+1

	if(pos.matrix[15,3]==1){
		pn = length( seq(pos.matrix[15,1], pos.matrix[15,2], by = 1) )
		for(i in 1:pn){
			pars[(nx+i), 1] = 0
			pars[(nx+i), 3] = 1
			nnx = paste("vxreg", i, sep="")
			sp = na.omit(match(fixed.names, nnx))
			if(length(sp)>0) pars[(nx+i), 2] = 1 else pars[(nx+i), 4] = 1
			pnames = c(pnames, nnx)
		}
	} else{
		pnames = c(pnames, "vxreg")

	}
	pidx[15,2] = nx+pn

	nx = nx + pn
	pn = 1
	pidx[16,1] = nx+1

	if(pos.matrix[16,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(substr(fixed.names, 1, 4)=="skew")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pidx[16,2] = nx+pn

	pnames = c(pnames, "skew")

	nx = nx + pn
	pn = 1
	pidx[17,1] = nx+1

	if(pos.matrix[17,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(substr(fixed.names, 1, 5)=="shape")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pnames = c(pnames, "shape")
	pidx[17,2] = nx+pn

	nx = nx + pn
	pn = 1
	pidx[18,1] = nx+1

	if(pos.matrix[18,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(substr(fixed.names, 1, 8)=="ghlambda")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pidx[18,2] = nx+pn
	pnames = c(pnames, "ghlambda")

	nx = nx + pn
	pn = 1
	pidx[19,1] = nx+1

	if(pos.matrix[19,3]==1){
		pars[nx+pn, 3] = 1
		pars[nx+pn, 1] = 0
		if(any(substr(fixed.names, 1, 2)=="xi")) pars[nx+pn,2] = 1 else pars[nx+pn,4] = 1
	}
	pidx[19,2] = nx+pn
	pnames = c(pnames, "xi")


	rownames(pars) = pnames

	zf = match(fixed.names, rownames(pars))
	if( length(zf)>0 ) pars[zf, 1] = unlist(fixed.pars)
	pars[,"LB"] = NA
	pars[,"UB"] = NA
	model = list(modelinc = modelinc, modeldesc = modeldesc, modeldata = modeldata, pars = pars,
			start.pars = start.pars, fixed.pars = fixed.pars, maxOrder = maxOrder,
			pos.matrix = pos.matrix, fmodel = fmodel, pidx = pidx)
	ans = new("uGARCHspec", model = model)

	return(ans)
}

setMethod(f = "ugarchspec", definition = .ugarchspec)

# extract spec from fit object
getspec = function(object)
{
	UseMethod("getspec")
}

.getspec = function(object)
{
	spec = ugarchspec(variance.model = list(model = object@model$modeldesc$vmodel, garchOrder = c(object@model$modelinc[8],object@model$modelinc[9]),
				submodel = object@model$modeldesc$vsubmodel, external.regressors = object@model$modeldata$vexdata),
		mean.model = list(armaOrder = c(object@model$modelinc[2],object@model$modelinc[3]),
				include.mean = object@model$modelinc[1],
				archm = ifelse(object@model$modelinc[5]>0,TRUE,FALSE), archpow = object@model$modelinc[5],
				arfima = object@model$modelinc[4], external.regressors = object@model$modeldata$mexdata,
				archex = object@model$modelinc[20]),
		distribution.model = object@model$modeldesc$distribution, start.pars  = object@model$start.pars,
		fixed.pars = object@model$fixed.pars)
	# should custom bounds be propagated?
	#idx = which(is.na(tmp@model$pars[,"LB"]))
	#tmp@model$pars[idx,"LB"] = object@model$pars[idx,"LB"]
	#idx = which(is.na(tmp@model$pars[,"UB"]))
	#tmp@model$pars[idx,"UB"] = object@model$pars[idx,"UB"]
	return(spec)
}

setMethod(f = "getspec", signature(object = "uGARCHfit"), definition = .getspec)


# internal function:
.model2spec = function(pars, model, type = "GARCH"){
	if(type == "GARCH"){
		ans = ugarchspec(variance.model = list(model = model$modeldesc$vmodel, garchOrder = c(model$modelinc[8],model$modelinc[9]),
						submodel = model$modeldesc$vsubmodel, external.regressors = model$modeldata$vexdata),
				mean.model = list(armaOrder = c(model$modelinc[2],model$modelinc[3]),
						include.mean = model$modelinc[1],
						archm = ifelse(model$modelinc[5]>0,TRUE,FALSE), archpow = model$modelinc[5],
						arfima = model$modelinc[4], external.regressors = model$modeldata$mexdata,
						archex = model$modelinc[20]),
				distribution.model = model$modeldesc$distribution)
		setfixed(ans)<-pars
	} else{
		ans = arfimaspec(
				mean.model = list(armaOrder = c(model$modelinc[2],model$modelinc[3]),
						include.mean = model$modelinc[1],
						arfima = model$modelinc[4], external.regressors = model$modeldata$mexdata),
				distribution.model = model$modeldesc$distribution)
		setfixed(ans)<-pars
	}
	return(ans)
}


# Set methods replace any existing fixed or starting parameters already present
setGeneric("setfixed<-", function(object, value){standardGeneric("setfixed<-")})

.setfixed = function(object, value){
	# get parameter values
	model = object@model
	ipars = model$pars
	pars = unlist(value)
	names(pars) = parnames = tolower(names(pars))
	# included parameters in model
	modelnames = rownames(ipars[which(ipars[,4]==1 | ipars[,2]==1), ,drop=FALSE])
	inc = NULL
	for(i in seq_along(parnames)){
		if(is.na(match(parnames[i], modelnames))){
			warning( (paste("Unrecognized Parameter in Fixed Values: ", parnames[i], "...Ignored", sep = "")))
		} else{
			inc = c(inc, i)
		}
	}
	fixed.pars = pars[inc]
	names(fixed.pars) = tolower(names(pars[inc]))
	# check for variance.targeting
	if(!is.na(model$modelinc[22])){
		vt = model$modelinc[22]
	} else{
		vt = as.logical(1-model$modelinc[7])
	}
	# set parameter values
	tmp = ugarchspec(variance.model = list(model = model$modeldesc$vmodel, garchOrder = c(model$modelinc[8], model$modelinc[9]),
					submodel = model$modeldesc$vsubmodel, external.regressors = model$modeldata$vexdata,
					variance.targeting = vt),
			mean.model = list(armaOrder = c(model$modelinc[2], model$modelinc[3]),
					include.mean = model$modelinc[1],
					archm = ifelse(model$modelinc[5]>0,TRUE,FALSE), archpow = model$modelinc[5],
					arfima = model$modelinc[4], external.regressors = model$modeldata$mexdata,
					archex = model$modelinc[20]),
			distribution.model = model$modeldesc$distribution, start.pars  = model$start.pars,
			fixed.pars = as.list(fixed.pars))
	# ToDo: Need to check that the parameters are not outside the bounds...
	idx = which(is.na(tmp@model$pars[,"LB"]))
	tmp@model$pars[idx,"LB"] = object@model$pars[idx,"LB"]
	idx = which(is.na(tmp@model$pars[,"UB"]))
	tmp@model$pars[idx,"UB"] = object@model$pars[idx,"UB"]
	return(tmp)
}
setReplaceMethod(f="setfixed", signature= c(object = "uGARCHspec", value = "vector"), definition = .setfixed)


setGeneric("setstart<-", function(object, value){standardGeneric("setstart<-")})

.setstart = function(object, value){
	# get parameter values
	model = object@model
	ipars = model$pars
	pars = unlist(value)
	names(pars) = parnames = tolower(names(pars))
	# included parameters in model
	modelnames = rownames(ipars[which(ipars[,4]==1 | ipars[,2]==1), , drop=FALSE])
	inc = NULL
	for(i in seq_along(parnames)){
		if(is.na(match(parnames[i], modelnames))){
			warning( (paste("Unrecognized Parameter in Start Values: ", parnames[i], "...Ignored", sep = "")))
		} else{
			inc = c(inc, i)
		}
	}
	start.pars = pars[inc]
	names(start.pars) = tolower(names(pars[inc]))
	# check for variance.targeting
	if(!is.na(model$modelinc[22])){
		vt = model$modelinc[22]
	} else{
		vt = as.logical(1-model$modelinc[7])
	}
	# set parameter values

	tmp = ugarchspec(variance.model = list(model = model$modeldesc$vmodel, garchOrder = c(model$modelinc[8], model$modelinc[9]),
					submodel = model$modeldesc$vsubmodel, external.regressors = model$modeldata$vexdata,
					variance.targeting = vt),
			mean.model = list(armaOrder = c(model$modelinc[2], model$modelinc[3]),
					include.mean = model$modelinc[1],
					archm = ifelse(model$modelinc[5]>0,TRUE,FALSE), archpow = model$modelinc[5],
					arfima = model$modelinc[4], external.regressors = model$modeldata$mexdata,
					archex = model$modelinc[20]),
			distribution.model = model$modeldesc$distribution, fixed.pars  = model$fixed.pars,
			start.pars = as.list(start.pars))
	# ToDo: Need to check that the parameters are not outside the bounds...
	idx = which(is.na(tmp@model$pars[,"LB"]))
	tmp@model$pars[idx,"LB"] = object@model$pars[idx,"LB"]
	idx = which(is.na(tmp@model$pars[,"UB"]))
	tmp@model$pars[idx,"UB"] = object@model$pars[idx,"UB"]
	return(tmp)
}

setReplaceMethod(f="setstart", signature= c(object = "uGARCHspec", value = "vector"), definition = .setstart)


.checkallfixed = function( spec ){
	# check that a given spec with fixed parameters
	model = spec@model
	pars = model$pars
	pnames = rownames(pars)
	estpars = pnames[as.logical(pars[,2] * pars[,3] + pars[,3] * pars[,4])]
	return( estpars )

}

setGeneric("setbounds<-", function(object, value){standardGeneric("setbounds<-")})

# Set the lower and upper bounds
# value is a list with names parameters taking 2 values (lower and upper)
# e.g. value  = list(alpha1 = c(0, 0.1), beta1 = c(0.9, 0.99))
.setbounds = function(object, value){
	model = object@model
	ipars = model$pars
	parnames = tolower(names(value))
	# included parameters in model
	modelnames = rownames(ipars[which(ipars[,4] == 1), ])
	sp = na.omit(match(parnames, modelnames))
	if(length(sp)>0){
		for(i in 1:length(sp)){
			#if(length(value[[modelnames[sp[i]]]])!=2)
			ipars[modelnames[sp[i]], 5] = as.numeric(value[[modelnames[sp[i]]]][1])
			ipars[modelnames[sp[i]], 6] = as.numeric(value[[modelnames[sp[i]]]][2])
		}
	}
	object@model$pars = ipars
	return(object)
}
setReplaceMethod(f="setbounds", signature= c(object = "uGARCHspec", value = "vector"), definition = .setbounds)

#----------------------------------------------------------------------------------
# univariate model dispatch methods
#----------------------------------------------------------------------------------

.ugarchfit = function(spec, data, out.sample = 0, solver = "solnp", solver.control = list(),
		fit.control = list(stationarity = 1, fixed.se = 0, scale = 0, rec.init = 'all'),
		numderiv.control = list(grad.eps=1e-4, grad.d=0.0001, grad.zero.tol=sqrt(.Machine$double.eps/7e-7),
				hess.eps=1e-4, hess.d=0.1, hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2),...)
{
	default.numd = list(grad.eps=1e-4, grad.d=0.0001, grad.zero.tol=sqrt(.Machine$double.eps/7e-7),
			hess.eps=1e-4, hess.d=0.1, hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2)
	idx1 = na.omit(match(names(numderiv.control), names(default.numd)))
	idx2 = na.omit(match(names(default.numd), names(numderiv.control)))
	if(length(idx1)>0){
		default.numd[idx1] = numderiv.control[idx2]
	}

	return( switch(spec@model$modeldesc$vmodel,
					sGARCH = .sgarchfit(spec = spec, data = data, out.sample = out.sample, solver = solver,
							solver.control = solver.control, fit.control = fit.control, numderiv.control = default.numd),
					iGARCH = .igarchfit(spec = spec, data = data, out.sample = out.sample, solver = solver,
							solver.control = solver.control, fit.control = fit.control, numderiv.control = default.numd),
					eGARCH = .egarchfit(spec = spec, data = data, out.sample = out.sample, solver = solver,
							solver.control = solver.control, fit.control = fit.control, numderiv.control = default.numd),
					gjrGARCH = .gjrgarchfit(spec = spec, data = data, out.sample = out.sample, solver = solver,
							solver.control = solver.control, fit.control = fit.control, numderiv.control = default.numd),
					apARCH = .aparchfit(spec = spec, data = data, out.sample = out.sample, solver = solver,
							solver.control = solver.control, fit.control = fit.control, numderiv.control = default.numd),
					fGARCH = .fgarchfit(spec = spec, data = data, out.sample = out.sample, solver = solver,
							solver.control = solver.control, fit.control = fit.control, numderiv.control = default.numd),
					csGARCH = .csgarchfit(spec = spec, data = data, out.sample = out.sample, solver = solver,
							solver.control = solver.control, fit.control = fit.control, numderiv.control = default.numd),
					mcsGARCH = .mcsgarchfit(spec = spec, data = data, out.sample = out.sample, solver = solver,
							solver.control = solver.control, fit.control = fit.control, numderiv.control = default.numd, ...),
					realGARCH = .realgarchfit(spec = spec, data = data, out.sample = out.sample, solver = solver,
							solver.control = solver.control, fit.control = fit.control, numderiv.control = default.numd, ...)) )
}

.ugarchfilter = function(spec, data, out.sample = 0, n.old = NULL, rec.init = 'all', ...)
{
	return( switch(spec@model$modeldesc$vmodel,
					sGARCH = .sgarchfilter(spec = spec, data = data, out.sample = out.sample, n.old = n.old, rec.init = rec.init),
					iGARCH = .igarchfilter(spec = spec, data = data, out.sample = out.sample, n.old = n.old, rec.init = rec.init),
					eGARCH = .egarchfilter(spec = spec, data = data, out.sample = out.sample, n.old = n.old, rec.init = rec.init),
					gjrGARCH = .gjrgarchfilter(spec = spec, data = data, out.sample = out.sample, n.old = n.old, rec.init = rec.init),
					apARCH = .aparchfilter(spec = spec, data = data, out.sample = out.sample, n.old = n.old, rec.init = rec.init),
					fGARCH = .fgarchfilter(spec = spec, data = data, out.sample = out.sample, n.old = n.old, rec.init = rec.init),
					csGARCH = .csgarchfilter(spec = spec, data = data, out.sample = out.sample, n.old = n.old, rec.init = rec.init),
					mcsGARCH = .mcsgarchfilter(spec = spec, data = data, out.sample = out.sample, n.old = n.old, rec.init = rec.init, ...),
					realGARCH = .realgarchfilter(spec = spec, data = data, out.sample = out.sample, n.old = n.old, rec.init = rec.init, ...)) )
}

.ugarchforecast1 = function(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0,
		external.forecasts = list(mregfor = NULL, vregfor = NULL), ...)
{
	return( switch(fitORspec@model$modeldesc$vmodel,
					sGARCH = .sgarchforecast(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll,
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					iGARCH = .igarchforecast(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll,
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					eGARCH = .egarchforecast(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll,
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					gjrGARCH = .gjrgarchforecast(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll,
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					apARCH = .aparchforecast(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll,
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					fGARCH = .fgarchforecast(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll,
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					csGARCH = .csgarchforecast(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll,
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					mcsGARCH = .mcsgarchforecast(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll,
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					realGARCH = .realgarchforecast(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll,
							out.sample = out.sample, external.forecasts = external.forecasts, ...)) )
}

.ugarchforecast2 = function(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0,
		external.forecasts = list(mregfor = NULL, vregfor = NULL), ...)
{
	return( switch(fitORspec@model$modeldesc$vmodel,
					sGARCH = .sgarchforecast2(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll,
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					iGARCH = .igarchforecast2(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll,
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					eGARCH = .egarchforecast2(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll,
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					gjrGARCH = .gjrgarchforecast2(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll,
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					apARCH = .aparchforecast2(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll,
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					fGARCH = .fgarchforecast2(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll,
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					csGARCH = .csgarchforecast2(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll,
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					mcsGARCH = .mcsgarchforecast2(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll,
							out.sample = out.sample, external.forecasts = external.forecasts, ...),
					realGARCH = .realgarchforecast2(fitORspec = fitORspec, data = data, n.ahead = n.ahead, n.roll = n.roll,
							out.sample = out.sample, external.forecasts = external.forecasts, ...)) )
}

.ugarchsim = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, startMethod = c("unconditional","sample"),
		presigma = NA, prereturns = NA, preresiduals = NA, rseed = NA,  custom.dist = list(name = NA, distfit = NA),
		mexsimdata = NULL, vexsimdata = NULL, ...)
{
	return( switch(fit@model$modeldesc$vmodel,
					sGARCH = .sgarchsim(fit = fit, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
							startMethod = startMethod, presigma = presigma, prereturns = prereturns,
							preresiduals = preresiduals, rseed = rseed,  custom.dist = custom.dist,
							mexsimdata = mexsimdata, vexsimdata = vexsimdata, ...),
					iGARCH = .igarchsim(fit = fit, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
							startMethod = startMethod, presigma = presigma, prereturns = prereturns,
							preresiduals = preresiduals, rseed = rseed,  custom.dist = custom.dist,
							mexsimdata = mexsimdata, vexsimdata = vexsimdata, ...),
					eGARCH = .egarchsim(fit = fit, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
							startMethod = startMethod, presigma = presigma, prereturns = prereturns,
							preresiduals = preresiduals, rseed = rseed,  custom.dist = custom.dist,
							mexsimdata = mexsimdata, vexsimdata = vexsimdata, ...),
					gjrGARCH = .gjrgarchsim(fit = fit, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
							startMethod = startMethod, presigma = presigma, prereturns = prereturns,
							preresiduals = preresiduals, rseed = rseed,  custom.dist = custom.dist,
							mexsimdata = mexsimdata, vexsimdata = vexsimdata, ...),
					apARCH = .aparchsim(fit = fit, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
							startMethod = startMethod, presigma = presigma, prereturns = prereturns,
							preresiduals = preresiduals, rseed = rseed,  custom.dist = custom.dist,
							mexsimdata = mexsimdata, vexsimdata = vexsimdata, ...),
					fGARCH = .fgarchsim(fit = fit, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
							startMethod = startMethod, presigma = presigma, prereturns = prereturns,
							preresiduals = preresiduals, rseed = rseed,  custom.dist = custom.dist,
							mexsimdata = mexsimdata, vexsimdata = vexsimdata, ...),
					csGARCH = .csgarchsim(fit = fit, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
							startMethod = startMethod, presigma = presigma, prereturns = prereturns,
							preresiduals = preresiduals, rseed = rseed,  custom.dist = custom.dist,
							mexsimdata = mexsimdata, vexsimdata = vexsimdata, ...),
					mcsGARCH = .mcsgarchsim(fit = fit, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
							startMethod = startMethod, presigma = presigma, prereturns = prereturns,
							preresiduals = preresiduals, rseed = rseed,  custom.dist = custom.dist,
							mexsimdata = mexsimdata, vexsimdata = vexsimdata, ...),
				   realGARCH = .realgarchsim(fit = fit, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
							startMethod = startMethod, presigma = presigma, prereturns = prereturns,
							preresiduals = preresiduals, rseed = rseed,  custom.dist = custom.dist,
							mexsimdata = mexsimdata, vexsimdata = vexsimdata, ...)) )
}

.ugarchpath = function(spec, n.sim = 1000, n.start = 0, m.sim = 1, presigma = NA, prereturns = NA, preresiduals = NA,
		rseed = NA, custom.dist = list(name = NA, distfit = NA), mexsimdata = NULL,  vexsimdata = NULL, ...)
{
	return( switch(spec@model$modeldesc$vmodel,
					sGARCH = .sgarchpath(spec = spec, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
							presigma = presigma, prereturns = prereturns, preresiduals = preresiduals,
							rseed = rseed,  custom.dist = custom.dist, mexsimdata = mexsimdata,
							vexsimdata = vexsimdata, ...),
					iGARCH = .igarchpath(spec = spec, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
							presigma = presigma, prereturns = prereturns, preresiduals = preresiduals,
							rseed = rseed,  custom.dist = custom.dist, mexsimdata = mexsimdata,
							vexsimdata = vexsimdata, ...),
					eGARCH = .egarchpath(spec = spec, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
							presigma = presigma, prereturns = prereturns, preresiduals = preresiduals,
							rseed = rseed,  custom.dist = custom.dist, mexsimdata = mexsimdata,
							vexsimdata = vexsimdata, ...),
					gjrGARCH = .gjrgarchpath(spec = spec, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
							presigma = presigma, prereturns = prereturns, preresiduals = preresiduals,
							rseed = rseed,  custom.dist = custom.dist, mexsimdata = mexsimdata,
							vexsimdata = vexsimdata, ...),
					apARCH = .aparchpath(spec = spec, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
							presigma = presigma, prereturns = prereturns, preresiduals = preresiduals,
							rseed = rseed,  custom.dist = custom.dist, mexsimdata = mexsimdata,
							vexsimdata = vexsimdata, ...),
					fGARCH = .fgarchpath(spec = spec, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
							presigma = presigma, prereturns = prereturns, preresiduals = preresiduals,
							rseed = rseed,  custom.dist = custom.dist, mexsimdata = mexsimdata,
							vexsimdata = vexsimdata, ...),
					csGARCH = .csgarchpath(spec = spec, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
							presigma = presigma, prereturns = prereturns, preresiduals = preresiduals,
							rseed = rseed,  custom.dist = custom.dist, mexsimdata = mexsimdata,
							vexsimdata = vexsimdata, ...),
				realGARCH = .realgarchpath(spec = spec, n.sim = n.sim, n.start = n.start, m.sim = m.sim,
							presigma = presigma, prereturns = prereturns, preresiduals = preresiduals,
							rseed = rseed,  custom.dist = custom.dist, mexsimdata = mexsimdata,
							vexsimdata = vexsimdata, ...)) )
}
#----------------------------------------------------------------------------------
# univariate filter method
#----------------------------------------------------------------------------------
ugarchfilter = function(spec, data, out.sample = 0, n.old = NULL, rec.init = 'all', ...)
{
	UseMethod("ugarchfilter")
}

setMethod("ugarchfilter", signature(spec = "uGARCHspec"), .ugarchfilter)
#----------------------------------------------------------------------------------
# univariate fit method
#----------------------------------------------------------------------------------
ugarchfit = function(spec, data, out.sample = 0, solver = "solnp", solver.control = list(),
		fit.control = list(stationarity = 1, fixed.se = 0, scale = 0, rec.init = 'all'),
		numderiv.control = list(grad.eps=1e-4, grad.d=0.0001, grad.zero.tol=sqrt(.Machine$double.eps/7e-7),
				hess.eps=1e-4, hess.d=0.1, hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2), ...)
{
	UseMethod("ugarchfit")
}

setMethod("ugarchfit", signature(spec = "uGARCHspec"), .ugarchfit)
#----------------------------------------------------------------------------------
# univariate forecast method
#----------------------------------------------------------------------------------
ugarchforecast = function(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0,
		external.forecasts = list(mregfor = NULL, vregfor = NULL), ...)
{
	UseMethod("ugarchforecast")
}

setMethod("ugarchforecast", signature(fitORspec = "uGARCHfit"), .ugarchforecast1)

#setMethod("ugarchforecast", signature(fitORspec = "anstGARCHfit"), .anstgarchforecast)
#alternative dispath method:
# we use the fitORspec rather than a method with fit, spec and data with missing
# methods as this requires implicit declaration of arguments

setMethod("ugarchforecast", signature(fitORspec = "uGARCHspec"), .ugarchforecast2)

#----------------------------------------------------------------------------------
# univariate simulation method
#----------------------------------------------------------------------------------

ugarchsim = function(fit, n.sim = 1000, n.start = 0, m.sim = 1,
		startMethod = c("unconditional","sample"),
		presigma = NA, prereturns = NA, preresiduals = NA, rseed = NA,
		custom.dist = list(name = NA, distfit = NA), mexsimdata = NULL,
		vexsimdata = NULL, ...)
{
	UseMethod("ugarchsim")
}


setMethod("ugarchsim", signature(fit = "uGARCHfit"), .ugarchsim)
#----------------------------------------------------------------------------------
# univariate path simulation method
#----------------------------------------------------------------------------------

ugarchpath = function(spec, n.sim = 1000, n.start = 0, m.sim = 1,
		presigma = NA, prereturns = NA, preresiduals = NA, rseed = NA,
		custom.dist = list(name = NA, distfit = NA), mexsimdata = NULL,
		vexsimdata = NULL, ...)
{
	UseMethod("ugarchpath")
}


setMethod("ugarchpath", signature(spec = "uGARCHspec"), .ugarchpath)
#----------------------------------------------------------------------------------
# resume method
#----------------------------------------------------------------------------------
resume = function(object, ...)
{
	UseMethod("resume")
}

# Need a special purpose resume function for the mcsGARCH method
.resume = function(object, ...)
{
	if(object@model$spec@model$modeldesc$vmodel=="mcsGARCH"){
		ans = .resumeroll.mcs(object, ...)
	} else if(object@model$spec@model$modeldesc$vmodel=="realGARCH"){
		ans = .resumeroll.real(object, ...)
	}else{
		ans = .resumeroll1(object, ...)
	}
	return( ans )
}
setMethod("resume", signature(object = "uGARCHroll"),  definition = .resume)

#----------------------------------------------------------------------------------
# univariate garch roll
#----------------------------------------------------------------------------------
# methods to recursively predict/filter/compare with refitting at every N points.
ugarchroll = function(spec, data, n.ahead = 1, forecast.length = 500,
		n.start = NULL, refit.every = 25, refit.window = c("recursive", "moving"),
		window.size = NULL, solver = "hybrid", fit.control = list(), solver.control = list(),
		calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), cluster = NULL,
		keep.coef = TRUE, ...)
{
	UseMethod("ugarchroll")
}

# Need a special purpose roll function for the mcsGARCH method

.ugarchroll = function(spec, data, n.ahead = 1, forecast.length = 500,
		n.start = NULL, refit.every = 25, refit.window = c("recursive", "moving"),
		window.size = NULL, solver = "hybrid", fit.control = list(), solver.control = list(),
		calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), cluster = NULL,
		keep.coef = TRUE, ...)
{
	if(spec@model$modeldesc$vmodel == "mcsGARCH"){
		if(is.null(list(...)$DailyVar)) stop("\nmcsGARCH models requires DailyVar forecast")
		DailyVar = list(...)$DailyVar
		ans = .rollfdensity.mcs(spec = spec, data = data, n.ahead = n.ahead,
				forecast.length = forecast.length, n.start = n.start,
				refit.every = refit.every, refit.window = refit.window[1],
				window.size = window.size, solver = solver,
				fit.control = fit.control, solver.control = solver.control,
				calculate.VaR = calculate.VaR, VaR.alpha = VaR.alpha,
				cluster = cluster, keep.coef = keep.coef, DailyVar = DailyVar)
	} else if(spec@model$modeldesc$vmodel == "realGARCH"){
		if(is.null(list(...)$realizedVol)) stop("\nrealGARCH model requires realizedVol")
		realizedVol = list(...)$realizedVol
		ans = .rollfdensity.real(spec = spec, data = data, n.ahead = n.ahead,
				forecast.length = forecast.length, n.start = n.start,
				refit.every = refit.every, refit.window = refit.window[1],
				window.size = window.size, solver = solver,
				fit.control = fit.control, solver.control = solver.control,
				calculate.VaR = calculate.VaR, VaR.alpha = VaR.alpha,
				cluster = cluster, keep.coef = keep.coef, realizedVol = realizedVol)
	} else{
		ans = .rollfdensity(spec = spec, data = data, n.ahead = n.ahead,
				forecast.length = forecast.length, n.start = n.start,
				refit.every = refit.every, refit.window = refit.window[1],
				window.size = window.size, solver = solver,
				fit.control = fit.control, solver.control = solver.control,
				calculate.VaR = calculate.VaR, VaR.alpha = VaR.alpha,
				cluster = cluster, keep.coef = keep.coef, ...)
	}
	return( ans )
}
setMethod("ugarchroll", signature(spec = "uGARCHspec"),  definition = .ugarchroll)

#----------------------------------------------------------------------------------
# univariate garch parameter distribution
#----------------------------------------------------------------------------------
ugarchdistribution = function(fitORspec, n.sim = 2000, n.start = 1,
		m.sim = 100, recursive = FALSE, recursive.length = 6000, recursive.window = 1000,
		presigma = NA, prereturns = NA, preresiduals = NA, rseed = NA,
		custom.dist = list(name = NA, distfit = NA), mexsimdata = NULL,
		vexsimdata = NULL, fit.control = list(), solver = "solnp",
		solver.control = list(), cluster = NULL, ...)
{
	UseMethod("ugarchdistribution")
}
setMethod("ugarchdistribution", signature(fitORspec = "uGARCHfit"), .ugarchdistribution)
setMethod("ugarchdistribution", signature(fitORspec = "uGARCHspec"), .ugarchdistribution)


#----------------------------------------------------------------------------------
# univariate garch bootstrap based forecast distribution
#----------------------------------------------------------------------------------
ugarchboot = function(fitORspec, data = NULL, method = c("Partial", "Full"),
		sampling = c("raw", "kernel", "spd"),
		spd.options = list(upper = 0.9, lower = 0.1, type = "pwm", kernel = "normal"),
		n.ahead = 10, n.bootfit = 100, n.bootpred = 500, out.sample = 0, rseed = NA, solver = "solnp",
		solver.control = list(), fit.control = list(), external.forecasts =  list(mregfor = NULL,
				vregfor = NULL), mexsimdata = NULL, vexsimdata = NULL, cluster = NULL,
		verbose = FALSE)
{
	UseMethod("ugarchboot")
}

setMethod("ugarchboot", signature(fitORspec = "uGARCHfit"), .ugarchbootfit)
setMethod("ugarchboot", signature(fitORspec = "uGARCHspec"), .ugarchbootspec)

#----------------------------------------------------------------------------------
# univariate plot method / seperate for fit,sim and forecast
#----------------------------------------------------------------------------------
setMethod(f = "plot", signature(x = "uGARCHfit", y = "missing"), .plotgarchfit)

setMethod(f = "plot", signature(x = "uGARCHfilter", y = "missing"), .plotgarchfilter)

setMethod(f = "plot", signature(x = "uGARCHsim", y = "missing"), .plotgarchsim)

setMethod(f = "plot", signature(x = "uGARCHforecast", y = "missing"), .plotgarchforecast)

setMethod(f = "plot", signature(x = "uGARCHpath", y = "missing"), .plotgarchpath)

setMethod(f = "plot", signature(x = "uGARCHroll", y = "missing"), .plotgarchroll)

setMethod(f = "plot", signature(x = "uGARCHdistribution", y = "missing"), .plotgarchdist)

setMethod(f = "plot", signature(x = "uGARCHboot", y = "missing"), .plotgarchboot)

#----------------------------------------------------------------------------------
# univariate show method / seperate for fit,sim and forecast
#----------------------------------------------------------------------------------

# spec show
setMethod("show",
		signature(object = "uGARCHspec"),
		function(object){
			model = object@model
			vmodel = model$modeldesc$vmodel
			modelinc = model$modelinc
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\n*       GARCH Model Spec          *", sep = ""))
			cat(paste("\n*---------------------------------*", sep = ""))
			cat("\n\nConditional Variance Dynamics \t")
			cat(paste("\n------------------------------------\n",sep=""))
			cat(paste("GARCH Model\t\t: ", vmodel, "(",modelinc[8],",",modelinc[9],")\n", sep = ""))
			if(vmodel == "fGARCH"){
				cat(paste("fGARCH Sub-Model\t: ", model$modeldesc$vsubmodel, "\n", sep = ""))
			}
			cat("Variance Targeting\t:", as.logical(1-modelinc[7]), "\n")
			if(modelinc[15]>0) cat(paste("Exogenous Regressor Dimension: ", modelinc[15], "\n",sep = ""))
			cat("\nConditional Mean Dynamics")
			cat(paste("\n------------------------------------\n",sep=""))
			cat("Mean Model\t\t: ARFIMA(", modelinc[2],",", ifelse(modelinc[4]>0, "d", 0),",",modelinc[3],")\n", sep = "")
			cat("Include Mean\t\t:", as.logical(modelinc[1]),"\n")
			cat("GARCH-in-Mean\t\t:", as.logical(modelinc[5]),"\n")
			if(modelinc[6]>0) cat(paste("Exogenous Regressor Dimension: ", modelinc[6],"\n",sep=""))
			cat("\nConditional Distribution")
			cat(paste("\n------------------------------------\n",sep=""))
			cat("Distribution\t: ", model$modeldesc$distribution,"\n")
			cat("Includes Skew\t: ", as.logical(modelinc[16]),"\n")
			cat("Includes Shape\t: ", as.logical(modelinc[17]),"\n")
			cat("Includes Lambda\t: ", as.logical(modelinc[18]),"\n\n")
			invisible(object)
		})

# fit show
setMethod("show",
		signature(object = "uGARCHfit"),
		function(object){
			vmodel = object@model$modeldesc$vmodel
			model = object@model
			modelinc = object@model$modelinc
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\n*          GARCH Model Fit        *", sep = ""))
			cat(paste("\n*---------------------------------*", sep = ""))
			cat("\n\nConditional Variance Dynamics \t")
			cat(paste("\n-----------------------------------", sep = ""))
			cat(paste("\nGARCH Model\t: ", vmodel,"(", modelinc[8], ",", modelinc[9], ")\n", sep=""))
			if(vmodel == "fGARCH"){
				cat(paste("fGARCH Sub-Model\t: ", model$modeldesc$vsubmodel, "\n", sep = ""))
			}
			cat("Mean Model\t: ARFIMA(", modelinc[2],",",ifelse(modelinc[4]>0, "d", 0),",",modelinc[3],")\n", sep = "")
			cat("Distribution\t:", model$modeldesc$distribution,"\n")
			if(object@fit$convergence == 0){
				cat("\nOptimal Parameters")
				cat(paste("\n------------------------------------\n",sep=""))
				print(round(object@fit$matcoef,6), digits = 5)
				cat("\nRobust Standard Errors:\n")
				print(round(object@fit$robust.matcoef,6), digits = 5)
				if(!is.null(object@fit$hessian.message)){
					cat(paste("\n", object@fit$hessian.message))
				}
				cat("\nLogLikelihood :", object@fit$LLH, "\n")
				stdresid = object@fit$residuals/object@fit$sigma
				itestm = infocriteria(object)
				cat("\nInformation Criteria")
				cat(paste("\n------------------------------------\n",sep=""))
				print(itestm,digits=5)
				cat("\nWeighted Ljung-Box Test on Standardized Residuals")
				cat(paste("\n------------------------------------\n",sep=""))
				tmp1 = .weightedBoxTest(stdresid, p = 1, df = sum(modelinc[2:3]))
				print(tmp1, digits = 4)
				cat(paste("d.o.f=", sum(modelinc[2:3]), sep=""))
				cat("\nH0 : No serial correlation\n")
				cat("\nWeighted Ljung-Box Test on Standardized Squared Residuals")
				cat(paste("\n------------------------------------\n",sep=""))
				tmp2 = .weightedBoxTest(stdresid, p = 2, df = sum(modelinc[8:9]))
				print(tmp2, digits = 4)
				cat(paste("d.o.f=", sum(modelinc[8:9]), sep=""))
				cat("\n\nWeighted ARCH LM Tests")
				cat(paste("\n------------------------------------\n",sep=""))
				gdf = sum(modelinc[8:9])
				L2 = .weightedarchlmtest(residuals(object), sigma(object), lags  = gdf+1, fitdf=gdf)
				L5 = .weightedarchlmtest(residuals(object), sigma(object), lags  = gdf+3, fitdf=gdf)
				L10 = .weightedarchlmtest(residuals(object), sigma(object), lags = gdf+5, fitdf=gdf)
				alm = matrix(0, ncol = 4, nrow = 3)
				alm[1,1:4] = as.numeric(c(L2$statistic, L2$parameter, L2$p.value))
				alm[2,1:4] = as.numeric(c(L5$statistic, L5$parameter, L5$p.value))
				alm[3,1:4] = as.numeric(c(L10$statistic, L10$parameter, L10$p.value))
				colnames(alm) = c("Statistic", "Shape", "Scale", "P-Value")
				rownames(alm) = c(paste("ARCH Lag[",gdf+1,"]",sep=""), paste("ARCH Lag[",gdf+3,"]",sep=""), paste("ARCH Lag[",gdf+5,"]",sep=""))
				print(alm,digits = 4)
				nyb = .nyblomTest(object)
				if(is.character(nyb$JointCritical)){
					colnames(nyb$IndividualStat)<-""
					cat("\nNyblom stability test")
					cat(paste("\n------------------------------------\n",sep=""))
					cat("Joint Statistic: ", "no.parameters>20 (not available)")
					cat("\nIndividual Statistics:")
					print(nyb$IndividualStat, digits = 4)
					cat("\nAsymptotic Critical Values (10% 5% 1%)")
					cat("\nIndividual Statistic:\t", round(nyb$IndividualCritical, 2))
					cat("\n\n")
				} else{
					colnames(nyb$IndividualStat)<-""
					cat("\nNyblom stability test")
					cat(paste("\n------------------------------------\n",sep=""))
					cat("Joint Statistic: ", round(nyb$JointStat,4))
					cat("\nIndividual Statistics:")
					print(nyb$IndividualStat, digits = 4)
					cat("\nAsymptotic Critical Values (10% 5% 1%)")
					cat("\nJoint Statistic:     \t", round(nyb$JointCritical, 3))
					cat("\nIndividual Statistic:\t", round(nyb$IndividualCritical, 2))
					cat("\n\n")
				}
				cat("Sign Bias Test")
				cat(paste("\n------------------------------------\n",sep=""))
				sgtest = signbias(object)
				print(sgtest, digits = 4)
				cat("\n")
				cat("\nAdjusted Pearson Goodness-of-Fit Test:")
				cat(paste("\n------------------------------------\n",sep=""))
				gofm = gof(object,c(20, 30, 40, 50))
				print(gofm, digits = 4)
				cat("\n")
				cat("\nElapsed time :", object@fit$timer,"\n\n")
			} else{
				cat("\nConvergence Problem:")
				cat("\nSolver Message:", object@fit$message,"\n\n")

			}
			invisible(object)
		})
# filter show
setMethod("show",
		signature(object = "uGARCHfilter"),
		function(object){
			vmodel = object@model$modeldesc$vmodel
			model = object@model
			modelinc = object@model$modelinc
			cat(paste("\n*------------------------------------*", sep = ""))
			cat(paste("\n*          GARCH Model Filter        *", sep = ""))
			cat(paste("\n*------------------------------------*", sep = ""))
			cat("\n\nConditional Variance Dynamics \t")
			cat(paste("\n--------------------------------------", sep = ""))
			cat(paste("\nGARCH Model\t: ", vmodel,"(", modelinc[8], ",", modelinc[9], ")\n", sep=""))
			if(vmodel == "fGARCH"){
				cat(paste("fGARCH Sub-Model\t: ", model$modeldesc$vsubmodel, "\n", sep = ""))
			}
			cat("Mean Model\t: ARFIMA(", modelinc[2],",",ifelse(modelinc[4]>0, "d", 0),",",modelinc[3],")\n", sep = "")
			cat("Distribution\t:", model$modeldesc$distribution,"\n")
			cat("\nFilter Parameters")
			cat(paste("\n---------------------------------------\n",sep=""))
			print(matrix(coef(object), ncol=1, dimnames = list(names(coef(object)), "")), digits = 5)
			cat("\nLogLikelihood :", object@filter$LLH, "\n")
			stdresid = object@filter$residuals/object@filter$sigma
			itestm = infocriteria(object)
			cat("\nInformation Criteria")
			cat(paste("\n---------------------------------------\n",sep=""))
			print(itestm,digits=5)
			cat("\nWeighted Ljung-Box Test on Standardized Residuals")
			cat(paste("\n---------------------------------------\n",sep=""))
			tmp1 = .weightedBoxTest(stdresid, p = 1, df = sum(modelinc[2:3]))
			print(tmp1, digits = 4)
			cat(paste("d.o.f=", sum(modelinc[2:3]), sep=""))
			cat("\nH0 : No serial correlation\n")
			cat("\nWeighted Ljung-Box Test on Standardized Squared Residuals")
			cat(paste("\n---------------------------------------\n",sep=""))
			tmp2 = .weightedBoxTest(stdresid, p = 2, df = sum(modelinc[8:9]))
			print(tmp2, digits = 4)
			cat(paste("d.o.f=", sum(modelinc[8:9]), sep=""))
			cat("\n\nWeighted ARCH LM Tests")
			cat(paste("\n---------------------------------------\n",sep=""))
			gdf = sum(modelinc[8:9])
			L2 = .weightedarchlmtest(residuals(object), sigma(object), lags  = gdf+1, fitdf=gdf)
			L5 = .weightedarchlmtest(residuals(object), sigma(object), lags  = gdf+3, fitdf=gdf)
			L10 = .weightedarchlmtest(residuals(object), sigma(object), lags = gdf+5, fitdf=gdf)
			alm = matrix(0,ncol = 4,nrow = 3)
			alm[1,1:4] = as.numeric(c(L2$statistic, L2$parameter, L2$p.value))
			alm[2,1:4] = as.numeric(c(L5$statistic, L5$parameter, L5$p.value))
			alm[3,1:4] = as.numeric(c(L10$statistic, L10$parameter, L10$p.value))
			colnames(alm) = c("Statistic", "Shape", "Scale", "P-Value")
			rownames(alm) = c(paste("ARCH Lag[",gdf+1,"]",sep=""), paste("ARCH Lag[",gdf+3,"]",sep=""), paste("ARCH Lag[",gdf+5,"]",sep=""))
			print(alm,digits = 4)
			cat("\n\n")
			cat("Sign Bias Test")
			cat(paste("\n---------------------------------------\n",sep=""))
			sgtest = signbias(object)
			print(sgtest, digits = 4)
			cat("\n")
			cat("\nAdjusted Pearson Goodness-of-Fit Test:")
			cat(paste("\n---------------------------------------\n",sep=""))
			gofm = gof(object,c(20, 30, 40, 50))
			print(gofm, digits = 4)
			cat("\n")
			invisible(object)
		})
# sim show
setMethod("show",
			signature(object = "uGARCHsim"),
			function(object){
			vmodel = object@model$modeldesc$vmodel
			model = object@model
			cat(paste("\n*------------------------------------*", sep = ""))
			cat(paste("\n*       GARCH Model Simulation       *", sep = ""))
			cat(paste("\n*------------------------------------*", sep = ""))
			cat(paste("\nModel : ", vmodel,sep=""))
			if(vmodel == "fGARCH"){
				cat(paste("fGARCH Sub-Model\t: ", model$modeldesc$vsubmodel, "\n", sep = ""))
			}
			sim = object@simulation
			sigma = sim$sigmaSim
			series = sim$seriesSim
			resids = sim$residSim
			m = dim(sigma)[2]
			N = dim(sigma)[1]
			cat(paste("\nHorizon: ",N))
			cat(paste("\nSimulations: ",m,"\n",sep=""))
			sd1 = apply(sigma^2, 2, FUN=function(x) mean(x))
			sd2 = apply(sigma^2, 2, FUN=function(x) range(x))
			rx1 = apply(series, 2, FUN=function(x) mean(x))
			rx2 = apply(series, 2, FUN=function(x) range(x))
			T = object@model$modeldata$T
			actual = c(0,mean(object@model$modeldata$sigma^2), min(object@model$modeldata$sigma^2),
					max(object@model$modeldata$sigma^2), mean(object@model$modeldata$data[1:T]),
					min(object@model$modeldata$data[1:T]), max(object@model$modeldata$data[1:T]))
			xspec = .model2spec(as.list(object@model$pars[object@model$pars[,3]==1,1]), object@model, type = "GARCH")
			setfixed(xspec)<-as.list(object@model$pars[which(object@model$pars[,3]==1),1])
			uv = uncvariance(xspec)
			um = uncmean(xspec)
			uncond = c(0, uv, NA, NA, um, NA, NA)
			dd = data.frame(Seed = object@seed, Sigma2.Mean = sd1, Sigma2.Min = sd2[1,],
					Sigma2.Max = sd2[2,], Series.Mean = rx1, Series.Min = rx2[1,],
					Series.Max = rx2[2,])
			meansim = apply(dd, 2, FUN = function(x) mean(x))
			meansim[1] = 0
			dd = rbind(dd, meansim, actual, uncond)
			rownames(dd) = c(paste("sim", 1:m, sep = ""), "Mean(All)", "Actual", "Unconditional")
			print(dd,digits = 3)
			cat("\n\n")
			})

# forecast show
setMethod("show",
		signature(object = "uGARCHforecast"),
		function(object){
			vmodel = object@model$modeldesc$vmodel
			model = object@model
			cat(paste("\n*------------------------------------*", sep = ""))
			cat(paste("\n*       GARCH Model Forecast         *", sep = ""))
			cat(paste("\n*------------------------------------*", sep = ""))
			cat(paste("\nModel: ", vmodel, sep = ""))
			if(vmodel == "fGARCH"){
				cat(paste("\nfGARCH Sub-Model: ", model$modeldesc$vsubmodel, "\n", sep = ""))
			}
			n.ahead = object@forecast$n.ahead
			cat(paste("\nHorizon: ", n.ahead, sep = ""))
			cat(paste("\nRoll Steps: ",object@forecast$n.roll, sep = ""))
			n.start = object@forecast$n.start
			if(n.start>0) infor = ifelse(n.ahead>n.start, n.start, n.ahead) else infor = 0
			cat(paste("\nOut of Sample: ", infor, "\n", sep = ""))
			cat(paste("\n0-roll forecast [T0=", as.character(object@model$modeldata$index[object@model$modeldata$T]), "]:\n", sep=""))
			if(vmodel=="csGARCH"){
				zz = cbind(object@forecast$seriesFor[,1], object@forecast$sigmaFor[,1], object@forecast$qFor[,1])
				colnames(zz) = c("Series", "Sigma[Transitory]", "Sigma[Permanent]")
				rownames(zz) = paste("T+",1:NROW(object@forecast$seriesFor), sep="")
			} else if(vmodel=="mcsGARCH"){
				zz = cbind(object@forecast$seriesFor[,1], object@forecast$sigmaFor[,1], object@forecast$qFor[,1])
				colnames(zz) = c("Series", "Sigma[Total]", "Sigma[Stochastic]")
				rownames(zz) = paste("T+",1:NROW(object@forecast$seriesFor), sep="")
			} else{
				zz = cbind(object@forecast$seriesFor[,1], object@forecast$sigmaFor[,1])
				colnames(zz) = c("Series", "Sigma")
				rownames(zz) = paste("T+",1:NROW(object@forecast$seriesFor), sep="")
			}
			print(zz, digits = 4)
			cat("\n\n")
		})

# path show
setMethod("show",
			signature(object = "uGARCHpath"),
			function(object){
			vmodel = object@model$modeldesc$vmodel
			model = object@model
			cat(paste("\n*------------------------------------*", sep = ""))
			cat(paste("\n*     GARCH Model Path Simulation    *", sep = ""))
			cat(paste("\n*------------------------------------*", sep = ""))
			cat(paste("\nModel: ", vmodel, sep = ""))
			if(vmodel == "fGARCH"){
				cat(paste("\nfGARCH Sub-Model: ", model$modeldesc$vsubmodel, "\n", sep = ""))
			}
			sim = object@path
			sigma = sim$sigmaSim
			series = sim$seriesSim
			resids = sim$residSim
			m = dim(sigma)[2]
			N = dim(sigma)[1]
			cat(paste("\nHorizon: ", N))
			cat(paste("\nSimulations: ", m, "\n", sep = ""))
			sd1 = apply(sigma^2, 2, FUN = function(x) mean(x))
			sd2 = apply(sigma^2, 2, FUN = function(x) range(x))
			rx1 = apply(series, 2, FUN = function(x) mean(x))
			rx2 = apply(series, 2, FUN = function(x) range(x))
			xspec = .model2spec(as.list(object@model$pars[object@model$pars[,3]==1,1]), object@model, type = "GARCH")
			setfixed(xspec)<-as.list(object@model$pars[which(object@model$pars[,3]==1),1])
			uv = uncvariance(xspec)
			um = uncmean(xspec)
			uncond = c(NA, uv, NA, NA, um, NA, NA)
			dd = data.frame(Seed = object@seed, Sigma2.Mean = sd1, Sigma2.Min = sd2[1,],
					Sigma2.Max = sd2[2,], Series.Mean = rx1, Series.Min = rx2[1,],
					Series.Max = rx2[2,])
			meansim = apply(dd, 2, FUN = function(x) mean(x))
			meansim[1] = 0
			dd = rbind(dd, meansim, uncond)
			rownames(dd) = c(paste("sim", 1:m, sep = ""), "Mean(All)", "Unconditional")
			print(dd, digits = 3)
			cat("\n\n")
			})

# distribution show
setMethod("show",
		signature(object = "uGARCHdistribution"),
		function(object){
			vmodel = object@model$modeldesc$vmodel
			vsubmodel = object@model$modeldesc$vsubmodel
			model = object@model
			cat(paste("\n*------------------------------------*", sep = ""))
			cat(paste("\n*    GARCH Parameter Distribution    *", sep = ""))
			cat(paste("\n*------------------------------------*", sep = ""))
			cat(paste("\nModel : ", vmodel, sep = ""))
			if(vmodel == "fGARCH"){
				cat(paste("\nfGARCH SubModel : ", vsubmodel, sep = ""))
			}
			cat(paste("\nNo. Paths (m.sim) : ", object@dist$details$m.sim, sep = ""))
			cat(paste("\nLength of Paths (n.sim) : ", object@dist$details$n.sim, sep = ""))
			cat(paste("\nRecursive : ", object@dist$details$recursive, sep = ""))
			if(object@dist$details$recursive){
				cat(paste("\nRecursive Length : ", object@dist$details$recursive.length, sep = ""))
				cat(paste("\nRecursive Window : ", object@dist$details$recursive.window, sep = ""))
			}
			cat("\n\n")
			cat("Coefficients: True vs Simulation Mean (Window-n)\n")
			nwindows = object@dist$details$nwindows
			nm = object@dist$details$n.sim + (0:(nwindows-1))*object@dist$details$recursive.window
			ns = matrix(0, ncol = dim(object@truecoef)[1], nrow = nwindows)
			for(i in 1:nwindows){
				ns[i,] = apply(as.data.frame(object, window = i), 2, FUN = function(x) mean(x, na.rm = T))
			}
			ns = rbind(object@truecoef[,1], ns)
			colnames(ns) = rownames(object@truecoef)
			rownames(ns) = c("true-coef",paste("window-", nm, sep=""))
			print(as.data.frame(ns), digits=5)
			for(i in 1:nwindows){
				if(any(object@dist[[i]]$convergence==1)) n = length(which(object@dist[[i]]$convergence==1)) else n = 0
				if(n>0) cat(paste("\nwindow-",nm[i]," no. of non-converged fits: ", n, "\n",sep=""))
			}
			cat("\n\n")
		})


# boot show
setMethod("show",
		signature(object = "uGARCHboot"),
		function(object){
			vmodel = object@model$modeldesc$vmodel
			vsubmodel = object@model$modeldesc$vsubmodel
			cat(paste("\n*-----------------------------------*", sep = ""))
			cat(paste("\n*     GARCH Bootstrap Forecast      *", sep = ""))
			cat(paste("\n*-----------------------------------*", sep = ""))
			cat(paste("\nModel : ", vmodel, sep = ""))
			if(vmodel == "fGARCH"){
				cat(paste("\nfGARCH SubModel : ", vsubmodel, sep = ""))
			}
			cat(paste("\nn.ahead : ", object@model$n.ahead, sep = ""))
			cat(paste("\nBootstrap method: ",object@model$type))
			cat("\nDate (T[0]):", as.character(object@model$indexT))
			sig = sigma(object@forc)
			ser = fitted(object@forc)
			zs = cbind(t(as.data.frame(object, which = "sigma", type = "summary")),  sig)
			zr = cbind(t(as.data.frame(object, which = "series", type = "summary")), ser)
			colnames(zr)[6] = colnames(zs)[6] = c("forecast[analytic]")
			hh = min(object@model$n.ahead, 10)
			cat("\n\nSeries (summary):\n")
			print(head(round((zr), 6), hh), digits = 5)
			cat(".....................\n")
			cat("\nSigma (summary):\n")
			print(head(round((zs), 6),hh), digits = 5)
			cat(".....................")
			cat("\n\n")
		})

setMethod("show",
		signature(object = "uGARCHroll"),
		function(object){
			if(!is.null(object@model$noncidx)){
				cat("\nObject contains non-converged estimation windows. Use resume method to re-estimate.\n")
				invisible(object)
			} else{
				cat(paste("\n*-------------------------------------*", sep = ""))
				cat(paste("\n*              GARCH Roll             *", sep = ""))
				cat(paste("\n*-------------------------------------*", sep = ""))
				N = object@model$n.refits
				model = object@model$spec@model
				modelinc = model$modelinc
				vmodel = model$modeldesc$vmodel
				cat("\nNo.Refits\t\t:", N)
				cat("\nRefit Horizon\t:", object@model$refit.every)
				cat("\nNo.Forecasts\t:", NROW(object@forecast$density))
				cat(paste("\nGARCH Model\t\t: ", vmodel, "(",modelinc[8],",",modelinc[9],")\n", sep = ""))
				if(vmodel == "fGARCH"){
					cat(paste("\nfGARCH SubModel\t: ", model$modeldesc$vsubmodel, "\n", sep = ""))
				}
				cat("Distribution\t:", model$modeldesc$distribution,"\n")
				cat("\nForecast Density:\n")
				print(round(head(object@forecast$density),4))
				cat("\n..........................\n")
				print(round(tail(object@forecast$density),4))
				cat("\nElapsed:", format(object@model$elapsed))
				cat("\n")
				invisible(object)
			}
		})
#-------------------------------------------------------------------------
# multi-methods
setMethod("show",
		signature(object = "uGARCHmultispec"),
		function(object){
			cat(paste("\n*-----------------------------*", sep = ""))
			cat(paste("\n*     GARCH Multi-Spec        *", sep = ""))
			cat(paste("\n*-----------------------------*", sep = ""))
			N = length(object@spec)
			cat(paste("\nMultiple Specifications\t: ", N, sep=""))
			cat(paste("\nMulti-Spec Type\t\t\t: ", object@type, sep=""))
			cat("\n")
			invisible(object)
		})

setMethod("show",
		signature(object = "uGARCHmultifit"),
		function(object){
			cat(paste("\n*----------------------------*", sep = ""))
			cat(paste("\n*     GARCH Multi-Fit        *", sep = ""))
			cat(paste("\n*----------------------------*", sep = ""))
			cat(paste("\nNo. Assets :", length(object@fit), sep=""))
			asset.names = object@desc$asset.names

			if(object@desc$type == "equal"){
				vmodel = object@fit[[1]]@model$modeldesc$vmodel
				cat(paste("\nGARCH Multi-Spec Type : Equal",sep=""))
				cat(paste("\nGARCH Model Spec",sep=""))
				cat(paste("\n--------------------------",sep=""))
				cat(paste("\nModel : ", vmodel,sep=""))
				if(vmodel == "fGARCH" ){
					cat(paste("\nfGARCH SubModel : ", object@fit[[1]]@model$modeldesc$vsubmodel, "\n", sep = ""))
				}
				if(object@fit[[1]]@model$modelinc[15]>0){
					cat("\nExogenous Regressors in variance equation: ", object@fit[[1]]@model$modelinc[15], "\n")
				} else{
					cat("\nExogenous Regressors in variance equation: none\n")
				}
				cat("\nMean Equation :")
				cat("\nInclude Mean : ", object@fit[[1]]@model$modelinc[1])
				cat(paste("\nAR(FI)MA Model : (",object@fit[[1]]@model$modelinc[2],",",
								ifelse(object@fit[[1]]@model$modelinc[4]>0, 1, "d"),
								",",object@fit[[1]]@model$modelinc[3],")",sep=""))
				cat("\nGARCH-in-Mean : ", as.logical(object@fit[[1]]@model$modelinc[5]))
				if(object@fit[[1]]@model$modelinc[6]>0){
					cat("\nExogenous Regressors in mean equation: ", object@fit[[1]]@model$modelinc[6])
				} else{
					cat("\nExogenous Regressors in mean equation: none")
				}
				cat("\nConditional Distribution: ",object@fit[[1]]@model$modeldesc$distribution,"\n")
				cv = sapply(object@fit, FUN = function(x) x@fit$convergence)
				if(any(cv != 0)){
					ncv = which(cv != 0)
					nncv = length(ncv)
					cat("\nNo. of non converged fits: ", ncv,"\n")
					if(ncv>0) cat("\nNon converged fits: ", nncv,"\n")

				} else{
				cat(paste("\nGARCH Model Fit", sep = ""))
				cat(paste("\n--------------------------", sep = ""))
				cat("\nOptimal Parameters:\n")
				ll = t(likelihood(object))
				rownames(ll) = "Log-Lik"
				cf = coef(object)
				colnames(cf) = asset.names
				print(round(rbind(cf, ll), digits = 5))
				cat("\n")
			}
		} else{
			cat(paste("\nGARCH Model Fit", sep = ""))
			cat(paste("\n--------------------------", sep = ""))
			cat("\nOptimal Parameters:\n")
			print(coef(object), digits = 5)
		}
		invisible(object)
	})

setMethod("show",
		signature(object = "uGARCHmultifilter"),
		function(object){
			asset.names = object@desc$asset.names
			cat(paste("\n*-------------------------------*", sep = ""))
			cat(paste("\n*     GARCH Multi-Filter        *", sep = ""))
			cat(paste("\n*-------------------------------*", sep = ""))
			cat(paste("\nNo. Assets :", length(object@filter), sep=""))
			if(object@desc$type == "equal"){
					cat(paste("\nGARCH Model Filter", sep = ""))
					cat(paste("\n--------------------------", sep = ""))
					cat("\nParameters:\n")
					cf = coef(object)
					colnames(cf) = asset.names
					print(round(cf, digits = 5))
			} else{
				cat(paste("\nGARCH Model Filter", sep = ""))
				cat(paste("\n--------------------------", sep = ""))
				cat("\nOptimal Parameters:\n")
				print(coef(object), digits = 5)
			}
			invisible(object)
		})

setMethod("show",
		signature(object = "uGARCHmultiforecast"),
		function(object){
			asset.names = object@desc$asset.names
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\n*       GARCH Multi-Forecast      *", sep = ""))
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\nNo. Assets :", length(object@forecast), sep=""))
			cat(paste("\nn.ahead    :", object@forecast[[1]]@forecast$n.ahead, sep=""))
			cat(paste("\nn.roll     :", object@forecast[[1]]@forecast$n.roll, sep=""))
			cat("\n")
			invisible(object)
		})
#----------------------------------------------------------------------------------
# report method
#----------------------------------------------------------------------------------
report = function(object, ...)
{
	UseMethod("report")
}

setMethod("report", signature(object = "uGARCHroll"), .ugarchrollreport)

#----------------------------------------------------------------------------------
# univariate fit extractors
#----------------------------------------------------------------------------------
# coef methods
.ugarchfitcoef = function(object)
{
	object@fit$coef
}

setMethod("coef", signature(object = "uGARCHfit"), .ugarchfitcoef)

.ugarchfiltercoef = function(object)
{
	cf = object@model$pars[object@model$pars[,2]==1, 1]
	names(cf) = rownames(object@model$pars[object@model$pars[,2]==1, 1,drop=FALSE])
	return(cf)
}

setMethod("coef", signature(object = "uGARCHfilter"), .ugarchfiltercoef)

# multi-fit and multi-filter coefficients:
.ugarchmultifitcoef = function(object)
{
	if(object@desc$type == "equal"){
		ans = sapply(object@fit, FUN = function(x) coef(x), simplify = TRUE)
	} else{
		ans = lapply(object@fit, FUN = function(x) coef(x))
	}
	return(ans)
}
setMethod("coef", signature(object = "uGARCHmultifit"), .ugarchmultifitcoef)

.ugarchmultifiltercoef = function(object)
{

	ans = sapply(object@filter, FUN = function(x) coef(x), simplify = TRUE)
	return(ans)
}

setMethod("coef", signature(object = "uGARCHmultifilter"), .ugarchmultifiltercoef)


.ugarchrollcoef = function(object)
{
	if(!is.null(object@model$noncidx)) stop("\nObject contains non-converged estimation windows.")
	return(object@model$coef)
}

setMethod("coef", signature(object = "uGARCHroll"), .ugarchrollcoef)
#----------------------------------------------------------------------------------
# as.data.frame method for distribution object
.ugarchdistdf = function(x, row.names = NULL, optional = FALSE, which = "coef", window = 1, ...)
{
	n = x@dist$details$nwindows
	if(window > n) stop("window size greater than actual available", call. = FALSE)

	if(which == "rmse"){
		ans = as.data.frame(t(x@dist[[window]]$rmse))
		colnames(ans) = rownames(x@truecoef)
	}

	if(which == "stats"){
		llh = x@dist[[window]]$likelist
		persist = x@dist[[window]]$persist
		uncvar = x@dist[[window]]$vlongrun
		uncmean = x@dist[[window]]$mlongrun
		maxret = x@dist[[window]]$simmaxdata[,1]
		minret = x@dist[[window]]$simmindata[,1]
		meanret = x@dist[[window]]$simmeandata[,1]
		kurtosis = x@dist[[window]]$simmomdata[,1]
		skewness = x@dist[[window]]$simmomdata[,2]
		maxsigma = x@dist[[window]]$simmaxdata[,3]
		minsigma = x@dist[[window]]$simmindata[,3]
		meansigma = x@dist[[window]]$simmeandata[,3]
		ans = data.frame(llh = llh, persist = persist, uncvar = uncvar, uncmean = uncmean,
				maxret = maxret, minret = minret, meanret = meanret, kurtosis = kurtosis,
				skewness  = skewness, maxsigma = maxsigma, minsigma = minsigma, meansigma = meansigma)
	}

	if(which == "coef"){
		cf = x@dist[[window]]$simcoef
		ans = data.frame(coef = cf)
		colnames(ans) = rownames(x@truecoef)
	}

	if(which == "coefse"){
		cfe = x@dist[[window]]$simcoefse
		ans = data.frame(coefse = cfe)
		colnames(ans) = rownames(x@truecoef)
	}

	ans
}

setMethod("as.data.frame", signature(x = "uGARCHdistribution"), .ugarchdistdf)
#----------------------------------------------------------------------------------
# as.data.frame method for bootstrap object
.ugarchbootdf = function(x, row.names = NULL, optional = FALSE, which = "sigma", type = "raw", qtile = c(0.01, 0.099))
{
	n.ahead = x@model$n.ahead
	if(which == "sigma")
	{
		if(type == "raw"){
			sigma = x@fsigma
			ans = data.frame(bootsigma = sigma)
			colnames(ans) = paste("t+", 1:n.ahead, sep="")
		}
		if(type == "q"){
			if(all(is.numeric(qtile)) && (all(qtile<1.0) && all(qtile >0.0))){
				sigma = x@fsigma
				ans = apply(sigma, 2, FUN = function(x) stats::quantile(x, qtile, na.rm=TRUE))
				ans = as.data.frame(ans)
				colnames(ans) = paste("t+", 1:n.ahead, sep="")
				rownames(ans) = paste("q", qtile, sep = "")
			} else{
				stop("\nfor type q, the qtile value must be numeric and between (>)0 and 1(<)\n", call.  = FALSE)
			}
		}
		if(type == "summary"){
			sigma = x@fsigma
			ans = apply(sigma, 2, FUN = function(x) c(min(x, na.rm=TRUE), stats::quantile(x, 0.25, na.rm=TRUE), mean(x, na.rm=TRUE), stats::quantile(x, 0.75, na.rm=TRUE), max(x, na.rm=TRUE) ))
			ans = as.data.frame(ans)
			colnames(ans) = paste("t+", 1:n.ahead, sep="")
			rownames(ans) = c("min", "q0.25", "mean", "q0.75", "max")
		}
	}
	if(which == "series")
	{
		if(type == "raw"){
			series = x@fseries
			ans = data.frame(bootseries = series)
			colnames(ans) = paste("t+", 1:n.ahead, sep="")
		}
		if(type == "q"){
			if(all(is.numeric(qtile)) && (all(qtile<1.0) && all(qtile >0.0))){
						series = x@fseries
						ans = apply(series, 2, FUN = function(x) stats::quantile(x, qtile, na.rm=TRUE))
						ans = as.data.frame(ans)
						colnames(ans) = paste("t+", 1:n.ahead, sep="")
						rownames(ans) = paste("q", qtile, sep = "")
					} else{
						stop("\nfor type q, the qtile value must be numeric and between (>)0 and 1(<)\n", call.  = FALSE)
					}
		}
		if(type == "summary"){
			series = x@fseries
			ans = apply(series, 2, FUN = function(x) c(min(x, na.rm=TRUE), stats::quantile(x, 0.25, na.rm=TRUE), mean(x, na.rm=TRUE), stats::quantile(x, 0.75, na.rm=TRUE), max(x, na.rm=TRUE) ))
			ans = as.data.frame(ans)
			colnames(ans) = paste("t+", 1:n.ahead, sep="")
			rownames(ans) = c("min", "q.25", "mean", "q.75", "max")
		}
	}
	ans
}
setMethod("as.data.frame", signature(x = "uGARCHboot"), .ugarchbootdf)


#----------------------------------------------------------------------------------
# as.data.frame method for roll object
# valid which = density, fpm, coefs
.ugarchrolldf = function(x, row.names = NULL, optional = FALSE, which = "density")
{
	if(!is.null(x@model$noncidx)) stop("\nObject contains non-converged estimation windows.")
	if(which == "density") ans =  x@forecast$density else ans = x@forecast$VaR
	return(ans)
}
setMethod("as.data.frame", signature(x = "uGARCHroll"), .ugarchrolldf)
#----------------------------------------------------------------------------------
# residuals method
.ugarchresids = function(object, standardize = FALSE)
{
	if(class(object)[1] == "uGARCHfit" | class(object)[1] == "uGARCHfilter"){
		D = object@model$modeldata$index[1:object@model$modeldata$T]
	}
	if(standardize){
		ans = switch(class(object)[1],
				uGARCHfit = xts(object@fit$residuals/object@fit$sigma, D),
				uGARCHfilter = xts(object@filter$residuals/object@filter$sigma, D),
				uGARCHmultifit = sapply(object@fit, FUN = function(x) residuals(x, standardize = TRUE), simplify = TRUE),
				uGARCHmultifilter = sapply(object@filter, FUN = function(x) residuals(x, standardize = TRUE), simplify = TRUE))
	} else{
		ans = switch(class(object)[1],
				uGARCHfit = xts(object@fit$residuals, D),
				uGARCHfilter = xts(object@filter$residuals, D),
				uGARCHmultifit = sapply(object@fit, FUN = function(x) residuals(x), simplify = TRUE),
				uGARCHmultifilter = sapply(object@filter, FUN = function(x) residuals(x), simplify = TRUE))
	}
	return(ans)
}

setMethod("residuals", signature(object = "uGARCHfit"), .ugarchresids)
setMethod("residuals", signature(object = "uGARCHfilter"), .ugarchresids)
setMethod("residuals", signature(object = "uGARCHmultifit"), .ugarchresids)
setMethod("residuals", signature(object = "uGARCHmultifilter"), .ugarchresids)
#----------------------------------------------------------------------------------
# sigma method
sigma = function(object, ...)
{
	UseMethod("sigma")
}

# xts is returned for fitted and filtered object. Does not make sense to
# returns anything but a matrix for the simulated and forecast.
.ugarchsigma = function(object)
{
	if(class(object)[1] == "uGARCHfit" | class(object)[1] == "uGARCHfilter"){
		D = object@model$modeldata$index[1:object@model$modeldata$T]
	}
	switch(class(object)[1],
			uGARCHfit = xts(object@fit$sigma, D),
			uGARCHfilter = xts(object@filter$sigma, D),
			uGARCHmultifit = sapply(object@fit, FUN = function(x) sigma(x), simplify = TRUE),
			uGARCHmultifilter = sapply(object@filter, FUN = function(x) sigma(x), simplify = TRUE),
			uGARCHsim = {
				ans = object@simulation$sigmaSim
				rownames(ans) = paste("T+",1:NROW(object@simulation$sigmaSim), sep="")
				return(ans)
			},
			uGARCHpath ={
				ans = object@path$sigmaSim
				rownames(ans) = paste("T+",1:NROW(object@path$sigmaSim), sep="")
				return(ans)
			},
			uGARCHforecast = object@forecast$sigmaFor
	)
}

setMethod("sigma", signature(object = "uGARCHfit"), .ugarchsigma)
setMethod("sigma", signature(object = "uGARCHfilter"), .ugarchsigma)
setMethod("sigma", signature(object = "uGARCHmultifit"), .ugarchsigma)
setMethod("sigma", signature(object = "uGARCHmultifilter"), .ugarchsigma)
setMethod("sigma", signature(object = "uGARCHsim"), .ugarchsigma)
setMethod("sigma", signature(object = "uGARCHpath"), .ugarchsigma)
setMethod("sigma", signature(object = "uGARCHforecast"), .ugarchsigma)

.ugarchsigmamf = function(object)
{
	n.assets = length(object@forecast)
	n.ahead = object@forecast[[1]]@forecast$n.ahead
	n.roll = object@forecast[[1]]@forecast$n.roll+1
	Z = array(NA, dim = c(n.ahead, n.roll, n.assets))
	for(i in 1:n.assets){
		Z[,,i] = sigma(object@forecast[[i]])
	}
	return(Z)
}
setMethod("sigma", signature(object = "uGARCHmultiforecast"), .ugarchsigmamf)

#----------------------------------------------------------------------------------
# nyblom method
nyblom = function(object)
{
	UseMethod("nyblom")
}

setMethod("nyblom", signature(object = "uGARCHfit"), .nyblomTest)
#----------------------------------------------------------------------------------
# signbias method
signbias = function(object)
{
	UseMethod("signbias")
}

setMethod("signbias", signature(object = "uGARCHfit"), .signbiasTest)
setMethod("signbias", signature(object = "uGARCHfilter"), .signbiasTest)
#----------------------------------------------------------------------------------
# goodness of fit method
gof = function(object,groups)
{
	UseMethod("gof")
}

setMethod("gof", signature(object = "uGARCHfit", groups = "numeric"), .gofTest)
setMethod("gof", signature(object = "uGARCHfilter", groups = "numeric"), .gofTest)

#----------------------------------------------------------------------------------
# Info Criteria method
infocriteria = function(object)
{
	UseMethod("infocriteria")
}
.ugarchinfocriteria = function(object)
{
	# indicator object@fit$ipars[,4] denotes the estimated parameters
	if(is(object, "uGARCHfilter")){
		# np = sum(object@filter$ipars[,2])
		# all parameters fixed
		np = 0
	} else{
		np = sum(object@fit$ipars[,4])
	}
	itest = .information.test(likelihood(object), nObs = length(fitted(object)),
			nPars = np)
	itestm = matrix(0, ncol = 1, nrow = 4)
	itestm[1,1] = itest$AIC
	itestm[2,1] = itest$BIC
	itestm[3,1] = itest$SIC
	itestm[4,1] = itest$HQIC
	colnames(itestm) = ""
	rownames(itestm) = c("Akaike", "Bayes", "Shibata", "Hannan-Quinn")
	return(itestm)
}

setMethod("infocriteria", signature(object = "uGARCHfit"), .ugarchinfocriteria)
setMethod("infocriteria", signature(object = "uGARCHfilter"), .ugarchinfocriteria)

#----------------------------------------------------------------------------------
# Likelihood method
likelihood = function(object)
{
	UseMethod("likelihood")
}
.ugarchLikelihood = function(object)
{
	switch(class(object)[1],
			uGARCHfit = object@fit$LLH,
			uGARCHfilter = object@filter$LLH,
			uGARCHmultifilter = sapply(object@filter, FUN = function(x) likelihood(x), simplify = TRUE),
			uGARCHmultifit = sapply(object@fit, FUN = function(x) likelihood(x), simplify = TRUE))
}

setMethod("likelihood", signature(object = "uGARCHfit"), .ugarchLikelihood)
setMethod("likelihood", signature(object = "uGARCHfilter"), .ugarchLikelihood)
setMethod("likelihood", signature(object = "uGARCHmultifilter"), .ugarchLikelihood)
setMethod("likelihood", signature(object = "uGARCHmultifit"), .ugarchLikelihood)
#----------------------------------------------------------------------------------
# Fitted method
.ugarchfitted = function(object)
{
	if(class(object)[1] == "uGARCHfit" | class(object)[1] == "uGARCHfilter"){
		D = object@model$modeldata$index[1:object@model$modeldata$T]
	}
	switch(class(object)[1],
			uGARCHfit = xts(object@fit$fitted.values, D),
			uGARCHfilter = xts(object@model$modeldata$data[1:object@model$modeldata$T] - object@filter$residuals, D),
			uGARCHmultifilter = sapply(object@filter, FUN = function(x) fitted(x), simplify = TRUE),
			uGARCHmultifit = sapply(object@fit, FUN = function(x) fitted(x), simplify = TRUE),
			uGARCHsim = {
				ans = object@simulation$seriesSim
				rownames(ans) = paste("T+",1:NROW(object@simulation$seriesSim), sep="")
				return(ans)
			},
			uGARCHpath ={
				ans = object@path$seriesSim
				rownames(ans) = paste("T+",1:NROW(object@path$seriesSim), sep="")
				return(ans)
			},
			uGARCHforecast = object@forecast$seriesFor
			)
}
setMethod("fitted", signature(object = "uGARCHfit"), .ugarchfitted)
setMethod("fitted", signature(object = "uGARCHfilter"), .ugarchfitted)
setMethod("fitted", signature(object = "uGARCHmultifit"), .ugarchfitted)
setMethod("fitted", signature(object = "uGARCHmultifilter"), .ugarchfitted)
setMethod("fitted", signature(object = "uGARCHsim"), .ugarchfitted)
setMethod("fitted", signature(object = "uGARCHpath"), .ugarchfitted)
setMethod("fitted", signature(object = "uGARCHforecast"), .ugarchfitted)


.ugarchfittedmf = function(object)
{
	n.assets = length(object@forecast)
	n.ahead = object@forecast[[1]]@forecast$n.ahead
	n.roll = object@forecast[[1]]@forecast$n.roll+1
	Z = array(NA, dim = c(n.ahead, n.roll, n.assets))
	for(i in 1:n.assets){
		Z[,,i] = fitted(object@forecast[[i]])
	}
	return(Z)
}
setMethod("fitted", signature(object = "uGARCHmultiforecast"), .ugarchfittedmf)
#----------------------------------------------------------------------------------
# reduce method
reduce = function(object, pvalue = 0.1, ...)
{
	UseMethod("reduce")
}

.reduce = function(object, pvalue = 0.1, use.robust=TRUE, ...){
	# does not yet contain a check for the case of all parameters eliminated
	if(use.robust) cf = object@fit$robust.matcoef else cf = object@fit$matcoef
	idx = which(cf[,4]>pvalue)
	if(length(idx)>0){
		pars = cf[idx,1]
		# zero the coefficients
		pars = pars*0
		names(pars) = rownames(cf)[idx]
		spec = getspec(object)
		setfixed(spec)<-as.list(pars)
		setstart(spec)<-as.list(cf[-idx,1])
		dat = xts(object@model$modeldata$data, object@model$modeldata$index)
		ns = object@model$n.start
		if(class(object)=="uGARCHfit"){
			fit = ugarchfit(spec, dat, out.sample = ns, ...)
		} else{
			fit = arfimafit(spec, dat, out.sample = ns, ...)
		}
		return(fit)
	} else{
		return(object)
	}
}

setMethod("reduce", signature(object = "uGARCHfit"), .reduce)
#----------------------------------------------------------------------------------
# quantile S4 method
.ugarchquantile = function(x, probs=c(0.01, 0.05))
{
	if(class(x)=="uGARCHroll"){
		d = x@model$spec@model$modeldesc$distribution
	} else{
		d = x@model$modeldesc$distribution
		s = sigma(x)
		m = fitted(x)
	}
	di = .DistributionBounds(d)
	if(di$include.skew)  skew  = x@model$pars["skew",1] else skew = 0
	if(di$include.shape) shape = x@model$pars["shape",1] else shape = 0
	if(di$include.ghlambda) ghlambda = x@model$pars["ghlambda",1] else ghlambda = 0

	if(class(x)=="uGARCHforecast"){
		ans = matrix(NA, dim(s)[1], dim(s)[2])
		if(length(probs)>1) stop("\nprobs must be a scalar for a uGARCHforecast object")
		for(i in 1:NCOL(ans)) ans[,i] = qdist(d, probs[1], mu = m[,i], sigma = s[,i],
					skew = skew, shape = shape, lambda = ghlambda)
		colnames(ans) = colnames(s)
		rownames(ans) = rownames(s)
	} else if(class(x)=="uGARCHsim" | class(x)=="uGARCHpath"){
		if(length(probs)>1) stop("\nprobs must be a scalar for a uGARCHsim and uGARCHpath objects")
		ans = matrix(NA, dim(s)[1], dim(s)[2])
		for(i in 1:NCOL(ans)) ans[,i] = qdist(d, probs[1], mu = m[,i], sigma = s[,i],
					skew = skew, shape = shape, lambda = ghlambda)
		colnames(ans) = colnames(s)
		rownames(ans) = rownames(s)
	} else if(class(x)=="uGARCHroll"){
		skew = x@forecast$density[,"Skew"]
		shape = x@forecast$density[,"Shape"]
		lambda = x@forecast$density[,"Shape(GIG)"]
		s = x@forecast$density[,"Sigma"]
		m = x@forecast$density[,"Mu"]
		ans = matrix(NA, ncol = length(probs), nrow = length(s))
		for(i in 1:length(probs)) ans[,i] = as.numeric(m) + qdist(d, probs[i], skew = skew, shape = shape,
					lambda = lambda) * as.numeric(s)
		colnames(ans) = paste("q[", probs,"]", sep="")
		ans = xts(ans, as.POSIXct(rownames(x@forecast$density)))
	} else{
		ans = matrix(NA, ncol = length(probs), nrow = NROW(s))
		for(i in 1:length(probs)) ans[,i] = as.numeric(m) + qdist(d, probs[i], skew = skew, shape = shape,
					lambda = ghlambda) * as.numeric(s)
		colnames(ans) = paste("q[", probs,"]", sep="")
		ans = xts(ans, index(s))
	}
	return(ans)
}
setMethod("quantile", signature(x = "uGARCHfit"), .ugarchquantile)
setMethod("quantile", signature(x = "uGARCHfilter"), .ugarchquantile)
setMethod("quantile", signature(x = "uGARCHforecast"), .ugarchquantile)
setMethod("quantile", signature(x = "uGARCHsim"), .ugarchquantile)
setMethod("quantile", signature(x = "uGARCHpath"), .ugarchquantile)
setMethod("quantile", signature(x = "uGARCHroll"), .ugarchquantile)

# pit S4 method
pit = function(object, ...)
{
	UseMethod("pit")
}

.ugarchpit = function(object)
{
	if(class(object)=="uGARCHroll"){
		d = object@model$spec@model$modeldesc$distribution
		skew = object@forecast$density[,"Skew"]
		shape = object@forecast$density[,"Shape"]
		lambda = object@forecast$density[,"Shape(GIG)"]
		s = object@forecast$density[,"Sigma"]
		m = object@forecast$density[,"Mu"]
		r = object@forecast$density[,"Realized"]
		ans =  pdist(d, q = r, mu = as.numeric(m), sigma = as.numeric(s),
				skew = skew, shape = shape, lambda = ghlambda)
		ans = xts(ans, as.POSIXct(rownames(object@forecast$density)))
		colnames(ans) = "pit"
	} else{
		d = object@model$modeldesc$distribution
		di = .DistributionBounds(d)
		if(di$include.skew)  skew  = object@model$pars["skew",1] else skew = 0
		if(di$include.shape) shape = object@model$pars["shape",1] else shape = 0
		if(di$include.ghlambda) ghlambda = object@model$pars["ghlambda",1] else ghlambda = 0
		s = sigma(object)
		m = fitted(object)
		r = object@model$modeldata$data[1:object@model$modeldata$T]
		ans =  pdist(d, q = r, mu = as.numeric(m), sigma = as.numeric(s),
				skew = skew, shape = shape, lambda = ghlambda)
		ans = xts(ans, index(s))
		colnames(ans) = "pit"
	}
	return(ans)
}
setMethod("pit", signature(object = "uGARCHfit"), .ugarchpit)
setMethod("pit", signature(object = "uGARCHfilter"), .ugarchpit)
setMethod("pit", signature(object = "uGARCHroll"), .ugarchpit)
# newsimpact curve method (not multiplied by unconditional sigma)
newsimpact = function(object, z = seq(-0.3, 0.3, length.out = 100))
{
	UseMethod("newsimpact")
}

.newsimpact = function(object, z = seq(-0.3, 0.3, length.out = 100))
{
	# For the mcsGARCH and csGARCH the stochastic and permanent components,
	# respectively, are shown, which are sGARCH
	vmodel = object@model$modeldesc$vmodel
	ans = switch(vmodel,
			sGARCH = .sgarchni(z, object),
			fGARCH = .fgarchni(z, object),
			gjrGARCH = .gjrgarchni(z, object),
			eGARCH = .egarchni(z, object),
			apARCH = .aparchni(z, object),
			iGARCH = .sgarchni(z, object),
			csGARCH = .sgarchni(z, object),
			mcsGARCH = .sgarchni(z, object),
			realGARCH = .realgarchni(z, object))
	return(ans)
}

setMethod("newsimpact", signature(object = "uGARCHfit"), .newsimpact)
setMethod("newsimpact", signature(object = "uGARCHfilter"), .newsimpact)

# the underlying news impact methods
.sgarchni = function(z, object)
{
	if(is.null(z)){
		if(is(object, "uGARCHfit")) mz = min(object@fit$residuals) else mz = min(object@filter$residuals)
		if(abs(mz) <1){
			zz = seq(round(mz, 2) , round(abs(mz), 2), length.out = 101)
		} else{
			zz = seq(round(mz, 0) , round(abs(mz), 0), length.out = 101)
		}
	} else{
		zz = z
	}
	if(is(object, "uGARCHfit")) ipars = object@fit$ipars else ipars = object@filter$ipars
	modelinc = object@model$modelinc
	idx = object@model$pidx
	omega = ipars[idx["omega",1],1]
	alpha = beta = 0
	if(modelinc[8]>0) alpha = ipars[idx["alpha",1],1]
	if(modelinc[9]>0) beta  = ipars[idx["beta",1],1]
	lrvar = rep(as.numeric(uncvariance(object)), length(zz))
	ans = omega + beta*lrvar + alpha*zz^2
	yexpr = expression(sigma[t]^2)
	xexpr = expression(epsilon[t-1])
	return(list(zy = as.numeric(ans), zx = zz, yexpr = yexpr, xexpr = xexpr))
}

.fgarchni = function(z, object)
{
	if(is.null(z)){
		if(is(object, "uGARCHfit")) mz = min(object@fit$residuals) else mz = min(object@filter$residuals)
		if(abs(mz) <1){
			zz = seq(round(mz, 2) , round(abs(mz), 2), length.out = 101)
		} else{
			zz = seq(round(mz, 0) , round(abs(mz), 0), length.out = 101)
		}
	} else{
		zz = z
	}
	if(is(object, "uGARCHfit")) ipars = object@fit$ipars else ipars = object@filter$ipars
	modelinc = object@model$modelinc
	idx = object@model$pidx
	omega = ipars[idx["omega",1],1]
	alpha = beta = lambda = delta = eta1 = eta2 = 0
	alpha = ipars[idx["alpha",1],1]
	beta = ipars[idx["beta",1],1]
	lambda = ipars[idx["lambda",1],1]
	delta = ipars[idx["delta",1],1]
	eta1 = ipars[idx["eta1",1],1]
	eta2 = ipars[idx["eta2",1],1]
	fk = object@model$fmodel$fpars$fk
	kdelta = delta + fk*lambda
	lrvar = rep(as.numeric(uncvariance(object))^(1/2), length(zz))
	ans = omega + alpha[1]*(lrvar^lambda)*(abs(zz/lrvar - eta2[1]) - eta1[1]*(zz/lrvar - eta2[1]))^kdelta + beta[1]*(lrvar^lambda)
	yexpr = expression(sigma[t]^2)
	xexpr = expression(epsilon[t-1])
	return(list(zy = ans^(2/lambda), zx = zz, yexpr = yexpr, xexpr = xexpr))
}

.egarchni = function(z, object)
{
	if(is.null(z)){
		if(is(object, "uGARCHfit")) mz = min(object@fit$residuals) else mz = min(object@filter$residuals)
		if(abs(mz) <1){
			zz = seq(round(mz, 2) , round(abs(mz), 2), length.out = 101)
		} else{
			zz = seq(round(mz, 0) , round(abs(mz), 0), length.out = 101)
		}
	} else{
		zz = z
	}
	if(is(object, "uGARCHfit")) ipars = object@fit$ipars else ipars = object@filter$ipars
	modelinc = object@model$modelinc
	idx = object@model$pidx
	omega = ipars[idx["omega",1],1]
	alpha = beta = gamma = skew = shape = ghlambda = 0
	if(modelinc[8]>0) alpha = ipars[idx["alpha",1],1]
	if(modelinc[9]>0) beta = ipars[idx["beta",1],1]
	if(modelinc[10]>0) gamma = ipars[idx["gamma",1],1]
	if(modelinc[16]>0) skew = ipars[idx["skew",1],1]
	if(modelinc[17]>0) shape = ipars[idx["shape",1],1]
	if(modelinc[18]>0) ghlambda = ipars[idx["ghlambda",1],1]
	k = egarchKappa(ghlambda, shape, skew, object@model$modeldesc$distribution)
	lrvar = rep(as.numeric(uncvariance(object)), length(zz))
	sqlr = sqrt(lrvar)
	ans=exp(omega + alpha[1]*zz/sqlr + gamma[1]*(abs(zz/sqlr)-k) + beta[1]*log(lrvar))
	yexpr = expression(sigma[t]^2)
	xexpr = expression(epsilon[t-1])
	return(list(zy = ans, zx = zz, yexpr = yexpr, xexpr = xexpr))
}

.gjrgarchni = function(z, object)
{
	if(is.null(z)){
		if(is(object, "uGARCHfit")) mz = min(object@fit$residuals) else mz = min(object@filter$residuals)
		if(abs(mz) <1){
			zz = seq(round(mz, 2) , round(abs(mz), 2), length.out = 101)
		} else{
			zz = seq(round(mz, 0) , round(abs(mz), 0), length.out = 101)
		}
	} else{
		zz = z
	}
	if(is(object, "uGARCHfit")) ipars = object@fit$ipars else ipars = object@filter$ipars
	modelinc = object@model$modelinc
	idx = object@model$pidx
	omega = ipars[idx["omega",1],1]
	alpha = beta = gamma = 0
	if(modelinc[8]>0) alpha = ipars[idx["alpha",1],1]
	if(modelinc[9]>0) beta = ipars[idx["beta",1],1]
	if(modelinc[10]>0) gamma = ipars[idx["gamma",1],1]
	lrvar = rep(as.numeric(uncvariance(object)), length(zz))
	ans = omega + alpha[1]*zz^2 + gamma[1]*(zz^2)*(zz<0) + beta[1]*lrvar
	yexpr = expression(sigma[t]^2)
	xexpr = expression(epsilon[t-1])
	return(list(zy = ans, zx = zz, yexpr = yexpr, xexpr = xexpr))
}

.aparchni = function(z, object)
{

	if(is.null(z)){
		if(is(object, "uGARCHfit")) mz = min(object@fit$residuals) else mz = min(object@filter$residuals)
		if(abs(mz) <1){
			zz = seq(round(mz, 2) , round(abs(mz), 2), length.out = 101)
		} else{
			zz = seq(round(mz, 0) , round(abs(mz), 0), length.out = 101)
		}
	} else{
		zz = z
	}
	if(is(object, "uGARCHfit")) ipars = object@fit$ipars else ipars = object@filter$ipars
	modelinc = object@model$modelinc
	idx = object@model$pidx
	omega = ipars[idx["omega",1],1]
	alpha = beta = gamma = 0
	delta = 2
	if(modelinc[8]>0) alpha = ipars[idx["alpha",1],1]
	if(modelinc[9]>0) beta = ipars[idx["beta",1],1]
	if(modelinc[10]>0) gamma = ipars[idx["gamma",1],1]
	if(modelinc[13]>0) delta = ipars[idx["delta",1],1]
	lrvar = rep(as.numeric(uncvariance(object))^(1/2), length(zz))
	ans = omega + alpha[1]*(abs(zz) - gamma[1]*(zz))^delta + beta[1]*(lrvar^delta)
	yexpr = expression(sigma[t]^2)
	xexpr = expression(epsilon[t-1])
	return(list(zy = ans^(2/delta), zx = zz, yexpr = yexpr, xexpr = xexpr))
}

.realgarchni = function(z, object)
{
	# alpha[1]*tau(z)
	# tau = function(z, eta1, eta2){ eta1*z + eta2*(z^2-1) }
	if(is.null(z)){
		if(is(object, "uGARCHfit")) mz = min(object@fit$residuals) else mz = min(object@filter$residuals)
		if(abs(mz) <1){
			zz = seq(round(mz, 2) , round(abs(mz), 2), length.out = 101)
		} else{
			zz = seq(round(mz, 0) , round(abs(mz), 0), length.out = 101)
		}
	} else{
		zz = z
	}
	if(is(object, "uGARCHfit")) ipars = object@fit$ipars else ipars = object@filter$ipars
	modelinc = object@model$modelinc
	idx = object@model$pidx
	alpha = eta1 = eta2 = 0
	alpha = ipars[idx["alpha",1],1]
	eta1 = ipars[idx["eta1",1],1]
	eta2 = ipars[idx["eta2",1],1]
	ans = alpha * (eta1*z + eta2*(z^2-1))
	yexpr = expression(paste(symbol('\045'),Delta(sigma[t])))
	xexpr = expression(z[t-1])
	return(list(zy = ans*100, zx = zz, yexpr = yexpr, xexpr = xexpr))
}
#----------------------------------------------------------------------------------
# Half-Life Method for the various garch models
# ln(0.5)/log(persistence)
halflife = function(object, pars, distribution = "norm", model = "sGARCH",
		submodel = "GARCH")
{
	UseMethod("halflife")
}

.halflife1<-function(object)
{
	ps = persistence(object)
	hlf = log(0.5)/log(ps)
	#names(hlf) = "Half-Life"
	return(hlf)
}

.halflife2 = function(pars, distribution = "norm", model = "sGARCH",
		submodel = "GARCH")
{
	ps = persistence(pars = pars, distribution = distribution, model = model, submodel = submodel)
	hlf = log(0.5)/log(ps)
	#names(hlf) = "Half-Life"
	return(hlf)
}
setMethod("halflife",signature(object = "uGARCHfilter", pars = "missing",
				distribution = "missing", model = "missing", submodel = "missing"),
		definition = .halflife1)

setMethod("halflife",signature(object = "uGARCHfit", pars = "missing",
				distribution = "missing", model = "missing", submodel = "missing"),
		definition = .halflife1)

setMethod("halflife",signature(object = "uGARCHspec", pars = "missing",
				distribution = "missing", model = "missing", submodel = "missing"),
		definition = .halflife1)

setMethod("halflife",signature(object = "missing", pars = "numeric",
				distribution = "character", model = "character", submodel = "ANY"),
		definition = .halflife2)
#----------------------------------------------------------------------------------
# Persistence
persistence = function(object, pars, distribution = "norm", model = "sGARCH",
		submodel = "GARCH")
{
	UseMethod("persistence")
}

# filter method
.filterpersistence = function(object)
{
	ans = object@filter$persistence
	#names(ans) = "persistence"
	return(ans)
}
setMethod("persistence", signature(object = "uGARCHfilter", pars = "missing",
				distribution = "missing", model = "missing", submodel = "missing"),
		definition = .filterpersistence)


# fit method
.persistence1 = function(object)
{
		pars = object@fit$ipars
		idx = object@model$pidx
		distribution = object@model$modeldesc$distribution
		vmodel = object@model$modeldesc$vmodel
		vsubmodel = object@model$modeldesc$vsubmodel
	ans = switch(vmodel,
			sGARCH = .persistsgarch1(pars, idx, distribution),
			eGARCH = .persistegarch1(pars, idx, distribution),
			gjrGARCH = .persistgjrgarch1(pars, idx, distribution),
			apARCH = .persistaparch1(pars, idx, distribution),
			fGARCH = .persistfgarch1(pars, idx, distribution, vsubmodel),
			csGARCH = .persistcsgarch1(pars, idx, distribution),
			mcsGARCH = .persistmcsgarch1(pars, idx, distribution),
			realGARCH = .persistrealgarch1(pars, idx, distribution),
			iGARCH = 1)
	#names(ans) = "persistence"
	return(ans)
}

.persistence2 = function(object)
{
	if(is.null(object@model$fixed.pars))
		stop("\nuncvariance with spec required fixed.pars list\n", call. = FALSE)
	# no other checks for now.
	model = object@model
	pars = unlist(model$fixed.pars)
	parnames = names(pars)
	modelnames = .checkallfixed(object)
	if(is.na(all(match(modelnames, parnames), 1:length(modelnames)))) {
		cat("\npersistence-->error: parameters names do not match specification\n")
		cat("Expected Parameters are: ")
		cat(paste(modelnames))
		cat("\n")
		stop("Exiting", call. = FALSE)
	}
	# once more into the spec
	setfixed(object)<-as.list(pars)
	model = object@model
	idx = model$pidx
	modelinc = model$modelinc
	pars = object@model$pars
	distribution = object@model$modeldesc$distribution
	vmodel = object@model$modeldesc$vmodel
	vsubmodel = object@model$modeldesc$vsubmodel
	ans = switch(vmodel,
			sGARCH = .persistsgarch1(pars, idx, distribution),
			eGARCH = .persistegarch1(pars, idx, distribution),
			gjrGARCH = .persistgjrgarch1(pars, idx, distribution),
			apARCH = .persistaparch1(pars, idx, distribution),
			fGARCH = .persistfgarch1(pars, idx, distribution, vsubmodel),
			csGARCH = .persistcsgarch1(pars, idx, distribution),
			mcsGARCH = .persistmcsgarch1(pars, idx, distribution),
			realGARCH = .persistrealgarch1(pars, idx, distribution),
			iGARCH = 1)
	#names(ans) = "persistence"
	return(ans)
}


.persistence3 = function(pars, distribution = "norm", model = "sGARCH",
		submodel = "GARCH")
{
	ans = switch(model,
			sGARCH = .persistsgarch2(pars, distribution),
			eGARCH = .persistegarch2(pars, distribution),
			gjrGARCH = .persistgjrgarch2(pars, distribution),
			apARCH = .persistaparch2(pars, distribution),
			fGARCH = .persistfgarch2(pars, distribution, submodel),
			csGARCH = .persistcsgarch2(pars, distribution),
			mcsGARCH = .persistmcsgarch2(pars, distribution),
			realGARCH = .persistrealgarch2(pars, distribution),
			iGARCH = 1)
	#names(ans) = "persistence"
	return(ans)
}


setMethod("persistence",signature(object = "uGARCHfit", pars = "missing",
				distribution = "missing", model = "missing", submodel = "missing"),
		definition = .persistence1)
setMethod("persistence",signature(object = "uGARCHspec", pars = "missing",
				distribution = "missing", model = "missing", submodel = "missing"),
		definition = .persistence2)
setMethod("persistence",signature(object = "missing", pars = "numeric",
				distribution = "character", model = "character", submodel = "ANY"),
		definition = .persistence3)

.persistsgarch1 = function(pars, idx, distribution = "norm"){
	ps = sum(pars[idx["alpha",1]:idx["alpha",2]]) + sum(pars[idx["beta",1]:idx["beta",2]])
	return(ps)
}

.persistsgarch2 = function(pars, idx, distribution = "norm"){
	Names = names(pars)
	if(any(substr(Names, 1, 5)=="alpha")){
		i = which(substr(Names, 1, 5)=="alpha")
		alpha = pars[i]
	} else{
		alpha = 0
	}
	if(any(substr(Names, 1, 4)=="beta")){
		i = which(substr(Names, 1, 4)=="beta")
		beta = pars[i]
	} else{
		beta = 0
	}
	ps = sum(alpha) + sum(beta)
	return(ps)
}


.persistmcsgarch1 = function(pars, idx, distribution = "norm"){
	ps = sum(pars[idx["alpha",1]:idx["alpha",2]]) + sum(pars[idx["beta",1]:idx["beta",2]])
	names(ps) = "Stochastic"
	return(ps)
}


.persistmcsgarch2 = function(pars, idx, distribution = "norm"){
	Names = names(pars)
	if(any(substr(Names, 1, 5)=="alpha")){
		i = which(substr(Names, 1, 5)=="alpha")
		alpha = pars[i]
	} else{
		alpha = 0
	}
	if(any(substr(Names, 1, 4)=="beta")){
		i = which(substr(Names, 1, 4)=="beta")
		beta = pars[i]
	} else{
		beta = 0
	}
	ps = sum(alpha) + sum(beta)
	names(ps) = "Stochastic"
	return(ps)
}

.persistcsgarch1 = function(pars, idx, distribution = "norm"){
	ps1 = sum(pars[idx["alpha",1]:idx["alpha",2]]) + sum(pars[idx["beta",1]:idx["beta",2]])
	ps2 = pars[idx["eta1",1]]
	names(ps1) = "Transitory"
	names(ps2) = "Permanent"
	return(c(ps1, ps2))
}

.persistcsgarch2 = function(pars, idx, distribution = "norm"){
	Names = names(pars)
	if(any(substr(Names, 1, 5)=="alpha")){
		i = which(substr(Names, 1, 5)=="alpha")
		alpha = pars[i]
	} else{
		alpha = 0
	}
	if(any(substr(Names, 1, 4)=="beta")){
		i = which(substr(Names, 1, 4)=="beta")
		beta = pars[i]
	} else{
		beta = 0
	}
	ps1 = sum(alpha) + sum(beta)
	if(any(substr(Names, 1, 5)=="eta11")){
		i = which(substr(Names, 1, 5)=="eta11")
		ps2 = pars[i]
	} else{
		ps2 = 0
	}
	names(ps1) = "Transitory"
	names(ps2) = "Permanent"
	return(c(ps1, ps2))
}

.persistgjrgarch1 = function(pars, idx, distribution = "norm"){
	alpha = pars[idx["alpha",1]:idx["alpha",2]]
	beta  = pars[idx["beta",1]:idx["beta",2]]
	gamma = pars[idx["gamma",1]:idx["gamma",2]]
	skew = pars[idx["skew",1]:idx["skew",2]]
	shape = pars[idx["shape",1]:idx["shape",2]]
	ghlambda = pars[idx["ghlambda",1]:idx["ghlambda",2]]

	ps = sum(beta)+ sum(alpha)+sum(apply(as.data.frame(gamma),1,FUN=function(x)
						x*pneg(ghlambda, shape, skew, distribution)))
	return(ps)
}

.persistgjrgarch2 = function(pars, distribution = "norm"){
	Names = names(pars)
	if(any(substr(Names, 1, 5)=="alpha")){
		i = which(substr(Names, 1, 5)=="alpha")
		alpha = pars[i]
	} else{
		alpha = 0
	}

	if(any(substr(Names, 1, 5)=="gamma")){
		i = which(substr(Names, 1, 5)=="gamma")
		gamma = pars[i]
	} else{
		gamma = 0
	}

	if(any(substr(Names, 1, 4)=="beta")){
		i = which(substr(Names, 1, 4)=="beta")
		beta = pars[i]
	} else{
		beta = 0
	}


	if(any(substr(Names, 1, 8)=="ghlambda")){
		i = which(substr(Names, 1, 8)=="ghlambda")
		ghlambda = pars[i]
	} else{
		ghlambda = 0
	}

	if(any(substr(Names, 1, 4)=="skew")){
		i = which(substr(Names, 1, 4)=="skew")
		skew = pars[i]
	} else{
		skew = 0
	}

	if(any(substr(Names, 1, 5)=="shape")){
		i = which(substr(Names, 1, 5)=="shape")
		shape = pars[i]
	} else{
		shape = 0
	}

	ps = sum(beta)+ sum(alpha)+sum(apply(as.data.frame(gamma),1,FUN=function(x)
						x*pneg(ghlambda, shape, skew, distribution)))
	return(ps)
}

.persistegarch1 = function(pars, idx, distribution = "norm"){
	beta  = pars[idx["beta",1]:idx["beta",2]]
	ps = sum(beta)
	return(ps)
}

.persistegarch2 = function(pars, distribution = "norm"){
	Names=names(pars)
	if(any(substr(Names, 1, 4)=="beta")){
		i=which(substr(Names, 1, 4)=="beta")
		beta=pars[i]
	} else{
		beta=0
	}
	ps=sum(beta)
	return(ps)
}

.persistaparch1 = function(pars, idx, distribution = "norm"){
	alpha = pars[idx["alpha",1]:idx["alpha",2]]
	beta  = pars[idx["beta",1]:idx["beta",2]]
	gamma = pars[idx["gamma",1]:idx["gamma",2]]
	delta = pars[idx["delta",1]:idx["delta",2]]
	skew = pars[idx["skew",1]:idx["skew",2]]
	shape = pars[idx["shape",1]:idx["shape",2]]
	ghlambda = pars[idx["ghlambda",1]:idx["ghlambda",2]]
	ps = sum(beta) + sum(apply(cbind(gamma,alpha), 1, FUN=function(x)
						x[2]*aparchKappa(x[1], delta, ghlambda, shape, skew,distribution)))
	return(ps)
}

.persistaparch2 = function(pars, distribution = "norm"){
	Names = names(pars)
	if(any(substr(Names, 1, 5)=="alpha")){
		i = which(substr(Names, 1, 5)=="alpha")
		alpha = pars[i]
	} else{
		alpha = 0
	}

	if(any(substr(Names, 1, 5)=="gamma")){
		i = which(substr(Names, 1, 5)=="gamma")
		gamma = pars[i]
	} else{
		gamma = rep(0, length(alpha))
	}

	if(any(substr(Names, 1, 5)=="delta")){
		i = which(substr(Names, 1, 5)=="delta")
		delta = pars[i]
	} else{
		delta=  2
	}

	if(any(substr(Names, 1, 4)=="beta")){
		i = which(substr(Names, 1, 4)=="beta")
		beta = pars[i]
	} else{
		beta = 0
	}
	if(any(substr(Names, 1, 8)=="ghlambda")){
		i = which(substr(Names, 1, 8)=="ghlambda")
		ghlambda = pars[i]
	} else{
		ghlambda = 0
	}

	if(any(substr(Names, 1, 4)=="skew")){
		i = which(substr(Names, 1, 4)=="skew")
		skew = pars[i]
	} else{
		skew = 0
	}

	if(any(substr(Names, 1, 5)=="shape")){
		i = which(substr(Names, 1, 5)=="shape")
		shape = pars[i]
	} else{
		shape = 0
	}
	ps = sum(beta) + sum(apply(cbind(gamma,alpha), 1, FUN=function(x)
						x[2]*aparchKappa(x[1], delta, ghlambda, shape, skew,distribution)))
	return(ps)
}

.persistfgarch1 = function(pars, idx, distribution = "norm", submodel){
	fm = .fgarchModel(submodel)
	alpha = pars[idx["alpha",1]:idx["alpha",2]]
	beta  = pars[idx["beta",1]:idx["beta",2]]
	eta1 = pars[idx["eta1",1]:idx["eta1",2]]
	eta2 = pars[idx["eta2",1]:idx["eta2",2]]
	lambda = pars[idx["lambda",1]:idx["lambda",2]]
	delta = pars[idx["delta",1]:idx["delta",2]]
	skew = pars[idx["skew",1]:idx["skew",2]]
	shape = pars[idx["shape",1]:idx["shape",2]]
	ghlambda = pars[idx["ghlambda",1]:idx["ghlambda",2]]
	fk = fm$parameters$fk
	ps = sum(beta) + sum(apply(cbind(alpha, eta1, eta2), 1, FUN=function(x)
						x[1]*fgarchKappa(lambda, delta, x[2], x[3], fk, ghlambda, shape, skew, distribution)))
	return(ps)
}

.persistfgarch2 = function(pars, distribution = "norm", submodel){
	fm = .fgarchModel(submodel)
	Names = names(pars)
	if(any(substr(Names, 1, 5)=="alpha")){
		i = which(substr(Names, 1, 5)=="alpha")
		alpha = pars[i]
	} else{
		alpha = 0
	}

	if(any(substr(Names, 1, 6)=="lambda")){
		i = which(substr(Names, 1, 6)=="lambda")
		lambda = pars[i]
	} else{
		lambda = fm$parameters$lambda
	}

	if(any(substr(Names, 1, 5)=="delta")){
		i = which(substr(Names, 1, 5)=="delta")
		delta = pars[i]
	} else{
		delta = fm$parameters$delta
	}

	if(any(substr(Names, 1, 4)=="eta1")){
		i = which(substr(Names, 1, 4)=="eta1")
		eta1 = pars[i]
	} else{
		eta1 = rep(0, length(alpha))
	}

	if(any(substr(Names, 1, 4)=="eta2")){
		i = which(substr(Names, 1, 4)=="eta2")
		eta2 = pars[i]
	} else{
		eta2 = rep(0, length(alpha))
	}

	if(any(substr(Names, 1, 4)=="beta")){
		i = which(substr(Names, 1, 4)=="beta")
		beta = pars[i]
	} else{
		beta = 0
	}

	if(any(substr(Names, 1, 8)=="ghlambda")){
		i = which(substr(Names, 1, 8)=="ghlambda")
		ghlambda = pars[i]
	} else{
		ghlambda = 0
	}

	if(any(substr(Names, 1, 4)=="skew")){
		i = which(substr(Names, 1, 4)=="skew")
		skew = pars[i]
	} else{
		skew = 0
	}

	if(any(substr(Names, 1, 5)=="shape")){
		i = which(substr(Names, 1, 5)=="shape")
		shape = pars[i]
	} else{
		shape = 0
	}
	fk = fm$parameters$fk
	ps = sum(beta) + sum(apply(cbind(alpha, eta1, eta2), 1, FUN=function(x)
						x[1]*fgarchKappa(lambda, delta, x[2], x[3], fk, ghlambda, shape, skew, distribution)))
	return(ps)
}

.persistanstgarch = function(pars, distribution = "norm", submodel){
	Names = names(pars)

	if(any(substr(Names, 1, 6)=="alpha1")){
		i = which(substr(Names, 1, 6)=="alpha1")
		alpha1 = pars[i]
	} else{
		alpha1 = 0
	}

	if(any(substr(Names, 1, 6)=="alpha2")){
		i = which(substr(Names, 1, 6)=="alpha2")
		alpha2 = pars[i]
	} else{
		alpha2 = 0
	}

		if(any(substr(Names, 1, 5)=="beta1")){
		i = which(substr(Names, 1, 5)=="beta1")
		beta1 = pars[i]
	} else{
		beta1 = 0
	}

	if(any(substr(Names, 1, 5)=="beta2")){
		i = which(substr(Names, 1, 5)=="beta2")
		beta2 = pars[i]
	} else{
		beta2 = 0
	}
	ps = 0.5 * (sum(alpha1) + sum(alpha2) + sum(beta1) + sum(beta2))

	return(ps)
}




.persistrealgarch1 = function(pars, idx, distribution = "norm"){
	alpha = pars[idx["alpha",1]:idx["alpha",2]]
	beta  = pars[idx["beta",1]:idx["beta",2]]
	delta = pars[idx["delta",1]:idx["delta",2]]
	ps = sum(beta) + delta*sum(alpha)
	return(ps)
}

.persistrealgarch2 = function(pars, distribution = "norm"){
	Names = names(pars)
	if(any(substr(Names, 1, 5)=="alpha")){
		i = which(substr(Names, 1, 5)=="alpha")
		alpha = pars[i]
	} else{
		alpha = 0
	}

	if(any(substr(Names, 1, 5)=="delta")){
		i = which(substr(Names, 1, 5)=="delta")
		delta = pars[i]
	} else{
		delta = 0
	}

	if(any(substr(Names, 1, 4)=="beta")){
		i = which(substr(Names, 1, 4)=="beta")
		beta = pars[i]
	} else{
		beta = 0
	}
	ps = sum(beta) + delta*sum(alpha)
	return(ps)
}


#----------------------------------------------------------------------------------
# Unconditional Variance
uncvariance = function(object, pars, distribution = "norm", model = "sGARCH",
		submodel = "GARCH", vexdata = NULL)
{
	UseMethod("uncvariance")
}

.unconditional1 = function(object)
{
	# special case check for eGARCH with variance targeting
	pars = object@fit$ipars[,1]
	idx = object@model$pidx
	distribution = object@model$modeldesc$distribution
	vmodel = object@model$modeldesc$vmodel
	vsubmodel = object@model$modeldesc$vsubmodel
	T = object@model$modeldata$T
	vexdata = object@model$modeldata$vexdata[1:T, ,drop = FALSE]
	ans = switch(vmodel,
			sGARCH = .uncsgarch1(pars, idx, distribution, vexdata),
			eGARCH = .uncegarch1(pars, idx, distribution, vexdata),
			gjrGARCH = .uncgjrgarch1(pars, idx, distribution, vexdata),
			apARCH = .uncaparch1(pars, idx, distribution, vexdata),
			fGARCH = .uncfgarch1(pars, idx, distribution, vsubmodel, vexdata),
			csGARCH = .unccsgarch1(pars, idx, distribution, vexdata),
			mcsGARCH = .uncmcsgarch1(pars, idx, distribution, vexdata),
			realGARCH = .uncrealgarch1(pars, idx, distribution, vexdata),
			iGARCH = Inf)
	#names(ans) = "unconditional"
	return(unname(ans))
}

.unconditional2 = function(object)
{
	if(is.null(object@model$fixed.pars))
		stop("\nuncvariance with spec required fixed.pars list\n", call. = FALSE)
	# no other checks for now.
	model = object@model
	pars = unlist(model$fixed.pars)
	parnames = names(pars)
	modelnames = .checkallfixed(object)
	if(is.na(all(match(modelnames, parnames), 1:length(modelnames)))) {
		cat("\nuncmean-->error: parameters names do not match specification\n")
		cat("Expected Parameters are: ")
		cat(paste(modelnames))
		cat("\n")
		stop("Exiting", call. = FALSE)
	}
	# once more into the spec
	setfixed(object)<-as.list(pars)
	model = object@model
	idx = model$pidx
	modelinc = model$modelinc
	pars = object@model$pars[,1]
	distribution = object@model$modeldesc$distribution
	vmodel = object@model$modeldesc$vmodel
	vsubmodel = object@model$modeldesc$vsubmodel
	#T = object@model$modeldata$T
	# assume no out of sample since this cannot be provided in the spec
	vexdata = object@model$modeldata$vexdata
	ans=switch(vmodel,
			sGARCH = 	.uncsgarch1(pars, idx, distribution, vexdata),
			eGARCH = 	.uncegarch1(pars, idx, distribution, vexdata),
			gjrGARCH = 	.uncgjrgarch1(pars, idx, distribution, vexdata),
			apARCH = 	.uncaparch1(pars, idx, distribution, vexdata),
			fGARCH = 	.uncfgarch1(pars, idx, distribution, vsubmodel, vexdata),
			csGARCH = 	.unccsgarch1(pars, idx, distribution, vexdata),
			mcsGARCH = 	.uncmcsgarch1(pars, idx, distribution, vexdata),
			realGARCH = .uncrealgarch1(pars, idx, distribution, vexdata),
			iGARCH = 	Inf)
	#names(ans) = "unconditional"
	return(unname(ans))
}

.unconditional3 = function(object)
{
	pars = object@filter$ipars[,1]
	distribution = object@model$modeldesc$distribution
	vmodel = object@model$modeldesc$vmodel
	vsubmodel = object@model$modeldesc$vsubmodel
	idx = object@model$pidx
	T = object@model$modeldata$T
	vexdata = object@model$modeldata$vexdata[1:T, ,drop = FALSE]
	ans=switch(vmodel,
			sGARCH = 	.uncsgarch1(pars, idx, distribution, vexdata),
			eGARCH = 	.uncegarch1(pars, idx, distribution, vexdata),
			gjrGARCH = 	.uncgjrgarch1(pars, idx, distribution, vexdata),
			apARCH = 	.uncaparch1(pars, idx, distribution, vexdata),
			fGARCH = 	.uncfgarch1(pars, idx, distribution, vsubmodel, vexdata),
			csGARCH = 	.unccsgarch1(pars, idx, distribution, vexdata),
			mcsGARCH = 	.uncmcsgarch1(pars, idx, distribution, vexdata),
			realGARCH = .uncrealgarch1(pars, idx, distribution, vexdata),
			iGARCH = 	Inf)
	#names(ans) = "unconditional"
	return(unname(ans))
}

.unconditional4 = function(pars, distribution = "norm", model = "sGARCH",
		submodel = "GARCH", vexdata = NULL)
{
	ans = switch(model,
			sGARCH = 	.uncsgarch2(pars, distribution, vexdata),
			eGARCH = 	.uncegarch2(pars, distribution, vexdata),
			gjrGARCH = 	.uncgjrgarch2(pars, distribution, vexdata),
			apARCH = 	.uncaparch2(pars, distribution, vexdata),
			fGARCH = 	.uncfgarch2(pars, distribution, submodel, vexdata),
			csGARCH = 	.unccsgarch2(pars, distribution, vexdata),
			mcsGARCH = 	.uncmcsgarch2(pars, distribution, vexdata),
			realGARCH = .uncrealgarch2(pars, distribution, vexdata),
			iGARCH = 	Inf)
	#names(ans) = "unconditional"
	return(unname(ans))
}


setMethod("uncvariance", signature(object = "uGARCHfit", pars = "missing",
				distribution = "missing", model = "missing", submodel = "missing", vexdata = "missing"),
		definition = .unconditional1)

setMethod("uncvariance", signature(object = "uGARCHspec", pars = "missing",
				distribution = "missing", model = "missing", submodel = "missing", vexdata = "missing"),
		definition = .unconditional2)

setMethod("uncvariance", signature(object = "uGARCHfilter", pars = "missing",
				distribution = "missing", model = "missing", submodel = "missing", vexdata = "missing"),
		definition = .unconditional3)

setMethod("uncvariance", signature(object = "missing", pars = "numeric",
				distribution = "character", model = "character", submodel = "ANY", vexdata = "ANY"),
		definition = .unconditional4)


.uncsgarch1 = function(pars, idx, distribution = "norm", vexdata = NULL){
	if(!is.null(vexdata)){
		vxreg = pars[idx["vxreg",1]:idx["vxreg",2]]
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	omega = pars[idx["omega",1]:idx["omega",2]]
	ps = sum(pars[idx["alpha",1]:idx["alpha",2]]) + sum(pars[idx["beta",1]:idx["beta",2]])
	uvol = (umeanvex+omega)/(1-ps)
	return(uvol)
}

.uncsgarch2 = function(pars, distribution = "norm", vexdata = NULL){
	Names = names(pars)
	if(!is.null(vexdata)){
		if(any(substr(Names, 1, 5) == "vxreg")){
			i = which(substr(Names, 1, 5) == "vxreg")
			vxreg = pars[i]
		} else{
			vxreg = 0
		}
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	if(any(substr(Names, 1, 5) == "omega")){
		i = which(substr(Names, 1, 5) == "omega")
		omega = pars[i]
	} else{
		omega = 0
	}
	ps = persistence(pars = pars, distribution = distribution, model = "sGARCH")
	uvol = (umeanvex+omega)/(1-ps)
	return(uvol)
}
.uncgjrgarch1 = function(pars, idx, distribution="norm", vexdata = NULL){
	if(!is.null(vexdata)){
		vxreg = pars[idx["vxreg",1]:idx["vxreg",2]]
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	omega = pars[idx["omega",1]:idx["omega",2]]
	ps = .persistgjrgarch1(pars, idx, distribution)
	uvol = (umeanvex+omega)/(1-ps)
	return(uvol)
}

.uncgjrgarch2 = function(pars, distribution = "norm", vexdata = NULL){
	Names = names(pars)
	if(!is.null(vexdata)){
		if(any(substr(Names, 1, 5) == "vxreg")){
			i = which(substr(Names, 1, 5) == "vxreg")
			vxreg = pars[i]
		} else{
			vxreg = 0
		}
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	if(any(substr(Names, 1, 5) == "omega")){
		i = which(substr(Names, 1, 5) == "omega")
		omega = pars[i]
	} else{
		omega = 0
	}
	ps = persistence(pars = pars, distribution = distribution, model = "gjrGARCH")
	uvol = (umeanvex+omega)/(1-ps)
	return(uvol)
}

.uncegarch1 = function(pars, idx, distribution="norm", vexdata = NULL){
	if(!is.null(vexdata)){
		vxreg = pars[idx["vxreg",1]:idx["vxreg",2]]
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	omega = pars[idx["omega",1]:idx["omega",2]]
	beta = pars[idx["beta",1]:idx["beta",2]]
	#gamma = pars[idx["gamma",1]:idx["gamma",2]]
	#kappa = egarchKappa(pars[idx["ghlambda",1]], pars[idx["shape",1]], ipars[idx["skew",1]], distribution)
	uvol = exp( (umeanvex+omega)/(1-sum(beta)) )
}

.uncegarch2 = function(pars,distribution="norm", vexdata = NULL){
	Names = names(pars)
	if(!is.null(vexdata)){
		if(any(substr(Names, 1, 5) == "vxreg")){
			i = which(substr(Names, 1, 5) == "vxreg")
			vxreg = pars[i]
		} else{
			vxreg = 0
		}
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	if(any(substr(Names, 1, 5)=="omega")){
		i=which(substr(Names, 1, 5)=="omega")
		omega=pars[i]
	} else{
		omega=0
	}
	if(any(substr(Names, 1, 4)=="beta")){
		i=which(substr(Names, 1, 4)=="beta")
		beta=pars[i]
	} else{
		beta=0
	}

	#if(any(substr(Names, 1, 5)=="gamma")){
	#	i=which(substr(Names, 1, 5)=="gamma")
	#	gamma=pars[i]
	#} else{
	#	gamma=0
	#}

	#if(any(substr(Names, 1, 8)=="ghlambda")){
	#	i=which(substr(Names, 1, 8)=="ghlambda")
	#	ghlambda=pars[i]
	#} else{
	#	ghlambda=0
	#}

	#if(any(substr(Names, 1, 5)=="shape")){
	#	i=which(substr(Names, 1, 5)=="shape")
	#	shape=pars[i]
	#} else{
	#	shape=0
	#}

	#if(any(substr(Names, 1, 4)=="skew")){
	#	i=which(substr(Names, 1, 4)=="skew")
	#	skew=pars[i]
	#} else{
	#	skew=0
	#}
	#kappa = egarchKappa(ghlambda, shape, skew, distribution)
	uvol = exp( (umeanvex + omega )/(1-sum(beta)) )
}

.uncaparch1 = function(pars, idx, distribution="norm", vexdata = NULL){
	if(!is.null(vexdata)){
		vxreg = pars[idx["vxreg",1]:idx["vxreg",2]]
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	omega = pars[idx["omega",1]:idx["omega",2]]
	delta = pars[idx["delta",1]:idx["delta",2]]
	ps = .persistaparch1(pars, idx, distribution)
	uvol = ((omega+umeanvex)/(1-ps))^(2/delta)
	return(uvol)
}


.uncaparch2 = function(pars,distribution="norm", vexdata = NULL){
	Names = names(pars)
	if(!is.null(vexdata)){
		if(any(substr(Names, 1, 5) == "vxreg")){
			i = which(substr(Names, 1, 5) == "vxreg")
			vxreg = pars[i]
		} else{
			vxreg = 0
		}
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	if(any(substr(Names, 1, 5)=="omega")){
		i=which(substr(Names, 1, 5)=="omega")
		omega=pars[i]
	} else{
		omega=0
	}
	if(any(substr(Names, 1, 5)=="delta")){
		i=which(substr(Names, 1, 5)=="delta")
		delta=pars[i]
	} else{
		delta=2
	}
	ps=persistence(pars=pars,distribution=distribution,model="apARCH")
	uvol = ((omega+umeanvex)/(1-ps))^(2/delta)
	return(uvol)
}

.uncfgarch1 = function(pars, idx, distribution="norm", submodel, vexdata = NULL){
	if(!is.null(vexdata)){
		vxreg = pars[idx["vxreg",1]:idx["vxreg",2]]
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	omega = pars[idx["omega",1]:idx["omega",2]]
	lambda = pars[idx["lambda",1]:idx["lambda",2]]
	ps = .persistfgarch1(pars, idx, distribution , submodel)
	uvol = ((omega+umeanvex)/(1-ps))^(2/lambda)
	return(uvol)
}

.uncfgarch2 = function(pars, distribution="norm", submodel, vexdata = NULL){
	Names = names(pars)
	if(!is.null(vexdata)){
		if(any(substr(Names, 1, 5) == "vxreg")){
			i = which(substr(Names, 1, 5) == "vxreg")
			vxreg = pars[i]
		} else{
			vxreg = 0
		}
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	fpars = .fgarchModel(submodel)$parameters
	if(any(substr(Names, 1, 6)=="lambda")){
		i=which(substr(Names, 1, 6)=="lambda")
		lambda=pars[i]
	} else{
		lambda = fpars$lambda
	}

	if(any(substr(Names, 1, 5)=="omega")){
		i=which(substr(Names, 1, 5)=="omega")
		omega=pars[i]
	} else{
		omega=0
	}

	ps = persistence(pars = pars, distribution = distribution, model = "fGARCH", submodel = submodel)
	uvol = ((omega+umeanvex)/(1-ps))^(2/lambda)
	return(uvol)
}




.unccsgarch1 = function(pars, idx, distribution = "norm", vexdata = NULL){
	if(!is.null(vexdata)){
		vxreg = pars[idx["vxreg",1]:idx["vxreg",2]]
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	omega = pars[idx["omega",1]:idx["omega",2]]
	ps = pars[idx["eta1",1]]
	uvol = (umeanvex+omega)/(1-ps)
	return(uvol)
}

.unccsgarch2 = function(pars, distribution = "norm", vexdata = NULL){
	Names = names(pars)
	if(!is.null(vexdata)){
		if(any(substr(Names, 1, 5) == "vxreg")){
			i = which(substr(Names, 1, 5) == "vxreg")
			vxreg = pars[i]
		} else{
			vxreg = 0
		}
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	if(any(substr(Names, 1, 5) == "omega")){
		i = which(substr(Names, 1, 5) == "omega")
		omega = pars[i]
	} else{
		omega = 0
	}
	if(any(substr(Names, 1, 5) == "eta11")){
		i = which(substr(Names, 1, 5) == "eta11")
		rho = pars[i]
	} else{
		rho = 1
	}
	uvol = (umeanvex+omega)/(1-rho)
	return(uvol)
}


.uncmcsgarch1 = function(pars, idx, distribution = "norm", vexdata = NULL){
	if(!is.null(vexdata)){
		vxreg = pars[idx["vxreg",1]:idx["vxreg",2]]
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	omega = pars[idx["omega",1]:idx["omega",2]]
	ps = sum(pars[idx["alpha",1]:idx["alpha",2]]) + sum(pars[idx["beta",1]:idx["beta",2]])
	uvol = (umeanvex+omega)/(1-ps)
	return(uvol)
}
.uncmcsgarch2 = function(pars, distribution = "norm", vexdata = NULL){
	Names = names(pars)
	if(!is.null(vexdata)){
		if(any(substr(Names, 1, 5) == "vxreg")){
			i = which(substr(Names, 1, 5) == "vxreg")
			vxreg = pars[i]
		} else{
			vxreg = 0
		}
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	if(any(substr(Names, 1, 5) == "omega")){
		i = which(substr(Names, 1, 5) == "omega")
		omega = pars[i]
	} else{
		omega = 0
	}
	ps = persistence(pars = pars, distribution = distribution, model = "mcsGARCH")
	uvol = (umeanvex+omega)/(1-ps)
	return(uvol)
}
# For the realGARCH model the derivation is based on Proposition 1 of the
# Hansen el al (2011) paper.
.uncrealgarch1 = function(pars, idx, distribution = "norm", vexdata = NULL){
	if(!is.null(vexdata)){
		vxreg = pars[idx["vxreg",1]:idx["vxreg",2]]
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	omega = pars[idx["omega",1]:idx["omega",2]]
	xi = pars[idx["xi",1]:idx["xi",2]]
	alpha = sum(pars[idx["alpha",1]:idx["alpha",2]])
	beta = sum(pars[idx["beta",1]:idx["beta",2]])
	delta = pars[idx["delta",1]:idx["delta",2]]
	ps = delta*alpha + beta
	uvol = (umeanvex+omega+alpha*xi)/(1-ps)
	return(exp(uvol))
}

.uncrealgarch2 = function(pars, distribution = "norm", vexdata = NULL){
	Names = names(pars)
	if(!is.null(vexdata)){
		if(any(substr(Names, 1, 5) == "vxreg")){
			i = which(substr(Names, 1, 5) == "vxreg")
			vxreg = pars[i]
		} else{
			vxreg = 0
		}
		meanvex = apply(vexdata, 2, "mean")
		umeanvex = sum(vxreg*meanvex)
	} else{
		umeanvex = 0
	}
	if(any(substr(Names, 1, 5) == "omega")){
		i = which(substr(Names, 1, 5) == "omega")
		omega = pars[i]
	} else{
		omega = 0
	}
	if(any(substr(Names, 1, 2) == "xi")){
		i = which(substr(Names, 1, 2) == "xi")
		xi = pars[i]
	} else{
		xi = 0
	}
	if(any(substr(Names, 1, 5) == "alpha")){
		i = which(substr(Names, 1, 5) == "alpha")
		alpha = sum(pars[i])
	} else{
		alpha = 0
	}
	if(any(substr(Names, 1, 4) == "beta")){
		i = which(substr(Names, 1, 4) == "beta")
		beta = sum(pars[i])
	} else{
		beta = 0
	}
	if(any(substr(Names, 1, 5) == "delta")){
		i = which(substr(Names, 1, 5) == "delta")
		delta = sum(pars[i])
	} else{
		delta = 0
	}
	uvol = (umeanvex+omega+alpha*xi)/(1-delta*alpha-beta)
	return(exp(uvol))
}

# realized vol (measurement equation) also has an unconditional value (same persistence
# as the vol process)
.uncrealgarchr = function(pars, idx, distribution = "norm"){
	omega = pars[idx["omega",1]:idx["omega",2]]
	xi = pars[idx["xi",1]:idx["xi",2]]
	alpha = sum(pars[idx["alpha",1]:idx["alpha",2]])
	beta = sum(pars[idx["beta",1]:idx["beta",2]])
	delta = pars[idx["delta",1]:idx["delta",2]]
	ps = delta*alpha + beta
	uvol = (delta*omega+(1-beta)*xi)/(1-ps)
	return(uvol)
}

#-------------------------------------------------------------------------
# Unconditional Mean
uncmean = function(object, method = c("analytical", "simulation"), n.sim = 20000, rseed = NA)
{
	UseMethod("uncmean")
}

.unconditionalmean1 = function(object, method = c("analytical", "simulation"), n.sim = 20000, rseed = NA)
{
	method = method[1]
	if(method == "analytical"){
		modelinc = object@model$modelinc
		idx = object@model$pidx
		if(is(object, "uGARCHfilter")){
			h = object@filter$sigma
			pars = object@filter$ipars[,1]
		} else{
			h = object@fit$sigma
			pars = object@fit$ipars[,1]
		}
		N = length(h)

		if(modelinc[6]>0){
			mxreg = matrix( pars[idx["mxreg",1]:idx["mxreg",2]], ncol = modelinc[6] )
			if(modelinc[20]==0){
				mexdata = matrix(object@model$modeldata$mexdata[1:N, ], ncol = modelinc[6])
				meanmex = apply(mexdata, 2, "mean")
				umeanmex = sum(mxreg*meanmex)
			} else{
				if(modelinc[20] == modelinc[6]){
					mexdata = matrix(object@model$modeldata$mexdata[1:N, ], ncol = modelinc[6])
					meanmex = apply(mexdata, 2, "mean")*(uncvariance(object)^(1/2))
					umeanmex = sum(mxreg*meanmex)
				} else{
					mexdata = matrix(object@model$modeldata$mexdata[1:N, ], ncol = modelinc[6])
					meanmex1 = apply(mexdata[,1:(modelinc[6]-modelinc[20]),drop=FALSE], 2, "mean")
					meanmex2 = apply(mexdata[,(modelinc[6]-modelinc[20]+1):modelinc[6],drop=FALSE], 2, "mean")*(uncvariance(object)^(1/2))
					umeanmex = sum(mxreg[,1:(modelinc[6]-modelinc[20])]*meanmex1)+sum(mxreg[,(modelinc[6]-modelinc[20]+1):modelinc[6]]*meanmex2)
				}
			}
		} else{
			umeanmex = 0
		}
		if(modelinc[5]>0){
			# this is obviously an approximation....
			if(modelinc[5] == 2){
				mh = uncvariance(object)
			} else{
				mh = uncvariance(object)^(1/2)
			}
			umeangim = mh * pars[idx["archm",1]]
		} else{
			umeangim = 0
		}
		if(modelinc[1]>0) mu = pars[idx["mu",1]] else mu=0
		umean = (mu + umeangim + umeanmex)
		return(unname(umean))
	 } else{
		 sim = ugarchsim(object, n.sim = n.sim, n.start = 1000, startMethod = "sample", rseed = rseed)
		 umean = mean(sim@simulation$seriesSim[,1])
		return(unname(umean))
	 }
}

.unconditionalmean2 = function(object, method = c("analytical", "simulation"), n.sim = 20000, rseed = NA)
{
	method = method[1]
	if(is.null(object@model$fixed.pars)) stop("uncmean with uGARCHspec requires fixed.pars list", call. = FALSE)
	if(method == "analytical"){
		model = object@model
		pars = unlist(model$fixed.pars)
		parnames = names(pars)
		modelnames = .checkallfixed(object)
		if(is.na(all(match(modelnames, parnames), 1:length(modelnames)))) {
			cat("\nuncmean-->error: parameters names do not match specification\n")
			cat("Expected Parameters are: ")
			cat(paste(modelnames))
			cat("\n")
			stop("Exiting", call. = FALSE)
		}
		# once more into the spec
		setfixed(object)<-as.list(pars)
		model = object@model
		idx = model$pidx
		modelinc = model$modelinc
		pars = object@model$pars[,1]


		if(modelinc[6]>0){
			mxreg = matrix( pars[idx["mxreg",1]:idx["mxreg",2]], ncol = modelinc[6] )
			if(modelinc[20]==0){
				mexdata = matrix(object@model$modeldata$mexdata, ncol = modelinc[6])
				meanmex = apply(mexdata, 2, "mean")
				umeanmex = sum(mxreg*meanmex)
			} else{
				if(modelinc[20] == modelinc[6]){
					mexdata = matrix(object@model$modeldata$mexdata, ncol = modelinc[6])
					meanmex = apply(mexdata, 2, "mean")*(uncvariance(object)^(1/2))
					umeanmex = sum(mxreg*meanmex)
				} else{
					mexdata = matrix(object@model$modeldata$mexdata, ncol = modelinc[6])
					meanmex1 = apply(mexdata[,1:(modelinc[6]-modelinc[20]),drop=FALSE], 2, "mean")
					meanmex2 = apply(mexdata[,(modelinc[6]-modelinc[20]+1):modelinc[6],drop=FALSE], 2, "mean")*(uncvariance(object)^(1/2))
					umeanmex = sum(mxreg[,1:(modelinc[6]-modelinc[20])]*meanmex1)+sum(mxreg[,(modelinc[6]-modelinc[20]+1):modelinc[6]]*meanmex2)
				}
			}
		} else{
			umeanmex = 0
		}

		if(modelinc[5]>0){
			# this is obviously an approximation....
			if(modelinc[5] == 2){
				mh = uncvariance(object)
			} else{
				mh = uncvariance(object)^(1/2)
			}
			umeangim = mh * pars[idx["archm",1]]
		} else{
			umeangim = 0
		}
		if(modelinc[1]>0) mu = pars[idx["mu",1]] else mu=0
		umean = (mu + umeangim + umeanmex)
		return(unname(umean))
	}  else{
		sim = ugarchpath(object, n.sim = n.sim, n.start = 1000, rseed = rseed)
		umean = mean(sim@path$seriesSim[,1])
		return(unname(umean))
	}
}
setMethod("uncmean", signature(object = "uGARCHfit"),    definition = .unconditionalmean1)
setMethod("uncmean", signature(object = "uGARCHfilter"), definition = .unconditionalmean1)
setMethod("uncmean", signature(object = "uGARCHspec"),   definition = .unconditionalmean2)

#----------------------------------------------------------------------------------
# The mult- methods
#----------------------------------------------------------------------------------
multispec = function( speclist )
{
	UseMethod("multispec")
}

setMethod("multispec", signature(speclist = "vector"),  definition = .multispecall)


multifit = function(multispec, data, out.sample = 0, solver = "solnp",
		solver.control = list(), fit.control = list(stationarity = 1, fixed.se = 0, scale = 0,
				rec.init = "all"), cluster = NULL, ...)
{
	UseMethod("multifit")
}

setMethod("multifit", signature(multispec = "uGARCHmultispec"),  definition = .multifitgarch)


multifilter = function(multifitORspec, data = NULL, out.sample = 0, n.old = NULL,
		rec.init = 'all', cluster = NULL, ...)
{
	UseMethod("multifilter")
}

setMethod("multifilter", signature(multifitORspec = "uGARCHmultifit"),  definition = .multifiltergarch1)
setMethod("multifilter", signature(multifitORspec = "uGARCHmultispec"),  definition = .multifiltergarch2)


multiforecast = function(multifitORspec, data = NULL, n.ahead = 1, n.roll = 0, out.sample = 0,
		external.forecasts = list(mregfor = NULL, vregfor = NULL),
		cluster = NULL, ...)
{
	UseMethod("multiforecast")
}

setMethod("multiforecast", signature(multifitORspec = "uGARCHmultifit"),  definition = .multiforecastgarch1)
setMethod("multiforecast", signature(multifitORspec = "uGARCHmultispec"),  definition = .multiforecastgarch2)

#----------------------------------------------------------------------------------
fpm = function( object, summary = TRUE, ...)
{
	UseMethod("fpm")
}


.fpm1 = function(object, summary = TRUE)
{
	n.ahead = object@forecast$n.ahead
	if(n.ahead == 1){
		n.roll = object@forecast$n.roll
		N = object@forecast$N
		ns = object@forecast$n.start
		if(n.roll == 0 | n.roll<4) stop("\nfpm-->error: Forecast Performance Measures require at least 5 out of sample data points (n.roll>3).")
		# get only the forecasts for which out.of.sample data is available
		if(n.roll >= ns){
			forecast = as.numeric(fitted(object))[1:(ns)]
			actual = object@model$modeldata$data[(N - ns + 1):(N - ns + n.roll)]
			rwn = as.character(object@model$modeldata$index[(N - ns + 1):(N - ns + n.roll)])
		} else{
			forecast = as.numeric(fitted(object))[1:(n.roll+1)]
			actual = object@model$modeldata$data[(N - ns + 1):(N - ns + n.roll+1)]
			rwn = as.character(object@model$modeldata$index[(N - ns + 1):(N - ns + n.roll+1)])
		}
		DAC = apply(cbind(actual, forecast), 1, FUN = function(x) as.integer(sign(x[1]) == sign(x[2])))
		if(summary){
			ans = data.frame(MSE = mean((forecast - actual)^2), MAE = mean(abs(forecast - actual)), DAC = mean(DAC))
		} else{
			ans = data.frame(SE = (forecast - actual)^2, AE = abs(forecast - actual), DAC = DAC)
			rownames(ans) = rwn
		}
	} else{
		if(n.ahead<4) stop("\nfpm-->error: Forecast Performance Measures require at least 5 out of sample data points (n.ahead>4).")
		T = object@model$modeldata$T
		Realized = object@model$modeldata$data
		N = NROW(Realized)
		Ser = fitted(object)
		T0 = as.POSIXct(colnames(Ser))
		index = as.POSIXct(as.character(object@model$modeldata$index))
		T1 = sapply(T0, function(x) min(which(x<index)))
		# if last value of T0>index[length(index)] we get Inf
		# but N - Inf + 1 = -Inf < 4 so excluded.
		L = N - T1 + 1
		inc =  which(L>4)
		if(summary){
			n.roll = object@forecast$n.roll
			actual = object@model$modeldata$data
			ans = matrix(NA, ncol = n.roll+1, nrow = 4)
			colnames(ans) = as.character(T0)
			rownames(ans) = c("MSE", "MAE", "DAC", "N")
			for(i in 1:length(inc)){
				ans[1,i] = mean((Ser[1:L[i],i] - Realized[T1[i]:N])^2)
				ans[2,i] = mean(abs((Ser[1:L[i],i] - Realized[T1[i]:N])))
				DAC = apply(cbind(Ser[1:L[i],i], Realized[T1[i]:N]), 1, FUN = function(x) as.integer(sign(x[1]) == sign(x[2])))
				ans[3,i] = mean(DAC)
				ans[4,i] = L[i]
			}
		} else{
			n.roll = object@forecast$n.roll
			actual = object@model$modeldata$data
			sol = vector(mode = "list", length = n.roll+1)
			names(sol) = paste("roll-", seq(0, n.roll), sep="")
			for(i in 1:(n.roll+1)){
				ans = matrix(NA, nrow = L[i], ncol = 3)
				colnames(ans) = c("SE", "AE", "DAC")
				ans[,1] = (Ser[1:L[i],i] - Realized[T1[i]:N])^2
				ans[,2] = abs(Ser[1:L[i],i] - Realized[T1[i]:N])
				ans[,3] = apply(cbind(Ser[1:L[i],i], Realized[T1[i]:N]), 1, FUN = function(x) as.integer(sign(x[1]) == sign(x[2])))
				ans = xts(ans, index[T1[i]:N])
				sol[[i]] = ans
			}
			ans = sol
		}
	}
	return( ans )
}

# VaR loss function used by Gonzalez-Rivera, Lee, and Mishra (2004)
# which can be used with the MCS test of Hansen, Lunde and Nason (2011) to compare
# models

VaRloss = function(alpha, actual,  VaR){
	actual = as.numeric(actual)
	VaR = as.numeric(VaR)
	alpha = as.numeric(alpha[1])
	loss = (alpha - (1 + exp(2500*(actual*100 - VaR*100)))^-1 ) * ( actual*100 - VaR*100)
	return(loss)
}

.fpm2 = function(object, summary = TRUE)
{
	if(!is.null(object@model$noncidx)) stop("\nObject contains non-converged estimation windows.")
	if(summary){
		forecast = object@forecast$density[,1]
		actual = object@forecast$density[,"Realized"]
		DAC = apply(cbind(actual, forecast), 1, FUN = function(x) as.integer(sign(x[1]) == sign(x[2])))
		tmp = c(mean( (forecast - actual)^2 ), mean(abs(forecast - actual)), mean(DAC))
		names(tmp) = c("MSE", "MAE", "DAC")
		tmp = data.frame(Stats = tmp)
	} else{
		forecast = as.numeric(object@forecast$density[,1])
		actual =  as.numeric(object@forecast$density[,"Realized"])
		DAC = apply(cbind(actual, forecast), 1, FUN = function(x) as.integer(sign(x[1]) == sign(x[2])))
		if(object@model$calculate.VaR){
			m = NCOL(object@forecast$VaR)-1
			V = NULL
			for(i in 1:m) V = cbind(V, VaRloss(object@model$VaR.alpha[i], actual, as.numeric(object@forecast$VaR[,i])))
			colnames(V) = paste("VaRLoss(",object@model$VaR.alpha,")",sep="")
			tmp = cbind( as.numeric( (forecast - actual)^2) , abs(forecast - actual), DAC, V)
			colnames(tmp)[1:3] = c("SE", "AE", "HIT")
			rownames(tmp) = rownames(object@forecast$density)
		} else{
			tmp = cbind( as.numeric( (forecast - actual)^2) , as.numeric(abs(forecast - actual)), DAC)
			colnames(tmp) = c("SE", "AE", "HIT")
			rownames(tmp) = rownames(object@forecast$density)
		}
	}
	return( tmp )
}
setMethod("fpm", signature(object = "uGARCHforecast"),  definition = .fpm1)
setMethod("fpm", signature(object = "uGARCHroll"),  definition = .fpm2)

convergence = function(object){
	UseMethod("convergence")
}
.convergence = function(object){
	return( object@fit$convergence )
}
setMethod("convergence", signature(object = "uGARCHfit"),  definition = .convergence)

.convergenceroll = function(object){
	nonc = object@model$noncidx
	if(is.null(nonc)){
		ans = 0
	} else{
		ans = 1
		attr(ans, 'nonconverged')<-nonc
	}
	return(ans)
}
setMethod("convergence", signature(object = "uGARCHroll"),  definition = .convergenceroll)

.vcov = function(object, robust = FALSE){
	if(robust){
		return( object@fit$robust.cvar )
	} else{
		return( object@fit$cvar)
	}
}

setMethod("vcov", signature(object = "uGARCHfit"),  definition = .vcov)


#-------------------------------------------------------------------------------------
.confint.garch <- function(object, parm, level = 0.95, robust=FALSE, ...)
{
	# extract the names of the estimated parameters
	pnames = rownames(object@model$pars[which(object@model$pars[,"Estimate"]==1),])
	cf <- coef(object)
	if(missing(parm)) parm <- pnames
	else if(is.numeric(parm)) parm <- pnames[parm]
	a <- (1 - level)/2
	a <- c(a, 1 - a)
	pct = paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%")
	fac <- qnorm(a)
	ci <- array(NA, dim = c(length(parm), 2L),
			dimnames = list(parm, pct))
	vc = vcov(object)
	colnames(vc) = rownames(vc) <- pnames
	ses <- sqrt(diag(vcov(object, robust)))
  	names(ses) = pnames
  	ses = ses[parm]
	ci[] <- cf[parm] + ses %o% fac
	return(ci)
}

setMethod("confint", signature(object = "uGARCHfit"),  definition = .confint.garch)
