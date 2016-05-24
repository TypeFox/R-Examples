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

#------------------------------------------------------------------------------
# specification
#------------------------------------------------------------------------------
.gogarchspec = function(
		mean.model = list(
				model = c("constant", "AR", "VAR"), robust = FALSE, lag = 1, 
				lag.max = NULL, lag.criterion = c("AIC", "HQ", "SC", "FPE"), 
				external.regressors = NULL, robust.control = list("gamma" = 0.25, 
						"delta" = 0.01, "nc" = 10, "ns" = 500)), 
		variance.model = list(model = "sGARCH", garchOrder = c(1,1), submodel = NULL, 
				variance.targeting = FALSE), distribution.model = c("mvnorm", "manig", "magh"), 
		ica = c("fastica", "radical"), ica.fix = list(A = NULL, K = NULL), ...)
{
	model = list()
	modeldata = list()
	modeldesc = list()
	if(is.null(distribution.model)) distribution.model = "mvnorm"
	vmodel = list(model = "sGARCH", garchOrder = c(1,1), submodel = NULL, variance.targeting = FALSE)
	mmodel = list(model = c("constant"), robust = FALSE, lag = 1, lag.max = NULL, 
			lag.criterion = c("AIC"), external.regressors = NULL, 
			robust.control = list("gamma" = 0.25, "delta" = 0.01, "nc" = 10, "ns" = 500))
	
	distribution = distribution.model[1]
	valid.distributions = c("mvnorm", "manig", "magh")
	if(!any(distribution == valid.distributions)) stop("\nInvalid Distribution Choice\n", call. = FALSE)
	marginal = switch(distribution,
			mvnorm = "norm",
			manig = "nig",
			magh = "ghyp")
	modeldesc$distribution = distribution
	
	mmatch = match(names(mean.model), c("model", "robust", "lag", "lag.max", 
					"lag.criterion", "external.regressors", "robust.control"))
	if(length(mmatch[!is.na(mmatch)]) > 0){
		mx = which(!is.na(mmatch))
		mmodel[mmatch[!is.na(mmatch)]] = mean.model[mx]
	}
	valid.model = c("constant", "AR", "VAR")
	if(!any(mmodel$model[1] == valid.model)) stop("\ngogarchspec-->error: Invalid demean choice\n")
	valid.criterion = c("AIC", "HQ", "SC", "FPE")
	if(!any(mmodel$lag.criterion[1] == valid.criterion)) stop("\ngogarchspec-->error: Invalid VAR criterion choice\n")
	
	vmatch = match(names(variance.model), c("model", "garchOrder", "submodel", "variance.targeting"))
	if(length(vmatch[!is.na(vmatch)]) > 0){
		vx = which(!is.na(vmatch))
		vmodel[vmatch[!is.na(vmatch)]] = variance.model[vx]
	}
	
	modelinc = rep(0, 7)
	names(modelinc) = c("constant", "ar", "var", "mvmxreg", "aux", "aux", "aux")
	
	modelinc[7] = which(valid.distributions == distribution)
	if(mmodel$model == "VAR"){
		if(is.null(mmodel$lag)) modelinc[3] = 1 else modelinc[3] = as.integer( mmodel$lag )
	}
	if(!is.null(mmodel$external.regressors)){
		if(!is.matrix(mmodel$external.regressors)) stop("\nexternal.regressors must be a matrix.")
		modelinc[4] = dim(mmodel$external.regressors)[2]
		modeldata$mexdata = mmodel$external.regressors
	} else{
		modeldata$mexdata = NULL
	}
	
	varmodel = list()
	if( modelinc[3]>0 ){
		varmodel$robust = mmodel$robust
		varmodel$lag.max = mmodel$lag.max
		varmodel$lag.criterion = mmodel$lag.criterion[1]
		varmodel$robust.control = mmodel$robust.control
		
	} else{
		if(mmodel$model == "constant"){
			varmodel$lag.max = 1
			varmodel$lag.criterion = "HQ"
			modelinc[1] = 1
		} else{
			varmodel$lag.max = 1
			varmodel$lag.criterion = "HQ"
			modelinc[2] = mmodel$lag
		}
	}
	
	umodel = list()
	umodel$vmodel = vmodel$model
	umodel$vsubmodel = vmodel$submodel
	umodel$garchOrder = vmodel$garchOrder
	umodel$variance.targeting = vmodel$variance.targeting
	umodel$distribution = marginal
	
	if(is.null(ica)){
		model$ica = "fastica"
	} else{
		valid.models = c("fastica", "radical")
		if(!any(tolower(ica[1]) == valid.models)) stop("\ngogarchspec-->error: Invalid ICA model chosen\n", call. = FALSE)
		model$ica = tolower(ica[1])
	}
	
	if(is.null(ica.fix)) model$ica.fix = list() else model$ica.fix = ica.fix
	
	model$modelinc = modelinc
	model$modeldata = modeldata
	model$varmodel = varmodel
	model$modeldesc = modeldesc
	
	ans = new("goGARCHspec",
			model = model,
			umodel = umodel)
	return(ans)
}

gogarchspec = function(mean.model = list(
				model = c("constant", "AR", "VAR"), robust = FALSE, lag = 1, 
				lag.max = NULL, lag.criterion = c("AIC", "HQ", "SC", "FPE"), 
				external.regressors = NULL, robust.control = list("gamma" = 0.25, 
						"delta" = 0.01, "nc" = 10, "ns" = 500)), 
		variance.model = list(model = "sGARCH", garchOrder = c(1,1), submodel = NULL, 
				variance.targeting = FALSE), 
		distribution.model = c("mvnorm", "manig", "magh"), 
		ica = c("fastica", "radical"), 
		ica.fix = list(A = NULL, K = NULL), ...)
{
	UseMethod("gogarchspec")
}

setMethod("gogarchspec",  definition = .gogarchspec)

#------------------------------------------------------------------------------
# estimation (fit)
#------------------------------------------------------------------------------
gogarchfit = function(spec, data, out.sample = 0, solver = "solnp", 
		fit.control = list(stationarity = 1), solver.control = list(), 
		cluster = NULL, VAR.fit = NULL, ARcoef = NULL, ...)
{
	UseMethod("gogarchfit")
}

setMethod("gogarchfit",  signature(spec = "goGARCHspec"), definition = .gogarchfit)


#------------------------------------------------------------------------------
# filter
#------------------------------------------------------------------------------
gogarchfilter = function(fit, data, out.sample = 0, n.old = NULL, 
		cluster = NULL, ...)
{
	UseMethod("gogarchfilter")
}

setMethod("gogarchfilter",  signature(fit = "goGARCHfit"), definition = .gogarchfilter)


#------------------------------------------------------------------------------
# forecast
#------------------------------------------------------------------------------
gogarchforecast = function(fit, n.ahead = 10, n.roll = 0, 
		external.forecasts = list(mregfor = NULL), cluster = NULL, ...)
{
	UseMethod("gogarchforecast")
}

setMethod("gogarchforecast",  signature(fit = "goGARCHfit"), definition = .gogarchforecast)

#------------------------------------------------------------------------------
# simulation
#------------------------------------------------------------------------------
gogarchsim = function(fit, n.sim = 1, n.start = 0, m.sim = 1, 
		startMethod = c("unconditional", "sample"), prereturns = NA, 
		preresiduals = NA, presigma = NA, mexsimdata = NULL, rseed = NULL, 
		cluster = NULL, ...)
{
	UseMethod("gogarchsim")
}

setMethod("gogarchsim",  signature(fit = "goGARCHfit"), definition = .gogarchsim)

#------------------------------------------------------------------------------
# rolling estimation
#------------------------------------------------------------------------------
gogarchroll = function(spec, data, n.ahead = 1, forecast.length = 50, 
		n.start = NULL, refit.every = 25, refit.window = c("recursive", "moving"), 
		window.size = NULL, solver = "solnp", solver.control = list(), 
		fit.control = list(), rseed = NULL,  cluster = NULL, save.fit = FALSE, 
		save.wdir = NULL, ...)
{
	UseMethod("gogarchroll")
}

setMethod("gogarchroll",  signature(spec = "goGARCHspec"), definition = .rollgogarch) 
#------------------------------------------------------------------------------
# show
#------------------------------------------------------------------------------
setMethod("show",
		signature(object = "goGARCHspec"),
		function(object){
			model = object@umodel$modeldesc$vmodel
			cat(paste("\n*------------------------------*", sep = ""))
			cat(paste("\n*       GO-GARCH Spec          *", sep = ""))
			cat(paste("\n*------------------------------*", sep = ""))
			cat(paste("\n\nMean Model\t\t: ", toupper(names(object@model$modelinc[1:3])[object@model$modelinc[1:3]>0]), sep = ""))
			if(sum(object@model$modelinc[2:3])>0){
				cat(paste("\n(Lag)\t\t\t: ", sum(object@model$modelinc[2:3]), sep = ""))
			}
			if(object@model$modelinc[3]>0){
				cat(paste("\n(Robust)\t\t: ", as.logical(object@model$varmodel$robust), sep = ""))
			}
			cat(paste("\nGARCH Model\t\t: ", object@umodel$vmodel, sep = ""))
			if(object@umodel$vmodel=="fGARCH"){
				cat(paste("\nfGARCH subModel\t\t: ", object@umodel$vsubmodel, sep = ""))
			}
			cat(paste("\nDistribution\t: ", object@model$modeldesc$distribution, sep = ""))
			cat(paste("\nICA Method\t\t: ", object@model$ica, sep = ""))
			cat("\n")
			invisible(object)
		})
		
setMethod("show",
		signature(object = "goGARCHfit"),
		function(object){
			vmodel = object@model$umodel$vmodel
			cat(paste("\n*------------------------------*", sep = ""))
			cat(paste("\n*        GO-GARCH Fit          *", sep = ""))
			cat(paste("\n*------------------------------*", sep = ""))
			cat(paste("\n\nMean Model\t\t: ", toupper(names(object@model$modelinc[1:3])[object@model$modelinc[1:3]>0]), sep = ""))
			if(sum(object@model$modelinc[2:3])>0){
				cat(paste("\n(Lag)\t\t\t: ", sum(object@model$modelinc[2:3]), sep = ""))
			}
			if(object@model$modelinc[3]>0){
				cat(paste("\n(Robust)\t\t: ", as.logical(object@model$varmodel$robust), sep = ""))
			}
			cat(paste("\nGARCH Model\t\t: ", object@model$umodel$vmodel, sep = ""))
			if(object@model$umodel$vmodel=="fGARCH"){
				cat(paste("\nfGARCH subModel\t: ", object@model$umodel$vsubmodel, sep = ""))
			}
			cat(paste("\nDistribution\t: ", object@model$modeldesc$distribution, sep = ""))
			cat(paste("\nICA Method\t\t: ", object@model$ica, sep = ""))
			cat(paste("\nNo. Factors\t\t: ", NCOL(object@mfit$Y), sep=""))
			cat(paste("\nNo. Periods\t\t: ", (object@model$modeldata$T), sep=""))
			cat(paste("\nLog-Likelihood\t: ", round(object@mfit$llh,2),sep=""))
			cat(paste("\n------------------------------------\n",sep=""))
			cat(paste("\nU (rotation matrix) : \n\n",sep=""))
			print(as.matrix(object, which = "U"), digits = 3)
			cat(paste("\nA (mixing matrix) : \n\n",sep=""))
			print(as.matrix(object, which = "A"), digits = 3)
			cat("\n")
			invisible(object)
		})


setMethod("show",
		signature(object = "goGARCHfilter"),
		function(object){
			vmodel = object@model$umodel$vmodel
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\n*        GO-GARCH Filter          *", sep = ""))
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\n\nMean Model\t\t: ", toupper(names(object@model$modelinc[1:3])[object@model$modelinc[1:3]>0]), sep = ""))
			if(sum(object@model$modelinc[2:3])>0){
				cat(paste("\n(Lag)\t\t\t: ", sum(object@model$modelinc[2:3]), sep = ""))
			}
			if(object@model$modelinc[3]>0){
				cat(paste("\n(Robust)\t\t: ", as.logical(object@model$varmodel$robust), sep = ""))
			}
			cat(paste("\nGARCH Model\t\t: ", object@model$umodel$vmodel, sep = ""))
			if(object@model$umodel$vmodel=="fGARCH"){
				cat(paste("\nfGARCH subModel\t\t: ", object@model$umodel$vsubmodel, sep = ""))
			}
			cat(paste("\nDistribution\t: ", object@model$modeldesc$distribution, sep = ""))
			cat(paste("\nICA Method\t\t: ", object@model$ica, sep = ""))
			cat(paste("\nNo. Factors\t\t: ", NCOL(object@mfilter$A), sep=""))
			cat(paste("\nNo. Assets\t\t: ", NROW(object@mfilter$A), sep=""))
			cat(paste("\nLog-Likelihood\t: ", round(likelihood(object),2),sep=""))
			cat(paste("\n------------------------------------\n",sep=""))
			cat(paste("\nU (rotation matrix) : \n\n",sep=""))
			print(as.matrix(object, which = "U"), digits = 3)
			cat(paste("\nA (mixing matrix) : \n\n",sep=""))
			print(as.matrix(object, which = "A"), digits = 3)
			cat("\n")
			invisible(object)
		})

setMethod("show",
		signature(object = "goGARCHforecast"),
		function(object){
			model = object@model
			cat(paste("\n*------------------------------*", sep = ""))
			cat(paste("\n*        GO-GARCH Forecast     *", sep = ""))
			cat(paste("\n*------------------------------*", sep = ""))
			cat(paste("\n\nMean Model\t\t: ", toupper(names(object@model$modelinc[1:3])[object@model$modelinc[1:3]>0]), sep = ""))
			if(sum(object@model$modelinc[2:3])>0){
				cat(paste("\n(Lag)\t\t\t: ", sum(object@model$modelinc[2:3]), sep = ""))
			}
			if(object@model$modelinc[3]>0){
				cat(paste("\n(Robust)\t\t: ", as.logical(object@model$varmodel$robust), sep = ""))
			}
			cat(paste("\nGARCH Model\t\t: ", object@model$umodel$vmodel, sep = ""))
			if(object@model$umodel$vmodel=="fGARCH"){
				cat(paste("\nfGARCH subModel\t\t: ", object@model$umodel$vsubmodel, sep = ""))
			}
			cat(paste("\nDistribution\t: ", object@model$modeldesc$distribution, sep = ""))
			cat(paste("\nICA Method\t\t: ", object@model$ica, sep = ""))
			cat(paste("\nNo. Factors\t\t: ", NCOL(object@mforecast$A), sep=""))
			cat(paste("\nNo. Assets\t\t: ", NROW(object@mforecast$A), sep=""))
			cat(paste("\nn.ahead\t\t\t: ", object@model$n.ahead, sep=""))
			cat(paste("\nn.roll\t\t\t: ", object@model$n.roll, sep=""))
			cat("\n\n")
			invisible(object)
		})


setMethod("show",
		signature(object = "goGARCHsim"),
		function(object){
			model = object@model
			cat(paste("\n*------------------------------*", sep = ""))
			cat(paste("\n*      GO-GARCH Simulation     *", sep = ""))
			cat(paste("\n*------------------------------*", sep = ""))
			cat(paste("\n\nMean Model\t\t: ", toupper(names(object@model$modelinc[1:3])[object@model$modelinc[1:3]>0]), sep = ""))
			if(sum(object@model$modelinc[2:3])>0){
				cat(paste("\n(Lag)\t\t\t: ", sum(object@model$modelinc[2:3]), sep = ""))
			}
			if(object@model$modelinc[3]>0){
				cat(paste("\n(Robust)\t\t: ", as.logical(object@model$varmodel$robust), sep = ""))
			}
			cat(paste("\nGARCH Model\t\t: ", object@model$umodel$vmodel, sep = ""))
			if(object@model$umodel$vmodel=="fGARCH"){
				cat(paste("\nfGARCH subModel\t\t: ", object@model$umodel$vsubmodel, sep = ""))
			}
			cat(paste("\nDistribution\t: ", object@model$modeldesc$distribution, sep = ""))
			cat(paste("\nICA Method\t\t: ", object@model$ica, sep = ""))
			cat(paste("\nNo. Factors\t\t: ", NCOL(object@msim$A), sep=""))
			cat(paste("\nNo. Assets\t\t: ", NROW(object@msim$A), sep=""))
			cat(paste("\nn.sim\t\t\t: ", object@msim$n.sim, sep=""))
			cat(paste("\nm.sim\t\t\t: ", object@msim$m.sim, sep=""))
			cat("\n\n")
			invisible(object)
		})

setMethod("show",
		signature(object = "goGARCHroll"),
		function(object){
			cat(paste("\n*---------------------------------------------------*", sep = ""))
			cat(paste("\n*                    GO-GARCH Roll                  *", sep = ""))
			cat(paste("\n*---------------------------------------------------*", sep = ""))
			cat("\n\nDistribution\t\t:", object@model$modeldesc$distribution)
			cat(paste("\nSimulation Horizon\t: ",  object@model$forecast.length, sep = ""))
			cat(paste("\nRefits\t\t\t\t: ",  object@model$n.refits, sep = ""))
			cat(paste("\nWindow\t\t\t\t: ",  object@model$refit.window, sep = ""))
			cat(paste("\nNo.Assets\t\t\t: ", length(object@model$asset.names), sep = ""))
			cat("\n")
			invisible(object)
		})
#------------------------------------------------------------------------------
# likelihood
#------------------------------------------------------------------------------
.gogarchllh = function(object)
{
	switch(class(object)[1],
			goGARCHfit = object@mfit$llh,
			goGARCHfilter = object@mfilter$llh)
}

setMethod("likelihood", signature(object = "goGARCHfit"), .gogarchllh)
setMethod("likelihood", signature(object = "goGARCHfilter"), .gogarchllh)
#------------------------------------------------------------------------------
# sigma
#------------------------------------------------------------------------------
.gogarchfitsig = function(object, factors = TRUE)
{
	T = object@model$modeldata$T
	D = object@model$modeldata$index[1:T]
	ans = switch(class(object)[1],
			goGARCHfit = object@mfit$factor.sigmas,
			goGARCHfilter = object@mfilter$factor.sigmas)
	if(factors){
		sol = xts(ans[1:T, ], D)
		colnames(sol) = paste("F_",1:NCOL(ans),sep="")
	} else{
		ans = ans[1:T,]^2
		A = as.matrix(object, which = "A")
		sol =  try(.Call("gogarchSigma", S = as.matrix(ans), 
						A = as.matrix(A), PACKAGE = "rmgarch"), silent=TRUE)
		if(inherits(sol, 'try-error')){
			sol = t(apply(ans, 1, function(x) sqrt(diag(A%*%diag(x)%*%t(A)))))
		} else{
			sol = sqrt(sol)
		}
		colnames(sol) = object@model$modeldata$asset.names
		sol = xts(sol[1:T,], D)
	}
	return(sol)
}

setMethod("sigma", signature(object = "goGARCHfit"), .gogarchfitsig)
setMethod("sigma", signature(object = "goGARCHfilter"), .gogarchfitsig)

.gogarchforcsig = function(object, factors = TRUE)
{
	T = object@model$modeldata$T
	m = NCOL(object@model$modeldata$data)
	n.roll = object@model$n.roll
	n.ahead = object@model$n.ahead
	T0 = object@model$modeldata$index[T:(T+n.roll)]
	Sx = object@mforecast$factor.sigmas
	if(factors){
		ans = Sx
		nams = paste("F_",1:dim(ans)[2], sep="")
		dimnames(ans) = list(paste("T+", 1:n.ahead,sep=""), nams, as.character(T0))
	} else{
		A = as.matrix(object, which = "A")
		ma = NCOL(A)
		ans = array(NA, dim=c(n.ahead, m, n.roll+1))
		for(i in 1:(n.roll+1)){
			ans[,,i] = sqrt(.Call("gogarchSigma", S = matrix(Sx[,,i]^2, ncol = ma), 
					A = A, PACKAGE = "rmgarch"))
		}
		nams = object@model$modeldata$asset.names
		dimnames(ans) = list(paste("T+", 1:n.ahead,sep=""), nams, as.character(T0))
	}
	return( ans )
}

setMethod("sigma", signature(object = "goGARCHforecast"), .gogarchforcsig)

.gogarchsimsig = function(object, sim = 1, factors = TRUE)
{
	T = object@model$modeldata$T
	m = NCOL(object@model$modeldata$data)
	n.sim = object@msim$n.sim
	m.sim = object@msim$m.sim
	if(sim > m.sim) stop("\nsim > m.sim!", call. = FALSE)
	ans = object@msim$factor.sigmaSim[[sim]]
	if(factors){
		colnames(ans) = paste("F_",1:dim(ans)[2], sep="")
	} else{
		A = as.matrix(object, which = "A")
		ma = NCOL(A)
		ans = sqrt(.Call("gogarchSigma", S = matrix(ans^2, ncol = ma), 
							A = A, PACKAGE = "rmgarch"))
		colnames(ans) = object@model$modeldata$asset.names
	}
	return( ans )
}

setMethod("sigma", signature(object = "goGARCHsim"), .gogarchsimsig)


# All roll objects will return the actual forecast date, not the T+0 date
# since n.ahead=1 and they are all inside the sample space
.gogarchsigroll = function(object, factors = TRUE)
{
	rind = object@model$rollind
	index = object@model$index
	s = object@model$out.sample
	n = length(object@forecast)
	if(factors) m = NROW(object@model$A[[1]]) else m = NCOL(object@model$A[[1]])
	X = NULL
	D = NULL
	for(i in 1:n){
		D = c(D, as.character(tail(index[rind[[i]]], s[i])))
		X = rbind(X, matrix(sigma(object@forecast[[i]], factors), ncol = m, byrow=TRUE))
	}
	colnames(X) = object@model$asset.names
	X = xts(X, as.POSIXct(D))
	return(X)
}
setMethod("sigma", signature(object = "goGARCHroll"), .gogarchsigroll)

#------------------------------------------------------------------------------
# coef
#------------------------------------------------------------------------------
.gogarchcoef = function(object)
{
	cf = switch(class(object)[1],
			goGARCHfit = object@mfit$garchcoef,
			goGARCHfilter = object@mfilter$garchcoef,
			goGARCHforecast = object@mforecast$garchcoef,
			goGARCHsim = object@msim$garchcoef)
	m = dim(cf)[2]
	colnames(cf) = paste("F_", 1:m, sep = "")
	return(cf)
}

setMethod("coef", signature(object = "goGARCHfit"), .gogarchcoef)
setMethod("coef", signature(object = "goGARCHfilter"), .gogarchcoef)
setMethod("coef", signature(object = "goGARCHforecast"), .gogarchcoef)
setMethod("coef", signature(object = "goGARCHsim"), .gogarchcoef)

.gogarchrollccoef = function(object)
{
	nf = length(object@forecast)
	cf = vector(mode="list", length = nf)
	for(i in 1:nf){
		cf[[i]] = coef(object@forecast[[i]])
	}
	return(cf)
}

setMethod("coef", signature(object = "goGARCHroll"), .gogarchrollccoef)
#------------------------------------------------------------------------------
# fitted
#------------------------------------------------------------------------------
.gogarchfitted = function(object)
{
	ans = switch(class(object)[1],
			goGARCHfit = object@mfit$mu,
			goGARCHfilter = object@mfilter$mu)
	T = object@model$modeldata$T
	D = object@model$modeldata$index[1:T]
	ans = xts(ans[1:T,], D)
	colnames(ans) = object@model$modeldata$asset.names
	return(ans)
}

setMethod("fitted", signature(object = "goGARCHfit"), .gogarchfitted)
setMethod("fitted", signature(object = "goGARCHfilter"), .gogarchfitted)

.gogarchforcfitted = function(object)
{
	T = object@model$modeldata$T
	m = NCOL(object@model$modeldata$data)
	n.roll = object@model$n.roll+1
	n.ahead = object@model$n.ahead
	T0 = object@model$modeldata$index[T:(T+n.roll-1)]
	ans = object@mforecast$mu
	dimnames(ans) = list(paste("T+", 1:n.ahead,sep=""), object@model$modeldata$asset.names, as.character(T0))
	return( ans )
}

setMethod("fitted", signature(object = "goGARCHforecast"), .gogarchforcfitted)

.gogarchfitroll = function(object)
{
	rind = object@model$rollind
	index = object@model$index
	s = object@model$out.sample
	n = length(object@forecast)
	m = NCOL(object@model$A[[1]])
	X = NULL
	D = NULL
	for(i in 1:n){
		D = c(D, as.character(tail(index[rind[[i]]], s[i])))
		X = rbind(X, matrix(fitted(object@forecast[[i]]), ncol = m, byrow=TRUE))
	}
	colnames(X) = object@model$asset.names
	X = xts(X, as.POSIXct(D))
	return(X)
}
setMethod("fitted", signature(object = "goGARCHroll"), .gogarchfitroll)

#------------------------------------------------------------------------------
# residuals
#------------------------------------------------------------------------------
.gogarchresids = function(object)
{
	T = object@model$modeldata$T
	D = object@model$modeldata$index[1:T]
	res = object@model$modeldata$data[1:T, ] - as.matrix(fitted(object))
	res = xts(res, D)
	colnames(res) = object@model$modeldata$asset.names
	return(res)
}

setMethod("residuals", signature(object = "goGARCHfit"), .gogarchresids)
setMethod("residuals", signature(object = "goGARCHfilter"), .gogarchresids)
#------------------------------------------------------------------------------
# ICA matrices
#------------------------------------------------------------------------------
as.matgofit = function(x, rownames.force = NA, which = "A")
{
	validmat = c("A", "W", "U", "K", "Kinv", "E", "D")
	if(!any(tolower(which) == tolower(validmat))) stop("Invalid type for gogarch matrix chosen. Valid choices as A, W, U, K, Kinv, E, and D.", call. = FALSE)
	Z = switch(class(x)[1], 
			goGARCHfit = x@mfit,
			goGARCHfilter = x@mfilter,
			goGARCHforecast = x@mforecast,
			goGARCHsim = x@msim)
	ans = switch(EXPR = which,
			A = Z$A,
			W = Z$W,
			U = Z$U,
			K = Z$K,
			Kinv = Z$Kinv,
			E = Z$E,
			D = Z$D)
	return( as.matrix( ans, rownames.force  = rownames.force ) )
}
setMethod("as.matrix", signature(x = "goGARCHfit"), as.matgofit)
setMethod("as.matrix", signature(x = "goGARCHfilter"), as.matgofit)
setMethod("as.matrix", signature(x = "goGARCHforecast"), as.matgofit)
setMethod("as.matrix", signature(x = "goGARCHsim"), as.matgofit)

#------------------------------------------------------------------------------
# covariance
#------------------------------------------------------------------------------
rcov = function(object, ...)
{
	UseMethod("rcov")
}

.rcovgarch = function(object, type = 1)
{
	A = as.matrix(object, which = "A")
	sig = as.matrix(sigma(object))^2
	n = NROW(sig)
	# for dimensionality reduced systems
	if(type == 1) m = NROW(A) else m = NCOL(A)
	V = array(data = NA, dim = c(m, m, n))
	T = object@model$modeldata$T
	sig = sig[1:T,]
	D = as.character(object@model$modeldata$index[1:T])
	if(type == 1){
		nam = object@model$modeldata$asset.names
		V = try(.Call("gogarchCov", S = sig, A = A, PACKAGE = "rmgarch"), silent=TRUE)
		if(inherits(V, 'try-error')){
			for(i in 1:n)
			{
				tmp =  diag(sig[i,])
				V[,,i] = A%*%tmp%*%t(A)
			}
		}
		dimnames(V) = list(nam, nam, D)
	} else{
		nam = paste("F_",1:NCOL(sig),sep="")
		for(i in 1:n)
		{
			V[,,i] =  diag(sig[i,])
		}
		dimnames(V) = list(nam, nam, D)
	}
	return( V )
}

setMethod("rcov", signature(object = "goGARCHfit"), .rcovgarch)
setMethod("rcov", signature(object = "goGARCHfilter"), .rcovgarch)

.rcovgarchf = function(object, type = 1)
{
	A = as.matrix(object, which = "A")
	sig = sigma(object)^2
	n.roll = object@model$n.roll
	n.ahead = object@model$n.ahead
	m = NCOL(A)
	H = vector(mode = "list", length = n.roll+1)
	T = object@model$modeldata$T
	D = as.character(object@model$modeldata$index[T:(T+n.roll)])
	names(H) = D
	if(type == 1){
		nam = object@model$modeldata$asset.names
		for(i in 1:(n.roll+1)){
			H[[i]] = try(.Call("gogarchCov", S = matrix(sig[,,i], ncol = m), A = A, PACKAGE = "rmgarch"), silent=TRUE)
			dimnames(H[[i]]) = list(nam, nam, paste("T+",1:n.ahead,sep=""))
		}
	} else{
		nam = paste("F_",1:m,sep="")
		for(i in 1:(n.roll+1)){
			V = array(data = NA, dim = c(m, m, n.ahead))
			for(j in 1:n.ahead)
			{
				V[,,j] =  diag(sig[j,,i])
			}
			H[[i]] = V
			dimnames(H[[i]]) = list(nam, nam, paste("T+",1:n.ahead,sep=""))
		}
	}
	#attr(H, "T0") = D
	return( H )
}
setMethod("rcov", signature(object = "goGARCHforecast"), .rcovgarchf)

.rcovgarchsim = function(object, type = 1, sim = 1)
{
	m.sim = object@msim$m.sim
	if(sim > m.sim) stop("\nsim > m.sim!", call. = FALSE)
	A = object@msim$A
	sig = sigma(object, sim = sim)^2
	if(!is.matrix(sig)) sig = matrix(sig, ncol = NCOL(A))	
	n = NROW(sig)
	# for dimensionality reduced systems
	if(type == 1) m = NROW(A) else m = NCOL(A)
	V = array(data = NA, dim = c(m, m, n))
	if(type == 1){
		nam = object@model$modeldata$asset.names
		V = try(.Call("gogarchCov", S = sig, A = A, PACKAGE = "rmgarch"), silent=TRUE)
		dimnames(V) = list(nam, nam, 1:n)
	} else{
		nam = paste("F_",1:NCOL(sig),sep="")
		for(i in 1:n)
		{
			V[,,i] =  diag(sig[i,])
		}
		dimnames(V) = list(nam, nam, 1:n)
	}
	return( V )
}

setMethod("rcov", signature(object = "goGARCHsim"), .rcovgarchsim)


.rcovroll = function(object, type = 1)
{
	rind = object@model$rollind
	index = object@model$index
	s = object@model$out.sample
	n = length(object@forecast)
	nam = object@model$asset.names
	X = NULL
	D = NULL
	X = rcov(object@forecast[[1]], type)
	dm = dim(X[[1]])
	X = array(sapply(X, function(x) x), dim = c(dm[1], dm[1], s[1]))
	D = as.character(tail(index[rind[[1]]], s[1]))
	if(n>1){ 
		for(i in 2:n){
		D = c(D, as.character(tail(index[rind[[i]]], s[i])))
		tmp = array(sapply(rcov(object@forecast[[i]], type), function(x) x),
				dim = c(dm[1], dm[1], s[i]))
		X = .abind(X, tmp)
		}
	}
	dimnames(X) = list(nam, nam, D)
	return(X)
}

setMethod("rcov", signature(object = "goGARCHroll"), .rcovroll)


#------------------------------------------------------------------------------
# correlation
#------------------------------------------------------------------------------
rcor = function(object, ...)
{
	UseMethod("rcor")
}

.rcorgarch = function(object)
{
	A = as.matrix(object, which = "A")
	sig = as.matrix(sigma(object))^2
	n = NROW(sig)
	m = NROW(A)
	T = object@model$modeldata$T
	sig = sig[1:T,]
	D = as.character(object@model$modeldata$index[1:T])
	nam = object@model$modeldata$asset.names
	R = try(.Call("gogarchCor", S = sig, A = A, PACKAGE = "rmgarch"), silent=TRUE)
	dimnames(R) = list(nam, nam, D)
	return( R )
}

setMethod("rcor", signature(object = "goGARCHfit"), .rcorgarch)
setMethod("rcor", signature(object = "goGARCHfilter"), .rcorgarch)

.rcorgarchf = function(object)
{
	A = as.matrix(object, which = "A")
	sig = sigma(object)^2
	n.roll = object@model$n.roll
	n.ahead = object@model$n.ahead
	m = NCOL(A)
	R = vector(mode = "list", length = n.roll+1)
	T = object@model$modeldata$T
	D = as.character(object@model$modeldata$index[T:(T+n.roll)])
	nam = object@model$modeldata$asset.names
	for(i in 1:(n.roll+1)){
		R[[i]] = try(.Call("gogarchCor", S = matrix(sig[,,i], ncol = m), A = A, PACKAGE = "rmgarch"), silent=TRUE)
		dimnames(R[[i]]) = list(nam, nam, paste("T+",1:n.ahead,sep=""))
	}
	# attr(R, "T0") = D
	names(R) = D
	return( R )
}
setMethod("rcor", signature(object = "goGARCHforecast"), .rcorgarchf)


.rcorgarchsim = function(object, sim = 1)
{
	m.sim = object@msim$m.sim
	if(sim > m.sim) stop("\nsim > m.sim!", call. = FALSE)
	A = object@msim$A
	sig = sigma(object, sim = sim)^2
	if(!is.matrix(sig)) sig = matrix(sig, ncol = NCOL(A))	
	n = NROW(sig)
	m = NROW(A)
	nam = object@model$modeldata$asset.names
	R = try(.Call("gogarchCor", S = sig, A = A, PACKAGE = "rmgarch"), silent=TRUE)
	dimnames(R) = list(nam, nam, 1:n)
	return( R )
}

setMethod("rcor", signature(object = "goGARCHsim"), .rcorgarchsim)

.rcorroll = function(object)
{
	rind = object@model$rollind
	index = object@model$index
	s = object@model$out.sample
	n = length(object@forecast)
	nam = object@model$asset.names
	X = NULL
	D = NULL
	X = rcor(object@forecast[[1]])
	dm = dim(X[[1]])
	X = array(sapply(X, function(x) x), dim = c(dm[1], dm[1], s[1]))
	D = as.character(tail(index[rind[[1]]], s[1]))
	if(n>1){ 
		for(i in 2:n){
			D = c(D, as.character(tail(index[rind[[i]]], s[i])))
			tmp = array(sapply(rcor(object@forecast[[i]]), function(x) x),
					dim = c(dm[1], dm[1], s[i]))
			X = .abind(X, tmp)
		}
	}
	dimnames(X) = list(nam, nam, D)
	return(X)
}

setMethod("rcor", signature(object = "goGARCHroll"), .rcorroll)

#------------------------------------------------------------------------------
# coskewness
#------------------------------------------------------------------------------
rcoskew = function(object, ...)
{
	UseMethod("rcoskew")
}
# roll = "all" will not be documented
.rcoskew = function(object, standardize = TRUE, from = 1, to = 1, roll = 0)
{
	A = as.matrix(object, which = "A")
	m = NROW(A)
	if(is.character(roll) && roll=="all" && class(object) == "goGARCHforecast"){
		from = 1
		to = 1
		roll = 1:(object@model$n.roll+1)
		n.roll = object@model$n.roll+1
		D3 = dimnames(sig<-sigma(object))[[3]][roll]
		sig = matrix(sig[from:to,,roll], ncol = NCOL(A), byrow=TRUE)
		mdist = object@model$modeldesc$distribution
		if(mdist == "mvnorm") stop("\nNo co-skewness for mvnorm distribution\n", call. = FALSE)
		xs = array(NA, dim = c(m , m^2, n.roll))
		skew = coef(object)["skew",]
		shape = coef(object)["shape",]
		if(mdist == "magh"){
			ghlambda = coef(object)["ghlambda",]
			S3 = dskewness("ghyp", skew = skew, shape = shape, lambda = ghlambda)
		} else{
			S3 = dskewness("nig", skew = skew, shape = shape)
		}
		# convert to moments since the standardized moments do not retain their 
		# geometric properties in transformation
		yskew = matrix(S3, ncol = NCOL(sig), nrow = NROW(sig), byrow = TRUE)*(sig^3)
		Z = as(t(A), "dgeMatrix")
		rhs = list(Matrix(t(A)), Z)
		for(i in 1:n.roll){
			K = .coskew.ind(yskew[i,])
			lhs = A%*%K
			xs[,,i] = fast_kron_M(rhs, lhs, m, p=2)
		}		
		if(standardize){
			for(i in 1:n.roll){
				SD = sqrt(diag(A%*%diag(sig[i,]^2)%*%t(A)))
				xs[,,i] = xs[,,i]/.coskew.sigma(SD)
			}
		}
		dimnames(xs) = list(NULL, NULL, D3)
	} else{
		if(class(object)=="goGARCHforecast"){
			n.roll = object@model$n.roll
			n.ahead = object@model$n.ahead
			roll = as.integer(roll[1])
			from = as.integer(from[1])
			to   = as.integer(to[1])
			if(roll>n.roll)  stop("\nroll>n.roll!")
			if(from>n.ahead) stop("\nfrom>n.ahead!")
			if(to>n.ahead)   stop("\nto>n.ahead!")
			D3 = dimnames(sig<-sigma(object))[[3]][roll+1]
			D1 = rownames(sig[from:to,,roll+1,drop=FALSE])
			sig = matrix(sig[from:to,,roll+1], ncol = NCOL(A), byrow=TRUE)
		} else{
			D1 = as.character(index(sig<-sigma(object)))[from:to]
			sig = matrix(sigma(object)[from:to,,], ncol = NCOL(A), byrow=TRUE)
		}
		mdist = object@model$modeldesc$distribution
		if(mdist == "mvnorm") stop("\nNo co-skewness for mvnorm distribution\n", call. = FALSE)
		n = length(seq(from, to, by = 1))
		if( n>100 ) stop("\nMaximum from/to is 100 points due to size restrictions\n", .call = FALSE)
		xs = array(NA, dim = c(m , m^2, n))
		skew = coef(object)["skew",]
		shape = coef(object)["shape",]
		if(mdist == "magh"){
			ghlambda = coef(object)["ghlambda",]
			S3 = dskewness("ghyp", skew = skew, shape = shape, lambda = ghlambda)
		} else{
			S3 = dskewness("nig", skew = skew, shape = shape)
		}
		# convert to moments since the standardized moments do not retain their 
		# geometric properties in transformation
		yskew = matrix(S3, ncol = NCOL(sig), nrow = NROW(sig), byrow = TRUE)*(sig^3)
		Z = as(t(A), "dgeMatrix")
		rhs = list(Matrix(t(A)), Z)
		for(i in 1:n){
			K = .coskew.ind(yskew[i,])
			lhs = A%*%K
			xs[,,i] = fast_kron_M(rhs, lhs, m, p=2)
		}
		if(standardize){
			for(i in 1:n){
				SD = sqrt(diag(A%*%diag(sig[i,]^2)%*%t(A)))
				xs[,,i] = xs[,,i]/.coskew.sigma(SD)
			}
		}
		dimnames(xs) = list(NULL, NULL, D1)
		if(class(object)=="goGARCHforecast"){
			attr(xs, "T+0") = D3
		}
	}
	return(xs)
}

setMethod("rcoskew", signature(object = "goGARCHfit"), .rcoskew)
setMethod("rcoskew", signature(object = "goGARCHfilter"), .rcoskew)
setMethod("rcoskew", signature(object = "goGARCHforecast"), .rcoskew)


.rcoskewroll = function(object, standardize = TRUE)
{
	rind = object@model$rollind
	index = object@model$index
	s = object@model$out.sample
	n = length(object@forecast)
	nam = object@model$asset.names
	X = NULL
	D = NULL
	X = rcoskew(object@forecast[[1]], standardize = standardize, roll = "all")
	D = as.character(tail(index[rind[[1]]], s[1]))
	if(n>1){
		for(i in 2:n){
			D = c(D, as.character(tail(index[rind[[i]]], s[i])))
			X = .abind(X, rcoskew(object@forecast[[i]], standardize = standardize, roll = "all"))
		}
	}
	dimnames(X) = list(NULL, NULL, D)
	return(X)
}

setMethod("rcoskew", signature(object = "goGARCHroll"), .rcoskewroll)


.rcoskewsim = function(object, standardize = TRUE, from = 1, to = 1, sim = 1)
{
	mdist = object@model$modeldesc$distribution
	if(mdist == "mvnorm") stop("\nNo co-skewness for mvnorm distribution\n", call. = FALSE)	
	m.sim = object@msim$m.sim
	if(sim > m.sim) stop("\nsim > m.sim!", call. = FALSE)
	n = length(seq(from, to, by = 1))
	if( n>100 ) stop("\nMaximum from/to is 100 points due to size restrictions\n", .call = FALSE)
	A = object@msim$A
	m = NROW(A)
	xs = array(NA, dim = c(m , m^2, n))
	sig = object@msim$factor.sigmaSim[[sim]]
	skew = coef(object)["skew",]
	shape = coef(object)["shape",]
	if(mdist == "magh"){
		ghlambda = coef(object)["ghlambda",]
		S3 = dskewness("ghyp", skew = skew, shape = shape, lambda = ghlambda)
	} else{
		S3 = dskewness("nig", skew = skew, shape = shape)
	}
	# convert to moments as since the standardized moments do not retain their geometric properties in transformation
	sig = sig[from:to,,drop = FALSE]
	yskew = matrix(S3, ncol = NCOL(sig), nrow = NROW(sig), byrow = TRUE)*(sig^3)
	Z = as(t(A), "dgeMatrix")
	rhs = list(Matrix(t(A)), Z)
	for(i in 1:n){
		K = .coskew.ind(yskew[i,])
		lhs = A%*%K
		xs[,,i] = fast_kron_M(rhs, lhs, m, p=2)
	}
	if(standardize){
		for(i in 1:n){
			SD = sqrt(diag(A%*%diag(sig[i,]^2)%*%t(A)))
			xs[,,i] = xs[,,i]/.coskew.sigma(SD)
		}
	}
	return(xs)
}

setMethod("rcoskew", signature(object = "goGARCHsim"), .rcoskewsim)

#------------------------------------------------------------------------------
# cokurtosis
#------------------------------------------------------------------------------
rcokurt = function(object, ...)
{
	UseMethod("rcokurt")
}
# roll="all" will not be documented
.rcokurt = function(object, standardize = TRUE, from = 1, to = 1, roll = 0)
{
	A = as.matrix(object, which = "A")
	m = NROW(A)
	if(is.character(roll) && roll=="all" && class(object) == "goGARCHforecast"){
		from = 1
		to = 1
		roll = 1:(object@model$n.roll+1)
		n.roll = object@model$n.roll+1
		D3 = dimnames(sig<-sigma(object))[[3]][roll]
		sig = matrix(sig[from:to,,roll], ncol = NCOL(A), byrow = TRUE)
		mdist = object@model$modeldesc$distribution
		if(mdist == "mvnorm") stop("\nNo co-skewness for mvnorm distribution\n", call. = FALSE)
		xs = array(NA, dim = c(m , m^3, n.roll))
		skew = coef(object)["skew",]
		shape = coef(object)["shape",]
		if(mdist == "magh"){
			ghlambda = coef(object)["ghlambda",]
			S4 = dkurtosis("ghyp", skew = skew, shape = shape, lambda = ghlambda)+3
		} else{
			S4 = dkurtosis("nig", skew = skew, shape = shape)+3
		}
		# convert to moments since the standardized moments do not retain their 
		# geometric properties in transformation
		ykurt = matrix(S4, ncol = NCOL(sig), nrow = NROW(sig), byrow = TRUE)*(sig^4)
		Z = as(kronecker(t(A), t(A)), "dgeMatrix")
		rhs = list(Matrix(t(A)), Z)
		for(i in 1:n.roll){
			K = .cokurt.ind(sig[i,]^2, ykurt[i,])
			lhs = A%*%K
			xs[,,i] = fast_kron_M(rhs, lhs, m, p=3)
		}
		if(standardize){
			for(i in 1:n.roll){
				SD = sqrt(diag(A%*%diag(sig[i,]^2)%*%t(A)))
				xs[,,i] = xs[,,i]/.cokurt.sigma(SD)
			}
		}
		dimnames(xs) = list(NULL, NULL, D3)
	} else{
		if(class(object)=="goGARCHforecast"){
			n.roll = object@model$n.roll
			n.ahead = object@model$n.ahead
			roll = as.integer(roll[1])
			from = as.integer(from[1])
			to   = as.integer(to[1])
			if(roll>n.roll)  stop("\nroll>n.roll!")
			if(from>n.ahead) stop("\nfrom>n.ahead!")
			if(to>n.ahead)   stop("\nto>n.ahead!")
			D3 = dimnames(sig<-sigma(object))[[3]][roll+1]
			D1 = rownames(sig[from:to,,roll+1,drop=FALSE])
			sig = matrix(sig[from:to,,roll+1], ncol = NCOL(A), byrow=TRUE)
		} else{
			D1 = as.character(index(sig<-sigma(object)))[from:to]
			sig = matrix(sig[from:to,,], ncol = NCOL(A), byrow=TRUE)
		}
		if(from>to)      stop("\nfrom>to!")
		mdist = object@model$modeldesc$distribution
		if(mdist == "mvnorm"){
			stop("\nNo co-kurtosis for mvnorm distribution\n", call. = FALSE)
		}
		n = length(seq(from, to, by = 1))
		if( n>100 ) stop("\nMaximum from/to is 100 points due to size restrictions\n", .call = FALSE)
	
		xs = array(NA, dim = c(m , m^3, n))
		skew = coef(object)["skew",]
		shape = coef(object)["shape",]
		if(mdist == "magh"){
			ghlambda = coef(object)["ghlambda",]
			S4 = dkurtosis("ghyp", skew = skew, shape = shape, lambda = ghlambda)+3
		} else{
			S4 = dkurtosis("nig", skew = skew, shape = shape)+3
		}
		# convert to moments as since the standardized moments do not retain their geometric properties in transformation
		ykurt = matrix(S4, ncol = NCOL(sig), nrow = NROW(sig), byrow = TRUE)*(sig^4)
		Z = as(kronecker(t(A), t(A)), "dgeMatrix")
		rhs = list(Matrix(t(A)), Z)
		for(i in 1:n){
			K = .cokurt.ind(sig[i,]^2, ykurt[i,])
			lhs = A%*%K
			xs[,,i] = fast_kron_M(rhs, lhs, m, p=3)
		}
		if(standardize){
			for(i in 1:n){
				SD = sqrt(diag(A%*%diag(sig[i,]^2)%*%t(A)))
				xs[,,i] = xs[,,i]/.cokurt.sigma(SD)
			}
		}
		dimnames(xs) = list(NULL, NULL, D1)
		if(class(object)=="goGARCHforecast"){
			attr(xs, "T+0") = D3
		}
	}
	return(xs)
}

setMethod("rcokurt", signature(object = "goGARCHfit"), .rcokurt)
setMethod("rcokurt", signature(object = "goGARCHfilter"), .rcokurt)
setMethod("rcokurt", signature(object = "goGARCHforecast"), .rcokurt)

.rcokurtroll = function(object, standardize = TRUE)
{
	rind = object@model$rollind
	index = object@model$index
	s = object@model$out.sample
	n = length(object@forecast)
	nam = object@model$asset.names
	X = NULL
	D = NULL
	X = rcokurt(object@forecast[[1]], standardize = standardize, roll="all")
	D = as.character(tail(index[rind[[1]]], s[1]))
	if(n>1){
		for(i in 2:n){
			D = c(D, as.character(tail(index[rind[[i]]], s[i])))
			X = .abind(X, rcokurt(object@forecast[[i]], standardize = standardize, roll="all"))
		}
	}
	dimnames(X) = list(NULL, NULL, D)
	return(X)
}

setMethod("rcokurt", signature(object = "goGARCHroll"), .rcokurtroll)

.rcokurtsim = function(object, standardize = TRUE, from = 1, to = 1, sim = 1)
{
	mdist = object@model$modeldesc$distribution
	if(mdist == "mvnorm"){
		stop("\nNo co-kurtosis for mvnorm distribution\n", call. = FALSE)
	}
	m.sim = object@msim$m.sim
	if(sim > m.sim) stop("\nsim > m.sim!", call. = FALSE)
	n = length(seq(from, to, by = 1))
	if( n>100 ) stop("\nMaximum from/to is 100 points due to size restrictions\n", .call = FALSE)
	A = object@msim$A
	m = NROW(A)
	xs = array(NA, dim = c(m , m^3, n))
	sig = object@msim$factor.sigmaSim[[sim]]
	skew = coef(object)["skew",]
	shape = coef(object)["shape",]
	if(mdist == "magh"){
		ghlambda = coef(object)["ghlambda",]
		S4 = dkurtosis("ghyp", skew = skew, shape = shape, lambda = ghlambda)+3
	} else{
		S4 = dkurtosis("nig", skew = skew, shape = shape)+3
	}
	# convert to moments as since the standardized moments do not retain their geometric properties in transformation
	sig = sig[from:to,,drop = FALSE]
	ykurt = matrix(S4, ncol = NCOL(sig), nrow = NROW(sig), byrow = TRUE)*(sig^4)
	Z = as(kronecker(t(A), t(A)), "dgeMatrix")
	rhs = list(Matrix(t(A)), Z)
	for(i in 1:n){
		K = .cokurt.ind(sig[i,]^2, ykurt[i,])
		lhs = A%*%K
		xs[,,i] = fast_kron_M(rhs, lhs, m, p=3)
	}	
	if(standardize){
		for(i in 1:n){
			SD = sqrt(diag(A%*%diag(sig[i,]^2)%*%t(A)))
			xs[,,i] = xs[,,i]/.cokurt.sigma(SD)
		}
	}
	return(xs)
}

setMethod("rcokurt", signature(object = "goGARCHsim"), .rcokurtsim)

###############################################################################
#------------------------------------------------------------------------------
# Geometric Portfolio Moments
# 
gportmoments = function(object, ...)
{
	UseMethod("gportmoments")
}

## filter out among the various distributions:
## TODO
.gportmoments = function(object, weights, debug = FALSE)
{
	A = as.matrix(object, which = "A")
	m = NROW(A)
	sig = matrix(sigma(object), ncol = m)
	n = NROW(sig)
	T = object@model$modeldata$T
	index = object@model$modeldata$index
	if(is.vector(weights)){
		if(length(weights) != m) stop("\nIncorrect weight vector length\n", call. = FALSE)
		weights = matrix(weights, ncol = m, nrow = n, byrow = TRUE)
	} else if(is.matrix(weights)){
		nn = NROW(weights)
		mm = NCOL(weights)
		if(mm != m) stop("\nIncorrect column dimension for weights matrix\n", call. = FALSE)
		if(nn != n) stop("\nIncorrect row dimension for weights matrix\n", call. = FALSE)
	}
	if(m < 100){
		pmom = matrix(NA, ncol = 4, nrow = n)
		colnames(pmom) = c("mu", "sigma", "skewness", "kurtosis")
	} else{
		pmom = matrix(NA, ncol = 3, nrow = n)
		colnames(pmom) = c("mu", "sigma", "skewness")
	}
	
	S = rcov(object)
	mu = matrix(fitted(object), ncol = m)
	# write C++ function
	for(i in 1:n){
		w = weights[i,,drop=FALSE]
		pmom[i, 1] = mu[i,]%*%t(w)
		pmom[i, 2] = sqrt(w%*%S[,,i]%*%t(w))
		pmom[i, 3] = (w%*%rcoskew(object, standardize = FALSE, from = i, to = i)[,,1]%*%kronecker(t(w), t(w)))/pmom[i, 2]^3
		if(m< 100) pmom[i, 4] = (w%*%rcokurt(object, standardize = FALSE, from = i, to = i)[,,1]%*%kronecker(t(w), kronecker(t(w), t(w)))/pmom[i, 2]^4)
		if(debug) print(i)
	}
	pmom = xts(pmom, index[1:T])
	return(pmom)
}

setMethod("gportmoments", signature(object = "goGARCHfit"), .gportmoments)
setMethod("gportmoments", signature(object = "goGARCHfilter"), .gportmoments)


.gportmomentsf = function(object, weights, debug = FALSE)
{
	A = as.matrix(object, which = "A")
	m = NROW(A)
	sig = sigma(object)
	n.roll = object@model$n.roll
	n.ahead = object@model$n.ahead
	n = max(n.roll+1, n.ahead)
	T = object@model$modeldata$T
	index = object@model$modeldata$index
	if(is.vector(weights)){
		if(length(weights) != m) stop("\nIncorrect weight vector length\n", call. = FALSE)
		weights = matrix(weights, ncol = m, nrow = n, byrow = TRUE)
	} else if(is.matrix(weights)){
		nn = NROW(weights)
		mm = NCOL(weights)
		if(mm != m) stop("\nIncorrect column dimension for weights matrix\n", call. = FALSE)
		if(nn != n) stop("\nIncorrect row dimension for weights matrix\n", call. = FALSE)
	}
	if(m < 100){
		pmom = array(NA, dim=c(n.ahead, 4, n.roll+1))
		dimnames(pmom) = list(paste("T+",1:n.ahead,sep=""), c("mu", "sigma", "skewness", "kurtosis"),
				as.character(index[T:(T+n.roll)]))
	} else{
		pmom = array(NA, dim=c(n.ahead, 3, n.roll+1))
		dimnames(pmom) = list(paste("T+",1:n.ahead,sep=""), c("mu", "sigma", "skewness"),
				as.character(index[T:(T+n.roll)]))
	}
	S = rcov(object)
	mu = fitted(object)
	# write C++ function
	if(n.ahead==1){
		for(i in 1:(n.roll+1)){
			w = weights[i,,drop=FALSE]
			pmom[1,1,i] = mu[,,i]%*%t(w)
			pmom[1,2,i] = sqrt(w%*%S[[i]][,,1]%*%t(w))
			pmom[1,3,i] = (w%*%rcoskew(object, standardize = FALSE, from = 1, to = 1, roll=i-1)[,,1]%*%kronecker(t(w), t(w)))/pmom[1,2,i]^3
			if(m< 100) pmom[1,4,i] = (w%*%rcokurt(object, standardize = FALSE, from = 1, to = 1, roll=i-1)[,,1]%*%kronecker(t(w), kronecker(t(w), t(w)))/pmom[1,2,i]^4)
			if(debug) print(i)
		}
		
	} else{
		for(i in 1:n.ahead){
			w = weights[i,,drop=FALSE]
			pmom[i,1,1] = mu[i,,1]%*%t(w)
			pmom[i,2,1] = sqrt(w%*%S[[1]][,,i]%*%t(w))
			pmom[i,3,1] = (w%*%rcoskew(object, standardize = FALSE, from = i, to = i, roll=0)[,,1]%*%kronecker(t(w), t(w)))/pmom[i,2,1]^3
			if(m< 100) pmom[i,4,1] = (w%*%rcokurt(object, standardize = FALSE, from = i, to = i, roll=0)[,,1]%*%kronecker(t(w), kronecker(t(w), t(w)))/pmom[i,2,1]^4)
			if(debug) print(i)
		}
	}
	return(pmom)
}
setMethod("gportmoments", signature(object = "goGARCHforecast"), .gportmomentsf)

.gportmoments.garchroll = function(object, weights, debug = FALSE)
{
	rind = object@model$rollind
	index = object@model$index
	s = object@model$out.sample
	n = length(object@forecast)
	nam = object@model$asset.names
	X = NULL
	D = NULL
	D = as.character(tail(index[rind[[1]]], s[1]))
	m = NCOL(object@model$A[[1]])
	T = sum(s)
	if(is.vector(weights)){
		if(length(weights) != m) stop("\nIncorrect weight vector length\n", call. = FALSE)
		weights = matrix(weights, ncol = m, nrow = T, byrow = TRUE)
	} else if(is.matrix(weights)){
		nn = NROW(weights)
		mm = NCOL(weights)
		if(mm != m) stop("\nIncorrect column dimension for weights matrix\n", call. = FALSE)
		if(nn != T) stop("\nIncorrect row dimension for weights matrix\n", call. = FALSE)
	}	
	pmom = matrix(gportmoments(object@forecast[[1]], weights = weights[1:s[1],,drop=FALSE]), ncol = 4, nrow = s[1], byrow=TRUE)
	if(n>1){
		for(i in 2:n){
			D = c(D, as.character(tail(index[rind[[i]]], s[i])))
			tmp = matrix(gportmoments(object@forecast[[i]], weights = weights[(sum(s[1:(i-1)])+1):(sum(s[1:(i)])),,drop=FALSE]), 
					ncol = 4, nrow = s[i], byrow=TRUE)
			pmom = rbind(pmom, tmp)
		}
	}
	dimnames(pmom) = list(D, c("mu", "sigma","skewness","kurtosis"))
	pmom = as.xts(pmom)
	return(pmom)
}

setMethod("gportmoments", signature(object = "goGARCHroll"), .gportmoments.garchroll)


.gportmoments.garchsim = function(object, weights, debug = FALSE, sim = 1)
{
	m.sim = object@msim$m.sim
	if(sim > m.sim) stop("\nsim > m.sim!", call. = FALSE)
	A = object@msim$A
	if(!is.matrix(object@msim$factor.sigmaSim[[sim]])) n = 1 else n = NROW(object@msim$factor.sigmaSim[[sim]])	
	m = NROW(A)
	if(is.vector(weights)){
		if(length(weights) != m) stop("\nIncorrect weight vector length\n", call. = FALSE)
		weights = matrix(weights, ncol = m, nrow = n, byrow = TRUE)
	} else if(is.matrix(weights)){
		nn = NROW(weights)
		mm = NCOL(weights)
		if(mm != m) stop("\nIncorrect column dimension for weights matrix\n", call. = FALSE)
		if(nn != n) stop("\nIncorrect row dimension for weights matrix\n", call. = FALSE)
	}
	if(m < 100) pmom = matrix(NA, ncol = 4, nrow = n) else pmom = matrix(NA, ncol = 3, nrow = n)
	
	S = rcov(object)
	for(i in 1:n){
		w = weights[i,,drop=FALSE]
		pmom[i, 1] = object@msim$seriesSim[[sim]][i,]%*%t(w)
		pmom[i, 2] = sqrt(w%*%S[,,i]%*%t(w))
		pmom[i, 3] = (w%*%rcoskew(object, standardize = FALSE, from = i, to = i)[,,1]%*%kronecker(t(w), t(w)))/pmom[i, 2]^3
		if(m< 100) pmom[i, 4] = (w%*%rcokurt(object, standardize = FALSE, from = i, to = i)[,,1]%*%kronecker(t(w), kronecker(t(w), t(w)))/pmom[i, 2]^4)
		if(debug) print(i)
	}
	return(pmom)
}

setMethod("gportmoments", signature(object = "goGARCHsim"), .gportmoments.garchsim)

###############################################################################
#------------------------------------------------------------------------------
# Numerical Portfolio Moments (FFT based density integration)
nportmoments = function(object, ...)
{
	UseMethod("nportmoments")
}


.nportmoments = function(object, subdivisions = 200, rel.tol = 1e-9, trace = 0, ...)
{
	md = object@dist$dist
	ans = switch(md,
			manig = .nmoments.fft(object, subdivisions = subdivisions, 
					rel.tol = rel.tol, trace = trace, ...),
			magh = .nmoments.fft(object, subdivisions = subdivisions, 
					rel.tol = rel.tol, trace = trace, ...),
			mvnorm = .nmoments.lin(object))
	return(ans)
}

setMethod("nportmoments", signature(object = "goGARCHfft"), .nportmoments)

.nmoments.fft = function(object, subdivisions = 400, rel.tol = 1e-9, trace = 0, ...)
{
	if( object@dist$support.method == "user" ){
		xpdf = seq(object@dist$fft.support[1], object@dist$fft.support[2], by = object@dist$fft.by)
		yd = object@dist$y
		N = dim(yd)[2]
		nmom = matrix(NA, ncol = 4, nrow = N)
		for(i in 1:N){
			dmad = .makeDNew(x = xpdf, dx = yd[,i], h = object@dist$fft.by, standM = "sum")
			nmom[i,1] = object@dist$wM1[i]
			m1 = nmom[i,1]
			if(is.na(m1)){
				nmom[i,2:4] = NA
			} else{
				fm1 = function(x) x * dmad(x)
				# its more accurate to evaluate m1
				m1 = integrate(fm1, -Inf, Inf, subdivisions = subdivisions, rel.tol = rel.tol, stop.on.error = FALSE)$value
				fm2 = function(x) ((x-m1)^2) * dmad(x)
				nmom[i,2] = sqrt(integrate(fm2, -Inf, Inf, subdivisions = subdivisions, rel.tol = rel.tol, stop.on.error = FALSE)$value)
				fm3 = function(x) ((x-m1)^3) * dmad(x)
				nmom[i,3] = (integrate(fm3, -Inf, Inf, subdivisions = subdivisions, rel.tol = rel.tol, stop.on.error = FALSE)$value)
				nmom[i,3] = nmom[i,3]/nmom[i,2]^3
				fm4 = function(x) ((x-m1)^4) * dmad(x)
				nmom[i,4] = (integrate(fm4, -Inf, Inf, subdivisions = subdivisions, rel.tol = rel.tol, stop.on.error = FALSE)$value)
				nmom[i,4] = (nmom[i,4]/nmom[i,2]^4)
			}
			if(trace) print(i)
		}
		colnames(nmom) = c("mu", "sigma", "skewness", "kurtosis")
	} else{
		yd = object@dist$y
		N = length(yd)
		nmom = matrix(NA, ncol = 4, nrow = N)
		
		for(i in 1:N){
			xpdf = seq(object@dist$support.user[i, 1], object@dist$support.user[i, 2], by = object@dist$fft.by)
			dmad = .makeDNew(x = xpdf, dx = yd[[i]], h = object@dist$fft.by, standM = "sum")
			nmom[i,1] = object@dist$wM1[i]
			m1 = nmom[i,1]
			if(is.na(m1)){
				nmom[i,2:4] = NA
			} else{
				fm2 = function(x) ((x-m1)^2) * dmad(x)
				nmom[i,2] = sqrt(integrate(fm2, -Inf, Inf, subdivisions = subdivisions, rel.tol = rel.tol, stop.on.error = FALSE)$value)
				fm3 = function(x) ((x-m1)^3) * dmad(x)
				nmom[i,3] = (integrate(fm3, -Inf, Inf, subdivisions = subdivisions, rel.tol = rel.tol, stop.on.error = FALSE)$value)
				nmom[i,3] = nmom[i,3]/nmom[i,2]^3
				fm4 = function(x) ((x-m1)^4) * dmad(x)
				nmom[i,4] = (integrate(fm4, -Inf, Inf, subdivisions = subdivisions, rel.tol = rel.tol, stop.on.error = FALSE)$value)
				nmom[i,4] = (nmom[i,4]/nmom[i,2]^4)
			}
			if(trace) print(i)
		}
		colnames(nmom) = c("mu", "sigma", "skewness", "kurtosis")
	}
	if(object@model$modtype=="goGARCHfit" | object@model$modtype=="goGARCHfilter"){
		index = object@model$modeldata$index
		T = object@model$modeldata$T
		nmom = xts(nmom, index[1:T])
	}
	if(object@model$modtype == "goGARCHforecast"){
		n.ahead = object@model$n.ahead
		n.roll = object@model$n.roll
		ans = array(NA, dim=c(n.ahead, dim(nmom)[2], n.roll+1))
		index = object@model$modeldata$index
		T = object@model$modeldata$T
		D = as.character(index[T:(T+n.roll)])
		if(n.ahead>1){
			ans[1:n.ahead, ,1] = nmom
		} else{
			for(i in 1:(n.roll+1)) ans[,,i] = nmom[i,]
		}
		dimnames(ans) = list(paste("T+",1:n.ahead,sep=""), colnames(nmom), D)
		nmom = ans
	}
	return(nmom)
}

.nmoments.lin = function(object)
{
	return( object@dist$y )
}

#------------------------------------------------------------------------------
###############################################################################


###############################################################################
#------------------------------------------------------------------------------
# Convolution
convolution = function(object, ...)
{
	UseMethod("convolution")
}


.convolution = function(object, weights, fft.step = 0.001, 
		fft.by = 0.0001, fft.support = c(-1, 1), support.method = c("user", "adaptive"), 
		use.ff = TRUE, cluster = NULL, trace = 0, ...)
{
	dist = object@model$modeldesc$distribution
	if(dist == "manig"){
		cf = cbind(coef(object)["skew", ], coef(object)["shape",])
		pars = .nigtransform(rho = cf[,1], zeta = cf[,2])
	} else if(dist == "magh"){
		cf = cbind(coef(object)["skew", ], coef(object)["shape",], coef(object)["ghlambda", ])
		pars = .ghyptransform(rho = cf[,1], zeta = cf[,2], lambda = cf[,3])
		# colnames(pars) = c("alpha", "beta", "delta", "mu", "lambda")
	} else{
		pars = NULL
	}
	M2 = rcov(object, type = 2)
	M1 = as.matrix(fitted(object))
	A  = as.matrix(object, which = "A")
	support.method = support.method[1]
	model = object@model
	modtype  = class(object)[1]
	model$modtype = modtype
	debug = as.logical(trace)
	ans = switch(dist,
			manig = .convolve.nig(nigpars = pars, M2 = M2, M1 = M1, A = A, 
					weights = weights, fft.step = fft.step, fft.by = fft.by, 
					fft.support = fft.support, use.ff = use.ff, 
					support.method = support.method, cluster = cluster, 
					model = model, debug = debug),
			magh = .convolve.gh(ghpars = pars, M2 = M2, M1 = M1, A = A, 
					weights = weights, fft.step = fft.step, fft.by = fft.by, 
					fft.support = fft.support, use.ff = use.ff, 
					support.method = support.method, cluster = cluster, 
					model = model, debug = debug),
			mvnorm = .convolve.norm(M2 = M2, M1 = M1, A = A, weights = weights,
					model = model))
	return(ans)
}
setMethod("convolution", signature(object = "goGARCHfit"), .convolution)
setMethod("convolution", signature(object = "goGARCHfilter"), .convolution)


.convolutionf = function(object, weights, fft.step = 0.001, 
		fft.by = 0.0001, fft.support = c(-1, 1), support.method = c("user", "adaptive"), 
		use.ff = TRUE, cluster = NULL, trace = 0, ...)
{
	dist = object@model$modeldesc$distribution
	if(dist == "manig"){
		cf = cbind(coef(object)["skew", ], coef(object)["shape",])
		pars = .nigtransform(rho = cf[,1], zeta = cf[,2])
	} else if(dist == "magh"){
		cf = cbind(coef(object)["skew", ], coef(object)["shape",], coef(object)["ghlambda", ])
		pars = .ghyptransform(rho = cf[,1], zeta = cf[,2], lambda = cf[,3])
		# colnames(pars) = c("alpha", "beta", "delta", "mu", "lambda")
	} else{
		pars = NULL
	}
	# either n.ahead == 1 || n.ahead>1
	# in either case flatten inputs and transform on exit to other functions
	n.ahead = object@model$n.ahead
	n.roll  = object@model$n.roll
	M2 = rcov(object, type = 2)
	M1 = fitted(object)
	A  = as.matrix(object, which = "A")
	if(n.ahead>1){
		M2 = M2[[1]]
		m  = dim(M2)[2]
		M1 = as.matrix(M1[,,1])
	} else{
		m = dim(M2[[1]])[2]
		M2 = array(sapply(M2, function(x) x[,,1]), dim=c(m, m, n.roll+1))
		M1 = matrix(M1, ncol = NROW(A), nrow = n.roll+1, byrow=TRUE)
	}
	support.method = support.method[1]
	model = object@model
	modtype  = class(object)[1]
	model$modtype = modtype
	debug = as.logical(trace)
	ans = switch(dist,
			manig = .convolve.nig(nigpars = pars, M2 = M2, M1 = M1, A = A, 
					weights = weights, fft.step = fft.step, fft.by = fft.by, 
					fft.support = fft.support, use.ff = use.ff, 
					support.method = support.method, cluster = cluster, 
					model = model, debug = debug),
			magh = .convolve.gh(ghpars = pars, M2 = M2, M1 = M1, A = A, 
					weights = weights, fft.step = fft.step, fft.by = fft.by, 
					fft.support = fft.support, use.ff = use.ff, 
					support.method = support.method, cluster = cluster, 
					model = model, debug = debug),
			mvnorm = .convolve.norm(M2 = M2, M1 = M1, A = A, weights = weights,
					model = model))
	return(ans)
}
setMethod("convolution", signature(object = "goGARCHforecast"), .convolutionf)


.convolution.gogarchsim = function(object, weights, fft.step = 0.001, 
		fft.by = 0.0001, fft.support = c(-1, 1), support.method = c("user", "adaptive"), 
		use.ff = TRUE, sim = 1, cluster = NULL, trace = 0, ...)
{
	dist = object@model$modeldesc$distribution
	if(dist == "manig"){
		cf = cbind(coef(object)["skew", ], coef(object)["shape",])
		pars = .nigtransform(rho = cf[,1], zeta = cf[,2])
	} else if(dist == "magh"){
		cf = cbind(coef(object)["skew", ], coef(object)["shape",], coef(object)["ghlambda", ])
		pars = .ghyptransform(rho = cf[,1], zeta = cf[,2], lambda = cf[,3])
		# colnames(pars) = c("alpha", "beta", "delta", "mu", "lambda")
	} else{
		pars = NULL
	}
	M2 = rcov(object, type = 2, sim = sim)
	M1	= object@msim$seriesSim[[sim]]
	A = object@msim$A
	support.method = support.method[1]
	model = object@model
	modtype  = class(object)[1]
	model$modtype = modtype
	debug = as.logical(trace)
	ans = switch(dist,
			manig = .convolve.nig(nigpars = pars, M2 = M2, M1 = M1, A = A, 
					weights = weights, fft.step = fft.step, fft.by = fft.by, 
					fft.support = fft.support, use.ff = use.ff, 
					support.method = support.method, cluster = cluster, 
					model = model, debug = debug),
			magh = .convolve.gh(ghpars = pars, M2 = M2, M1 = M1, A = A, 
					weights = weights, fft.step = fft.step, fft.by = fft.by, 
					fft.support = fft.support, use.ff = use.ff, 
					support.method = support.method, cluster = cluster, 
					model = model, debug = debug),
			mvnorm = .convolve.norm(M2 = M2, M1 = M1, A = A, weights = weights, model = model))
	return(ans)
}

setMethod("convolution", signature(object = "goGARCHsim"), .convolution.gogarchsim)

.convolution.gogarchroll = function(object, weights, fft.step = 0.001, 
		fft.by = 0.0001, fft.support = c(-1, 1), support.method = c("user", "adaptive"), 
		use.ff = TRUE, cluster = NULL, trace = 0, ...)
{
	dist = object@model$modeldesc$distribution
	support.method = support.method[1]
	model = object@model
	modtype  = class(object)[1]
	model$modtype = modtype
	model$n.ahead = 1
	debug = as.logical(trace)
	ans = switch(dist,
			mvnorm = .convolve.roll.mn(object, weights, model, ...),
			manig = .convolve.roll.nig(object, weights, fft.step, fft.by, 
					fft.support, support.method, use.ff, cluster = cluster, 
					model, trace, ...),
			magh = .convolve.roll.nig(object, weights, fft.step, fft.by, 
					fft.support, support.method, use.ff, 
					cluster = cluster, model, trace, ...))
	return(ans)
}

.convolve.roll.mn = function(object, weights, model, ...)
{
	dist = object@model$modeldesc$distribution
	fseq = object@model$out.sample
	forecast.length = sum(fseq)
	m = NROW(object@forecast[[1]]@mforecast$A)
	
	if(is.matrix(weights)){
		if(dim(weights)[2] != m) stop("\nconvolution-->error: Weights column dimension not equal to data column dimension!")
		if(dim(weights)[1] != forecast.length) stop("\nconvolution-->error: Weights row dimension must equal floor(forecast.length/refit.every) * refit.every")
	} else{
		weights = matrix(weights, ncol = m, nrow = forecast.length, byrow = TRUE)
	}
	nf = length(object@forecast)
	idx = matrix(NA, ncol = 2, nrow = nf)
	cfseq = cumsum(fseq)
	idx[,1] = c(1, cfseq[-length(cfseq)]+1)
	idx[,2] = cfseq
	rpdf = matrix(NA, ncol = 2, nrow = forecast.length)
	
	for(i in 1:nf){
		tmp = .convolution(object@forecast[[i]], weights = weights[idx[i,1]:idx[i,2],,drop = FALSE])
		rpdf[idx[i,1]:idx[i,2], ] = tmp@dist$y
	}
	wM1 = rep(NA, NROW(weights))
	fx = as.matrix(fitted(object))
	for(i in 1:NROW(weights)){
		wM1[i] = fx[i, ] %*% weights[i, ]
	}
	sol = list()
	sol$wM1 = wM1
	sol$dist = "mvnorm"
	sol$y = rpdf
	sol$fft.support = NULL
	sol$fft.step = NULL
	sol$fft.by = NULL
	# we want to pass the model spec but without too much added weight
	model$modeldata$data = NULL
	model$modeldata$mexdata = NULL
	model$modeldata$vexdata = NULL
	ans = new("goGARCHfft",
			dist = sol,
			model = model)
	return(ans)
}

.convolve.roll.nig = function(object, weights, fft.step, fft.by, 
		fft.support, support.method = c("user", "adaptive"), 
		use.ff = TRUE, cluster = NULL, model, trace = 0, ...)
{
	dist = object@model$modeldesc$distribution
	fseq = object@model$out.sample
	forecast.length = sum(fseq)
	m = NROW(object@forecast[[1]]@mforecast$A)
	
	if(is.matrix(weights)){
		if(dim(weights)[2] != m) stop("\nconvolution-->error: Weights column dimension not equal to data column dimension!")
		if(dim(weights)[1] != forecast.length) stop("\nconvolution-->error: Weights row dimension must equal floor(forecast.length/refit.every) * refit.every")
	} else{
		weights = matrix(weights, ncol = m, nrow = forecast.length, byrow = TRUE)
	}
	nf = length(object@forecast)
	
	idx = matrix(NA, ncol = 2, nrow = nf)
	cfseq = cumsum(fseq)
	idx[,1] = c(1, cfseq[-length(cfseq)]+1)
	idx[,2] = cfseq
	
	support.user = NULL
	
	if( support.method == "user" ){
		zz = seq(fft.support[1], fft.support[2], by = fft.by)
		if(use.ff){
			rpdf = ff(vmode="double", dim=c(length(zz), forecast.length)) 
		} else{
			rpdf = matrix(NA, nrow = length(zz), ncol = forecast.length)
		}
		
		for(i in 1:nf){
			tmp = .convolutionf(object@forecast[[i]], 
					weights = weights[idx[i,1]:idx[i,2],,drop = FALSE], 
					fft.step = fft.step, fft.by = fft.by, fft.support = fft.support, 
					support.method = support.method, use.ff = use.ff, 
					cluster = cluster, trace = trace)
			rpdf[1:length(zz), idx[i,1]:idx[i,2]] = tmp@dist$y[1:length(zz), 1:fseq[i]]
		}
	} else{
		rpdf = NULL
		for(i in 1:nf){
			tmp = .convolutionf(object@forecast[[i]], 
					weights =  weights[idx[i,1]:idx[i,2],,drop = FALSE], 
					fft.step = fft.step, fft.by = fft.by, fft.support = fft.support, 
					support.method = support.method, use.ff = use.ff, 
					cluster = cluster, trace = trace)
			rpdf = c(rpdf, tmp@dist$y)
			support.user = rbind(support.user, tmp@dist$support.user)
		}
	}
	wM1 = rep(NA, NROW(weights))
	fx = as.matrix(fitted(object))
	for(i in 1:NROW(weights)){
		wM1[i] = fx[i, ] %*% weights[i, ]
	}
	sol = list()
	sol$wM1 = wM1
	sol$dist = dist
	sol$y = rpdf
	sol$support.method = support.method
	sol$support.user = support.user
	sol$fft.support = fft.support
	sol$fft.step = fft.step
	sol$fft.by = fft.by
	model$modeldata$data = NULL
	model$modeldata$mexdata = NULL
	model$modeldata$vexdata = NULL
	ans = new("goGARCHfft",
			dist = sol,
			model = model)
	return(ans)
}

setMethod("convolution", signature(object = "goGARCHroll"), .convolution.gogarchroll)


# maNIG convolution of independent densities
.convolve.nig = function(nigpars, M2, M1, A, weights, fft.step = 0.001, 
		fft.by = 0.0001, fft.support = c(-1, 1), support.method = c("user", "adaptive"), 
		use.ff = TRUE, cluster = NULL, model, debug = 0)
{
	sol = list()
	n = dim(M2)[3]
	m = NCOL(A)
	if(!is.matrix(weights)) weights = matrix(weights[1:NROW(A)], ncol = NROW(A), nrow = n, byrow = TRUE)
	if(NCOL(weights) != NROW(A))
		stop("\nconvolution-->error: wrong no. of columns for weights matrix.", call. = FALSE)
	if(NROW(weights) != n) 
		stop("\nconvolution-->error: wrong no. of rows for weights matrix.", call. = FALSE)
	# The weighting vector for the distribution
	w.hat = matrix(NA, ncol = m, nrow = n)
	# also return the weighted first moment
	wM1 = rep(NA, n)
	for(i in 1:n){
		dS = sqrt(M2[,,i])
		w.hat[i,] = weights[i,] %*% (A %*% dS)
		wM1[i] = M1[i, ] %*% weights[i, ]
	}
	# weights[i,] * (A %*% diag(dS))
	# zres = sapply(fit@mfit$ufit@fit, FUN = function(x) residuals(x, TRUE))
	# port = rep(0,n)
	# port2 = rep(0, n)
	# XX = tail(fit@model$modeldata$data[1:fit@model$modeldata$T,],n)
	# for(i in 1:n)  port[i] = weights[i,]%*%M1[i,] + weights[i,]%*%(A %*% sqrt(M2[,,i]) %*% (zres[i,]))
	# for(i in 1:n)  port2[i] = weights[i,]%*%XX[i,]
	# all.equal(port, port2)
	wpars = array(data = NA, dim = c(m, 4, n))
	# Scaling of parameters (Blaesild)
	for(j in 1:m){
		tmp = matrix(nigpars[j,1:4], ncol = 4, nrow = n, byrow = TRUE)
		tmp = tmp*cbind(1/abs(w.hat[1:n,j]),1/w.hat[1:n,j],abs(w.hat[1:n,j]),w.hat[1:n,j]) + cbind(matrix(0, nrow = n, ncol = 3), M1[,j] * weights[,j])
		wpars[j,,] = as.array(t(tmp), dim = c(1, 4, n))
	}
	if( support.method == "user" ){
		support.user = NULL
		zz = seq(fft.support[1], fft.support[2], by = fft.by)
		zzLength = length(zz)
		if(use.ff){
			nigpdf = ff(vmode="double", dim = c(zzLength, n))
		} else{
			nigpdf = matrix(NA, nrow = zzLength, ncol = n)
		}
		if( !is.null(cluster) ){
			clusterEvalQ(cluster, require(rmgarch))
			clusterExport(cluster, c("zz", "fft.step", "wpars","cfinv","nigmvcf"), 
					envir = environment())
			xpdf = parLapply(cluster, as.list(1:n), fun = function(i){
						cfinv(z = zz, f = nigmvcf, 
								step = fft.step, alpha = wpars[,1,i], 
								beta = wpars[,2,i], delta = wpars[,3,i], 
								mu = wpars[,4,i])
					})
			for(i in 1:n) nigpdf[,i] = xpdf[[i]]
			rm(xpdf)
		} else{	
			for(i in 1:n){
				nigpdf[,i] = cfinv(z = zz, f = nigmvcf, step = fft.step, 
						alpha = wpars[,1,i], beta = wpars[,2,i], 
						delta = wpars[,3,i], mu = wpars[,4,i])
				if(debug) print(i)
			}
		}
	} else{
		nigpdf = vector(mode = "list", length = n)
		support.user = matrix(NA, ncol = 2, nrow = n)
		if( !is.null(cluster) ){
			clusterEvalQ(cluster, require(rmgarch))
			clusterExport(cluster, c("fft.step", "fft.by", "wpars","cfinv","nigmvcf"), envir = environment())
			xsol = parLapply(cluster, as.list(1:n), fun = function(i){
						xmin = min(apply(cbind(wpars[,1,i], wpars[,2,i], wpars[,3,i], wpars[,4,i]), 1, 
										FUN = function(x) rugarch:::qnig(0.0000001, 
													alpha = x[1], beta = x[2], 
													delta = x[3], mu = x[4])))
						xmax = max(apply(cbind(wpars[,1,i], wpars[,2,i], wpars[,3,i], wpars[,4,i]), 1, 
										FUN = function(x) rugarch:::qnig(1-0.0000001, 
													alpha = x[1], beta = x[2], 
													delta = x[3], mu = x[4])))
						zz = seq(xmin, xmax, by = fft.by)
						ans = cfinv(z = zz, f = nigmvcf, 
								step = fft.step, alpha = wpars[,1,i], 
								beta = wpars[,2,i], delta = wpars[,3,i], 
								mu = wpars[,4,i])
						sol = list(d = ans, support = c(xmin, xmax))
						return(sol)
					})
			for(i in 1:n){
				nigpdf[[i]] = xsol[[i]]$d
				support.user[i, ] = xsol[[i]]$support
			}
		} else{
			for(i in 1:n){
				xmin = min(apply(cbind(wpars[,1,i], wpars[,2,i], wpars[,3,i], wpars[,4,i]), 1, 
								FUN = function(x) rugarch:::qnig(0.0000001, alpha = x[1], 
											beta = x[2], delta = x[3], mu = x[4])))
				xmax = max(apply(cbind(wpars[,1,i], wpars[,2,i], wpars[,3,i], wpars[,4,i]), 1, 
								FUN = function(x) rugarch:::qnig(1-0.0000001, alpha = x[1], 
											beta = x[2], delta = x[3], mu = x[4])))
				zz = seq(xmin, xmax, by = fft.by)
				nigpdf[[i]] = cfinv(z = zz, f = nigmvcf, step = fft.step, 
						alpha = wpars[,1,i], beta = wpars[,2,i], 
						delta = wpars[,3,i], mu = wpars[,4,i])
				support.user[i, ] = c(xmin, xmax)
				if(debug) print(i)
			}
		}
	}
	# i.e. for each data point(n points) there is a density of size zzLength
	#nigpdf = big.matrix(nrow = zzLength, ncol = n, type = "double", init = NA, dimnames = NULL)
	sol$dist = "manig"
	sol$y = nigpdf
	sol$wM1 = wM1
	sol$support.method = support.method
	sol$support.user = support.user
	sol$fft.support = fft.support
	sol$fft.step = fft.step
	sol$fft.by = fft.by
	model$modeldata$data = NULL
	model$modeldata$mexdata = NULL
	model$modeldata$vexdata = NULL
	ans = new("goGARCHfft",
			dist = sol,
			model = model)
	return(ans)
}


.convolve.gh = function(ghpars, M2, M1, A, weights, fft.step = 0.001, fft.by = 0.0001, 
		fft.support = c(-1, 1), support.method = c("user", "adaptive"), 
		use.ff = TRUE, cluster = NULL, model, debug = 0)
{
	sol = list()
	n = dim(M2)[3]
	m = NCOL(A)
	if(!is.matrix(weights)) weights = matrix(weights, ncol = NROW(A), nrow = n, byrow = TRUE)
	if(NCOL(weights) != NROW(A)) 
		stop("\nconvolution-->error: wrong no. of columns for weights matrix.", call. = FALSE)
	if(NROW(weights) != n) 
		stop("\nconvolution-->error: wrong no. of rows for weights matrix.", call. = FALSE)
	# The weighting vector for the distribution
	w.hat = matrix(NA, ncol = m, nrow = n)
	# also return the weighted first moment
	wM1 = rep(NA, n)
	for(i in 1:n){
		dS = sqrt(M2[,,i])
		w.hat[i,] = weights[i,] %*% (A %*% dS)
		wM1[i] = M1[i, ] %*% weights[i, ]
	}
	# weights[i,] * (A %*% diag(dS))
	# port = rep(0,n)
	# for(i in 1:n)  port[i] = weights[i,]%*%Mu[i,] + weights[i,]%*%(t(A) %*% sqrt(Sigma[,,i]) %*% (zres[i,]))
	# Mu[i,] + (A %*% dS %*% (zres[i,]))
	wpars = array(data = NA, dim = c(m, 5, n))
	
	# Scaling of parameters (Blaesild proof)
	for(j in 1:m){
		tmp = matrix(ghpars[j,1:5], ncol = 5, nrow = n, byrow = TRUE)
		tmp = tmp*cbind(1/abs(w.hat[1:n,j]), 1/w.hat[1:n,j], abs(w.hat[1:n,j]), w.hat[1:n,j], rep(1,n)) + 
				cbind(matrix(0, ncol = 3, nrow = n), M1[,j] * weights[,j], rep(0, n))
		wpars[j,,] = as.array(t(tmp), dim = c(1, 5, n))
	}
	
	
	if( support.method == "user" ){
		support.user = NULL
		zz = seq(fft.support[1], fft.support[2], by = fft.by)
		zzLength = length(zz)
		if(use.ff){
			ghpdf = ff(vmode="double", dim = c(zzLength, n))
		} else{
			ghpdf = matrix(NA, nrow = zzLength, ncol = n)
		}
		if( !is.null(cluster) ){
			clusterEvalQ(cluster, require(rmgarch))
			clusterExport(cluster, c("zz", "fft.step", "wpars","cfinv","ghypmvcf"), 
					envir = environment())
			xpdf = parLapply(cluster, 1:n, fun = function(i){
						cfinv(z = zz, f = ghypmvcf, 
								step = fft.step, lambda = wpars[,5,i], 
								alpha = wpars[,1,i], beta = wpars[,2,i], 
								delta = wpars[,3,i], mu = wpars[,4,i])
					})
			for(i in 1:n) ghpdf[,i] = xpdf[[i]]
			rm(xpdf)
		} else{
			for(i in 1:n){
				ghpdf[,i] = cfinv(z = zz, f = ghypmvcf, step = fft.step, lambda = wpars[,5,i], 
						alpha = wpars[,1,i], beta = wpars[,2,i], delta = wpars[,3,i], mu = wpars[,4,i])
				if(debug) print(i)
			}
		}	
	} else{
		ghpdf = vector(mode = "list", length = n)
		support.user = matrix(NA, ncol = 2, nrow = n)
		if( !is.null(cluster) ){
			clusterEvalQ(cluster, require(rmgarch))
			clusterExport(cluster, c("zz", "fft.step", "fft.by", 
							"wpars","cfinv","ghypmvcf"), envir = environment())
			xsol = parLapply(cluster, as.list(1:n), fun = function(i){
						xmin = min(apply(cbind(wpars[,1,i], wpars[,2,i], wpars[,3,i], wpars[,4,i], wpars[,5,i]), 1, 
										FUN = function(x) rugarch:::qgh(0.0000001, 
													alpha = x[1], beta = x[2], 
													delta = x[3], mu = x[4], 
													lambda = x[5])))
						xmax = max(apply(cbind(wpars[,1,i], wpars[,2,i], wpars[,3,i], wpars[,4,i], wpars[,5,i]), 1, 
										FUN = function(x) rugarch:::qgh(1-0.0000001, 
													alpha = x[1], beta = x[2], 
													delta = x[3], mu = x[4], 
													lambda = x[5])))
						zz = seq(xmin, xmax, by = fft.by)
						ans = cfinv(z = zz, f = ghypmvcf, 
								step = fft.step, lambda = wpars[,5,i], 
								alpha = wpars[,1,i], beta = wpars[,2,i], 
								delta = wpars[,3,i], mu = wpars[,4,i])
							sol = list(d = ans, support = c(xmin, xmax))
							return(sol)
						})
			for(i in 1:n){
				ghpdf[[i]] = xsol[[i]]$d
				support.user[i, ] = xsol[[i]]$support
			}
		} else{
			for(i in 1:n){
				xmin = min(apply(cbind(wpars[,1,i], wpars[,2,i], wpars[,3,i], wpars[,4,i], wpars[,5,i]), 1, 
								FUN = function(x) rugarch:::qgh(0.0000001, 
											alpha = x[1], beta = x[2], delta = x[3], 
											mu = x[4], lambda = x[5])))
				xmax = max(apply(cbind(wpars[,1,i], wpars[,2,i], wpars[,3,i], wpars[,4,i], wpars[,5,i]), 1, 
								FUN = function(x) rugarch:::qgh(1-0.0000001, 
											alpha = x[1], beta = x[2], delta = x[3], 
											mu = x[4], lambda = x[5])))
				zz = seq(xmin, xmax, by = fft.by)
				ghpdf[[i]] = cfinv(z = zz, f = ghypmvcf, step = fft.step, 
						lambda = wpars[,5,i], alpha = wpars[,1,i], 
						beta = wpars[,2,i], delta = wpars[,3,i], mu = wpars[,4,i])
				support.user[i, ] = c(xmin, xmax)
				if(debug) print(i)
			}
		}
	}
	sol$dist = "magh"
	sol$y = ghpdf
	sol$wM1 = wM1
	sol$fft.support = fft.support
	sol$support.method = support.method
	sol$support.user = support.user
	sol$fft.step = fft.step
	sol$fft.by = fft.by
	model$modeldata$data = NULL
	model$modeldata$mexdata = NULL
	model$modeldata$vexdata = NULL
	ans = new("goGARCHfft",
			dist = sol,
			model = model)
	return(ans)
}

.convolve.norm = function(M2, M1, A, weights, model, debug = FALSE)
{
	sol = list()
	n = dim(M2)[3]
	m = NCOL(A)
	if(!is.matrix(weights)) weights = matrix(weights, ncol = NROW(A), nrow = n, byrow = TRUE)
	if(NCOL(weights) != NROW(A)) stop("\nconvolution-->error: wrong no. of columns for weights matrix.", call. = FALSE)
	if(NROW(weights) != n) stop("\nconvolution-->error: wrong no. of rows for weights matrix.", call. = FALSE)
	
	# A matrix with the weighte portfolio sigma and mu
	normpdf = matrix(NA, ncol = 2, nrow = n)
	for(i in 1:n){
		w = weights[i, , drop = FALSE]
		normpdf[i, 1] = M1[i,] %*% t(w)
		normpdf[i, 2] = sqrt(w%*%(A%*%M2[,,i]%*%t(A))%*%t(w))
		if(debug) print(i)
	}
	colnames(normpdf) = c("mu", "sigma")
	sol$dist = "mvnorm"
	sol$y = normpdf
	sol$fft.support = NULL
	sol$fft.step = NULL
	sol$fft.by = NULL
	model$modeldata$data = NULL
	model$modeldata$mexdata = NULL
	model$modeldata$vexdata = NULL
	ans = new("goGARCHfft",
			dist = sol,
			model = model)
	return(ans)
}
#------------------------------------------------------------------------------
# interpolation of d, p, q portfolio density
dfft = function(object, index = 1)
{
	UseMethod("dfft")
}

setMethod("dfft", signature(object = "goGARCHfft"), .dfft)

pfft = function(object, index = 1)
{
	UseMethod("pfft")
}

setMethod("pfft", signature(object = "goGARCHfft"), .pfft)


qfft = function(object, index = 1)
{
	UseMethod("qfft")
}

setMethod("qfft", signature(object = "goGARCHfft"), .qfft)

#------------------------------------------------------------------------------
###############################################################################

nisurface = function(object, ...)
{
	UseMethod("nisurface")
}

.newsimpact.gogarch = function(object, type = "cov", pair = c(1,2), factor = c(1,2), 
		plot = TRUE, plot.type = c("surface", "contour"))
{
	ans = switch(type,
			cov = .newsimpact.gogarch.cov(object = object, pair = pair, 
					factor = factor, plot = plot, type = plot.type),
			cor = .newsimpact.gogarch.cor(object = object, pair = pair, 
					factor = factor, plot = plot, type = plot.type)
			)
	return(ans)
}

setMethod("nisurface", signature(object = "goGARCHfit"), .newsimpact.gogarch)
setMethod("nisurface", signature(object = "goGARCHfilter"), .newsimpact.gogarch)

.newsimpact.gogarch.cov = function(object, pair = c(1,2), factor = c(1,2), 
		plot = TRUE, type = c("surface", "contour"))
{
	if(is(object, "goGARCHfit")){
		garchcoef = object@mfit$garchcoef
		Y = object@mfit$Y
		A = object@mfit$A
	} else{
		garchcoef = object@mfilter$garchcoef
		Y = object@mfilter$Y
		A = object@mfilter$A
	}
	# factors
	m = NCOL(A)
	nilist = vector(mode = "list", length = m)
	filt = vector(mode = "list", length = m)
	cnames = object@model$modeldata$asset.names
	for(i in 1:m){
		specx = ugarchspec(mean.model = list(include.mean = FALSE, armaOrder = c(0,0)),
				variance.model = list(model = object@model$umodel$vmodel,
						submodel = object@model$umodel$vsubmodel, 
						variance.targeting = object@model$umodel$variance.targeting),
				distribution.model = object@model$umodel$distribution, 
				fixed.pars = as.list(garchcoef[,i]))
		filt[[i]] = ugarchfilter(specx, data = Y[,i], out.sample = object@model$modeldata$n.start)
	}
	nilist[[1]] = newsimpact(z = NULL, object = filt[[1]])
	# we want a common epsilon
	zeps = nilist[[1]]$zx
	for(i in 2:m) nilist[[i]]  = newsimpact(z = zeps, object = filt[[i]])
	sigmas = sapply(nilist, FUN = function(x) x$zy)
	N = dim(sigmas)[1]
	H = array(NA, dim = c(NROW(A), NROW(A), N))
	ni = matrix(NA, ncol = N, nrow = N)
	zr = as.numeric(which(round(nilist[[1]]$zx,12) == 0)[1])
	# to get the off diagonals we need to keep 1 fixed and
	# revolve the epsilon of the other
	s2 = sapply(nilist, FUN = function(x) x$zy)
	rownames(s2) = round(as.numeric(zeps),4)
	np = 1:m
	for(j in 1:N){
		for(i in 1:N){	
			sigmas = s2[zr,]
			sigmas1 = as.numeric(s2[j, factor[1]])
			sigmas2 = as.numeric(s2[i, factor[2]])
			sigmas[factor[1]] = sigmas1
			sigmas[factor[2]] = sigmas2
			H = A%*%diag(sigmas)%*%t(A)
			ni[j,i] = H[pair[1], pair[2]]
		}
	}
	if(plot){
		if(tolower(type[1]) == "surface"){
			x1 = shape::drapecol(ni, col = shape::femmecol(100), NAcol = "white")
			persp(  x = zeps,
					y = zeps,
					z = ni,  col = x1, theta = 45, phi = 25, expand = 0.5,
					ltheta = 120, shade = 0.75, ticktype = "detailed", xlab = paste("shock[f", factor[1], "]",sep=""),
					ylab = paste("shock[f", factor[2], "]",sep=""), zlab = "cov",
					cex.axis = 0.7,  cex.main = 0.8, main = paste("GOGARCH News Impact Covariance Surface\n", 
							cnames[pair[1]], "-", cnames[pair[2]], sep = ""))
		} else{
			tcol <- terrain.colors(12)
			contour(  x = zeps,
					y = zeps,
					z = ni,  col = tcol[2], lty = "solid",  cex.axis = 0.7,  cex.main = 0.8,
					main = paste("GOGARCH News Impact Covariance Contour\n", 
							cnames[pair[1]], "-", cnames[pair[2]], sep = ""))
		}
	}
	return(list(nisurface = ni, axis = zeps))
}


.newsimpact.gogarch.cor = function(object, pair = c(1,2), factor = c(1,2), 
		plot = TRUE, type = c("surface", "contour"))
{
	if(is(object, "goGARCHfit")){
		garchcoef = object@mfit$garchcoef
		Y = object@mfit$Y
		A = object@mfit$A
	} else{
		garchcoef = object@mfilter$garchcoef
		Y = object@mfilter$Y
		A = object@mfilter$A
	}
	cnames = object@model$modeldata$asset.names
	m = NCOL(A)
	nilist = vector(mode = "list", length = m)
	filt = vector(mode = "list", length = m)
	for(i in 1:m){
		specx = ugarchspec(mean.model = list(include.mean = FALSE, armaOrder = c(0,0)),
				variance.model = list(model = object@model$umodel$vmodel,
						submodel = object@model$umodel$vsubmodel, 
						variance.targeting = object@model$umodel$variance.targeting),
				distribution.model = object@model$umodel$distribution, 
				fixed.pars = as.list(garchcoef[,i]))
		filt[[i]] = ugarchfilter(specx, data = Y[,i], out.sample = object@model$modeldata$n.start)
	}
	nilist[[1]] = newsimpact(z = NULL, object = filt[[1]])	
	zeps = nilist[[1]]$zx
	for(i in 2:m) nilist[[i]]  = newsimpact(z = zeps, object = filt[[i]])
	sigmas = sapply(nilist, FUN = function(x) x$zy)
	N = dim(sigmas)[1]
	H = array(NA, dim = c(NROW(A), NROW(A), N))
	ni = matrix(NA, ncol = N, nrow = N)
	zr = as.numeric(which(round(nilist[[1]]$zx,12) == 0)[1])
	s2 = sapply(nilist, FUN = function(x) x$zy)
	rownames(s2) = round(as.numeric(zeps),4)
	np = 1:m
	for(j in 1:N){
		for(i in 1:N){	
			sigmas = s2[zr,]
			sigmas1 = as.numeric(s2[j, factor[1]])
			sigmas2 = as.numeric(s2[i, factor[2]])
			sigmas[factor[1]] = sigmas1
			sigmas[factor[2]] = sigmas2
			H = cov2cor(A%*%diag(sigmas)%*%t(A))
			ni[j,i] = H[pair[1], pair[2]]
		}
	}
	if(plot){
		if(tolower(type[1]) == "surface"){
			x1 = shape::drapecol(ni, col = shape::femmecol(100), NAcol = "white")
			persp(  x = zeps,
					y = zeps,
					z = ni,  col = x1, theta = 45, phi = 25, expand = 0.5,
					ltheta = 120, shade = 0.75, ticktype = "detailed", xlab = paste("shock[f", factor[1], "]",sep=""),
					ylab = paste("shock[f", factor[2], "]",sep=""), zlab = "cor",
					cex.axis = 0.7,  cex.main = 0.8, main = paste("GOGARCH News Impact Correlation Surface\n", 
							cnames[pair[1]], "-", cnames[pair[2]], sep = ""))
		} else{
			tcol <- terrain.colors(12)
			contour(  x = zeps,
					y = zeps,
					z = ni,  col = tcol[2], lty = "solid",  cex.axis = 0.7,  cex.main = 0.8,
					main = paste("GOGARCH News Impact Correlation Contour\n", 
							cnames[pair[1]], "-", cnames[pair[2]], sep = ""))
		}
	}
	return(list(nisurface = ni, axis = zeps))
}

#--------------------------------------------------------------------------------

betacovar = function(object, ...)
{
	UseMethod("betacovar")
}

.betacovar.fit = function(object, weights, asset = 1)
{
	A = as.matrix(object, which = "A")
	m = NROW(A)
	sig = as.matrix(sigma(object))
	n = NROW(sig)
	if(is.vector(weights)){
		if(length(weights) != m) stop("\nIncorrect weight vector length\n", call. = FALSE)
		weights = matrix(weights, ncol = m, nrow = n, byrow = TRUE)
	} else if(is.matrix(weights)){
		nn = NROW(weights)
		mm = NCOL(weights)
		if(mm != m) stop("\nIncorrect column dimension for weights matrix\n", call. = FALSE)
		if(nn != n) stop("\nIncorrect row dimension for weights matrix\n", call. = FALSE)
	}
	V = rcov(object)
	d = c(dim(V), asset-1)
	b = .Call("tvbetacovar", wi = weights, V_ = V, di = as.integer(d), PACKAGE="rmgarch")
	b = xts(b, object@model$modeldata$index[1:object@model$modeldata$T])
	colnames(b) = object@model$modeldata$asset.names[asset]
	return(b)
}

# TODO: methods for forecast and simulation

setMethod("betacovar", signature(object = "goGARCHfit"), .betacovar.fit)
setMethod("betacovar", signature(object = "goGARCHfilter"), .betacovar.fit)

.betacovar.forecast = function(object, weights, asset = 1)
{
	A = as.matrix(object, which = "A")
	m = NROW(A)
	sig = sigma(object)
	if(object@model$n.roll>0){
		n = object@model$n.roll+1
	} else{
		n = object@model$n.ahead
	}
	if(is.vector(weights)){
		if(length(weights) != m) stop("\nIncorrect weight vector length\n", call. = FALSE)
		weights = matrix(weights, ncol = m, nrow = n, byrow = TRUE)
	} else if(is.matrix(weights)){
		nn = NROW(weights)
		mm = NCOL(weights)
		if(mm != m) stop("\nIncorrect column dimension for weights matrix\n", call. = FALSE)
		if(nn != n) stop("\nIncorrect row dimension for weights matrix\n", call. = FALSE)
	}
	V = array(unlist(rcov(object)), dim = c(m,m,n))
	d = c(dim(V), asset-1)
	b = .Call("tvbetacovar", wi = weights, V_ = V, di = as.integer(d), PACKAGE="rmgarch")
	if(object@model$n.roll>0){
		b = xts(b, as.POSIXct(dimnames(fitted(object))[[3]]))
		colnames(b)<-paste(object@model$modeldata$asset.names[asset],"[T+1]",sep="")
	} else{
		b = matrix(b, ncol = 1)
		rownames(b) = dimnames(fitted(object))[[1]]
		colnames(b) = paste(object@model$modeldata$asset.names[asset],"[T0=",dimnames(fitted(object))[[3]],"]",sep="")
	}
	return(b)
}
setMethod("betacovar", signature(object = "goGARCHforecast"), .betacovar.forecast)
#--------------------------------------------------------------------------------
betacoskew = function(object, ...)
{
	UseMethod("betacoskew")
}


.betacoskew.fit = function(object, weights, asset = 1, cluster = NULL)
{
	A = as.matrix(object, which = "A")
	m = NROW(A)
	sig = as.matrix(sigma(object))
	n = NROW(sig)
	if(is.vector(weights)){
		if(length(weights) != m) stop("\nIncorrect weight vector length\n", call. = FALSE)
		weights = matrix(weights, ncol = m, nrow = n, byrow = TRUE)
	} else if(is.matrix(weights)){
		nn = NROW(weights)
		mm = NCOL(weights)
		if(mm != m) stop("\nIncorrect column dimension for weights matrix\n", call. = FALSE)
		if(nn != n) stop("\nIncorrect row dimension for weights matrix\n", call. = FALSE)
	}
	if(n<100){
		idx = cbind(1, n)
		k = 1
	} else{
		idx1 = seq(1, n, by = 50)
		if(idx1[length(idx1)]!=n) idx1 = c(idx1, n)			
		idx1 = idx1[-1]
		idx2 = c(1, idx1[-length(idx1)]+1)
		idx = cbind(idx2, idx1)
		k = NROW(idx)
	}
	b = rep(0, n)
	if(!is.null(cluster)){
		clusterExport(cluster, c("object","idx","asset","weights"), envir = environment())
		clusterEvalQ(cluster, loadNamespace('rmgarch'))
		ans = parLapply(cluster, 1:k, function(i){
			S = rmgarch::rcoskew(object, from = idx[i,1], to = idx[i,2], standardize = FALSE)
			d = c(dim(S), asset-1)
			sol = as.numeric( .Call("tvbetacoskew", wi = weights[idx[i,1]:idx[i,2],,drop=FALSE], Si = S, di = as.integer(d), PACKAGE="rmgarch") )
			return(sol)
			})
		for(i in 1:k) b[idx[i,1]:idx[i,2]] = ans[[i]]
	} else{
		for(i in 1:k){
			S = rcoskew(object, from = idx[i,1], to = idx[i,2], standardize = FALSE)
			d = c(dim(S), asset-1)
			b[idx[i,1]:idx[i,2]] = as.numeric( .Call("tvbetacoskew", wi = weights[idx[i,1]:idx[i,2],,drop=FALSE], Si = S, di = as.integer(d), PACKAGE="rmgarch") )
		}
	}
	b = xts(b, object@model$modeldata$index[1:object@model$modeldata$T])
	colnames(b) = object@model$modeldata$asset.names[asset]
	return(b)
}

setMethod("betacoskew", signature(object = "goGARCHfit"), .betacoskew.fit)
setMethod("betacoskew", signature(object = "goGARCHfilter"), .betacoskew.fit)

.betacoskew.forecast = function(object, weights, asset = 1)
{
	A = as.matrix(object, which = "A")
	m = NROW(A)
	sig = sigma(object)
	if(object@model$n.roll>0){
		n = object@model$n.roll+1
		useroll = TRUE
	} else{
		n = object@model$n.ahead
		useroll = FALSE
		if(n<100){
			idx = cbind(1, n)
			k = 1
		} else{
			idx1 = seq(1, n, by = 50)
			if(idx1[length(idx1)]!=n) idx1 = c(idx1, n)			
			idx1 = idx1[-1]
			idx2 = c(1, idx1[-length(idx1)]+1)
			idx = cbind(idx2, idx1)
			k = NROW(idx)
		}
	}
	if(is.vector(weights)){
		if(length(weights) != m) stop("\nIncorrect weight vector length\n", call. = FALSE)
		weights = matrix(weights, ncol = m, nrow = n, byrow = TRUE)
	} else if(is.matrix(weights)){
		nn = NROW(weights)
		mm = NCOL(weights)
		if(mm != m) stop("\nIncorrect column dimension for weights matrix\n", call. = FALSE)
		if(nn != n) stop("\nIncorrect row dimension for weights matrix\n", call. = FALSE)
	}
	b = rep(0, n)
	if(useroll){
		if(n<100){
			S = rcoskew(object, from = 1, to = 1, roll = "all", standardize = FALSE)
			d = c(dim(S), asset-1)
			b = as.numeric( .Call("tvbetacoskew", wi = weights, Si = S, di = as.integer(d), PACKAGE="rmgarch") )
		} else{
			for(i in 1:n){
				S = rcoskew(object, from = 1, to = 1, roll = i-1, standardize = FALSE)
				d = c(dim(S), asset-1)
				b[i] = as.numeric( .Call("tvbetacoskew", wi = weights[i,,drop=FALSE], Si = S, di = as.integer(d), PACKAGE="rmgarch") )
			}
		}
	} else{
		for(i in 1:k){
			S = rcoskew(object, from = idx[i,1], to = idx[i,2], roll=0, standardize = FALSE)
			d = c(dim(S), asset-1)
			b[idx[i,1]:idx[i,2]] = as.numeric( .Call("tvbetacoskew", wi = weights[idx[i,1]:idx[i,2],,drop=FALSE], Si = S, di = as.integer(d), PACKAGE="rmgarch") )
		}
	}
	if(object@model$n.roll>0){
		b = xts(b, as.POSIXct(dimnames(fitted(object))[[3]]))
		colnames(b)<-paste(object@model$modeldata$asset.names[asset],"[T+1]",sep="")
	} else{
		b = matrix(b, ncol = 1)
		rownames(b) = dimnames(fitted(object))[[1]]
		colnames(b) = paste(object@model$modeldata$asset.names[asset],"[T0=",dimnames(fitted(object))[[3]],"]",sep="")
	}
	return(b)
}

setMethod("betacoskew", signature(object = "goGARCHforecast"), .betacoskew.forecast)
#--------------------------------------------------------------------------------
betacokurt = function(object, ...)
{
	UseMethod("betacokurt")
}
.betacokurt.fit = function(object, weights, asset = 1, cluster = NULL)
{
	A = as.matrix(object, which = "A")
	m = NROW(A)
	sig = as.matrix(sigma(object))
	n = NROW(sig)
	if(is.vector(weights)){
		if(length(weights) != m) stop("\nIncorrect weight vector length\n", call. = FALSE)
		weights = matrix(weights, ncol = m, nrow = n, byrow = TRUE)
	} else if(is.matrix(weights)){
		nn = NROW(weights)
		mm = NCOL(weights)
		if(mm != m) stop("\nIncorrect column dimension for weights matrix\n", call. = FALSE)
		if(nn != n) stop("\nIncorrect row dimension for weights matrix\n", call. = FALSE)
	}
	if(n<100){
		idx = cbind(1, n)
		k = 1
	} else{
		idx1 = seq(1, n, by = 50)
		if(idx1[length(idx1)]!=n) idx1 = c(idx1, n)			
		idx1 = idx1[-1]
		idx2 = c(1, idx1[-length(idx1)]+1)
		idx = cbind(idx2, idx1)
		k = NROW(idx)
	}
	b = rep(0, n)
	
	if(!is.null(cluster)){
		clusterExport(cluster, c("object","idx","asset","weights"), envir = environment())
		clusterEvalQ(cluster, loadNamespace('rmgarch'))
		ans = parLapply(cluster, 1:k, function(i){
					K = rmgarch::rcokurt(object, from = idx[i,1], to = idx[i,2], standardize = FALSE)
					d = c(dim(K), asset-1)
					sol = as.numeric( .Call("tvbetacokurt", wi = weights[idx[i,1]:idx[i,2],,drop=FALSE], Ki = K, di = as.integer(d), PACKAGE="rmgarch") )
					return(sol)
				})
		for(i in 1:k) b[idx[i,1]:idx[i,2]] = ans[[i]]
	} else{
		for(i in 1:k){
			K = rcokurt(object, from = idx[i,1], to = idx[i,2], standardize = FALSE)
			d = c(dim(K), asset-1)
			b[idx[i,1]:idx[i,2]] = as.numeric( .Call("tvbetacokurt", wi = weights[idx[i,1]:idx[i,2],,drop=FALSE], Ki = K, di = as.integer(d), PACKAGE="rmgarch") )
		}
	}
	b = xts(b, object@model$modeldata$index[1:object@model$modeldata$T])
	colnames(b) = object@model$modeldata$asset.names[asset]
	return(b)
}

setMethod("betacokurt", signature(object = "goGARCHfit"), .betacokurt.fit)
setMethod("betacokurt", signature(object = "goGARCHfilter"), .betacokurt.fit)


.betacokurt.forecast = function(object, weights, asset = 1)
{
	A = as.matrix(object, which = "A")
	m = NROW(A)
	sig = sigma(object)
	if(object@model$n.roll>0){
		n = object@model$n.roll+1
		useroll = TRUE
	} else{
		n = object@model$n.ahead
		useroll = FALSE
		if(n<100){
			idx = cbind(1, n)
			k = 1
		} else{
			idx1 = seq(1, n, by = 50)
			if(idx1[length(idx1)]!=n) idx1 = c(idx1, n)			
			idx1 = idx1[-1]
			idx2 = c(1, idx1[-length(idx1)]+1)
			idx = cbind(idx2, idx1)
			k = NROW(idx)
		}
	}
	if(is.vector(weights)){
		if(length(weights) != m) stop("\nIncorrect weight vector length\n", call. = FALSE)
		weights = matrix(weights, ncol = m, nrow = n, byrow = TRUE)
	} else if(is.matrix(weights)){
		nn = NROW(weights)
		mm = NCOL(weights)
		if(mm != m) stop("\nIncorrect column dimension for weights matrix\n", call. = FALSE)
		if(nn != n) stop("\nIncorrect row dimension for weights matrix\n", call. = FALSE)
	}
	b = rep(0, n)
	if(useroll){
		if(n<50){
			K = rcokurt(object, from = 1, to = 1, roll = "all", standardize = FALSE)
			d = c(dim(K), asset-1)
			b = as.numeric( .Call("tvbetacokurt", wi = weights, Ki = K, di = as.integer(d), PACKAGE="rmgarch") )
		} else{
			for(i in 1:n){
				K = rcokurt(object, from = 1, to = 1, roll = i-1, standardize = FALSE)
				d = c(dim(K), asset-1)
				b[i] = as.numeric( .Call("tvbetacokurt", wi = weights[i,,drop=FALSE], Ki = K, di = as.integer(d), PACKAGE="rmgarch") )
			}
		}
	} else{
		for(i in 1:k){
			K = rcokurt(object, from = idx[i,1], to = idx[i,2], roll=0, standardize = FALSE)
			d = c(dim(K), asset-1)
			b[idx[i,1]:idx[i,2]] = as.numeric( .Call("tvbetacokurt", wi = weights[idx[i,1]:idx[i,2],,drop=FALSE], Ki = K, di = as.integer(d), PACKAGE="rmgarch") )
		}
	}
	if(object@model$n.roll>0){
		b = xts(b, as.POSIXct(dimnames(fitted(object))[[3]]))
		colnames(b)<-paste(object@model$modeldata$asset.names[asset],"[T+1]",sep="")
	} else{
		b = matrix(b, ncol = 1)
		rownames(b) = dimnames(fitted(object))[[1]]
		colnames(b) = paste(object@model$modeldata$asset.names[asset],"[T0=",dimnames(fitted(object))[[3]],"]",sep="")
	}
	return(b)
}
setMethod("betacokurt", signature(object = "goGARCHforecast"), .betacokurt.forecast)

#-----------------------------------------------------------------------------
.abind = function (x, y) 
{
	m  = dim(x)[1]
	mm = dim(x)[2]
	n1 = dim(x)[3]
	n2 = dim(y)[3]
	dimnames(x)<-NULL
	dimnames(y)<-NULL
	nw = array(NA, dim = c(m, mm, n1 + n2))
	nw[, , 1:n1] = x[, , 1:n1]
	nw[, , (n1 + 1):(n1 + n2)] = y[, , 1:n2]
	return(nw)
}