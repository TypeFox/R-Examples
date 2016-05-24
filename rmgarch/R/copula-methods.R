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
cgarchspec = function(uspec, VAR = FALSE, robust = FALSE, lag = 1, lag.max = NULL, 
		lag.criterion = c("AIC", "HQ", "SC", "FPE"), external.regressors = NULL,
		robust.control = list("gamma" = 0.25, "delta" = 0.01, "nc" = 10, "ns" = 500), 
		dccOrder = c(1,1), asymmetric = FALSE, distribution.model = list(copula = c("mvnorm", "mvt"), 
				method = c("Kendall", "ML"), time.varying = FALSE, 
				transformation = c("parametric", "empirical", "spd")),
		start.pars = list(), fixed.pars = list())
{
	UseMethod("cgarchspec")	
}

setMethod("cgarchspec", signature(uspec = "uGARCHmultispec"), .cgarchspec)

.setfixedcg = function(object, value){
	model = object@model
	umodel = object@umodel
	ipars = model$pars
	pars = unlist(value)
	names(pars) = parnames = names(pars)
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
	names(fixed.pars) = names(pars[inc])
	
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, umodel$vt)
	# set parameter values
	tmp = cgarchspec(uspec = mspec, 
			VAR = ifelse(model$modelinc[1]>0, TRUE, FALSE), 
			robust = ifelse(!is.null(model$varmodel$robust), model$varmodel$robust, FALSE), 
			lag = model$modelinc[1], lag.max = model$varmodel$lag.max, 
			lag.criterion = model$varmodel$lag.criterion, 
			external.regressors = if(model$modelinc[2]>0) model$modeldata$mexdata else NULL, 
			robust.control = if(!is.null(model$varmodel$robust.control)) model$varmodel$robust.control else NULL, 
			dccOrder = model$modelinc[4:5], asymmetric = ifelse(model$modelinc[6]>0, TRUE, FALSE), 
			distribution.model = list(copula = model$modeldesc$distribution, 
					method = model$modeldesc$cor.method, 
					time.varying = model$modeldesc$timecopula, 
					transformation = model$modeldesc$transformation), 
			start.pars = if(is.null(model$start.pars)) model$start.pars else list(), 
			fixed.pars = as.list(fixed.pars))
	return(tmp)
}

setReplaceMethod(f="setfixed", signature= c(object = "cGARCHspec", value = "vector"), definition = .setfixedcg)


.setstartcg = function(object, value){
	model = object@model
	umodel = object@umodel
	ipars = model$pars
	pars = unlist(value)
	names(pars) = parnames = names(pars)
	# included parameters in model
	modelnames = rownames(ipars[which(ipars[,4]==1 | ipars[,2]==1), ,drop=FALSE])
	inc = NULL
	for(i in seq_along(parnames)){
		if(is.na(match(parnames[i], modelnames))){
			warning( (paste("Unrecognized Parameter in Start Values: ", parnames[i], "...Ignored", sep = "")))
		} else{
			inc = c(inc, i)
		}
	}
	start.pars = pars[inc]
	names(start.pars) = names(pars[inc])
	
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, umodel$vt)
	
	# set parameter values
	tmp = cgarchspec(uspec = mspec, VAR = ifelse(model$modelinc[1]>0, TRUE, FALSE), 
			robust = ifelse(!is.null(model$varmodel$robust), model$varmodel$robust, FALSE), 
			lag = model$modelinc[1], lag.max = model$varmodel$lag.max, 
			lag.criterion = model$varmodel$lag.criterion, 
			external.regressors = if(model$modelinc[2]>0) model$modeldata$mexdata else NULL, 
			robust.control = if(!is.null(model$varmodel$robust.control)) model$varmodel$robust.control else NULL, 
			dccOrder = model$modelinc[4:5], asymmetric = ifelse(model$modelinc[6]>0, TRUE, FALSE), 
			distribution.model = list(copula = model$modeldesc$distribution, 
					method = model$modeldesc$cor.method, 
					time.varying = model$modeldesc$timecopula, 
					transformation = model$modeldesc$transformation), 
			start.pars = as.list(start.pars),  fixed.pars =  model$fixed.pars)
	return(tmp)
}

setReplaceMethod(f="setstart", signature= c(object = "cGARCHspec", value = "vector"), definition = .setstartcg)


cgarchfit = function(spec, data, spd.control = list(lower = 0.1, upper = 0.9, 
				type = "pwm", kernel = "epanech"), 
		fit.control = list(eval.se = TRUE, stationarity = TRUE, scale = FALSE), 
		solver = "solnp", solver.control = list(), out.sample = 0, cluster = NULL, 
		fit = NULL, VAR.fit = NULL, realizedVol = NULL, ...)
{
	UseMethod("cgarchfit")
}

setMethod("cgarchfit", signature(spec = "cGARCHspec"), .cgarchfit)


cgarchfilter = function(spec, data, out.sample = 0, filter.control = list(n.old = NULL), 
		spd.control = list(lower = 0.1, upper = 0.9, type = "pwm", kernel = "epanech"), 
		cluster = NULL, varcoef = NULL, realizedVol = NULL, ...)
{
	UseMethod("cgarchfilter")
}

setMethod("cgarchfilter", signature(spec = "cGARCHspec"), .cgarchfilter)

cgarchsim = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, 
		startMethod = c("unconditional", "sample"), presigma = NULL, 
		preresiduals = NULL, prereturns = NULL, preR = NULL, preQ = NULL,
		preZ = NULL, rseed = NULL, mexsimdata = NULL, vexsimdata = NULL, 
		cluster = NULL, only.density= FALSE, prerealized = NULL, ...)
{
	UseMethod("cgarchsim")
}

setMethod("cgarchsim", signature(fit = "cGARCHfit"), .cgarchsim)

#----------------------------------------------------------------------------------
# show methods
#----------------------------------------------------------------------------------
setMethod("show",
		signature(object = "cGARCHspec"),
		function(object){
			m = dim(object@umodel$modelinc)[2]
			dccpars = sum(object@model$modelinc[4:8])
			mlpars = sum(object@model$modelinc[3])
			garchpars = sum( object@umodel$modelinc[1:18,] )
			# VAR = mxm x lags + lags*constant + lags*mxreg
			varpars = object@model$modelinc[1]*(m*m + m + object@model$modelinc[2])
			NPx = mlpars + dccpars + garchpars + varpars + ( (m^2 - m)/2 )
			cat(paste("\n*--------------------------------*", sep = ""))
			cat(paste("\n*       Copula GARCH Spec        *", sep = ""))
			cat(paste("\n*--------------------------------*", sep = ""))
			cat("\n\nDistribution\t\t: ", object@model$modeldesc$distribution)
			if(object@model$modeldesc$timecopula){
				cat("\nModel\t\t\t\t: ", paste(object@model$modeldesc$dccmodel, "(", object@model$modelinc[4], ",", object@model$modelinc[5],")", sep=""))
			}
			cat("\nTransformation\t\t: ", object@model$modeldesc$transformation)
			if(!object@model$modeldesc$timecopula) cat("\nCorrelation Estimate: ", object@model$modeldesc$cor.method)
			cat("\nNo. of Parameters\t: ", NPx)
			if(!object@model$modeldesc$timecopula){
				NP = paste("[",varpars, "+", garchpars,"+", (m^2 - m)/2,"]", sep="")
				cat("\n[VAR GARCH COV]\t\t:", NP)
			} else{
				NP = paste("[",varpars, "+", garchpars,"+", dccpars, "+",(m^2 - m)/2,"]", sep="")
				cat("\n[VAR GARCH DCC UncQ]:", NP)
			}
			cat("\nNo. of Series\t\t: ", m)
			cat("\n\n")
			invisible(object)
		})

# fit show
setMethod("show",
		signature(object = "cGARCHfit"),
		function(object){
			m = dim(object@model$umodel$modelinc)[2]
			cat(paste("\n*-------------------------------------------------*", sep = ""))
			cat(paste("\n*                  Copula GARCH Fit               *", sep = ""))
			cat(paste("\n*-------------------------------------------------*", sep = ""))	
			cat("\n\nDistribution\t\t: ", object@model$modeldesc$distribution)
			if(object@model$modelinc[1]>0){
				npvar = dim(object@model$varcoef)[1] * dim(object@model$varcoef)[2]
			} else{
				npvar = 0
			}
			if(object@model$modeldesc$timecopula){
				cat("\nDCC Order\t\t\t: ", object@model$modelinc[4:5])
				cat("\nAsymmetric\t\t\t: ", ifelse(object@model$modelinc[6]>0, TRUE, FALSE))
				NP = paste("[",npvar, "+", length(object@mfit$garchnames),"+",length(object@mfit$dccnames), "+",(m^2 - m)/2,"]", sep="")
				cat("\nNo. of Parameters\t: ", npvar+length(object@mfit$garchnames) + length(object@mfit$dccnames) + ( (m^2 - m)/2 ))
				cat("\n[VAR GARCH DCC UncQ]:", NP)
			} else{
				NP = paste("[",npvar, "+", length(object@mfit$garchnames),"+",length(object@mfit$dccnames),"]", sep="")
				cat("\nNo. of Parameters\t: ", npvar+length(object@mfit$garchnames) + length(object@mfit$dccnames))
				cat("\n[VAR GARCH CC]\t\t:", NP)
			}			
			cat("\nNo. of Series\t\t: ", m)
			cat("\nNo. of Observations\t: ", object@model$modeldata$T)
			cat("\nLog-Likelihood\t\t: ", object@mfit$llh)
			cat("\nAv.Log-Likelihood\t: ", round(object@mfit$llh/object@model$modeldata$T,3), "\n")
			cat("\nOptimal Parameters")
			cat(paste("\n---------------------------------------------------\n", sep = ""))
			print(round(object@mfit$matcoef,6), digits = 5)
			itest = rugarch:::.information.test(object@mfit$llh, nObs = object@model$modeldata$T, nPars = length(object@mfit$matcoef[,1]))
			itestm = matrix(0, ncol = 1, nrow = 4)
			itestm[1,1] = itest$AIC
			itestm[2,1] = itest$BIC
			itestm[3,1] = itest$SIC
			itestm[4,1] = itest$HQIC
			colnames(itestm) = ""
			rownames(itestm) = c("Akaike", "Bayes", "Shibata", "Hannan-Quinn")
			cat("\nInformation Criteria")
			cat(paste("\n---------------------\n", sep = ""))
			print(itestm,digits=5)
			cat("\n")
			cat("\nElapsed time :", object@mfit$timer,"\n\n")
			invisible(object)
})

setMethod("show",
		signature(object = "cGARCHfilter"),
		function(object){
			m = dim(object@model$umodel$modelinc)[2]
			cat(paste("\n*-------------------------------------------------*", sep = ""))
			cat(paste("\n*               Copula GARCH Filter               *", sep = ""))
			cat(paste("\n*-------------------------------------------------*", sep = ""))	
			cat("\n\nDistribution\t\t: ", object@model$modeldesc$distribution)
			if(object@model$modelinc[1]>0){
				npvar = dim(object@model$varcoef)[1] * dim(object@model$varcoef)[2]
			} else{
				npvar = 0
			}
			if(object@model$modeldesc$timecopula){
				cat("\nDCC Order\t\t\t: ", object@model$modelinc[4:5])
				cat("\nAsymmetric\t\t\t: ", ifelse(object@model$modelinc[6]>0, TRUE, FALSE))
				NP = paste("[",npvar, "+", length(object@mfilter$garchnames),"+",length(object@mfilter$dccnames), "+",(m^2 - m)/2,"]", sep="")
				cat("\nNo. of Parameters\t: ", npvar+length(object@mfilter$garchnames) + length(object@mfilter$dccnames) + ( (m^2 - m)/2 ))
				cat("\n[VAR GARCH DCC UncQ]:", NP)
			} else{
				NP = paste("[",npvar, "+", length(object@mfilter$garchnames),"+",length(object@mfilter$dccnames),"]", sep="")
				cat("\nNo. of Parameters\t: ", npvar+length(object@mfilter$garchnames) + length(object@mfilter$dccnames))
				cat("\n[VAR GARCH CC]\t\t:", NP)
			}			
			cat("\nNo. of Series\t\t: ", m)
			cat("\nNo. of Observations\t: ", object@model$modeldata$T)
			cat("\nLog-Likelihood\t\t: ", object@mfilter$llh)
			cat("\nAv.Log-Likelihood\t: ", round(object@mfilter$llh/object@model$modeldata$T,3), "\n")
			cat("\nOptimal Parameters")
			cat(paste("\n---------------------------------------------------\n", sep = ""))
			cf = data.frame(object@mfilter$coef, row.names = names(object@mfilter$coef))
			colnames(cf) = "Value"
			print(round(cf, 4), digits = 5)
			itest = rugarch:::.information.test(object@mfilter$llh, nObs = object@model$modeldata$T, nPars = length(object@mfilter$coef))
			itestm = matrix(0, ncol = 1, nrow = 4)
			itestm[1,1] = itest$AIC
			itestm[2,1] = itest$BIC
			itestm[3,1] = itest$SIC
			itestm[4,1] = itest$HQIC
			colnames(itestm) = ""
			rownames(itestm) = c("Akaike", "Bayes", "Shibata", "Hannan-Quinn")
			cat("\nInformation Criteria")
			cat(paste("\n---------------------\n", sep = ""))
			print(itestm,digits=5)
			cat("\n")
			cat("\nElapsed time :", object@mfilter$timer,"\n\n")
			invisible(object)
		})


setMethod("show",
		signature(object = "cGARCHsim"),
		function(object){
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\n*      Copula GARCH Simulation    *", sep = ""))
			cat(paste("\n*---------------------------------*", sep = ""))
			cat("\n\nDistribution\t\t:", object@model$modeldesc$distribution)
			cat("\nTransformation\t\t:", object@model$modeldesc$transformation)
			cat("\nTime-Varying\t\t:", object@model$modeldesc$timecopula)
			cat(paste("\nSimulation Horizon\t: ",  object@model$n.sim, sep = ""))
			cat(paste("\nBurn In\t\t\t\t: ",  object@model$n.start, sep = ""))
			cat(paste("\nNo. of Simulations\t: ",object@model$m.sim, sep = ""))
			cat("\n\n")
			invisible(object)
		})
#----------------------------------------------------------------------------------
# fitted
#----------------------------------------------------------------------------------
.fitted.cgarchfit = function(object)
{
	T = object@model$modeldata$T
	ans = xts(object@model$mu[1:T,,drop=FALSE], object@model$modeldata$index[1:T])
	colnames(ans) = object@model$modeldata$asset.names
	return( ans )
}

setMethod("fitted", signature(object = "cGARCHfit"), .fitted.cgarchfit)
setMethod("fitted", signature(object = "cGARCHfilter"), .fitted.cgarchfit)


.fitted.cgarchsim = function(object, sim = 1)
{
	n = object@model$m.sim
	m.sim = as.integer(sim)
	if( m.sim > n | m.sim < 1 ) stop("\rmgarch-->error: fitted (simulation) sim index out of bounds!")
	ans = object@msim$simX[[m.sim]]
	colnames(ans) = object@model$modeldata$asset.names
	rownames(ans) = NULL
	return( ans )
}

setMethod("fitted", signature(object = "cGARCHsim"), .fitted.cgarchsim)
#----------------------------------------------------------------------------------
# residuals
#----------------------------------------------------------------------------------
.residuals.cgarchfit = function(object)
{
	T = object@model$modeldata$T
	ans = xts(object@model$residuals[1:T,,drop=FALSE], object@model$modeldata$index[1:T])
	colnames(ans) = object@model$modeldata$asset.names
	return( ans )
}

setMethod("residuals", signature(object = "cGARCHfit"), .residuals.cgarchfit)
setMethod("residuals", signature(object = "cGARCHfilter"), .residuals.cgarchfit)
#----------------------------------------------------------------------------------
# sigma
#----------------------------------------------------------------------------------
.sigma.cgarchfit = function(object)
{
	T = object@model$modeldata$T
	H = rcov(object)[,,1:T]
	m = dim(H)[2]
	sig = sqrt(.Call("ArrayDiag", H, c(m,m,T), PACKAGE="rmgarch"))
	# sig = sqrt(t(apply(H, 3, FUN = function(x) diag(x))))
	sig = xts(sig[1:T,,drop=FALSE], object@model$modeldata$index[1:T])
	colnames(sig) = object@model$modeldata$asset.names
	return( sig )
}

setMethod("sigma", signature(object = "cGARCHfit"), .sigma.cgarchfit)
setMethod("sigma", signature(object = "cGARCHfilter"), .sigma.cgarchfit)

.sigma.cgarchsim = function(object, sim = 1)
{
	n = object@model$m.sim
	m.sim = as.integer(sim)
	if( m.sim > n | m.sim < 1 ) stop("\rmgarch-->error: fitted (simulation) sim index out of bounds!")
	H = rcov(object, m.sim)
	dm = dim(H)
	ans = sqrt(.Call("ArrayDiag", H, c(dm[1],dm[2],dm[3]), PACKAGE="rmgarch"))
	colnames(ans) = object@model$modeldata$asset.names
	rownames(ans) = NULL
	return( ans )
}

setMethod("sigma", signature(object = "cGARCHsim"), .sigma.cgarchsim)

#----------------------------------------------------------------------------------
# rcov
#----------------------------------------------------------------------------------
.rcov.cgarchfit = function(object)
{
	ans = switch(class(object)[1],
			cGARCHfit = object@mfit$H,
			cGARCHfilter = object@mfilter$H)
	nam = object@model$modeldata$asset.names
	D = as.character(object@model$modeldata$index[1:object@model$modeldata$T])
	dimnames(ans)<-list(nam, nam, D)
	return(ans)
}

setMethod("rcov", signature(object = "cGARCHfit"), .rcov.cgarchfit)
setMethod("rcov", signature(object = "cGARCHfilter"), .rcov.cgarchfit)


.rcov.cgarchsim = function(object, sim = 1)
{
	n = object@model$m.sim
	m.sim = as.integer(sim)
	if( m.sim > n | m.sim < 1 ) stop("\rmgarch-->error: rcor sim index out of bounds!")	
	ans = object@msim$simH[[sim]]
	nam = object@model$modeldata$asset.names
	dimnames(ans) = list(nam, nam, 1:dim(ans)[3])
	return( ans )
}

setMethod("rcov", signature(object = "cGARCHsim"), .rcov.cgarchsim)

.rcor.cgarchfit = function(object)
{
	if(object@model$modeldesc$timecopula){
		if(class(object)[1]=="cGARCHfit") R = object@mfit$Rt else R = object@mfilter$Rt
		n = length(R)
		m = dim(R[[1]])[1]
		ans = array(NA, dim = c(m,m,n))
		ans[,,1:n] = sapply(R, FUN = function(x) x)
		nam = object@model$modeldata$asset.names
		D = as.character(object@model$modeldata$index[1:object@model$modeldata$T])
		dimnames(ans)<-list(nam, nam, D)
	} else{
		if(class(object)[1]=="cGARCHfit") ans = object@mfit$Rt else ans = object@mfilter$Rt
		nam = object@model$modeldata$asset.names
		colnames(ans) = rownames(ans) = nam
	}
	return(ans)
}
setMethod("rcor", signature(object = "cGARCHfit"), .rcor.cgarchfit)
setMethod("rcor", signature(object = "cGARCHfilter"), .rcor.cgarchfit)

#----------------------------------------------------------------------------------
# rcor
#----------------------------------------------------------------------------------
.rcor.cgarchsim = function(object, sim = 1)
{
	n = length(object@msim$simH)
	m.sim = as.integer(sim)
	if( m.sim > n | m.sim < 1 ) stop("\rmgarch-->error: rcor sim index out of bounds!")
	
	if(object@model$modeldesc$timecopula){
		ans = object@msim$simR[[sim]]
		nam = object@model$modeldata$asset.names
		dimnames(ans) = list(nam, nam, 1:dim(ans)[3])
	} else{
		# fixed (no uncertainty)
		ans = cov2cor(object@msim$simH[[m.sim]][,,1])
		nam = object@model$modeldata$asset.names
		colnames(ans) = rownames(ans) = nam
	}
	return( ans )
}

setMethod("rcor", signature(object = "cGARCHsim"), .rcor.cgarchsim)

#----------------------------------------------------------------------------------
# coef
#----------------------------------------------------------------------------------
.coef.cgarchfit = function(object, type = "all")
{
	mpars = object@model$mpars
	m = dim(mpars)[2]-1
	if( type == "all" ){
		cf = object@mfit$coef
	} else if( type == "garch" ){
		cf = mpars[which(object@model$eidx[,1:m]==1)]
		names(cf) = object@mfit$garchnames
	} else{
		cf = mpars[which(object@model$eidx[,m+1]==1), m+1]
		names(cf) = object@mfit$dccnames
	}
	return( cf )
}

setMethod("coef", signature(object = "cGARCHfit"), .coef.cgarchfit)

.coef.cgarchfilter = function(object, type = "all")
{
	m = dim(object@model$umodel$modelinc)[2]
	mpars = object@model$mpars
	if( type == "all" ){
		cf = object@mfilter$matcoef[,1]
	} else if( type == "garch" ){
		cf = mpars[which(object@model$midx[,1:m]==1)]
		names(cf) = object@mfilter$garchnames
	} else{
		if(any(object@model$midx[,m+1]>0)){
			cf = mpars[which(object@model$midx[,m+1]==1),m+1]
			names(cf) = object@mfilter$dccnames
		} else{
			cf = NULL
		}
	}
	return( cf )
}

setMethod("coef", signature(object = "cGARCHfilter"), .coef.cgarchfilter)
#----------------------------------------------------------------------------------
# likelihood
#----------------------------------------------------------------------------------
.likelihood.cgarchfit = function(object)
{
	switch(class(object)[1],
			cGARCHfit = object@mfit$llh,
			cGARCHfilter = object@mfilter$llh)
}

setMethod("likelihood", signature(object = "cGARCHfit"), .likelihood.cgarchfit)
setMethod("likelihood", signature(object = "cGARCHfilter"), .likelihood.cgarchfit)


.shape.cgarchfit = function(object)
{
	sh = NA
	if( object@model$modelinc[7]>0 ){
		cf = object@model$mpars[,dim(object@model$mpars)[2]]
		nx = which( substr(names(cf), 1, 6) == "mshape" )
		sh = cf[nx]
	}
	return( sh )
}

setMethod("rshape", signature(object = "cGARCHfit"), .shape.cgarchfit)
setMethod("rshape", signature(object = "cGARCHfilter"), .shape.cgarchfit)



.skew.cgarchfit = function(object)
{
	sk = NA
	if( object@model$modelinc[8]>0 ){
		cf = object@model$mpars[,dim(object@model$mpars)[2]]
		nx = which( substr(names(cf), 1, 5) == "mskew" )
		sk = cf[nx]
	}
	return( sk )
}

setMethod("rskew", signature(object = "cGARCHfit"), .skew.cgarchfit)
setMethod("rskew", signature(object = "cGARCHfilter"), .skew.cgarchfit)