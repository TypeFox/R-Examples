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
arfimaspec = function( mean.model = list(armaOrder = c(1,1), include.mean = TRUE, arfima = FALSE, 
				external.regressors = NULL), 
		distribution.model = "norm", start.pars = list(), fixed.pars = list(), ...)
{
	UseMethod("arfimaspec")
}


.xarfimaspec = function( mean.model = list(armaOrder = c(1,1), include.mean = TRUE, arfima = FALSE, 
				external.regressors = NULL), 
		distribution.model = "norm", start.pars = list(), fixed.pars = list(), ...)
{
	mmodel = mean.model
	dmodel = distribution.model
	# make the temporary subsitution so that it will be accepted by ugarchspec
	if(!is.null(fixed.pars$sigma)){
		fixed.pars$omega = fixed.pars$sigma
		fixed.pars$sigma = NULL
	}
	if(!is.null(start.pars$sigma)){
		start.pars$omega = start.pars$sigma
		start.pars$sigma = NULL
	}
	mm = match(names(mean.model), c("armaOrder", "include.mean", "arfima", "external.regressors"))
	if(any(is.na(mm))){
		idx = which(is.na(mm))
		enx = NULL
		for(i in 1:length(idx)) enx = c(enx, names(mean.model)[idx[i]])
		warning(paste(c("unidentified option(s) in mean.model:\n", enx), sep="", collapse=" "), call. = FALSE, domain = NULL)
	}
	ans = ugarchspec(mean.model = list(armaOrder = mmodel$armaOrder, include.mean = mmodel$include.mean,
					arfima = mmodel$arfima, external.regressors = mmodel$external.regressors, archm = FALSE,
	archpow = 1), variance.model = list(garchOrder = c(0,0), model = "sGARCH"), distribution.model = dmodel,
	start.pars = start.pars, fixed.pars = fixed.pars)

	if(!is.null(ans@model$fixed.pars$omega)){
		ans@model$fixed.pars$sigma = ans@model$fixed.pars$omega
		ans@model$fixed.pars$omega = NULL
	}
	if(!is.null(ans@model$start.pars$omega)){
		ans@model$start.pars$sigma = ans@model$start.pars$omega
		ans@model$start.pars$omega = NULL
	}
	model = ans@model
	# change the omega name to sigma
	names(model$modelinc)[7] = "sigma"
	model$modeldesc$vmodel = "constant"
    idx = which(rownames(model$pars) == "omega")
	rownames(model$pars)[idx] = "sigma"
	rownames(model$pos.matrix)[7] = "sigma"
	rownames(model$pidx)[7] = "sigma"
	sol = new("ARFIMAspec", model = model)
	return(sol)
}

setMethod(f = "arfimaspec", definition = .xarfimaspec)

##############################################################################
# custom arfima for non consecutive arma orders
.zarfimaspec = function( arOrder = c(1,1), maOrder = c(1,1), include.mean = TRUE, 
		arfima = FALSE, external.regressors = NULL, distribution.model = "norm", 
		start.pars = list(), fixed.pars = list(), ...){
	
	fx = list()
	if(sum(arOrder)>0){
		arm = max(which(arOrder==1))
		arOrder = arOrder[1:arm]
	} else{
		arm = 0
	}
	if(sum(maOrder)>0){
		mam = max(which(maOrder==1)) 
		maOrder = maOrder[1:mam]
	} else{
		mam = 0
	}
	if(any(arOrder==0)){
		idx = which(arOrder==0)
		for(i in 1:length(idx)){
			eval(parse(text=paste("fx$ar", idx[i],"=0",sep="")))
		}
	}
	if(any(maOrder==0)){
		idx = which(maOrder==0)
		for(i in 1:length(idx)){
			eval(parse(text=paste("fx$ma", idx[i],"=0",sep="")))
		}
	}
	fx = c(fx, fixed.pars)
	spec = arfimaspec(mean.model = list(armaOrder = c(arm, mam),
					include.mean = include.mean, arfima = arfima, 
					external.regressors = external.regressors), 
					distribution.model = distribution.model, 
					start.pars = start.pars, fixed.pars = fx)
	return(spec)
}

.getarfimaspec = function(object)
{
	spec = arfimaspec(mean.model = list(armaOrder = c(object@model$modelinc[2], object@model$modelinc[3]),
					include.mean = object@model$modelinc[1], 
					arfima = object@model$modelinc[4], external.regressors = object@model$modeldata$mexdata), 
			distribution.model = object@model$modeldesc$distribution, start.pars  = object@model$start.pars, 
			fixed.pars = object@model$fixed.pars)
	return(spec)
}
setMethod(f = "getspec", signature(object = "ARFIMAfit"), definition = .getarfimaspec)

.setfixedarfima = function(object, value){
	# get parameter values
	model = object@model
	ipars = model$pars
	pars = unlist(value)
	names(pars) = parnames = tolower(names(pars))
	# included parameters in model
	modelnames = rownames(ipars[which(ipars[,4]==1 | ipars[,2]==1), ])
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
	# set parameter values
	tmp = arfimaspec(mean.model = list(armaOrder = c(model$modelinc[2], model$modelinc[3]), 
					include.mean = model$modelinc[1], 
					arfima = model$modelinc[4], external.regressors = model$modeldata$mexdata), 
			distribution.model = model$modeldesc$distribution, start.pars  = model$start.pars, 
			fixed.pars = as.list(fixed.pars))
	# ToDo: Need to check that the parameters are not outside the bounds...
	idx = which(is.na(tmp@model$pars[,"LB"]))
	tmp@model$pars[idx,"LB"] = object@model$pars[idx,"LB"]
	idx = which(is.na(tmp@model$pars[,"UB"]))
	tmp@model$pars[idx,"UB"] = object@model$pars[idx,"UB"]
	return(tmp)
}
setReplaceMethod(f="setfixed", signature= c(object = "ARFIMAspec", value = "vector"), definition = .setfixedarfima)

.setstartarfima = function(object, value){
	# get parameter values
	model = object@model
	ipars = model$pars
	pars = unlist(value)
	names(pars) = parnames = tolower(names(pars))
	# included parameters in model
	modelnames = rownames(ipars[which(ipars[,4]==1 | ipars[,2]==1), ])
	inc = NULL
	for(i in seq_along(parnames)){
		if(is.na(match(parnames[i], modelnames))){
			warning( (paste("Unrecognized Parameter in Fixed Values: ", parnames[i], "...Ignored", sep = "")))
		} else{
			inc = c(inc, i)
		}
	}
	start.pars = pars[inc]
	names(start.pars) = tolower(names(pars[inc]))
	# set parameter values
	tmp = arfimaspec(mean.model = list(armaOrder = c(model$modelinc[2], model$modelinc[3]), 
					include.mean = model$modelinc[1], 
					arfima = model$modelinc[4], external.regressors = model$modeldata$mexdata), 
			distribution.model = model$modeldesc$distribution, fixed.pars  = model$fixed.pars, 
			start.pars = as.list(start.pars))
	# ToDo: Need to check that the parameters are not outside the bounds...
	idx = which(is.na(tmp@model$pars[,"LB"]))
	tmp@model$pars[idx,"LB"] = object@model$pars[idx,"LB"]
	idx = which(is.na(tmp@model$pars[,"UB"]))
	tmp@model$pars[idx,"UB"] = object@model$pars[idx,"UB"]
	return(tmp)
}

setReplaceMethod(f="setstart", signature= c(object = "ARFIMAspec", value = "vector"), definition = .setstartarfima)

setReplaceMethod(f="setbounds", signature= c(object = "ARFIMAspec", value = "vector"), definition = .setbounds)

arfimafit = function(spec, data, out.sample = 0, solver = "solnp", solver.control = list(), 
		fit.control = list(fixed.se = 0, scale = 0), 
		numderiv.control = list(grad.eps=1e-4, grad.d=0.0001, grad.zero.tol=sqrt(.Machine$double.eps/7e-7),
				hess.eps=1e-4, hess.d=0.1, hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2), ...)
{
	UseMethod("arfimafit")
}

setMethod("arfimafit", signature(spec = "ARFIMAspec"), .arfimafit)

arfimafilter = function(spec, data, out.sample = 0, n.old = NULL, ...)
{
	UseMethod("arfimafilter")
}

setMethod("arfimafilter", signature(spec = "ARFIMAspec"), .arfimafilter)


arfimaforecast = function(fitORspec, data = NULL, n.ahead = 10, n.roll = 0, out.sample = 0, 
		external.forecasts = list(mregfor = NULL), ...)
{
	UseMethod("arfimaforecast")
}

setMethod("arfimaforecast", signature(fitORspec = "ARFIMAfit"), .arfimaforecast)

setMethod("arfimaforecast", signature(fitORspec = "ARFIMAspec"), .arfimaforecast2)

arfimasim = function(fit, n.sim = 1000, n.start = 0, m.sim = 1, startMethod = c("unconditional","sample"), 
		prereturns = NA, preresiduals = NA, rseed = NA, custom.dist = list(name = NA, distfit = NA, type = "z"), 
		mexsimdata = NULL, ...)
{
	UseMethod("arfimasim")
}


setMethod("arfimasim", signature(fit = "ARFIMAfit"), .arfimasim)

arfimapath = function(spec, n.sim = 1000, n.start = 0, m.sim = 1, prereturns = NA, preresiduals = NA, 
		rseed = NA,  custom.dist = list(name = NA, distfit = NA, type = "z"), mexsimdata = NULL, ...)
{
	UseMethod("arfimapath")
}


setMethod("arfimapath", signature(spec = "ARFIMAspec"), .arfimapath)


arfimaroll = function(spec, data, n.ahead = 1, forecast.length = 500, 
		n.start = NULL, refit.every = 25, refit.window = c("recursive", "moving"), 
		window.size = NULL, solver = "hybrid", fit.control = list(), 
		solver.control = list(), calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05), 
		cluster = NULL, keep.coef = TRUE, ...)
{
	setMethod("arfimaroll")
}

setMethod("arfimaroll", signature(spec = "ARFIMAspec"),  definition = .arfimaroll)

setMethod("resume", signature(object = "ARFIMAroll"),  definition = .resumeroll2)


arfimadistribution = function(fitORspec, n.sim = 2000, n.start = 1, 
		m.sim = 100, recursive = FALSE, recursive.length = 6000, recursive.window = 1000, 
		prereturns = NA, preresiduals = NA, rseed = NA, 
		custom.dist = list(name = NA, distfit = NA, type = "z"), 
		mexsimdata = NULL, fit.control = list(), solver = "solnp", 
		solver.control = list(), cluster = NULL, ...)
{
	setMethod("arfimadistribution")
}
setMethod("arfimadistribution", signature(fitORspec = "ARFIMAfit"), .arfimadistribution)
setMethod("arfimadistribution", signature(fitORspec = "ARFIMAspec"), .arfimadistribution)

setMethod("show",
		signature(object = "ARFIMAspec"),
		function(object){
			model = object@model
			modelinc = model$modelinc
			cat(paste("\n*----------------------------------*", sep = ""))
			cat(paste("\n*       ARFIMA Model Spec          *", sep = ""))
			cat(paste("\n*----------------------------------*", sep = ""))
			cat("\nConditional Mean Dynamics")
			cat(paste("\n------------------------------------\n",sep=""))
			cat("Mean Model\t\t\t: ARFIMA(", modelinc[2],",", ifelse(modelinc[4]>0, "d", 0),",",modelinc[3],")\n", sep = "")
			cat("Include Mean\t\t:", as.logical(modelinc[1]),"\n")
			if(modelinc[6]>0) cat(paste("Exogenous Regressor Dimension: ", modelinc[6],"\n",sep=""))
			cat("\nConditional Distribution")
			cat(paste("\n------------------------------------\n",sep=""))
			cat("Distribution\t: ", model$modeldesc$distribution,"\n")
			cat("Includes Skew\t: ", as.logical(modelinc[16]),"\n")
			cat("Includes Shape\t: ", as.logical(modelinc[17]),"\n")
			cat("Includes Lambda\t: ", as.logical(modelinc[18]),"\n\n")
			return(invisible(object))
		})

# fit show
# fit show
setMethod("show",
		signature(object = "ARFIMAfit"),
		function(object){
			model = object@model
			modelinc = model$modelinc
			cat(paste("\n*----------------------------------*", sep = ""))
			cat(paste("\n*          ARFIMA Model Fit        *", sep = ""))
			cat(paste("\n*----------------------------------*", sep = ""))
			cat("\nMean Model\t: ARFIMA(", modelinc[2],",",ifelse(modelinc[4]>0, "d", 0),",",modelinc[3],")\n", sep = "")
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
				stdresid = object@fit$residuals/coef(object)["sigma"]
				itestm = infocriteria(object)
				cat("\nInformation Criteria")
				cat(paste("\n------------------------------------\n",sep=""))
				print(itestm,digits=5)
				cat("\nWeighted Ljung-Box Test on Standardized Residuals")
				cat(paste("\n------------------------------------\n",sep=""))
				tmp1 = .weightedBoxTest(stdresid, p = 1, df = sum(modelinc[2:3]))
				print(tmp1, digits = 4)
				cat("\nH0 : No serial correlation\n")
				cat("\nWeighted Ljung-Box Test on Standardized Squared Residuals")
				cat(paste("\n------------------------------------\n",sep=""))
				tmp2 = .weightedBoxTest(stdresid, p = 2, df = sum(modelinc[8:9]))
				print(tmp2, digits = 4)
				cat("\n\nARCH LM Tests")
				cat(paste("\n------------------------------------\n",sep=""))
				L2 = .archlmtest(stdresid, lags = 2)
				L5 = .archlmtest(stdresid, lags = 5)
				L10 = .archlmtest(stdresid, lags = 10)
				alm = matrix(0,ncol = 3,nrow = 3)
				alm[1,1:3] = c(L2$statistic, L2$parameter, L2$p.value)
				alm[2,1:3] = c(L5$statistic, L5$parameter, L5$p.value)
				alm[3,1:3] = c(L10$statistic, L10$parameter, L10$p.value)
				colnames(alm) = c("Statistic", "DoF", "P-Value")
				rownames(alm) = c("ARCH Lag[2]", "ARCH Lag[5]", "ARCH Lag[10]")
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
				cat("\nElapsed time :", object@fit$timer,"\n\n")
			} else{
				cat("\nConvergence Problem:")
				cat("\nSolver Message:", object@fit$message,"\n\n")
				
			}
			return(invisible(object))
		})
# filter show
setMethod("show",
		signature(object = "ARFIMAfilter"),
		function(object){
			model = object@model
			modelinc = object@model$modelinc
			cat(paste("\n*-------------------------------------*", sep = ""))
			cat(paste("\n*          ARFIMA Model Filter        *", sep = ""))
			cat(paste("\n*-------------------------------------*", sep = ""))
			cat("\nMean Model\t: ARFIMA(", modelinc[2],",",ifelse(modelinc[4]>0, "d", 0),",",modelinc[3],")\n", sep = "")
			cat("Distribution\t:", model$modeldesc$distribution,"\n")
			cat("\nFilter Parameters")
			cat(paste("\n---------------------------------------\n",sep=""))
			print(matrix(coef(object), ncol=1, dimnames = list(names(coef(object)), "")), digits = 5)
			cat("\nLogLikelihood :", object@filter$LLH, "\n")
			stdresid = object@filter$residuals/object@model$pars["sigma", 1]
			itestm = infocriteria(object)
			cat("\nInformation Criteria")
			cat(paste("\n---------------------------------------\n",sep=""))
			print(itestm,digits=5)
			cat("\nQ-Statistics on Standardized Residuals")
			cat(paste("\n---------------------------------------\n",sep=""))
			tmp1 = .box.test(stdresid, p = 1, df = sum(modelinc[2:3]))
			print(tmp1, digits = 4)
			cat("\nH0 : No serial correlation\n")
			cat("\nQ-Statistics on Standardized Squared Residuals")
			cat(paste("\n---------------------------------------\n",sep=""))
			tmp2 = .box.test(stdresid, p = 2, df = sum(modelinc[2:3]))
			print(tmp2, digits = 4)
			cat("\nARCH LM Tests")
			cat(paste("\n---------------------------------------\n",sep=""))
			L2 = .archlmtest(stdresid, lags = 2)
			L5 = .archlmtest(stdresid, lags = 5)
			L10 = .archlmtest(stdresid, lags = 10)
			alm = matrix(0,ncol = 3,nrow = 3)
			alm[1,1:3] = c(L2$statistic, L2$parameter, L2$p.value)
			alm[2,1:3] = c(L5$statistic, L5$parameter, L5$p.value)
			alm[3,1:3] = c(L10$statistic, L10$parameter, L10$p.value)
			colnames(alm) = c("Statistic", "DoF", "P-Value")
			rownames(alm) = c("ARCH Lag[2]", "ARCH Lag[5]", "ARCH Lag[10]")
			print(alm,digits = 4)
			cat("\n\n")
			return(invisible(object))
		})
# sim show
setMethod("show",
		signature(object = "ARFIMAsim"),
		function(object){
			model = object@model
			modelinc = model$modelinc
			cat(paste("\n*-------------------------------------*", sep = ""))
			cat(paste("\n*       ARFIMA Model Simulation       *", sep = ""))
			cat(paste("\n*-------------------------------------*", sep = ""))
			sim = object@simulation
			dates = object@model$dates
			series = sim$seriesSim
			resids = sim$residSim
			m = dim(series)[2]
			N = dim(series)[1]
			cat(paste("\nHorizon: ",N))
			cat(paste("\nSimulations: ",m,"\n",sep=""))
			rx1 = apply(series, 2, FUN=function(x) mean(x))
			rx2 = apply(series, 2, FUN=function(x) range(x))
			T = object@model$modeldata$T
			xspec = .model2spec(as.list(object@model$pars[object@model$pars[,3]==1,1]), object@model, type = "ARFIMA")
			actual = c(0, mean(object@model$modeldata$data[1:T]), 
					min(object@model$modeldata$data[1:T]), max(object@model$modeldata$data[1:T]))
			uncond = c(0, uncmean(xspec), NA, NA)
			dd = data.frame(Seed = object@seed, Series.Mean = rx1, Series.Min = rx2[1,],
					Series.Max = rx2[2,])
			meansim = apply(dd, 2, FUN = function(x) mean(x))
			meansim[1] = 0
			dd = rbind(dd, meansim, actual, uncond)
			rownames(dd) = c(paste("sim", 1:m, sep = ""), "Mean(All)", "Actual", "Unconditional")
			print(dd,digits = 3)
			cat("\n\n")
			return(invisible(object))
		})
		
# forecast show
setMethod("show",
		signature(object = "ARFIMAforecast"),
		function(object){
			cat(paste("\n*----------------------------------*", sep = ""))
			cat(paste("\n*        ARFIMA Model Forecast     *", sep = ""))
			cat(paste("\n*----------------------------------*", sep = ""))
			n.ahead = object@forecast$n.ahead
			cat(paste("\n\nHorizon: ", n.ahead, sep = ""))
			cat(paste("\nRoll Steps: ",object@forecast$n.roll, sep = ""))
			n.start = object@forecast$n.start
			if(n.start>0) infor = ifelse(n.ahead>n.start, n.start, n.ahead) else infor = 0
			cat(paste("\nOut of Sample: ", infor, "\n", sep = ""))
			cat("\n0-roll forecast: \n")
			zz = object@forecast$seriesFor[,1]
			print(zz, digits = 4)
			cat("\n\n")
			return(invisible(object))
		})

# path show
setMethod("show",
		signature(object = "ARFIMApath"),
		function(object){
			model = object@model
			modelinc = model$modelinc
			cat(paste("\n*--------.........------------------------*", sep = ""))
			cat(paste("\n*        ARFIMA Model Path Simulation     *", sep = ""))
			cat(paste("\n*-----------------------------------------*", sep = ""))
			sim = object@path
			series = sim$seriesSim
			resids = sim$residSim
			m = dim(series)[2]
			N = dim(series)[1]
			cat(paste("\n\nHorizon: ", N))
			cat(paste("\nSimulations: ", m, "\n", sep = ""))
			T = object@model$modeldata$T
			xspec = .model2spec(as.list(object@model$pars[object@model$pars[,3]==1,1]), object@model, type = "ARFIMA")
			uncond = c(0, uncmean(xspec), NA, NA)
			rx1 = apply(series, 2, FUN = function(x) mean(x))
			rx2 = apply(series, 2, FUN = function(x) range(x))
			dd = data.frame(Seed = object@seed, Series.Mean = rx1, Series.Min = rx2[1,], 
					Series.Max = rx2[2,])
			meansim = apply(dd, 2, FUN = function(x) mean(x))
			meansim[1] = 0			
			dd = rbind(dd, meansim, uncond)
			rownames(dd) = c(paste("sim", 1:m, sep = ""), "Mean(All)", "Unconditional")			
			print(dd, digits = 3)
			cat("\n\n")
			return(invisible(object))
		})

setMethod("show",
		signature(object = "ARFIMAroll"),
		function(object){
			if(!is.null(object@model$noncidx)){
				cat("\nObject containts non-converged estimation windows. Use resume method to re-estimate.\n")
				return(invisible(object))
			} else{
				cat(paste("\n*--------------------------------------*", sep = ""))
				cat(paste("\n*              ARFIMA Roll             *", sep = ""))
				cat(paste("\n*--------------------------------------*", sep = ""))
				N = object@model$n.refits
				model = object@model$spec@model
				modelinc = model$modelinc
				cat("\nNo.Refits\t\t:", N)
				cat("\nRefit Horizon\t:", object@model$refit.every)
				cat("\nNo.Forecasts\t:", NROW(object@forecast$density))
				cat("\nDistribution\t:", model$modeldesc$distribution,"\n")
				cat("\nForecast Density:\n")
				print(round(head(object@forecast$density),4))
				cat("\n..........................\n")
				print(round(tail(object@forecast$density),4))
				cat("\n\n")
				return(invisible(object))
			}
		})
# distribution show
# distribution show
setMethod("show",
		signature(object = "ARFIMAdistribution"),
		function(object){
			cat(paste("\n*-------------------------------------*", sep = ""))
			cat(paste("\n*    ARFIMA Parameter Distribution    *", sep = ""))
			cat(paste("\n*-------------------------------------*", sep = ""))
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
			return(invisible(object))
		})

#-------------------------------------------------------------------------
# multi-methods
setMethod("show",
		signature(object = "ARFIMAmultispec"),
		function(object){
			cat(paste("\n*---------------------------------*", sep = ""))
			cat(paste("\n*        ARFIMA Multi-Spec        *", sep = ""))
			cat(paste("\n*---------------------------------*", sep = ""))
			N = length(object@spec)
			cat(paste("\n\nMultiple Specifications\t: ", N, sep=""))
			cat(paste("\nMulti-Spec Type\t\t\t: ", object@type, sep=""))
			cat("\n\n")
			return(invisible(object))
		})		

setMethod("show",
		signature(object = "ARFIMAmultifit"),
		function(object){
			cat(paste("\n*--------------------------------*", sep = ""))
			cat(paste("\n*       ARFIMA Multi-Fit         *", sep = ""))
			cat(paste("\n*--------------------------------*", sep = ""))
			cat(paste("\n\nNo. Assets :", length(object@fit), sep=""))
			asset.names = object@desc$asset.names
			if(object@desc$type == "equal"){
				cat(paste("\nMulti-Spec Type : Equal",sep=""))
				cat(paste("\n\nModel Spec",sep=""))
				cat(paste("\n-------------------------------\n",sep=""))
				cat("\nInclude Mean\t: ", as.logical( object@fit[[1]]@model$modelinc[1] ) )
				cat(paste("\nAR(FI)MA Model  : (",object@fit[[1]]@model$modelinc[2],",",
								ifelse(object@fit[[1]]@model$modelinc[4]>0, 1, "d"),
								",",object@fit[[1]]@model$modelinc[3],")",sep=""))
				if(object@fit[[1]]@model$modelinc[6]>0){
					cat("\nExogenous Regressors : ", object@fit[[1]]@model$modelinc[6])
				} else{
					cat("\nExogenous Regressors : none")
				}
				cat(paste("\nConditional Distribution: ",object@fit[[1]]@model$modeldesc$distribution,"\n", sep=""))
				cv = sapply(object@fit, FUN = function(x) x@fit$convergence)
				if(any(cv != 0)){
					ncv = which(cv != 0)
					nncv = length(ncv)
					cat("\nNo. of non converged fits: ", ncv,"\n")
					if(ncv>0) cat("\nNon converged fits: ", nncv,"\n")
					
				} else{
					cat(paste("\nModel Fit", sep = ""))
					cat(paste("\n-------------------------------\n",sep=""))
					cat("\n")
					ll = t(likelihood(object))
					rownames(ll) = "Log-Lik"
					cf = coef(object)
					colnames(cf) = asset.names
					print(round(rbind(cf, ll), digits = 5))
					cat("\n")
				}
			} else{
				cat(paste("\nARFIMA Model Fit", sep = ""))
				cat(paste("\n-------------------------------\n",sep=""))
				cat(paste("\nModel Fit", sep = ""))
				cat(paste("\n-------------------------------\n",sep=""))
				cat("\n")
				print(coef(object), digits = 5)
			}
			cat("\n\n")
			return(invisible(object))
		})

setMethod("show",
		signature(object = "ARFIMAmultifilter"),
		function(object){
			cat(paste("\n*--------------------------------*", sep = ""))
			cat(paste("\n*       ARFIMA Multi-Filter         *", sep = ""))
			cat(paste("\n*--------------------------------*", sep = ""))
			cat(paste("\n\nNo. Assets :", length(object@filter), sep=""))
			asset.names = object@desc$asset.names
			if(object@desc$type == "equal"){
				cat(paste("\nMulti-Spec Type : Equal",sep=""))
				cat(paste("\n\nModel Spec",sep=""))
				cat(paste("\n-------------------------------\n",sep=""))
				cat("\nInclude Mean\t: ", as.logical( object@filter[[1]]@model$modelinc[1] ) )
				cat(paste("\nAR(FI)MA Model : (",object@filter[[1]]@model$modelinc[2],",",
								ifelse(object@filter[[1]]@model$modelinc[4]>0, 1, "d"),
								",",object@filter[[1]]@model$modelinc[3],")",sep=""))
				if(object@filter[[1]]@model$modelinc[6]>0){
					cat("\nExogenous Regressors in mean equation: ", object@filter[[1]]@model$modelinc[6])
				} else{
					cat("\nExogenous Regressors in mean equation: none")
				}
				cat("\nConditional Distribution: ",object@filter[[1]]@model$modeldesc$distribution,"\n")			
				cat(paste("\nModel Filter", sep = ""))
				cat(paste("\n-------------------------------\n",sep=""))
				cat("\n")
				ll = t(likelihood(object))
				rownames(ll) = "Log-Lik"
				cf = coef(object)
				colnames(cf) = asset.names
				print(round(rbind(cf, ll), digits = 5))
				cat("\n")
			} else{
				cat(paste("\nARFIMA Model Filter", sep = ""))
				cat(paste("\n-------------------------------\n",sep=""))
				cat(paste("\nModel Fit", sep = ""))
				cat(paste("\n-------------------------------\n",sep=""))
				cat("\n")
				print(coef(object), digits = 5)
			}
			cat("\n\n")
			return(invisible(object))
		})

setMethod("show",
		signature(object = "ARFIMAmultiforecast"),
		function(object){
			asset.names = object@desc$asset.names
			cat(paste("\n*----------------------------------*", sep = ""))
			cat(paste("\n*       ARFIMA Multi-Forecast      *", sep = ""))
			cat(paste("\n*----------------------------------*", sep = ""))
			cat(paste("\n\nNo. Assets :", length(object@forecast), sep=""))
			cat(paste("\n--------------------------\n",sep=""))
			cat("\n\n")
			return(invisible(object))
		})

#----------------------------------------------------------------------------------
# univariate fit extractors
#----------------------------------------------------------------------------------
# coef methods
.arfimacoef = function(object)
{
	switch(class(object)[1],
			ARFIMAfit = object@fit$coef,
			ARFIMAfilter = object@model$pars[object@model$pars[,2]==1, 1],
			ARFIMAmultifit = {
				if(object@desc$type == "equal"){
					sapply(object@fit, FUN = function(x) coef(x), simplify = TRUE)
				} else{
					lapply(object@fit, FUN = function(x) coef(x))
				}
			},
			ARFIMAmultifilter = {
				if(object@desc$type == "equal"){
					sapply(object@filter, FUN = function(x) coef(x), simplify = TRUE)
				} else{
					lapply(object@filter, FUN = function(x) coef(x))
				}
			},
			ARFIMAroll = {
				if(!is.null(object@model$noncidx)) stop("\nObject containts non-converged estimation windows.")
				object@model$coef
			})
}

setMethod("coef", signature(object = "ARFIMAfit"), .arfimacoef)
setMethod("coef", signature(object = "ARFIMAfilter"), .arfimacoef)
setMethod("coef", signature(object = "ARFIMAmultifit"), .arfimacoef)
setMethod("coef", signature(object = "ARFIMAmultifilter"), .arfimacoef)
setMethod("coef", signature(object = "ARFIMAroll"), .arfimacoef)
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# Fitted method
.arfimafitted = function(object)
{
	if(class(object)[1] == "ARFIMAfit" | class(object)[1] == "ARFIMAfilter"){
		D = object@model$modeldata$index[1:object@model$modeldata$T]
	}
	switch(class(object)[1],
			ARFIMAfit = xts(object@fit$fitted.values, D), 
			ARFIMAfilter = xts(object@model$modeldata$data[1:object@model$modeldata$T] - object@filter$residuals, D),
			ARFIMAmultifilter = sapply(object@filter, FUN = function(x) fitted(x), simplify = TRUE),
			ARFIMAmultifit = sapply(object@fit, FUN = function(x) fitted(x), simplify = TRUE),
			ARFIMAsim = {
				ans = object@simulation$seriesSim
				rownames(ans) = paste("T+",1:NROW(object@simulation$seriesSim), sep="")
				return(ans)
			},
			ARFIMApath ={
				ans = object@path$seriesSim
				rownames(ans) = paste("T+",1:NROW(object@path$seriesSim), sep="")
				return(ans)
			},
			ARFIMAforecast = object@forecast$seriesFor
	)
}
setMethod("fitted", signature(object = "ARFIMAfit"), .arfimafitted)
setMethod("fitted", signature(object = "ARFIMAfilter"), .arfimafitted)
setMethod("fitted", signature(object = "ARFIMAmultifit"), .arfimafitted)
setMethod("fitted", signature(object = "ARFIMAmultifilter"), .arfimafitted)
setMethod("fitted", signature(object = "ARFIMAsim"), .arfimafitted)
setMethod("fitted", signature(object = "ARFIMApath"), .arfimafitted)
setMethod("fitted", signature(object = "ARFIMAforecast"), .arfimafitted)

.arfimafittedmf = function(object)
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
setMethod("fitted", signature(object = "ARFIMAmultiforecast"), .arfimafittedmf)

#----------------------------------------------------------------------------------
# as.data.frame method for distribution object
.arfimadistdf = function(x, row.names = NULL, optional = FALSE, which = "coef", window = 1, ...)
{
	n = x@dist$details$nwindows
	if(window > n) stop("window size greater than actual available", call. = FALSE)
	
	if(which == "rmse"){
		ans = as.data.frame(t(x@dist[[window]]$rmse))
		colnames(ans) = rownames(x@truecoef)
	}
	
	if(which == "stats"){	
		llh = x@dist[[window]]$likelist
		uncmean = x@dist[[window]]$mlongrun
		maxret = x@dist[[window]]$simmaxdata[,1]
		minret = x@dist[[window]]$simmindata[,1]
		meanret = x@dist[[window]]$simmeandata[,1]
		kurtosis = x@dist[[window]]$simmomdata[,1]
		skewness = x@dist[[window]]$simmomdata[,2]
		ans = data.frame(llh = llh, uncmean = uncmean, maxret = maxret, minret = minret, 
				meanret = meanret, kurtosis = kurtosis, skewness  = skewness)
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

setMethod("as.data.frame", signature(x = "ARFIMAdistribution"), .arfimadistdf)
#----------------------------------------------------------------------------------
# as.data.frame method for bootstrap object
.arfimabootdf = function(x, row.names = NULL, optional = FALSE, type = "raw", qtile = c(0.01, 0.099))
{
	n.ahead = x@model$n.ahead
	if(type == "raw"){
		series = x@fseries
		ans = data.frame(bootseries = series)
		colnames(ans) = paste("t+", 1:n.ahead, sep="")
	}
	if(type == "q"){
		if(all(is.numeric(qtile)) && (all(qtile<1.0) && all(qtile >0.0))){
			series = x@fseries
			ans = apply(series, 2, FUN = function(x) quantile(x, qtile))
			ans = as.data.frame(ans)
			colnames(ans) = paste("t+", 1:n.ahead, sep="")
			rownames(ans) = paste("q", qtile, sep = "")
		} else{
			stop("\nfor type q, the qtile value must be numeric and between (>)0 and 1(<)\n", call.  = FALSE)
		} 
	}
	if(type == "summary"){
		series = x@fseries
		ans = apply(series, 2, FUN = function(x) c(min(x), quantile(x, 0.25), mean(x), quantile(x, 0.75), max(x) ))
		ans = as.data.frame(ans)
		colnames(ans) = paste("t+", 1:n.ahead, sep="")
		rownames(ans) = c("min", "q.25", "mean", "q.75", "max")
	}
	ans
}
#setMethod("as.data.frame", signature(x = "ARFIMAboot"), .arfimabootdf)


#----------------------------------------------------------------------------------
# as.data.frame method for roll object
.arfimarolldf = function(x, row.names = NULL, optional = FALSE, which = "density")
{
	if(!is.null(x@model$noncidx)) stop("\nObject containts non-converged estimation windows.")
	if(which == "density") ans =  x@forecast$density else ans = x@forecast$VaR
	return(ans)
}
setMethod("as.data.frame", signature(x = "ARFIMAroll"), .arfimarolldf)

#----------------------------------------------------------------------------------
# residuals method
.arfimaresids = function(object, standardize = FALSE)
{
	if(class(object)[1] == "ARFIMAfit" | class(object)[1] == "ARFIMAfilter"){
		D = object@model$modeldata$index[1:object@model$modeldata$T]
		s = object@model$pars["sigma",1]
	}
	if(standardize){
		ans = switch(class(object)[1],
				ARFIMAfit = xts(object@fit$residuals/s, D),
				ARFIMAfilter = xts(object@filter$residuals/s, D),
				ARFIMAmultifit = sapply(object@fit, FUN = function(x) residuals(x, standardize = TRUE), simplify = TRUE),
				ARFIMAmultifilter = sapply(object@filter, FUN = function(x) residuals(x, standardize = TRUE), simplify = TRUE))
	} else{
		ans = switch(class(object)[1],
				ARFIMAfit = xts(object@fit$residuals, D),
				ARFIMAfilter = xts(object@filter$residuals, D),
				ARFIMAmultifit = sapply(object@fit, FUN = function(x) residuals(x), simplify = TRUE),
				ARFIMAmultifilter = sapply(object@filter, FUN = function(x) residuals(x), simplify = TRUE))
	}
	return(ans)
}

setMethod("residuals", signature(object = "ARFIMAfit"), .arfimaresids)
setMethod("residuals", signature(object = "ARFIMAfilter"), .arfimaresids)
setMethod("residuals", signature(object = "ARFIMAmultifit"), .arfimaresids)
setMethod("residuals", signature(object = "ARFIMAmultifilter"), .arfimaresids)

#----------------------------------------------------------------------------------
# Likelihood method
.arfimaLikelihood = function(object)
{
	switch(class(object)[1],
			ARFIMAfit = object@fit$LLH,
			ARFIMAfilter = object@filter$LLH,
			ARFIMAmultifilter = sapply(object@filter, FUN = function(x) likelihood(x), simplify = TRUE),
			ARFIMAmultifit = sapply(object@fit, FUN = function(x) likelihood(x), simplify = TRUE))
}

setMethod("likelihood", signature(object = "ARFIMAfit"), .arfimaLikelihood)
setMethod("likelihood", signature(object = "ARFIMAfilter"), .arfimaLikelihood)
setMethod("likelihood", signature(object = "ARFIMAmultifilter"), .arfimaLikelihood)
setMethod("likelihood", signature(object = "ARFIMAmultifit"), .arfimaLikelihood)
#----------------------------------------------------------------------------------
# Info Criteria method
.arfimainfocriteria = function(object)
{
	if(is(object, "ARFIMAfit")){
		np = sum(object@fit$ipars[,4])
	} else{
		np = length(coef(object))
	}
	itest = .information.test(likelihood(object), nObs = NROW(fitted(object)), 
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

setMethod("infocriteria", signature(object = "ARFIMAfit"), .arfimainfocriteria)
setMethod("infocriteria", signature(object = "ARFIMAfilter"), .arfimainfocriteria)


#----------------------------------------------------------------------------------
# The mult- methods
#----------------------------------------------------------------------------------
.multispecarfima = function( speclist )
{
	# first create a spec which goes through validation process
	tp = 1
	ans = new("ARFIMAmultispec", 
			spec = speclist,
			type = "equal")
	# then check type
	n = length(speclist)
	for(i in 2:n){
		modelnames1 = rownames( speclist[[i]]@model$pars[speclist[[i]]@model$pars[,3]==1, ] )
		modelnames2 = rownames( speclist[[i-1]]@model$pars[speclist[[i-1]]@model$pars[,3]==1, ] )
		if(length(modelnames1) != length(modelnames2)){
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
	ans = new("ARFIMAmultispec",
			spec = speclist,
			type = type)
	return(ans)
}

setMethod("multifit", signature(multispec = "ARFIMAmultispec"),  definition = .multifitarfima)
setMethod("multifilter", signature(multifitORspec = "ARFIMAmultifit"),  definition = .multifilterarfima1)
setMethod("multifilter", signature(multifitORspec = "ARFIMAmultispec"),  definition = .multifilterarfima2)
setMethod("multiforecast", signature(multifitORspec = "ARFIMAmultifit"),  definition = .multiforecastarfima1)
setMethod("multiforecast", signature(multifitORspec = "ARFIMAmultispec"),  definition = .multiforecastarfima2)
# Unconditional Mean
.unconditionalmean11 = function(object, method = c("analytical", "simulation"), n.sim = 20000, rseed = NA)
{
	method = method[1]
	N = object@model$modeldata$T
	if(is(object, "ARFIMAfit")) pars = object@fit$ipars[,1] else pars = object@filter$ipars[,1]
	if(method == "analytical"){
		modelinc = object@model$modelinc
		idx = object@model$pidx
		if(modelinc[6]>0){
			mxreg = pars[idx["mxreg",1]:idx["mxreg",2]]
			mexdata = matrix(object@model$modeldata$mexdata[1:N, ], ncol = modelinc[6])
			meanmex = apply(mexdata, 2, "mean")
			umeanmex = sum(mxreg*meanmex)
		} else{
			umeanmex = 0
		}
		if(modelinc[1]>0) mu = pars[idx["mu",1]] else mu=0
		umean = (mu + umeanmex)
		return(umean)
	} else{
		if(is(object, "ARFIMAfit")){
			sim = arfimasim(object, n.sim = n.sim, n.start = 1000, startMethod = "sample", rseed = rseed)
			umean = mean(as.vector(sim@simulation$seriesSim))
			return(umean)
		} else{
			stop("\nuncmean by simulation not available for ARFIMAfilter class object (used spec instead).")
		}
	}
}

.unconditionalmean21 = function(object, method = c("analytical", "simulation"), n.sim = 20000, rseed = NA)
{
	method = method[1]
	if(is.null(object@model$fixed.pars)) stop("uncmean with ARFIMAspec requires fixed.pars list", call. = FALSE)
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
			mxreg = pars[idx["mxreg",1]:idx["mxreg",2]]
			meanmex = apply(object@model$modeldata$mexdata, 2, "mean")
			umeanmex = sum(mxreg*meanmex)
		} else{
			umeanmex = 0
		}
		if(modelinc[1]>0) mu = pars[idx["mu",1]] else mu=0
		umean = (mu + umeanmex)
		return(umean)
	}  else{
		sim = arfimapath(object, n.sim = n.sim, n.start = 1000, rseed = rseed)
		umean = mean(as.vector(sim@path$seriesSim))
		return(umean)
	}
}
setMethod("uncmean", signature(object = "ARFIMAfit"),    definition = .unconditionalmean11)
setMethod("uncmean", signature(object = "ARFIMAfilter"), definition = .unconditionalmean11)
setMethod("uncmean", signature(object = "ARFIMAspec"),   definition = .unconditionalmean21)

# forecast performance measures
.arfimarollreport = function(object, type = "VaR", VaR.alpha = 0.01, conf.level = 0.95)
{
	if(!is.null(object@model$noncidx)) stop("\nObject containts non-converged estimation windows.")
	switch(type,
			VaR = .rollVaRreport2(object, VaR.alpha, conf.level),
			fpm = .rollfpmreport2(object))
	invisible(object)
}

setMethod("report", signature(object = "ARFIMAroll"), .arfimarollreport)

.rollfpmreport2 = function(object)
{
	vmodel = object@model$spec@model$modeldesc$vmodel
	vsubmodel = object@model$spec@model$modeldesc$vsubmodel
	cat(paste("\nARFIMA Roll Mean Forecast Performance Measures", sep = ""))
	cat(paste("\n---------------------------------------------", sep = ""))
	cat(paste("\nNo.Refits\t: ", object@model$n.refits, sep = ""))
	cat(paste("\nNo.Forecasts: ", NROW(object@forecast$density), sep = ""))
	cat("\n\n")
	tmp = fpm(object)
	print(signif(tmp, 4))
	cat("\n\n")
}

.rollVaRreport2 = function(object, VaR.alpha = 0.01, conf.level = 0.95)
{
	VaR.alpha = VaR.alpha[1]
	v.a = object@model$VaR.alpha
	if(!object@model$calculate.VaR) stop("\nplot-->error: VaR was not calculated for this object\n", call.=FALSE)
	if(!any(v.a==VaR.alpha[1])) stop("\nplot-->error: VaR.alpha chosen is invalid for the object\n", call.=FALSE)
	dvar = object@forecast$VaR
	m = NROW(dvar)
	idx = match(VaR.alpha, v.a)
	.VaRreport(if(is.null(object@model$datanames)) "" else object@model$datanames, "ARFIMA", 
			object@model$spec@model$modeldesc$distribution, 
			p = VaR.alpha, actual = as.numeric(dvar[,"realized"]), 
			VaR = dvar[, idx], 
			conf.level = conf.level)
	invisible(object)
}

.getspec2 = function(object)
{
	spec = arfimaspec(
			mean.model = list(armaOrder = c(object@model$modelinc[2],object@model$modelinc[3]),
					include.mean = object@model$modelinc[1],
					arfima = object@model$modelinc[4], external.regressors = object@model$modeldata$mexdata),
			distribution.model = object@model$modeldesc$distribution, fixed.pars = object@model$fixed.pars)
	# should custom bounds be propagated?
	#idx = which(is.na(tmp@model$pars[,"LB"]))
	#tmp@model$pars[idx,"LB"] = object@model$pars[idx,"LB"]
	#idx = which(is.na(tmp@model$pars[,"UB"]))
	#tmp@model$pars[idx,"UB"] = object@model$pars[idx,"UB"]
	return(spec)
}

setMethod(f = "getspec", signature(object = "ARFIMAfit"), definition = .getspec2)

#######################
setMethod("convergence", signature(object = "ARFIMAfit"),  definition = .convergence)
setMethod("vcov", signature(object = "ARFIMAfit"),  definition = .vcov)

setMethod("fpm", signature(object = "ARFIMAforecast"),  definition = .fpm1)
setMethod("fpm", signature(object = "ARFIMAroll"),  definition = .fpm2)

setMethod("reduce", signature(object = "ARFIMAfit"), .reduce)
