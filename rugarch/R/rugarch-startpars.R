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

# The starting parameters for the mean equation should be estimated first
# The specification will depend on which aspects of the mean equation are
# chosen i.e. arfima/inmean/arma/
# For the time being, we ignore fixed parameters in the initialization process
# and treat all as non-fixed (although arima function provides for a fixed pars
# option, I am not sure what to do with the lm (I think the 'offset' options does 
# something...)

# Care should be taken with using arima0 and arima with regards to
# 1: parameter names (arima returns xreg while arima0 returns mexdata),
# i.e. one returns the dataframe name and the other the colname.
# 2: using fixed parameters and the subsequent warnings (not used right now)

.arfimastart = function(pars, arglist)
{
	.eps = .Machine$double.eps
	data = arglist$data
	dscale = arglist$dscale
	model = arglist$model
	start.pars = model$start.pars
	start.names = names(start.pars)
	T = model$modeldata$T
	fixed.pars = model$fixed.pars
	fixed.names = names(fixed.pars)
	idx = model$pidx
	modelinc = model$modelinc
	mexdata = model$modeldata$mexdata[1:T, , drop=FALSE]
	if(!is.null(mexdata)) colnames(mexdata)  = NULL
	# account for archex
	if(modelinc[20]>0 && modelinc[6]>0){
		mexdata[,(modelinc[6]-modelinc[20]+1):modelinc[6]] = mexdata[,(modelinc[6]-modelinc[20]+1):modelinc[6]]*sd(data)
	}
	if( modelinc[2]!=0 | modelinc[3]!=0 ){
		ow <- options("warn")
		options(warn = -1)
		ttemp = arima(data, order = c(modelinc[2], 0, modelinc[3]), include.mean = modelinc[1], xreg = mexdata, method = "CSS")
		options(ow)
		fit.mean = ttemp$coef
		#res=ttemp$residuals
		if(modelinc[1]>0){
			pars[idx["mu",1]:idx["mu",2],1] = fit.mean["intercept"]
		}
		if(modelinc[2]!=0) pars[idx["ar",1]:idx["ar",2],1] = fit.mean[c(paste("ar",1:modelinc[2],sep=""))]
		if(modelinc[3]!=0) pars[idx["ma",1]:idx["ma",2],1] = fit.mean[c(paste("ma",1:modelinc[3],sep=""))]
		if(modelinc[6]>0){
			i = which(substr(names(fit.mean), 1, 7) == "mexdata")
			pars[idx["mxreg",1]:idx["mxreg",2],1] = fit.mean[i]
		}
	} else{
		y = data
		if( modelinc[6]>0 ){
			if(modelinc[1]>0){
				fit.mean = lm(y~mexdata)
				i = which(substr(names(fit.mean$coef), 1, 7) == "mexdata")
				pars[idx["mu",1]:idx["mu",2],1] = fit.mean$coef["(Intercept)"]
				pars[idx["mxreg",1]:idx["mxreg",2],1] = fit.mean$coef[i]
			} else{
				fit.mean = lm(y~mexdata-1)
				i = which(substr(names(fit.mean$coef), 1, 7) == "mexdata")
				pars[idx["mxreg",1]:idx["mxreg",2],1] = fit.mean$coef[i]
			}
		} else{
			pars[idx["mu",1]:idx["mu",2],1] = 0
			res = data
		}
	}
	
	if(modelinc[1]>0){
		# Need to control for this special case which sometimes occurs
		if(mean(data)==0){
			if(is.na(pars[idx["mu", 1]:idx["mu", 2], 5])) pars[idx["mu", 1]:idx["mu", 2], 5] = -2*sd(data)
			if(is.na(pars[idx["mu", 1]:idx["mu", 2], 6])) pars[idx["mu", 1]:idx["mu", 2], 6] =  2*sd(data)
		} else{
			if(is.na(pars[idx["mu", 1]:idx["mu", 2], 5])) pars[idx["mu", 1]:idx["mu", 2], 5] = -100*abs(mean(data))
			if(is.na(pars[idx["mu", 1]:idx["mu", 2], 6])) pars[idx["mu", 1]:idx["mu", 2], 6] =  100*abs(mean(data))
		}
		if(!is.null(start.pars$mu)) pars[idx["mu", 1]:idx["mu", 2], 1] = start.pars$mu[1]/dscale
		if(any(substr(fixed.names, 1, 2)=="mu")){
			pars[idx["mu", 1]:idx["mu", 2], 1] = as.numeric(fixed.pars$mu)
			pars[idx["mu", 1]:idx["mu", 2], 5] = fixed.pars$mu
			pars[idx["mu", 1]:idx["mu", 2], 6] = fixed.pars$mu
		}
	}
	
	# ar (we changed the naming of darfima into arfima which creates some extra problems
	# to be caught when using the substr function
	if(modelinc[2]>0){
		arnames=paste("ar",1:modelinc[2],sep="")
		# ar1 stationarity bounds else...
		if(modelinc[2]==1){
			pxd = which(is.na(pars[idx["ar", 1]:idx["ar", 2], 5]))
			if(length(pxd)>0) pars[(idx["ar", 1]:idx["ar", 2])[pxd], 5] = -1+TinY
			pxd = which(is.na(pars[idx["ar", 1]:idx["ar", 2], 6]))
			if(length(pxd)>0) pars[(idx["ar", 1]:idx["ar", 2])[pxd], 6] =  1+TinY
		} else{
			pxd = which(is.na(pars[idx["ar", 1]:idx["ar", 2], 5]))
			if(length(pxd)>0) pars[(idx["ar", 1]:idx["ar", 2])[pxd], 5] = -4+TinY
			pxd = which(is.na(pars[idx["ar", 1]:idx["ar", 2], 6]))
			if(length(pxd)>0) pars[(idx["ar", 1]:idx["ar", 2])[pxd], 6] =  4+TinY
		}
		sp = na.omit(match(start.names, arnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[arnames[sp[i]], 1] = as.numeric(start.pars[arnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, arnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[arnames[sp[i]], 1] = as.numeric(fixed.pars[arnames[sp[i]]])
				pars[arnames[sp[i]], 5] = as.numeric(fixed.pars[arnames[sp[i]]])
				pars[arnames[sp[i]], 6] = as.numeric(fixed.pars[arnames[sp[i]]])
			}
		}
	}

	# ma
	if(modelinc[3]>0){
		manames = paste("ma",1:modelinc[3],sep="")
		if(modelinc[3]==1){
			pxd = which(is.na(pars[idx["ma", 1]:idx["ma", 2], 5]))
			if(length(pxd)>0) pars[(idx["ma", 1]:idx["ma", 2])[pxd], 5] = -1+TinY
			pxd = which(is.na(pars[idx["ma", 1]:idx["ma", 2], 6]))
			if(length(pxd)>0) pars[(idx["ma", 1]:idx["ma", 2])[pxd], 6] =  1+TinY
		} else{
			pxd = which(is.na(pars[idx["ma", 1]:idx["ma", 2], 5]))
			if(length(pxd)>0) pars[(idx["ma", 1]:idx["ma", 2])[pxd], 5] = -4+TinY
			pxd = which(is.na(pars[idx["ma", 1]:idx["ma", 2], 6]))
			if(length(pxd)>0) pars[(idx["ma", 1]:idx["ma", 2])[pxd], 6] =  4+TinY
		}
		sp = na.omit(match(start.names, manames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[manames[sp[i]], 1] = as.numeric(start.pars[manames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, manames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[manames[sp[i]], 1] = as.numeric(fixed.pars[manames[sp[i]]])
				pars[manames[sp[i]], 5] = as.numeric(fixed.pars[manames[sp[i]]])
				pars[manames[sp[i]], 6] = as.numeric(fixed.pars[manames[sp[i]]])
			}
		}
	}

	# garch in mean
	if(modelinc[5]>0){
		pars[idx["archm", 1]:idx["archm", 2], 5] = -10
		pars[idx["archm", 1]:idx["archm", 2], 6] =  10
		if(!is.null(start.pars$archm)) pars[idx["archm", 1]:idx["archm", 2], 1] = start.pars$archm[1]
		if(any(substr(fixed.names, 1, 5)=="archm")){
			pars[idx["archm", 1]:idx["archm", 2], 1] = as.numeric(fixed.pars$archm)
			pars[idx["archm", 1]:idx["archm", 2], 5] = fixed.pars$archm
			pars[idx["archm", 1]:idx["archm", 2], 6] = fixed.pars$archm
		}
	}
	
	# arfima
	if(modelinc[4]>0){
		if(is.na(pars[idx["arfima", 1]:idx["arfima", 2], 5])) pars[idx["arfima", 1]:idx["arfima", 2], 5] = 1e-8
		if(is.na(pars[idx["arfima", 1]:idx["arfima", 2], 6])) pars[idx["arfima", 1]:idx["arfima", 2], 6] = 0.5
		if(is.null(start.pars$arfima)) pars[idx["arfima", 1]:idx["arfima", 2], 1] = .absHmom(data)-0.5 else pars[idx["arfima", 1]:idx["arfima", 2], 1] = start.pars$arfima[1]
		if(any(substr(fixed.names, 1, 6)=="arfima")){
			pars[idx["arfima", 1]:idx["arfima", 2], 1] = as.numeric(fixed.pars$arfima)
			pars[idx["arfima", 1]:idx["arfima", 2], 5] = as.numeric(fixed.pars$arfima)
			pars[idx["arfima", 1]:idx["arfima", 2], 6] = as.numeric(fixed.pars$arfima)
		}
	}

	# exogenous regressors
	if(modelinc[6]>0){
		mxnames = paste("mxreg",1:modelinc[6],sep="")
		pxd = which(is.na(pars[idx["mxreg", 1]:idx["mxreg", 2], 5]))
		if(length(pxd)>0) pars[(idx["mxreg", 1]:idx["mxreg", 2])[pxd], 5] = -100
		pxd = which(is.na(pars[idx["mxreg", 1]:idx["mxreg", 2], 6]))
		if(length(pxd)>0) pars[(idx["mxreg", 1]:idx["mxreg", 2])[pxd], 6] =  100
		sp = na.omit(match(start.names, mxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[mxnames[sp[i]], 1] = as.numeric(start.pars[mxnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, mxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[mxnames[sp[i]], 1] = as.numeric(fixed.pars[mxnames[sp[i]]])
				pars[mxnames[sp[i]], 5] = as.numeric(fixed.pars[mxnames[sp[i]]])
				pars[mxnames[sp[i]], 6] = as.numeric(fixed.pars[mxnames[sp[i]]])
			}
		}
	}
	if(modelinc[7]>0){
		.sd = sd(data)
		# enforce positivity of lower bound
		if(is.na(pars[idx["sigma", 1]:idx["sigma", 2], 5])) pars[idx["sigma", 1]:idx["sigma", 2], 5] = .eps else pars[idx["sigma", 1]:idx["sigma", 2], 5] = abs(pars[idx["sigma", 1]:idx["sigma", 2], 5])
		if(is.na(pars[idx["sigma", 1]:idx["sigma", 2], 6])) pars[idx["sigma", 1]:idx["sigma", 2], 6] = 100*.sd
		if(is.null(start.pars$omega)) pars[idx["sigma", 1]:idx["sigma", 2], 1] = .sd else pars[idx["sigma", 1]:idx["sigma", 2], 1] = abs(start.pars$sigma[1])/dscale
		if(any(substr(fixed.names, 1, 5) == "sigma")){
			pars[idx["sigma", 1]:idx["sigma", 2], 1] = as.numeric(fixed.pars$sigma)
			pars[idx["sigma", 1]:idx["sigma", 2], 5] = fixed.pars$sigma
			pars[idx["sigma", 1]:idx["sigma", 2], 6] = fixed.pars$sigma
		}
	}
	
	dbounds = .DistributionBounds(distribution = model$modeldesc$distribution)
	if(modelinc[16]>0){
		if(is.na(pars[idx["skew", 1]:idx["skew", 2], 5])) pars[idx["skew", 1]:idx["skew", 2], 5] = dbounds$skew.LB
		if(is.na(pars[idx["skew", 1]:idx["skew", 2], 6])) pars[idx["skew", 1]:idx["skew", 2], 6] = dbounds$skew.UB
		if(is.null(start.pars$skew)) pars[idx["skew", 1]:idx["skew", 2], 1] = dbounds$skew else pars[idx["skew", 1]:idx["skew", 2], 1] = start.pars$skew[1]
		if(any(substr(fixed.names, 1, 4) == "skew")){
			pars[idx["skew", 1]:idx["skew", 2], 1] = as.numeric(fixed.pars$skew)
			pars[idx["skew", 1]:idx["skew", 2], 5] = as.numeric(fixed.pars$skew)
			pars[idx["skew", 1]:idx["skew", 2], 6] = as.numeric(fixed.pars$skew)
		}
	}
	if(modelinc[17]>0){
		if(is.na(pars[idx["shape", 1]:idx["shape", 2], 5])) pars[idx["shape", 1]:idx["shape", 2], 5] = dbounds$shape.LB
		if(is.na(pars[idx["shape", 1]:idx["shape", 2], 6])) pars[idx["shape", 1]:idx["shape", 2], 6] = dbounds$shape.UB
		if(is.null(start.pars$shape)) pars[idx["shape", 1]:idx["shape", 2], 1] = dbounds$shape else pars[idx["shape", 1]:idx["shape", 2], 1] = start.pars$shape[1]
		if(any(substr(fixed.names, 1, 5) == "shape")){
			pars[idx["shape", 1]:idx["shape", 2], 1] = as.numeric(fixed.pars$shape)
			pars[idx["shape", 1]:idx["shape", 2], 5] = as.numeric(fixed.pars$shape)
			pars[idx["shape", 1]:idx["shape", 2], 6] = as.numeric(fixed.pars$shape)
		}
	}
	if(modelinc[18]>0){
		if(is.na(pars[idx["ghlambda", 1]:idx["ghlambda", 2], 5])) pars[idx["ghlambda", 1]:idx["ghlambda", 2], 5] = dbounds$ghlambda.LB
		if(is.na(pars[idx["ghlambda", 1]:idx["ghlambda", 2], 6])) pars[idx["ghlambda", 1]:idx["ghlambda", 2], 6] = dbounds$ghlambda.UB
		if(is.null(start.pars$ghlambda)) pars[idx["ghlambda", 1]:idx["ghlambda", 2], 1] = dbounds$ghlambda else pars[idx["ghlambda", 1]:idx["ghlambda", 2], 1] = start.pars$ghlambda[1]
		if(any(substr(fixed.names, 1, 8) == "ghlambda")){
			pars[idx["ghlambda", 1]:idx["ghlambda", 2], 1] = as.numeric(fixed.pars$ghlambda)
			pars[idx["ghlambda", 1]:idx["ghlambda", 2], 5] = as.numeric(fixed.pars$ghlambda)
			pars[idx["ghlambda", 1]:idx["ghlambda", 2], 6] = as.numeric(fixed.pars$ghlambda)
		}
	}
	return( pars )
}

# [mu ar ma arfima im mxreg omega alpha beta gamma gamma11 gamma21 delta lambda vxreg skew shape dlamda aux aux aux aux]

.meqstartpars = function(pars, arglist)
{
	# Case 1 garchInMean yes:
	data = arglist$data
	model = arglist$model
	N = length(as.numeric(unlist(data)))
	modelinc = model$modelinc
	modeldesc = model$modeldesc
	idx = model$pidx
	mxreg = numeric()
	if(modelinc[6] > 0){
		mexdata = model$modeldata$mexdata[1:N,, drop = FALSE]
		# account for archex
		if(modelinc[20]>0){
			mexdata[,(modelinc[6]-modelinc[20]+1):modelinc[6]] = mexdata[,(modelinc[6]-modelinc[20]+1):modelinc[6]]*sd(data)
		}
		# easier to NULL the names of the data rather than search them by name later
		# (search now is 'mexdata'..1:mxn)
		colnames(mexdata) = NULL
	} else{
		mexdata = NULL
	}
	tmph = 0
	# get the sigma vector should it be needed (garchInMean)
	if(modelinc[5] > 0 && is.null(model$start.pars$archm)){
		ctrl = list(eval.max = 1000, trace = 0, iter.max = 1000)
		tempspec = ugarchspec(variance.model = list(model = modeldesc$vmodel, garchOrder = c(modelinc[8], modelinc[9]), 
						submodel = modeldesc$vsubmodel), mean.model = list(armaOrder = c(modelinc[2], modelinc[3]),
						include.mean = modelinc[1], archm = FALSE), distribution.model = "norm")
		# we try first with stationarity conditions enforced
		tempfit = ugarchfit(data = data, spec = tempspec, solver="solnp", solver.control = ctrl, 
				fit.control = list(stationarity = 1))
		if(tempfit@fit$convergence!=0){
			tempfit = ugarchfit(data = data, spec = tempspec, solver = "solnp", solver.control = ctrl, fit.control = list(stationarity = 0))
			if(tempfit@fit$convergence!=0) stop("\nugarchfit-->error: could not find appropriate starting values for recursion\n")
		}
		tmph = tempfit@fit$sigma
	}
	# arima without garchInMean
	if( (modelinc[5] == 0 ||  !is.null(model$start.pars$archm)) && (modelinc[2]>0 | modelinc[3]>0)){
		ttemp = arima(data, order = c(modelinc[2], 0, modelinc[3]), include.mean = modelinc[1], 
				xreg = mexdata, method = "CSS")
		fit.mean = ttemp$coef
		#res=ttemp$residuals
		if(modelinc[1]>0){
			pars[idx["mu", 1]:idx["mu", 2], 1] = fit.mean["intercept"]
		}
		if(modelinc[2]>0) pars[idx["ar", 1]:idx["ar", 2], 1] = fit.mean[c(paste("ar",1:modelinc[2],sep=""))]
		if(modelinc[3]>0) pars[idx["ma", 1]:idx["ma", 2], 1] = fit.mean[c(paste("ma",1:modelinc[3],sep=""))]		
		if(modelinc[6]>0){
			i = which(substr(names(fit.mean), 1, 7) == "mexdata")
			pars[idx["mxreg", 1]:idx["mxreg", 2], 1] = fit.mean[i]
		}
	}
	
	# arima with garchInMean	
	if((modelinc[5]>0) && (modelinc[2]+modelinc[3])>0 && is.null(model$start.pars$archm)){
			mexdata = cbind(mexdata, tmph^modelinc[5])
			mxn = modelinc[6]+1
			colnames(mexdata) = paste("xreg", 1:mxn,sep="")
			ttemp = arima0(data, order = c(modelinc[2], 0, modelinc[3]), include.mean = modelinc[1], xreg = mexdata)
			fit.mean = ttemp$coef
			if(modelinc[1]>0){
				pars[idx["mu", 1]:idx["mu", 2], 1] = fit.mean["intercept"]
			}
			if(modelinc[2]>0) pars[idx["ar", 1]:idx["ar", 2], 1] = fit.mean[c(paste("ar",1:modelinc[2],sep=""))]
			if(modelinc[3]>0) pars[idx["ma", 1]:idx["ma", 2], 1] = fit.mean[c(paste("ma",1:modelinc[3],sep=""))]	
			i = which(substr(names(fit.mean), 1, 4) == "xreg")
			z = length(i) # at a minimum it is 2 ex+inmean
			if(modelinc[6]>0){
				pars[idx["mxreg", 1]:idx["mxreg", 2], 1] = fit.mean[i[1:(z-1)]]
			}
			pars[idx["archm", 1]:idx["archm", 2], 1] = fit.mean[i[z]]
	}
	
	# lm for garchInMean without arma
	if((modelinc[5]>0) && (modelinc[2] == 0 && modelinc[3] == 0) && is.null(model$start.pars$archm)){
			mexdata = cbind(mexdata, tmph^modelinc[5])
			mxn = modelinc[6]
			y = data
			if(modelinc[1]>0){
				fit.mean = lm(y~mexdata)
				pars[idx["mu", 1]:idx["mu", 2], 1] = fit.mean$coef["(Intercept)"]				
				i = which(substr(names(fit.mean$coef), 1, 7) == "mexdata")
				z = length(i) # at a minimum it is 2 ex+inmean
				if(modelinc[6]>0){
					pars[idx["mxreg", 1]:idx["mxreg", 2], 1] = fit.mean$coef[i[1:(z-1)]]
				}
				pars[idx["archm", 1]:idx["archm", 2], 1] = fit.mean$coef[i[z]]				
				#res=as.numeric(fit.mean$residuals)
			} else{
				fit.mean = lm(y~mexdata-1)
				if(modelinc[6]>0) pars[idx["mxreg", 1]:idx["mxreg", 2], 1] = fit.mean$coef[c(paste("mexdata",1:(mxn-1),sep=""))]			
				pars[idx["archm", 1]:idx["archm", 2], 1] = fit.mean$coef[c(paste("mexdata", mxn,sep=""))]
				#res=as.numeric(fit.mean$residuals)
			}
	}
	
	if(modelinc[5]==0 && modelinc[2] == 0 && modelinc[3] == 0){
		y = data
		if(modelinc[6]>0){
			mxn = modelinc[6]
			fit.mean = lm(y~mexdata)
			i = which(substr(names(fit.mean$coef), 1, 7) == "mexdata")
			if(modelinc[1]>0){
				pars[idx["mu", 1]:idx["mu", 2], 1] = fit.mean$coef["(Intercept)"]
			}
			pars[idx["mxreg", 1]:idx["mxreg", 2], 1] = fit.mean$coef[i]
			#res=as.numeric(fit.mean$residuals)
		} else{
			if(modelinc[1]>0){
				pars[idx["mu", 1]:idx["mu", 2], 1] = mean(y)
			}
			#res=(data-mu)
		}
	}
	# assign("tmph", tmph, envir = garchenv)
	return(list(pars = pars, tmph = tmph))
}

# common to all specifications is the mean equation:
.meqstart = function(pars, arglist)
{
	data = arglist$data
	dscale = arglist$dscale
	model = arglist$model
	start.pars = model$start.pars
	start.names = names(start.pars)
	
	fixed.pars = model$fixed.pars
	fixed.names = names(fixed.pars)
	
	idx = model$pidx
	modelinc = model$modelinc
	
	# this is where we fix the bounds for the fixed parameters
	# fill and then the fixed.names to overwrite it...also for the
	# bounds
	
	if(modelinc[1]>0){
		# Need to control for this special case which sometimes occurs
		if(mean(data)==0){
			if(is.na(pars[idx["mu", 1]:idx["mu", 2], 5])) pars[idx["mu", 1]:idx["mu", 2], 5] = -0.5
			if(is.na(pars[idx["mu", 1]:idx["mu", 2], 6])) pars[idx["mu", 1]:idx["mu", 2], 6] =  0.5
		} else{			
			if(is.na(pars[idx["mu", 1]:idx["mu", 2], 5])) pars[idx["mu", 1]:idx["mu", 2], 5] = -100*abs(mean(data))
			if(is.na(pars[idx["mu", 1]:idx["mu", 2], 6])) pars[idx["mu", 1]:idx["mu", 2], 6] =  100*abs(mean(data))
		}
			if(!is.null(start.pars$mu)) pars[idx["mu", 1]:idx["mu", 2], 1] = start.pars$mu[1]/dscale
			if(any(substr(fixed.names, 1, 2)=="mu")){
				pars[idx["mu", 1]:idx["mu", 2], 1] = as.numeric(fixed.pars$mu)
				pars[idx["mu", 1]:idx["mu", 2], 5] = fixed.pars$mu
				pars[idx["mu", 1]:idx["mu", 2], 6] = fixed.pars$mu
			}
	}
	
	# ar (we changed the naming of darfima into arfima which creates some extra problems
	# to be caught when using the substr function
	if(modelinc[2]>0){
		arnames=paste("ar",1:modelinc[2],sep="")
		# ar1 stationarity bounds else...
		if(modelinc[2]==1){
			pxd = which(is.na(pars[idx["ar", 1]:idx["ar", 2], 5]))
			if(length(pxd)>0) pars[(idx["ar", 1]:idx["ar", 2])[pxd], 5] = -1+TinY
			pxd = which(is.na(pars[idx["ar", 1]:idx["ar", 2], 6]))
			if(length(pxd)>0) pars[(idx["ar", 1]:idx["ar", 2])[pxd], 6] =  1+TinY
		} else{
			pxd = which(is.na(pars[idx["ar", 1]:idx["ar", 2], 5]))
			if(length(pxd)>0) pars[(idx["ar", 1]:idx["ar", 2])[pxd], 5] = -4+TinY
			pxd = which(is.na(pars[idx["ar", 1]:idx["ar", 2], 6]))
			if(length(pxd)>0) pars[(idx["ar", 1]:idx["ar", 2])[pxd], 6] =  4+TinY
		}
		sp = na.omit(match(start.names, arnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[arnames[sp[i]], 1] = as.numeric(start.pars[arnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, arnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[arnames[sp[i]], 1] = as.numeric(fixed.pars[arnames[sp[i]]])
				pars[arnames[sp[i]], 5] = as.numeric(fixed.pars[arnames[sp[i]]])
				pars[arnames[sp[i]], 6] = as.numeric(fixed.pars[arnames[sp[i]]])
			}
		}
	}
	# ma
	if(modelinc[3]>0){
		manames = paste("ma",1:modelinc[3],sep="")
		if(modelinc[3]==1){
			pxd = which(is.na(pars[idx["ma", 1]:idx["ma", 2], 5]))
			if(length(pxd)>0) pars[(idx["ma", 1]:idx["ma", 2])[pxd], 5] = -1+TinY
			pxd = which(is.na(pars[idx["ma", 1]:idx["ma", 2], 6]))
			if(length(pxd)>0) pars[(idx["ma", 1]:idx["ma", 2])[pxd], 6] =  1+TinY
		} else{
			pxd = which(is.na(pars[idx["ma", 1]:idx["ma", 2], 5]))
			if(length(pxd)>0) pars[(idx["ma", 1]:idx["ma", 2])[pxd], 5] = -4+TinY
			pxd = which(is.na(pars[idx["ma", 1]:idx["ma", 2], 6]))
			if(length(pxd)>0) pars[(idx["ma", 1]:idx["ma", 2])[pxd], 6] =  4+TinY
		}
		sp = na.omit(match(start.names, manames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[manames[sp[i]], 1] = as.numeric(start.pars[manames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, manames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[manames[sp[i]], 1] = as.numeric(fixed.pars[manames[sp[i]]])
				pars[manames[sp[i]], 5] = as.numeric(fixed.pars[manames[sp[i]]])
				pars[manames[sp[i]], 6] = as.numeric(fixed.pars[manames[sp[i]]])
			}
		}
	}
	
	# garch in mean
	if(modelinc[5]>0){
		if(is.na(pars[idx["archm", 1]:idx["archm", 2], 5])) pars[idx["archm", 1]:idx["archm", 2], 5] = -10
		if(is.na(pars[idx["archm", 1]:idx["archm", 2], 6])) pars[idx["archm", 1]:idx["archm", 2], 6] =  10
		if(!is.null(start.pars$archm)) pars[idx["archm", 1]:idx["archm", 2], 1] = start.pars$archm[1]
		if(any(substr(fixed.names, 1, 5)=="archm")){
			pars[idx["archm", 1]:idx["archm", 2], 1] = as.numeric(fixed.pars$archm)
			pars[idx["archm", 1]:idx["archm", 2], 5] = fixed.pars$archm
			pars[idx["archm", 1]:idx["archm", 2], 6] = fixed.pars$archm
		}
	}
	
	# arfima
	if(modelinc[4]>0){
		if(is.na(pars[idx["arfima", 1]:idx["arfima", 2], 5])) pars[idx["arfima", 1]:idx["arfima", 2], 5] = 1e-8
		if(is.na(pars[idx["arfima", 1]:idx["arfima", 2], 6])) pars[idx["arfima", 1]:idx["arfima", 2], 6] = 0.5
		if(is.null(start.pars$arfima)) pars[idx["arfima", 1]:idx["arfima", 2], 1] = .absHmom(data)-0.5 else pars[idx["arfima", 1]:idx["arfima", 2], 1] = start.pars$arfima[1]
		if(any(substr(fixed.names, 1, 6)=="arfima")){
			pars[idx["arfima", 1]:idx["arfima", 2], 1] = as.numeric(fixed.pars$arfima)
			pars[idx["arfima", 1]:idx["arfima", 2], 5] = as.numeric(fixed.pars$arfima)
			pars[idx["arfima", 1]:idx["arfima", 2], 6] = as.numeric(fixed.pars$arfima)
		}
	}
	
	# exogenous regressors
	if(modelinc[6]>0){
		mxnames = paste("mxreg",1:modelinc[6],sep="")
		pxd = which(is.na(pars[idx["mxreg", 1]:idx["mxreg", 2], 5]))
		if(length(pxd)>0) pars[(idx["mxreg", 1]:idx["mxreg", 2])[pxd], 5] = -100
		pxd = which(is.na(pars[idx["mxreg", 1]:idx["mxreg", 2], 6]))
		if(length(pxd)>0) pars[(idx["mxreg", 1]:idx["mxreg", 2])[pxd], 6] =  100
		sp = na.omit(match(start.names, mxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[mxnames[sp[i]], 1] = as.numeric(start.pars[mxnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, mxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[mxnames[sp[i]], 1] = as.numeric(fixed.pars[mxnames[sp[i]]])
				pars[mxnames[sp[i]], 5] = as.numeric(fixed.pars[mxnames[sp[i]]])
				pars[mxnames[sp[i]], 6] = as.numeric(fixed.pars[mxnames[sp[i]]])
			}
		}
	}
	return( pars )
}

# starting parameters s.t. specification
.garchstart = function(pars, arglist)
{
	tmp = .meqstartpars(pars, arglist)
	pars = tmp$pars
	tmph = tmp$temph
	ans = switch(arglist$model$modeldesc$vmodel,
			sGARCH = .aparchstart(pars, arglist),
			fGARCH = .fgarchstart(pars, arglist),
			gjrGARCH = .aparchstart(pars, arglist),
			apARCH = .aparchstart(pars, arglist),
			eGARCH = .egarchstart(pars, arglist),
			iGARCH = .igarchstart(pars, arglist),
			csGARCH = .csgarchstart(pars, arglist),
			hyGARCH = .hygarchstart(pars, arglist),
			mcsGARCH = .mcaparchstart(pars, arglist),
			realGARCH = .realgarchstart(pars, arglist))
	#anstGARCH = .anstgarchstart(arglist))
	return(list(pars = ans, tmph = tmph))
}
# [mu ar ma arfima im mxreg omega alpha beta gamma gamma11 gamma21 delta lambda vxreg skew shape dlamda aux aux aux aux]

.csgarchstart = function(pars, arglist)
{
	data = arglist$data
	model = arglist$model
	dscale = arglist$dscale
	modelinc = model$modelinc
	start.pars = model$start.pars
	fixed.pars = model$fixed.pars
	idx = model$pidx
	fixed.names = names(fixed.pars)
	start.names = names(start.pars)
	pars = .meqstart(pars, arglist)
	if(modelinc[7]>0){
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 5])) pars[idx["omega", 1]:idx["omega", 2], 5] = 1e-12
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 6])) pars[idx["omega", 1]:idx["omega", 2], 6] = 1
		if(is.null(start.pars$omega)) pars[idx["omega", 1]:idx["omega", 2], 1] = (var(data, na.rm = TRUE))/1000 else 
			pars[idx["omega", 1]:idx["omega", 2], 1] = start.pars$omega[1]/dscale
		if(any(substr(fixed.names, 1, 5) == "omega")){
			pars[idx["omega", 1]:idx["omega", 2], 1] = as.numeric(fixed.pars$omega)
			pars[idx["omega", 1]:idx["omega", 2], 5] = fixed.pars$omega
			pars[idx["omega", 1]:idx["omega", 2], 6] = fixed.pars$omega
		}
	}
	if(modelinc[8]>0){
		gpnames = paste("alpha",1:modelinc[8],sep="")
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 5]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 6]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 6] =  1-TinY
		pars[idx["alpha", 1]:idx["alpha", 2], 1] = rep(0.05/modelinc[8], modelinc[8])
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	if(modelinc[9] > 0){
		gqnames = paste("beta",1:modelinc[9],sep="")
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 5]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 6]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 6] =  1-TinY
		pars[idx["beta", 1]:idx["beta", 2], 1] = rep(0.7/modelinc[9], modelinc[9])
		sp = na.omit(match(start.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gqnames[sp[i]], 1] = as.numeric(start.pars[gqnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gqnames[sp[i]], 1] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 5] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 6] = as.numeric(fixed.pars[gqnames[sp[i]]])
			}
		}
	}
	# \rho in the paper notation
	if(modelinc[11] > 0){
		gqnames = "eta11"
		if(is.na(pars[idx["eta1", 1]:idx["eta1", 2], 5])) pars[idx["eta1", 1]:idx["eta1", 2], 5] = TinY
		if(is.na(pars[idx["eta1", 1]:idx["eta1", 2], 6])) pars[idx["eta1", 1]:idx["eta1", 2], 6] = 1-TinY
		pars[idx["eta1", 1]:idx["eta1", 2], 1] = 0.98
		if(any(substr(start.names, 1, 4) == "eta1")){
			j = which(substr(start.names, 1, 4) == "eta1")
			gqmatch = charmatch(start.names[j],gqnames)
			pars[gqnames[gqmatch], 1] = as.numeric(start.pars[j])
		}
		if(any(substr(fixed.names, 1, 4) == "eta1")){
			j = which(substr(fixed.names, 1, 4) == "eta1")
			gqmatch = charmatch(fixed.names[j],gqnames)
			pars[gqnames[gqmatch], 1] = as.numeric(fixed.pars[j])
			pars[gqnames[gqmatch], 5] = as.numeric(fixed.pars[j])
			pars[gqnames[gqmatch], 6] = as.numeric(fixed.pars[j])
		}
	}
	# \phi in the paper notation
	if(modelinc[12] > 0){
		gqnames = "eta21"
		if(is.na(pars[idx["eta2", 1]:idx["eta2", 2], 5])) pars[idx["eta2", 1]:idx["eta2", 2], 5] = TinY
		if(is.na(pars[idx["eta2", 1]:idx["eta2", 2], 6])) pars[idx["eta2", 1]:idx["eta2", 2], 6] = 1-TinY
		pars[idx["eta2", 1]:idx["eta2", 2], 1] = 0.05
		if(any(substr(start.names, 1, 4) == "eta2")){
			j = which(substr(start.names, 1, 4) == "eta2")
			gqmatch = charmatch(start.names[j],gqnames)
			pars[gqnames[gqmatch], 1] = as.numeric(start.pars[j])
		}
		if(any(substr(fixed.names, 1, 4) == "eta2")){
			j = which(substr(fixed.names, 1, 4) == "eta2")
			gqmatch = charmatch(fixed.names[j],gqnames)
			pars[gqnames[gqmatch], 1] = as.numeric(fixed.pars[j])
			pars[gqnames[gqmatch], 5] = as.numeric(fixed.pars[j])
			pars[gqnames[gqmatch], 6] = as.numeric(fixed.pars[j])
		}
	}
	if(modelinc[15]>0){
		vxnames = paste("vxreg",1:modelinc[15],sep="")
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 5]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 6]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 6] = 100		
		pars[idx["vxreg", 1]:idx["vxreg", 2], 1] = rep(TinY, modelinc[15])
		sp = na.omit(match(start.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[vxnames[sp[i]], 1] = as.numeric(start.pars[vxnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[vxnames[sp[i]], 1] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 5] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 6] = as.numeric(fixed.pars[vxnames[sp[i]]])
			}
		}
	}
	pars = .distributionstart(pars, model, start.pars, fixed.pars)
	return( pars )
}


.realgarchstart = function(pars, arglist)
{
	data = arglist$data
	realized = arglist$realized
	model = arglist$model
	dscale = arglist$dscale
	modelinc = model$modelinc
	start.pars = model$start.pars
	fixed.pars = model$fixed.pars
	idx = model$pidx
	fixed.names = names(fixed.pars)
	start.names = names(start.pars)
	pars = .meqstart(pars, arglist)
	
	tmp = cbind(log(realized[1:length(data)]), log(data^2))
	tmp[!is.finite(tmp[,1]),1] = NA
	tmp[!is.finite(tmp[,2]),2] = NA
	tmp = na.omit(tmp)
	colnames(tmp)<-c("Real","Sqr")
	mod = lm(Real~Sqr, data = as.data.frame(tmp))
	xistart = coef(mod)[1]
	lambdastart = sd(mod$residuals)/2
	
	if(modelinc[7]>0){
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 5])) pars[idx["omega", 1]:idx["omega", 2], 5] = -9
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 6])) pars[idx["omega", 1]:idx["omega", 2], 6] = 5
		if(is.null(start.pars$omega)) pars[idx["omega", 1]:idx["omega", 2], 1] = -1 else 
			pars[idx["omega", 1]:idx["omega", 2], 1] = start.pars$omega[1]/dscale
		if(any(substr(fixed.names, 1, 5) == "omega")){
			pars[idx["omega", 1]:idx["omega", 2], 1] = as.numeric(fixed.pars$omega)
			pars[idx["omega", 1]:idx["omega", 2], 5] = fixed.pars$omega
			pars[idx["omega", 1]:idx["omega", 2], 6] = fixed.pars$omega
		}
	}
	if(modelinc[8]>0){
		gpnames = paste("alpha",1:modelinc[8],sep="")
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 5]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 6]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 6] =  1
		pars[idx["alpha", 1]:idx["alpha", 2], 1] = rep(0.05/modelinc[8], modelinc[8])
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	if(modelinc[9] > 0){
		gqnames = paste("beta",1:modelinc[9],sep="")
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 5]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 6]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 6] = 1
		pars[idx["beta", 1]:idx["beta", 2], 1] = rep(0.7/modelinc[9], modelinc[9])
		sp = na.omit(match(start.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gqnames[sp[i]], 1] = as.numeric(start.pars[gqnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gqnames[sp[i]], 1] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 5] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 6] = as.numeric(fixed.pars[gqnames[sp[i]]])
			}
		}
	}
	# \rho in the paper notation
	if(modelinc[11] > 0){
		gqnames = "eta11"
		if(is.na(pars[idx["eta1", 1]:idx["eta1", 2], 5])) pars[idx["eta1", 1]:idx["eta1", 2], 5] = -0.5
		if(is.na(pars[idx["eta1", 1]:idx["eta1", 2], 6])) pars[idx["eta1", 1]:idx["eta1", 2], 6] = 0.5
		pars[idx["eta1", 1]:idx["eta1", 2], 1] = 0.1
		if(any(substr(start.names, 1, 4) == "eta1")){
			j = which(substr(start.names, 1, 4) == "eta1")
			gqmatch = charmatch(start.names[j],gqnames)
			pars[gqnames[gqmatch], 1] = as.numeric(start.pars[j])
		}
		if(any(substr(fixed.names, 1, 4) == "eta1")){
			j = which(substr(fixed.names, 1, 4) == "eta1")
			gqmatch = charmatch(fixed.names[j],gqnames)
			pars[gqnames[gqmatch], 1] = as.numeric(fixed.pars[j])
			pars[gqnames[gqmatch], 5] = as.numeric(fixed.pars[j])
			pars[gqnames[gqmatch], 6] = as.numeric(fixed.pars[j])
		}
	}
	# \phi in the paper notation
	if(modelinc[12] > 0){
		gqnames = "eta21"
		if(is.na(pars[idx["eta2", 1]:idx["eta2", 2], 5])) pars[idx["eta2", 1]:idx["eta2", 2], 5] = -0.5
		if(is.na(pars[idx["eta2", 1]:idx["eta2", 2], 6])) pars[idx["eta2", 1]:idx["eta2", 2], 6] = 0.5
		pars[idx["eta2", 1]:idx["eta2", 2], 1] = 0.05
		if(any(substr(start.names, 1, 4) == "eta2")){
			j = which(substr(start.names, 1, 4) == "eta2")
			gqmatch = charmatch(start.names[j],gqnames)
			pars[gqnames[gqmatch], 1] = as.numeric(start.pars[j])
		}
		if(any(substr(fixed.names, 1, 4) == "eta2")){
			j = which(substr(fixed.names, 1, 4) == "eta2")
			gqmatch = charmatch(fixed.names[j],gqnames)
			pars[gqnames[gqmatch], 1] = as.numeric(fixed.pars[j])
			pars[gqnames[gqmatch], 5] = as.numeric(fixed.pars[j])
			pars[gqnames[gqmatch], 6] = as.numeric(fixed.pars[j])
		}
	}
	
	if(modelinc[13] > 0){
		gqnames = "delta"
		if(is.na(pars[idx["delta", 1]:idx["delta", 2], 5])) pars[idx["delta", 1]:idx["delta", 2], 5] = 0.1
		if(is.na(pars[idx["delta", 1]:idx["delta", 2], 6])) pars[idx["delta", 1]:idx["delta", 2], 6] = 5
		pars[idx["delta", 1]:idx["delta", 2], 1] = 1
		if(any(substr(start.names, 1, 5) == "delta")){
			j = which(substr(start.names, 1, 5) == "delta")
			gqmatch = charmatch(start.names[j],gqnames)
			pars[gqnames[gqmatch], 1] = as.numeric(start.pars[j])
		}
		if(any(substr(fixed.names, 1, 5) == "delta")){
			j = which(substr(fixed.names, 1, 5) == "delta")
			gqmatch = charmatch(fixed.names[j],gqnames)
			pars[gqnames[gqmatch], 1] = as.numeric(fixed.pars[j])
			pars[gqnames[gqmatch], 5] = as.numeric(fixed.pars[j])
			pars[gqnames[gqmatch], 6] = as.numeric(fixed.pars[j])
		}
	}
	
	if(modelinc[14] > 0){
		gqnames = "lambda"
		if(is.na(pars[idx["lambda", 1]:idx["lambda", 2], 5])) pars[idx["lambda", 1]:idx["lambda", 2], 5] = TinY
		if(is.na(pars[idx["lambda", 1]:idx["lambda", 2], 6])) pars[idx["lambda", 1]:idx["lambda", 2], 6] = 5
		pars[idx["lambda", 1]:idx["lambda", 2], 1] = lambdastart
		if(any(substr(start.names, 1, 6) == "lambda")){
			j = which(substr(start.names, 1, 6) == "lambda")
			gqmatch = charmatch(start.names[j],gqnames)
			pars[gqnames[gqmatch], 1] = as.numeric(start.pars[j])
		}
		if(any(substr(fixed.names, 1, 6) == "lambda")){
			j = which(substr(fixed.names, 1, 6) == "lambda")
			gqmatch = charmatch(fixed.names[j],gqnames)
			pars[gqnames[gqmatch], 1] = as.numeric(fixed.pars[j])
			pars[gqnames[gqmatch], 5] = as.numeric(fixed.pars[j])
			pars[gqnames[gqmatch], 6] = as.numeric(fixed.pars[j])
		}
	}
	
	
	if(modelinc[15]>0){
		vxnames = paste("vxreg",1:modelinc[15],sep="")
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 5]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 6]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 6] = 100		
		pars[idx["vxreg", 1]:idx["vxreg", 2], 1] = rep(TinY, modelinc[15])
		sp = na.omit(match(start.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[vxnames[sp[i]], 1] = as.numeric(start.pars[vxnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[vxnames[sp[i]], 1] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 5] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 6] = as.numeric(fixed.pars[vxnames[sp[i]]])
			}
		}
	}
	
	if(modelinc[19] > 0){
		gqnames = "xi"
		if(is.na(pars[idx["xi", 1]:idx["xi", 2], 5])) pars[idx["xi", 1]:idx["xi", 2], 5] = -10
		if(is.na(pars[idx["xi", 1]:idx["xi", 2], 6])) pars[idx["xi", 1]:idx["xi", 2], 6] = 10
		pars[idx["xi", 1]:idx["xi", 2], 1] = xistart
		if(any(substr(start.names, 1, 2) == "xi")){
			j = which(substr(start.names, 1, 2) == "xi")
			gqmatch = charmatch(start.names[j],gqnames)
			pars[gqnames[gqmatch], 1] = as.numeric(start.pars[j])
		}
		if(any(substr(fixed.names, 1, 2) == "xi")){
			j = which(substr(fixed.names, 1, 2) == "xi")
			gqmatch = charmatch(fixed.names[j],gqnames)
			pars[gqnames[gqmatch], 1] = as.numeric(fixed.pars[j])
			pars[gqnames[gqmatch], 5] = as.numeric(fixed.pars[j])
			pars[gqnames[gqmatch], 6] = as.numeric(fixed.pars[j])
		}
	}
	pars = .distributionstart(pars, model, start.pars, fixed.pars)
	return( pars )
}

.hygarchstart = function(pars, arglist)
{
	data = arglist$data
	model = arglist$model
	dscale = arglist$dscale
	modelinc = model$modelinc
	start.pars = model$start.pars
	fixed.pars = model$fixed.pars
	idx = model$pidx
	fixed.names = names(fixed.pars)
	start.names = names(start.pars)
	pars = .meqstart(pars, arglist)
	if(modelinc[7]>0){
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 5])) pars[idx["omega", 1]:idx["omega", 2], 5] = var(data)/100000
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 6])) pars[idx["omega", 1]:idx["omega", 2], 6] = var(data)*100000
		if(is.null(start.pars$omega)) pars[idx["omega", 1]:idx["omega", 2], 1] = (var(data, na.rm = TRUE))/1000 else 
			pars[idx["omega", 1]:idx["omega", 2], 1] = start.pars$omega[1]/dscale
		if(any(substr(fixed.names, 1, 5) == "omega")){
			pars[idx["omega", 1]:idx["omega", 2], 1] = as.numeric(fixed.pars$omega)
			pars[idx["omega", 1]:idx["omega", 2], 5] = fixed.pars$omega
			pars[idx["omega", 1]:idx["omega", 2], 6] = fixed.pars$omega
		}
	}
	if(modelinc[8]>0){
		gpnames = paste("alpha",1:modelinc[8],sep="")
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 5]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 5] = -5
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 6]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 6] =  1-TinY
		pars[idx["alpha", 1]:idx["alpha", 2], 1] = rep(0.05/modelinc[8], modelinc[8])
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	if(modelinc[9] > 0){
		gqnames = paste("beta",1:modelinc[9],sep="")
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 5]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 5] = -5
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 6]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 6] =  1-TinY
		pars[idx["beta", 1]:idx["beta", 2], 1] = rep(0.9/modelinc[9], modelinc[9])
		sp = na.omit(match(start.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gqnames[sp[i]], 1] = as.numeric(start.pars[gqnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gqnames[sp[i]], 1] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 5] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 6] = as.numeric(fixed.pars[gqnames[sp[i]]])
			}
		}
	}
	if(modelinc[13]>0){
		if(is.na(pars[idx["delta", 1]:idx["delta", 2], 5])) pars[idx["delta", 1]:idx["delta", 2], 5] = 0
		if(is.na(pars[idx["delta", 1]:idx["delta", 2], 6])) pars[idx["delta", 1]:idx["delta", 2], 6] = 1
		if(is.null(start.pars$delta)) pars[idx["delta", 1]:idx["delta", 2], 1] = 0 else pars[idx["delta", 1]:idx["delta", 2], 1] = start.pars$delta[1]
		if(any(substr(fixed.names, 1, 5) == "delta")){
			pars[idx["delta", 1]:idx["delta", 2], 1] = as.numeric(fixed.pars$delta)
			pars[idx["delta", 1]:idx["delta", 2], 5] = as.numeric(fixed.pars$delta)
			pars[idx["delta", 1]:idx["delta", 2], 6] = as.numeric(fixed.pars$delta)
		}
	}
	if(modelinc[14]>0){
		if(is.na(pars[idx["lambda", 1]:idx["lambda", 2], 5])) pars[idx["lambda", 1]:idx["lambda", 2], 5] = 0
		if(is.na(pars[idx["lambda", 1]:idx["lambda", 2], 6])) pars[idx["lambda", 1]:idx["lambda", 2], 6] = 100
		if(is.null(start.pars$lambda)) pars[idx["lambda", 1]:idx["lambda", 2], 1] = 0 else pars[idx["lambda", 1]:idx["lambda", 2], 1] = start.pars$lambda[1]
		if(any(substr(fixed.names, 1, 6) == "lambda")){
			pars[idx["lambda", 1]:idx["lambda", 2], 1] = as.numeric(fixed.pars$lambda)
			pars[idx["lambda", 1]:idx["lambda", 2], 5] = as.numeric(fixed.pars$lambda)
			pars[idx["lambda", 1]:idx["lambda", 2], 6] = as.numeric(fixed.pars$lambda)
		}
	}
	if(modelinc[15]>0){
		vxnames = paste("vxreg",1:modelinc[15],sep="")
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 5]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 6]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 6] = 100
		pars[idx["vxreg", 1]:idx["vxreg", 2], 1] = rep(TinY, modelinc[15])
		sp = na.omit(match(start.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[vxnames[sp[i]], 1] = as.numeric(start.pars[vxnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[vxnames[sp[i]], 1] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 5] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 6] = as.numeric(fixed.pars[vxnames[sp[i]]])
			}
		}
	}
	pars = .distributionstart(pars, model, start.pars, fixed.pars)
	
	return( pars )
}

# Hentchel Family Model Initialization:
# model specification
# modified Hentschel model to include indicator power for switching
.fgarchModel = function(model)
{
	ans = switch(model,
			GARCH    = list(parameters = 
							list(lambda = 2, delta = 2, eta1 = 0, eta2 = 0, fk = 0),
							indicator = c(0, 0, 0, 0, 0), garchtype = 1),
			TGARCH   = list(parameters = 
							list(lambda = 1, delta = 1, eta1 = 0.05, eta2 = 0, fk = 0),
							indicator = c(0, 0, 1, 0, 0), garchtype = 2),
			AVGARCH  = list(parameters = 
							list(lambda = 1, delta = 1, eta1 = 0.02, eta2 = 0.05, fk = 0), 
							indicator = c(0, 0, 1, 1, 0), garchtype = 3),
			NGARCH   = list(parameters = 
							list(lambda = 2, delta = 0, eta1 = 0, eta2 = 0, fk = 1), 
							indicator = c(1, 0, 0, 0, 0), garchtype = 4),
			NAGARCH  = list(parameters = 
							list(lambda = 2, delta = 2, eta1 = 0, eta2 = 0.05, fk = 0), 
							indicator = c(0, 0, 0, 1, 0), garchtype = 5),
			APARCH   = list(parameters = 
							list(lambda = 1, delta = 0, eta1 = 0.05, eta2 = 0, fk = 1), 
							indicator = c(1, 0, 1, 0, 0), garchtype = 6),
			ALLGARCH = list(parameters = 
							list(lambda = 2, delta = 0, eta1 = 0.05, eta2 = 0.05, fk = 1), 
							indicator = c(1, 0, 1, 1, 0), garchtype = 7),
			GJRGARCH = list(parameters = 
							list(lambda = 2, delta = 2, eta1 = 0.05, eta2 = 0, fk = 0), 
							indicator = c(0, 0, 1, 0, 0), garchtype = 8),
			EGARCH   = list(parameters = 
							list(lambda = eps, delta = 1, eta1 = 0.05, eta2 = 0, fk = 0), 
							indicator = c(0, 0, 1, 0, 0), garchtype = 0))
	return(ans)
}

.fmodelBounds = function(model)
{
	# lambda, delta, eta1, eta2, fk
	ans = switch(model,
			GARCH    = list( LB = c(    2, 2,       0,   0,  0 ), UB = c(   2, 2,      0,  0, 0 ) ),
			TGARCH   = list( LB = c(    1, 1, -1+TinY,   0,  0 ), UB = c(   1, 1, 1-TinY,  0, 0 ) ),
			AVGARCH  = list( LB = c(    1, 1, -1+TinY, -10,  0 ), UB = c(   1, 1, 1-TinY, 10, 0 ) ),
			NGARCH   = list( LB = c( 0.01, 0,       0,   0,  1 ), UB = c(   4, 0,      0,  0, 1 ) ),
			NAGARCH  = list( LB = c(    2, 2,       0, -10,  0 ), UB = c(   2, 2,      0, 10, 0 ) ),
			APARCH   = list( LB = c( 0.01, 0, -1+TinY,   0,  1 ), UB = c(   4, 0, 1-TinY,  0, 1 ) ),
			ALLGARCH = list( LB = c( 0.01, 0, -1+TinY, -10,  1 ), UB = c(   4, 0, 1-TinY, 10, 1 ) ),
			GJRGARCH = list( LB = c(    2, 2, -1+TinY,   0,  0 ), UB = c(   2, 2, 1-TinY,  0, 0 ) ),
			EGARCH   = list( LB = c( 0.01, 1,      -1,   0,  0 ), UB = c( eps, 1,      1,  0, 0 ) ))
	return(ans)
}


# fGARCH model start parameters
.fgarchstart = function(pars, arglist)
{
	data = arglist$data
	model = arglist$model
	dscale = arglist$dscale
	modelinc = model$modelinc
	start.pars = model$start.pars
	fixed.pars = model$fixed.pars
	idx = model$pidx
	fixed.names = names(fixed.pars)
	start.names = names(start.pars)
	pars = .meqstart(pars, arglist)
	fmodel = model$fmodel
	
	if(modelinc[7]>0){
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 5])) pars[idx["omega", 1]:idx["omega", 2], 5] = var(data)/100000
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 6])) pars[idx["omega", 1]:idx["omega", 2], 6] = var(data)*100000
		if(is.null(start.pars$omega)) pars[idx["omega", 1]:idx["omega", 2], 1] = (var(data, na.rm = TRUE))/1000 else 
			pars[idx["omega", 1]:idx["omega", 2], 1] = start.pars$omega[1]/dscale
		if(any(substr(fixed.names, 1, 5) == "omega")){
			pars[idx["omega", 1]:idx["omega", 2], 1] = as.numeric(fixed.pars$omega)
			pars[idx["omega", 1]:idx["omega", 2], 5] = fixed.pars$omega
			pars[idx["omega", 1]:idx["omega", 2], 6] = fixed.pars$omega
		}
	}
	
	if(modelinc[8]>0){
		gpnames = paste("alpha",1:modelinc[8],sep="")
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 5]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 6]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 6] =  1-TinY
		pars[idx["alpha", 1]:idx["alpha", 2], 1] = rep(0.05/modelinc[8], modelinc[8])
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	# Upper and Lower Custom Bounds on some of the fGARCH model parameters is
	# not allowed
	if(modelinc[11] > 0){
		gqnames = paste("eta1",1:modelinc[11],sep="")
		pars[idx["eta1", 1]:idx["eta1", 2], 5] = fmodel$fbounds$LB[3]
		pars[idx["eta1", 1]:idx["eta1", 2], 6] = fmodel$fbounds$UB[3]
		pars[idx["eta1", 1]:idx["eta1", 2], 1] = rep(fmodel$fpars$eta1/modelinc[11], modelinc[11])
		sp = na.omit(match(start.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gqnames[sp[i]], 1] = as.numeric(start.pars[gqnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gqnames[sp[i]], 1] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 5] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 6] = as.numeric(fixed.pars[gqnames[sp[i]]])
			}
		}
	}
	
	if(modelinc[12] > 0){
		gqnames = paste("eta2",1:modelinc[12],sep="")
		pars[idx["eta2", 1]:idx["eta2", 2], 5] = fmodel$fbounds$LB[4]
		pars[idx["eta2", 1]:idx["eta2", 2], 6] = fmodel$fbounds$UB[4]
		pars[idx["eta2", 1]:idx["eta2", 2], 1] = rep(fmodel$fpars$eta2/modelinc[12], modelinc[12])
		sp = na.omit(match(start.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gqnames[sp[i]], 1] = as.numeric(start.pars[gqnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gqnames[sp[i]], 1] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 5] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 6] = as.numeric(fixed.pars[gqnames[sp[i]]])
			}
		}
	}
	if(modelinc[9] > 0){
		gqnames = paste("beta",1:modelinc[9],sep="")
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 5]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 6]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 6] =  1-TinY
		pars[idx["beta", 1]:idx["beta", 2], 1] = rep(0.9/modelinc[9], modelinc[9])
		sp = na.omit(match(start.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gqnames[sp[i]], 1] = as.numeric(start.pars[gqnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gqnames[sp[i]], 1] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 5] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 6] = as.numeric(fixed.pars[gqnames[sp[i]]])
			}
		}
	}
	
	if(modelinc[14]>0){
		pars[idx["lambda", 1]:idx["lambda", 2], 5] = fmodel$fbounds$LB[1]
		pars[idx["lambda", 1]:idx["lambda", 2], 6] = fmodel$fbounds$UB[1]
		if(is.null(start.pars$lambda)) pars[idx["lambda", 1]:idx["lambda", 2], 1] = fmodel$fpars$lambda else pars[idx["lambda", 1]:idx["lambda", 2], 1] = start.pars$lambda[1]
		if(any(substr(fixed.names, 1, 6) == "lambda")){
			pars[idx["lambda", 1]:idx["lambda", 2], 1] = as.numeric(fixed.pars$lambda)
			pars[idx["lambda", 1]:idx["lambda", 2], 5] = as.numeric(fixed.pars$lambda)
			pars[idx["lambda", 1]:idx["lambda", 2], 6] = as.numeric(fixed.pars$lambda)
		}
	}
		
	
	if(modelinc[13]>0){
		pars[idx["delta", 1]:idx["delta", 2], 5] = fmodel$fbounds$LB[2]
		pars[idx["delta", 1]:idx["delta", 2], 6] = fmodel$fbounds$UB[2]
		if(is.null(start.pars$delta)) pars[idx["delta", 1]:idx["delta", 2], 1] = fmodel$fpars$delta else pars[idx["delta", 1]:idx["delta", 2], 1] = start.pars$delta[1]
		if(any(substr(fixed.names, 1, 5) == "delta")){
			pars[idx["delta", 1]:idx["delta", 2], 1] = as.numeric(fixed.pars$delta)
			pars[idx["delta", 1]:idx["delta", 2], 5] = as.numeric(fixed.pars$delta)
			pars[idx["delta", 1]:idx["delta", 2], 6] = as.numeric(fixed.pars$delta)
		}
	}
	
	if(modelinc[15]>0){
		vxnames = paste("vxreg",1:modelinc[15],sep="")
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 5]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 6]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 6] = 100
		pars[idx["vxreg", 1]:idx["vxreg", 2], 1] = rep(TinY, modelinc[15])
		sp = na.omit(match(start.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[vxnames[sp[i]], 1] = as.numeric(start.pars[vxnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[vxnames[sp[i]], 1] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 5] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 6] = as.numeric(fixed.pars[vxnames[sp[i]]])
			}
		}
	}
	pars = .distributionstart(pars, model, start.pars, fixed.pars)
	
	return( pars )
}

.egarchstart = function(pars, arglist)
{
	data = arglist$data
	model = arglist$model
	dscale = arglist$dscale
	modelinc = model$modelinc
	start.pars = model$start.pars
	fixed.pars = model$fixed.pars
	idx = model$pidx
	fixed.names = names(fixed.pars)
	start.names = names(start.pars)
	pars = .meqstart(pars, arglist)
	if(modelinc[7]>0){
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 5])) pars[idx["omega", 1]:idx["omega", 2], 5] = -10
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 6])) pars[idx["omega", 1]:idx["omega", 2], 6] = 10
		if(is.null(start.pars$omega)) pars[idx["omega", 1]:idx["omega", 2], 1] = (var(data, na.rm = TRUE)) else pars[idx["omega", 1]:idx["omega", 2], 1] = start.pars$omega[1]/dscale
		if(any(substr(fixed.names, 1, 5) == "omega")){
			pars[idx["omega", 1]:idx["omega", 2], 1] = as.numeric(fixed.pars$omega)
			pars[idx["omega", 1]:idx["omega", 2], 5] = fixed.pars$omega
			pars[idx["omega", 1]:idx["omega", 2], 6] = fixed.pars$omega
		}
	}
	if(modelinc[8]>0){
		gpnames = paste("alpha",1:modelinc[8],sep="")
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 5]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 5] = -10
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 6]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 6] =  10
		pars[idx["alpha", 1]:idx["alpha", 2], 1] = rep(0.05/modelinc[8], modelinc[8])
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	
	if(modelinc[10] > 0){
		gqnames = paste("gamma",1:modelinc[10],sep="")
		pxd = which(is.na(pars[idx["gamma", 1]:idx["gamma", 2], 5]))
		if(length(pxd)>0) pars[(idx["gamma", 1]:idx["gamma", 2])[pxd], 5] = -10
		pxd = which(is.na(pars[idx["gamma", 1]:idx["gamma", 2], 6]))
		if(length(pxd)>0) pars[(idx["gamma", 1]:idx["gamma", 2])[pxd], 6] =  10
		pars[idx["gamma", 1]:idx["gamma", 2], 1] = rep(0.1/modelinc[10],modelinc[10])
		sp = na.omit(match(start.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gqnames[sp[i]], 1] = as.numeric(start.pars[gqnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gqnames[sp[i]], 1] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 5] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 6] = as.numeric(fixed.pars[gqnames[sp[i]]])
			}
		}
	}
	
	if(modelinc[9] > 0){
		gqnames = paste("beta",1:modelinc[9],sep="")
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 5]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 5] = -1+TinY
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 6]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 6] =  1-TinY
		pars[idx["beta", 1]:idx["beta", 2], 1] = rep(0.9/modelinc[9], modelinc[9])
		sp = na.omit(match(start.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gqnames[sp[i]], 1] = as.numeric(start.pars[gqnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gqnames[sp[i]], 1] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 5] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 6] = as.numeric(fixed.pars[gqnames[sp[i]]])
			}
		}
	}
	if(modelinc[15]>0){
		vxnames = paste("vxreg",1:modelinc[15],sep="")
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 5]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 5] = -100
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 6]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 6] = 100		
		pars[idx["vxreg", 1]:idx["vxreg", 2], 1] = rep(TinY, modelinc[15])
		sp = na.omit(match(start.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[vxnames[sp[i]], 1] = as.numeric(start.pars[vxnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[vxnames[sp[i]], 1] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 5] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 6] = as.numeric(fixed.pars[vxnames[sp[i]]])
			}
		}
	}
	pars = .distributionstart(pars, model, start.pars, fixed.pars)
	
	return( pars )
}

# apARCH model start parameters
.aparchstart = function(pars, arglist)
{
	data = arglist$data
	model = arglist$model
	dscale = arglist$dscale
	modelinc = model$modelinc
	start.pars = model$start.pars
	fixed.pars = model$fixed.pars
	idx = model$pidx
	fixed.names = names(fixed.pars)
	start.names = names(start.pars)
	pars = .meqstart(pars, arglist)
	if(modelinc[7]>0){
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 5])) pars[idx["omega", 1]:idx["omega", 2], 5] = eps
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 6])) pars[idx["omega", 1]:idx["omega", 2], 6] = var(data)*1000
		if(is.null(start.pars$omega)) pars[idx["omega", 1]:idx["omega", 2], 1] = (var(data, na.rm = TRUE))/1000 else 
			pars[idx["omega", 1]:idx["omega", 2], 1] = start.pars$omega[1]/dscale
		if(any(substr(fixed.names, 1, 5) == "omega")){
			pars[idx["omega", 1]:idx["omega", 2], 1] = as.numeric(fixed.pars$omega)
			pars[idx["omega", 1]:idx["omega", 2], 5] = fixed.pars$omega
			pars[idx["omega", 1]:idx["omega", 2], 6] = fixed.pars$omega
		}
	}
	if(modelinc[8]>0){
		gpnames = paste("alpha",1:modelinc[8],sep="")
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 5]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 6]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 6] =  1-TinY
		pars[idx["alpha", 1]:idx["alpha", 2], 1] = rep(0.05/modelinc[8], modelinc[8])
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	
	if(modelinc[10] > 0){
		gqnames = paste("gamma",1:modelinc[10],sep="")
		pxd = which(is.na(pars[idx["gamma", 1]:idx["gamma", 2], 5]))
		if(length(pxd)>0) pars[(idx["gamma", 1]:idx["gamma", 2])[pxd], 5] = -1+TinY
		pxd = which(is.na(pars[idx["gamma", 1]:idx["gamma", 2], 6]))
		if(length(pxd)>0) pars[(idx["gamma", 1]:idx["gamma", 2])[pxd], 6] =  1-TinY
		pars[idx["gamma", 1]:idx["gamma", 2], 1] = rep(0.05/modelinc[10],modelinc[10])
		sp = na.omit(match(start.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gqnames[sp[i]], 1] = as.numeric(start.pars[gqnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gqnames[sp[i]], 1] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 5] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 6] = as.numeric(fixed.pars[gqnames[sp[i]]])
			}
		}
	}
	
	if(modelinc[13]>0){
		if(is.na(pars[idx["delta", 1]:idx["delta", 2], 5])) pars[idx["delta", 1]:idx["delta", 2], 5] = 0.01
		if(is.na(pars[idx["delta", 1]:idx["delta", 2], 6])) pars[idx["delta", 1]:idx["delta", 2], 6] = 3.5
		if(is.null(start.pars$delta)) pars[idx["delta", 1]:idx["delta", 2], 1] = 2 else pars[idx["delta", 1]:idx["delta", 2], 1] = start.pars$delta[1]
		if(any(substr(fixed.names, 1, 5) == "delta")){
			pars[idx["delta", 1]:idx["delta", 2], 1] = as.numeric(fixed.pars$delta)
			pars[idx["delta", 1]:idx["delta", 2], 5] = as.numeric(fixed.pars$delta)
			pars[idx["delta", 1]:idx["delta", 2], 6] = as.numeric(fixed.pars$delta)
		}
	}
	
	if(modelinc[9] > 0){
		gqnames = paste("beta",1:modelinc[9],sep="")
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 5]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 6]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 6] =  1-TinY
		pars[idx["beta", 1]:idx["beta", 2], 1] = rep(0.9/modelinc[9], modelinc[9])
		sp = na.omit(match(start.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gqnames[sp[i]], 1] = as.numeric(start.pars[gqnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gqnames[sp[i]], 1] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 5] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 6] = as.numeric(fixed.pars[gqnames[sp[i]]])
			}
		}
	}
	if(modelinc[15]>0){
		vxnames = paste("vxreg",1:modelinc[15],sep="")
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 5]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 6]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 6] = 100
		pars[idx["vxreg", 1]:idx["vxreg", 2], 1] = rep(TinY, modelinc[15])
		sp = na.omit(match(start.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[vxnames[sp[i]], 1] = as.numeric(start.pars[vxnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[vxnames[sp[i]], 1] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 5] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 6] = as.numeric(fixed.pars[vxnames[sp[i]]])
			}
		}
	}
	pars = .distributionstart(pars, model, start.pars, fixed.pars)
	
	return( pars )
}

# multiplicative component apARCH model start parameters
.mcaparchstart = function(pars, arglist)
{
	data = arglist$data
	model = arglist$model
	dscale = arglist$dscale
	modelinc = model$modelinc
	start.pars = model$start.pars
	fixed.pars = model$fixed.pars
	idx = model$pidx
	fixed.names = names(fixed.pars)
	start.names = names(start.pars)
	pars = .meqstart(pars, arglist)
	if(modelinc[7]>0){
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 5])) pars[idx["omega", 1]:idx["omega", 2], 5] = eps
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 6])) pars[idx["omega", 1]:idx["omega", 2], 6] = 5
		if(is.null(start.pars$omega)) pars[idx["omega", 1]:idx["omega", 2], 1] = 0.05 else 
			pars[idx["omega", 1]:idx["omega", 2], 1] = start.pars$omega[1]/dscale
		if(any(substr(fixed.names, 1, 5) == "omega")){
			pars[idx["omega", 1]:idx["omega", 2], 1] = as.numeric(fixed.pars$omega)
			pars[idx["omega", 1]:idx["omega", 2], 5] = fixed.pars$omega
			pars[idx["omega", 1]:idx["omega", 2], 6] = fixed.pars$omega
		}
	}
	if(modelinc[8]>0){
		gpnames = paste("alpha",1:modelinc[8],sep="")
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 5]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 6]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 6] =  1-TinY
		pars[idx["alpha", 1]:idx["alpha", 2], 1] = rep(0.05/modelinc[8], modelinc[8])
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
	}
	
	if(modelinc[10] > 0){
		gqnames = paste("gamma",1:modelinc[10],sep="")
		pxd = which(is.na(pars[idx["gamma", 1]:idx["gamma", 2], 5]))
		if(length(pxd)>0) pars[(idx["gamma", 1]:idx["gamma", 2])[pxd], 5] = -1+TinY
		pxd = which(is.na(pars[idx["gamma", 1]:idx["gamma", 2], 6]))
		if(length(pxd)>0) pars[(idx["gamma", 1]:idx["gamma", 2])[pxd], 6] =  1-TinY
		pars[idx["gamma", 1]:idx["gamma", 2], 1] = rep(0.05/modelinc[10],modelinc[10])
		sp = na.omit(match(start.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gqnames[sp[i]], 1] = as.numeric(start.pars[gqnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gqnames[sp[i]], 1] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 5] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 6] = as.numeric(fixed.pars[gqnames[sp[i]]])
			}
		}
	}
	
	if(modelinc[13]>0){
		if(is.na(pars[idx["delta", 1]:idx["delta", 2], 5])) pars[idx["delta", 1]:idx["delta", 2], 5] = 0.01
		if(is.na(pars[idx["delta", 1]:idx["delta", 2], 6])) pars[idx["delta", 1]:idx["delta", 2], 6] = 3.5
		if(is.null(start.pars$delta)) pars[idx["delta", 1]:idx["delta", 2], 1] = 2 else pars[idx["delta", 1]:idx["delta", 2], 1] = start.pars$delta[1]
		if(any(substr(fixed.names, 1, 5) == "delta")){
			pars[idx["delta", 1]:idx["delta", 2], 1] = as.numeric(fixed.pars$delta)
			pars[idx["delta", 1]:idx["delta", 2], 5] = as.numeric(fixed.pars$delta)
			pars[idx["delta", 1]:idx["delta", 2], 6] = as.numeric(fixed.pars$delta)
		}
	}
	
	if(modelinc[9] > 0){
		gqnames = paste("beta",1:modelinc[9],sep="")
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 5]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["beta", 1]:idx["beta", 2], 6]))
		if(length(pxd)>0) pars[(idx["beta", 1]:idx["beta", 2])[pxd], 6] =  1-TinY
		pars[idx["beta", 1]:idx["beta", 2], 1] = rep(0.9/modelinc[9], modelinc[9])
		sp = na.omit(match(start.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gqnames[sp[i]], 1] = as.numeric(start.pars[gqnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gqnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gqnames[sp[i]], 1] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 5] = as.numeric(fixed.pars[gqnames[sp[i]]])
				pars[gqnames[sp[i]], 6] = as.numeric(fixed.pars[gqnames[sp[i]]])
			}
		}
	}
	if(modelinc[15]>0){
		vxnames = paste("vxreg",1:modelinc[15],sep="")
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 5]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 6]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 6] = 100
		pars[idx["vxreg", 1]:idx["vxreg", 2], 1] = rep(TinY, modelinc[15])
		sp = na.omit(match(start.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[vxnames[sp[i]], 1] = as.numeric(start.pars[vxnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[vxnames[sp[i]], 1] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 5] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 6] = as.numeric(fixed.pars[vxnames[sp[i]]])
			}
		}
	}
	pars = .distributionstart(pars, model, start.pars, fixed.pars)
	
	return( pars )
}

# iGARCH model start parameters
.igarchstart = function(pars, arglist)
{
	data = arglist$data
	model = arglist$model
	dscale = arglist$dscale
	modelinc = model$modelinc
	start.pars = model$start.pars
	fixed.pars = model$fixed.pars
	idx = model$pidx
	fixed.names = names(fixed.pars)
	start.names = names(start.pars)
	pars = .meqstart(pars, arglist)
	if(modelinc[7]>0){
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 5])) pars[idx["omega", 1]:idx["omega", 2], 5] = var(data)/100000
		if(is.na(pars[idx["omega", 1]:idx["omega", 2], 6])) pars[idx["omega", 1]:idx["omega", 2], 6] = var(data)*100000
		if(is.null(start.pars$omega)) pars[idx["omega", 1]:idx["omega", 2], 1] = (var(data, na.rm = TRUE))/1000 else 
			pars[idx["omega", 1]:idx["omega", 2], 1] = start.pars$omega[1]/dscale
		if(any(substr(fixed.names, 1, 5) == "omega")){
			pars[idx["omega", 1]:idx["omega", 2], 1] = as.numeric(fixed.pars$omega)
			pars[idx["omega", 1]:idx["omega", 2], 5] = fixed.pars$omega
			pars[idx["omega", 1]:idx["omega", 2], 6] = fixed.pars$omega
		}
	}
	# we add the sum of alpha in order to calculate the iGARCH beta
	sumalpha = 0
	if(modelinc[8]>0){
		gpnames = paste("alpha",1:modelinc[8],sep="")
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 5]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["alpha", 1]:idx["alpha", 2], 6]))
		if(length(pxd)>0) pars[(idx["alpha", 1]:idx["alpha", 2])[pxd], 6] =  1-TinY		
		pars[idx["alpha", 1]:idx["alpha", 2], 1] = rep(0.05/modelinc[8], modelinc[8])
		sp = na.omit(match(start.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[gpnames[sp[i]], 1] = as.numeric(start.pars[gpnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, gpnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[gpnames[sp[i]], 1] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 5] = as.numeric(fixed.pars[gpnames[sp[i]]])
				pars[gpnames[sp[i]], 6] = as.numeric(fixed.pars[gpnames[sp[i]]])
			}
		}
		sumalpha = sum(pars[idx["alpha", 1]:idx["alpha", 2], 1])
	}
	# Custom bounds for beta on iGARCH not allowed
	if(modelinc[9] > 0){
		if(modelinc[9] == 1){
			pars[idx["beta", 1]:idx["beta", 2], 5] = 0
			pars[idx["beta", 1]:idx["beta", 2], 6] = 1-TinY
			pars[idx["beta", 1]:idx["beta", 2], 1] = 1 - sumalpha
			pars[idx["beta", 1]:idx["beta", 2], 2] = 1
			pars[idx["beta", 1]:idx["beta", 2], 4] = 0
		} else{
			# this now excludes the last beta
			gqnames = paste("beta",1:(modelinc[9]-1),sep="")
			pars[idx["beta", 1]:(idx["beta", 2]-1), 5] = 0
			pars[idx["beta", 1]:(idx["beta", 2]-1), 6] = 1-TinY
			pars[idx["beta", 1]:(idx["beta", 2]-1), 1] = rep(0.6/modelinc[9], modelinc[9]-1)
			sumbeta = sum(rep(0.6/modelinc[9], modelinc[9]-1))
			pars[idx["beta", 2], 5] = 0
			pars[idx["beta", 2], 6] = 1-TinY
			pars[idx["beta", 2], 1] = 1 - sumalpha - sumbeta
			pars[idx["beta", 2], 2] = 1
			pars[idx["beta", 2], 4] = 0
			sp = na.omit(match(start.names, gqnames))
			if(length(sp)>0){
				for(i in 1:length(sp)) pars[gqnames[sp[i]], 1] = as.numeric(start.pars[gqnames[sp[i]]])
			}
			sp = na.omit(match(fixed.names, gqnames))
			if(length(sp)>0){
				for(i in 1:length(sp)){
					pars[gqnames[sp[i]], 1] = as.numeric(fixed.pars[gqnames[sp[i]]])
					pars[gqnames[sp[i]], 5] = as.numeric(fixed.pars[gqnames[sp[i]]])
					pars[gqnames[sp[i]], 6] = as.numeric(fixed.pars[gqnames[sp[i]]])
				}
			}
		}
	}
	
	if(modelinc[15]>0){
		vxnames = paste("vxreg",1:modelinc[15],sep="")
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 5]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 5] = 0
		pxd = which(is.na(pars[idx["vxreg", 1]:idx["vxreg", 2], 6]))
		if(length(pxd)>0) pars[(idx["vxreg", 1]:idx["vxreg", 2])[pxd], 6] = 100
		pars[idx["vxreg", 1]:idx["vxreg", 2], 1] = rep(TinY, modelinc[15])
		sp = na.omit(match(start.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)) pars[vxnames[sp[i]], 1] = as.numeric(start.pars[vxnames[sp[i]]])
		}
		sp = na.omit(match(fixed.names, vxnames))
		if(length(sp)>0){
			for(i in 1:length(sp)){
				pars[vxnames[sp[i]], 1] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 5] = as.numeric(fixed.pars[vxnames[sp[i]]])
				pars[vxnames[sp[i]], 6] = as.numeric(fixed.pars[vxnames[sp[i]]])
			}
		}
	}
	pars = .distributionstart(pars, model, start.pars, fixed.pars)
	
	return( pars )
}

################################################################################
# modified rsFit/absvalFit method from the fArma package:
.rsfit = function(x, levels = 50, minnpts = 3, cut.off = 10^c(0.7, 2.5))
{  
	# A functions implemented by Diethelm Wuertz
	
	# Description:
	#   R/S Statistic Method - [Taqqu 3.6]
	
	# Arguments:
	#   x - a numeric vector, a 'timeSeries' object, or any other
	#       object which can be transformed into a vector by the
	#       function 'as.vector'.
	#   levels - the number of aggregation levels
	#   minnpts - the minimum block size
	#   cut.off - the lower and upper cut off for the fit
	
	# FUNCTION:
	
	# Settings:
	call = match.call()
	data = list(x = x)
	x = as.vector(x)
	n = length(x)
	increment = (log10(n/minnpts))/levels
	M = floor(10^((1:levels)*increment))
	M = M[M > 1]
	
	# R/S Method:
	Y = cumsum(x)
	Y2 = cumsum(x*x)
	RS = NULL
	for (m in M) {
		S = sqrt(Y2[m]/m - (Y[m]/m)^2)
		Z = Y[1:m]-(1:m)*Y[m]/m
		STATS = (max(Z) - min(Z))/S
		RS = c(RS, STATS)
	}
	# Fit:
	wt = trunc((sign((M-cut.off[1])*(cut.off[2]-M))+1)/2)
	fit = lsfit(log10(M), log10(RS), wt)
	fitH = fit
	fitH$wt = NULL
	diag = as.data.frame(ls.print(fitH, print.it = FALSE)[[2]][[1]])
	beta = fit$coef[[2]]
	return(beta)
}

.absHmom = function(x, levels = 50, minnpts = 3, cut.off = 10^c(0.7, 2.5), moment = 1){
	x = as.vector(x)
	n = length(x)
	increment = (log10(n/minnpts))/levels
	M = floor(10^((1:levels) * increment))
	ABSVAL = NULL
	for(m in M) {
		nCols = n%/%m
		X = matrix(x[1:(m * nCols)], byrow = FALSE, ncol = nCols)
		Y = colMeans(X)
		MEAN = mean(Y)
		STATS = sum((abs(Y - MEAN))^moment)/(length(Y) - 1)
		ABSVAL = c(ABSVAL, STATS)
	}
	wt = trunc((sign((M - cut.off[1]) * (cut.off[2] - M)) + 1)/2)
	fit = lsfit(log10(M), log10(ABSVAL), wt)
	fitH = lsfit(log10(M), log10(ABSVAL * M^moment)/moment, wt)
	fitH$wt = NULL
	diag = as.data.frame(ls.print(fitH, print.it = FALSE)[[2]][[1]])
	beta = fit$coef[[2]]
	H = beta/moment + 1
	return(H)
}
################################################################################

.distributionstart = function(pars, model, start.pars, fixed.pars){
	modelinc = model$modelinc
	idx = model$pidx
	fixed.names = names(fixed.pars)
	start.names = names(start.pars)
	dbounds = .DistributionBounds(distribution = model$modeldesc$distribution)
	if(modelinc[16]>0){		
		if(is.na(pars[idx["skew", 1]:idx["skew", 2], 5])) pars[idx["skew", 1]:idx["skew", 2], 5] = dbounds$skew.LB
		if(is.na(pars[idx["skew", 1]:idx["skew", 2], 6])) pars[idx["skew", 1]:idx["skew", 2], 6] = dbounds$skew.UB
		if(is.null(start.pars$skew)) pars[idx["skew", 1]:idx["skew", 2], 1] = dbounds$skew else pars[idx["skew", 1]:idx["skew", 2], 1] = start.pars$skew[1]
		if(any(substr(fixed.names, 1, 4) == "skew")){
			pars[idx["skew", 1]:idx["skew", 2], 1] = as.numeric(fixed.pars$skew)
			pars[idx["skew", 1]:idx["skew", 2], 5] = as.numeric(fixed.pars$skew)
			pars[idx["skew", 1]:idx["skew", 2], 6] = as.numeric(fixed.pars$skew)
		}
	}
	if(modelinc[17]>0){
		if(is.na(pars[idx["shape", 1]:idx["shape", 2], 5])) pars[idx["shape", 1]:idx["shape", 2], 5] = dbounds$shape.LB
		if(is.na(pars[idx["shape", 1]:idx["shape", 2], 6])) pars[idx["shape", 1]:idx["shape", 2], 6] = dbounds$shape.UB
		if(is.null(start.pars$shape)) pars[idx["shape", 1]:idx["shape", 2], 1] = dbounds$shape else pars[idx["shape", 1]:idx["shape", 2], 1] = start.pars$shape[1]
		if(any(substr(fixed.names, 1, 5) == "shape")){
			pars[idx["shape", 1]:idx["shape", 2], 1] = as.numeric(fixed.pars$shape)
			pars[idx["shape", 1]:idx["shape", 2], 5] = as.numeric(fixed.pars$shape)
			pars[idx["shape", 1]:idx["shape", 2], 6] = as.numeric(fixed.pars$shape)
		}
	}
	if(modelinc[18]>0){
		if(is.na(pars[idx["ghlambda", 1]:idx["ghlambda", 2], 5])) pars[idx["ghlambda", 1]:idx["ghlambda", 2], 5] = dbounds$ghlambda.LB
		if(is.na(pars[idx["ghlambda", 1]:idx["ghlambda", 2], 6])) pars[idx["ghlambda", 1]:idx["ghlambda", 2], 6] = dbounds$ghlambda.UB
		if(is.null(start.pars$ghlambda)) pars[idx["ghlambda", 1]:idx["ghlambda", 2], 1] = dbounds$ghlambda else pars[idx["ghlambda", 1]:idx["ghlambda", 2], 1] = start.pars$ghlambda[1]
		if(any(substr(fixed.names, 1, 8) == "ghlambda")){
			pars[idx["ghlambda", 1]:idx["ghlambda", 2], 1] = as.numeric(fixed.pars$ghlambda)
			pars[idx["ghlambda", 1]:idx["ghlambda", 2], 5] = as.numeric(fixed.pars$ghlambda)
			pars[idx["ghlambda", 1]:idx["ghlambda", 2], 6] = as.numeric(fixed.pars$ghlambda)
		}
	}
	return(pars)
}