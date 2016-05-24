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


#####################################################################################
# Transformation Functions
#------------------------------------------------------------------------------------

.pparametric = function(mfit, zres)
{
	m = dim(zres)[2]
	n = dim(zres)[1] 
	ures = matrix(NA, ncol = m, nrow = n)
	
	for(i in 1:m){
		gdist = mfit@fit[[i]]@model$modeldesc$distribution
		lambda = ifelse(mfit@fit[[i]]@model$modelinc[18]>0, mfit@fit[[i]]@fit$ipars["ghlambda",1], 0)
		skew = ifelse(mfit@fit[[i]]@model$modelinc[16]>0, mfit@fit[[i]]@fit$ipars["skew",1],  0)
		shape = ifelse(mfit@fit[[i]]@model$modelinc[17]>0, mfit@fit[[i]]@fit$ipars["shape",1],  0)
		ures[,i] = pdist(gdist, zres[,i], mu = 0, sigma = 1, lambda = lambda, skew = skew, shape = shape)
	}
	return(ures)
}

.pparametric.filter = function(mflt, zres)
{
	m = dim(zres)[2]
	n = dim(zres)[1] 
	ures = matrix(NA, ncol = m, nrow = n)
	
	for(i in 1:m){
		gdist = mflt@filter[[i]]@model$modeldesc$distribution
		lambda = ifelse(mflt@filter[[i]]@model$modelinc[18]>0, mflt@filter[[i]]@model$pars["ghlambda",1], 0)
		skew = ifelse(mflt@filter[[i]]@model$modelinc[16]>0, mflt@filter[[i]]@model$pars["skew",1],  0)
		shape = ifelse(mflt@filter[[i]]@model$modelinc[17]>0, mflt@filter[[i]]@model$pars["shape",1],  0)		
		ures[,i] = pdist(gdist, zres[,i], mu = 0, sigma = 1, lambda = lambda, skew = skew, shape = shape)
	}
	return(ures)
}

.pempirical = function(zres)
{
	m = dim(zres)[2]
	n = dim(zres)[1] 
	ures = matrix(NA, ncol = m, nrow = n)
	for(i in 1:m){
		fn = ecdf(sort(zres[,i]))
		ures[,i] = fn(zres[,i])
	}
	return(ures)
}


.pempirical.filter = function(zres, dcc.old)
{
	m = dim(zres)[2]
	n = dim(zres)[1] 
	ures = matrix(NA, ncol = m, nrow = n)
	for(i in 1:m){
		fn = ecdf(sort(zres[1:dcc.old,i]))
		ures[,i] = fn(zres[,i])
	}
	return(ures)
}


.pspd = function(zres, spd.control)
{
	m = dim(zres)[2]
	n = dim(zres)[1] 
	ures = matrix(NA, ncol = m, nrow = n)
	sfit = vector(mode = "list", length = m)
	sfit = lapply(as.list(1:m), function(i) spdfit(zres[,i], upper = spd.control$upper, lower = spd.control$lower, 
						tailfit = "GPD", type = spd.control$type, kernelfit = spd.control$kernel, information = "observed"))
	for(i in 1:m){
		ures[,i] = pspd(zres[,i], sfit[[i]])
	}
	return(list(ures = ures, sfit = sfit))
}

.pspd.filter = function(zres, spd.control, dcc.old)
{
	m = dim(zres)[2]
	n = dim(zres)[1] 
	ures = matrix(NA, ncol = m, nrow = n)
	sfit = vector(mode = "list", length = m)
	sfit = lapply(as.list(1:m), function(i) spdfit(zres[1:dcc.old,i], upper = spd.control$upper, lower = spd.control$lower, 
						tailfit = "GPD", type = spd.control$type, kernelfit = spd.control$kernel, information = "observed"))
	# temporarary workaround for problem with spd
	for(i in 1:m){
			ures[,i] = pspd(zres[,i], sfit[[i]])
	}

	return(list(ures = ures, sfit = sfit))
}


#------------------------------------------------------------------------------------
.qparametric = function(ures, modelinc, pars)
{
	m = dim(ures)[2]
	zres = matrix(NA, ncol = m, nrow = dim(ures)[1])
	distn = c("norm", "snorm", "std", "sstd","ged", "sged", "nig", "ghyp", "jsu")
	for(i in 1:m){
		gdist = distn[modelinc[21,i]]
		lambda = ifelse(modelinc[18,i]>0, pars["ghlambda",i], 0)
		skew =   ifelse(modelinc[16,i]>0, pars["skew",i],     0)
		shape =  ifelse(modelinc[17,i]>0, pars["shape",i],    0)		
		zres[,i] = qdist(gdist, ures[,i], mu = 0, sigma = 1, lambda = lambda, skew = skew, shape = shape)
	}
	return(zres)
}

.qempirical = function(ures, oldz)
{
	zres = matrix(NA, ncol = dim(ures)[2], nrow = dim(ures)[1])	
	for(i in 1:dim(ures)[2]){
		zres[,i] = quantile(oldz[,i], ures[,i], type = 1)
	}
	return(zres)
}

.qspd = function(ures, sfit)
{
	zres = matrix(NA, ncol = dim(ures)[2], nrow = dim(ures)[1])	
	for(i in 1:dim(ures)[2]){
		zres[,i] = qspd(ures[,i], sfit[[i]])
	}
	return(zres)
}
#------------------------------------------------------------------------------------
