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


copula.tvnormalLLH1 = function(pars, arglist)
{
	mgarchenv = arglist$mgarchenv
	model = arglist$model
	modelinc = model$modelinc
	estidx = arglist$estidx
	idx = model$pidx
	ipars = arglist$ipars
	ipars[estidx, 1] = pars
	trace = arglist$trace
	m = arglist$m
	mx = model$maxdccOrder
	fit.control = arglist$fit.control
	sumdcc = sum(ipars[idx["dcca",1]:idx["dcca",2],1])+sum(ipars[idx["dccb",1]:idx["dccb",2],1])
	sumdccg = sum(ipars[idx["dccg",1]:idx["dccg",2],1])
	udata = arglist$ures
	T = dim(udata)[1]
	Z = qnorm(udata)
	Qbar = cov(Z)
	# Take care of the Asymmetry Matrices
	if(modelinc[6]>0){
		Ibar = .asymI(Z)
		aZ = Ibar*Z
		Nbar = cov(aZ)
	} else{
		Ibar = .asymI(Z)
		# zero what is not needed
		aZ = Ibar*Z*0
		Nbar = matrix(0, m, m)
	}
	
	Z = rbind( matrix(0, nrow = mx, ncol = m), Z )
	aZ = rbind( matrix(0, nrow = mx, ncol = m), aZ )
	
	arglist$Nbar = Nbar
	arglist$Qbar = Qbar

	if( fit.control$stationarity ){
		if(modelinc[6]>0){
			persist = .adcccon(pars, arglist)
		} else{
			persist = .dcccon(pars, arglist)
		}
		if( !is.na(persist) && persist >= 1 ) return(llh = get("rmgarch_llh", mgarchenv) + 0.1*(abs(get("rmgarch_llh", mgarchenv))))
	}
	res = .Call("copulaNormalC2", model = as.integer(modelinc), 
			pars = as.numeric(ipars[,1]), idx = as.integer(idx[,1]-1), Qbar = Qbar, 
			Nbar = Nbar, Z = Z, N = aZ, epars = c(sumdcc, sumdccg, mx), 
			PACKAGE = "rmgarch")
	
	llh = res[[3]]
	
	if(is.finite(llh) | !is.na(llh) | !is.nan(llh)){
		assign("rmgarch_llh", llh, envir = mgarchenv)
	} else {
		llh = (get("rmgarch_llh", mgarchenv) + 0.1*(abs(get("rmgarch_llh", mgarchenv))))
	}
		
	ans = switch(arglist$returnType,
			lik = res[[2]][(mx+1):(T+mx)],
			llh = llh,
			ALL = list(lik = res[[2]][(mx+1):(T+mx)], llh = res[[3]], Rt = res[[4]][(mx+1):(T+mx)], Qt = res[[1]][(mx+1):(T+mx)]))
	return(ans)
}

copula.tvnormalLLH2 = function(pars, arglist)
{
	.eps = 	.Machine$double.eps
	model = arglist$model
	umodel = arglist$umodel
	modelinc = model$modelinc
	estidx = arglist$estidx
	eidx = arglist$eidx
	midx = arglist$midx
	mpars = arglist$mpars
	idx = model$pidx	
	ipars = arglist$ipars
	m = arglist$m
	T = arglist$T
	mpars[which(eidx==1, arr.ind = TRUE)] = pars
	# assign the pars to the ipars (used in the DCC diltering)
	ipars[estidx,1] = mpars[which(eidx[,m+1]==1),m+1]
	trace = arglist$trace
	data = arglist$data
	n.start = model$modeldata$n.start
	mx = model$maxdccOrder
	sumdcc = sum(ipars[idx["dcca",1]:idx["dcca",2],1])+sum(ipars[idx["dccb",1]:idx["dccb",2],1])
	sumdccg = sum(ipars[idx["dccg",1]:idx["dccg",2],1])
	garch.llhvec = resids = H = matrix(0, nrow = T, ncol = m)
	
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, umodel$vt)
	specx = vector(mode = "list", length = m)
	for(i in 1:m){
		specx[[i]] = mspec@spec[[i]]
		setfixed(specx[[i]]) = as.list(mpars[which(midx[,i]==1), i])
	}
	flt = multifilter(multifitORspec = multispec(specx), data = xts(data, arglist$index[1:nrow(data)]), 
			out.sample = n.start, realizedVol = arglist$realizedVol[1:nrow(data),])
	garch.llhvec = sapply(flt@filter, FUN = function(x) x@filter$log.likelihoods)
	H = sigma(flt)
	resids = residuals(flt)
	stdres = resids/H
	T = dim(stdres)[1]
	transformation = model$modeldesc$transformation
	spd.control = model$spd.control
	ans = switch(transformation,
			parametric = .pparametric.filter(flt, stdres),
			empirical = .pempirical(stdres),
			spd = .pspd(stdres, spd.control))
	if(transformation == "spd") {
		ures = ans$ures
		sfit = ans$sfit
	} else{
		ures = ans
		sfit = NULL
	}
	# make tail adjustments in order to avoid problem with the quantile functions in
	# optimization
	if(any(ures > 0.99999)){
		xn = which(ures > 0.99999)
		ures[xn] = 0.99999
	}
	if(any(ures < .eps)){
		xn = which(ures < (1.5*.eps))
		ures[xn] = .eps
	}
	
	Z = qnorm(ures)
	Qbar = cov(Z)
	# Take care of the Asymmetry Matrices
	if(modelinc[6]>0){
		Ibar = .asymI(Z)
		aZ = Ibar*Z
		Nbar = cov(aZ)
	} else{
		Ibar = .asymI(Z)
		# zero what is not needed
		aZ = Ibar*Z*0
		Nbar = matrix(0, m, m)
	}
	Z = rbind( matrix(0, nrow = mx, ncol = m), Z )
	aZ = rbind( matrix(0, nrow = mx, ncol = m), aZ )
	
	sumdcc = sum(ipars[idx["dcca",1]:idx["dcca",2],1])+sum(ipars[idx["dccb",1]:idx["dccb",2],1])
	sumdccg = sum(ipars[idx["dccg",1]:idx["dccg",2],1])	

	
	res = .Call("copulaNormalC2", model = as.integer(modelinc), 
			pars = as.numeric(ipars[,1]), idx = as.integer(idx[,1]-1), 
			Qbar = Qbar, Nbar = Nbar, Z = Z, N = aZ, 
			epars = c(sumdcc, sumdccg, mx), PACKAGE = "rmgarch")
	
	copula.llhvec = res[[2]][(mx+1):(T+mx)]
	full.llhvec = apply(cbind(garch.llhvec, -copula.llhvec), 1, "sum")
	#ttmp = sum(-copula.llhvec) - sum(apply(garch.llhvec, 1, "sum"))
	#.ttmp <<-ttmp
	sol = list()
	sol$llhvec = full.llhvec
	sol$llh = sum(full.llhvec)
	
	ans = switch(arglist$returnType,
			lik = sol$llhvec,
			llh = sum(full.llhvec),
			ALL = list(lik = sol$llhvec, llh = sum(full.llhvec), Rt = res[[4]][(mx+1):(T+mx)], Qt = res[[1]][(mx+1):(T+mx)],
					Z = Z, Nbar = Nbar, Qbar = Qbar))
	return(ans)
}


copula.tvnormalLLH3 = function(arglist){
	.eps = 	.Machine$double.eps
	model = arglist$model
	umodel = arglist$umodel
	modelinc = model$modelinc
	estidx = arglist$estidx
	eidx = arglist$eidx
	midx = arglist$midx
	mpars = arglist$mpars
	idx = model$pidx
	ipars = arglist$ipars
	m = arglist$m
	T = arglist$T
	filter.control = arglist$filter.control
	spd.control = arglist$spd.control
	n.old = arglist$n.old
	flt = arglist$filterlist
	trace = arglist$trace
	data = arglist$data
	n.start = model$modeldata$n.start
	mx = model$maxdccOrder
	garch.llhvec = resids = H = matrix(0, nrow = T, ncol = m)
	garch.llhvec = sapply(flt@filter, FUN = function(x) x@filter$log.likelihoods)
	H = sigma(flt)
	resids = residuals(flt)
	stdres = resids/H
	if(is.null(filter.control$n.old)) dcc.old = dim(stdres)[1] else dcc.old = n.old
	
	T = dim(stdres)[1]
	transformation = model$modeldesc$transformation
	ans = switch(transformation,
			parametric = .pparametric.filter(flt, stdres),
			empirical = .pempirical.filter(stdres, dcc.old),
			spd = .pspd.filter(stdres, spd.control, dcc.old))
	if(transformation == "spd") {
		ures = ans$ures
		sfit = ans$sfit
	} else{
		ures = ans
		sfit = NULL
	}
	# make tail adjustments in order to avoid problem with the quantile functions in
	# optimization
	if(any(ures > 0.99999)){
		xn = which(ures > 0.99999)
		ures[xn] = 0.99999
	}
	if(any(ures < .eps)){
		xn = which(ures < (1.5*.eps))
		ures[xn] = .eps
	}
	
	Z = qnorm(ures)
	Qbar = cov(Z[1:dcc.old, , drop = FALSE])
	# Take care of the Asymmetry Matrices
	if(modelinc[6]>0){
		Ibar = .asymI(Z)
		aZ = Ibar*Z
		Nbar = cov(aZ[1:dcc.old, , drop = FALSE])
	} else{
		Ibar = .asymI(Z)
		# zero what is not needed
		aZ = Ibar*Z*0
		Nbar = matrix(0, m, m)
	}
	Z = rbind( matrix(0, nrow = mx, ncol = m), Z )
	aZ = rbind( matrix(0, nrow = mx, ncol = m), aZ )
	
	
	sumdcc = sum(ipars[idx["dcca",1]:idx["dcca",2],1])+sum(ipars[idx["dccb",1]:idx["dccb",2],1])
	sumdccg = sum(ipars[idx["dccg",1]:idx["dccg",2],1])	
	
	
	res = .Call("copulaNormalC2", model = as.integer(modelinc), 
			pars = as.numeric(ipars[,1]), idx = as.integer(idx[,1]-1), 
			Qbar = Qbar, Nbar = Nbar, Z = Z, N = aZ, 
			epars = c(sumdcc, sumdccg, mx), PACKAGE = "rmgarch")
	
	copula.llhvec = res[[2]][(mx+1):(T+mx)]
	full.llhvec = apply(cbind(garch.llhvec, -copula.llhvec), 1, "sum")
	#ttmp = sum(-copula.llhvec) - sum(apply(garch.llhvec, 1, "sum"))
	#.ttmp <<-ttmp
	sol = list()
	sol$llhvec = full.llhvec
	sol$llh = sum(full.llhvec)
	
	ans = switch(arglist$returnType,
			lik = sol$llhvec,
			llh = sum(full.llhvec),
			ALL = list(lik = sol$llhvec, llh = sum(full.llhvec), Rt = res[[4]][(mx+1):(T+mx)], Qt = res[[1]][(mx+1):(T+mx)],
					Z = Z, Nbar = Nbar, Qbar = Qbar))
	return(ans)
}

copula.normalLLH1 = function(pars, arglist)
{
	mgarchenv = arglist$mgarchenv
	model = arglist$model
	modelinc = model$modelinc
	estidx = arglist$estidx
	idx = model$pidx
	ipars = arglist$ipars
	ipars[estidx, 1] = pars
	trace = arglist$trace
	m = arglist$m
	fit.control = arglist$fit.control
	udata = arglist$ures
	T = dim(udata)[1]
	Z = qnorm(udata)
	if(modelinc[3]>0){
		Rbar = .Pconstruct(ipars[idx["C", 1]:idx["C", 2],1])
	} else{
		Rbar = cor(Z)
	}
	diag(Rbar) = 1
	
	eigR = eigen(Rbar)$values
	#print(eigR)
	if(any(is.complex(eigR)) || min(eigR) < 0){
		warning("\nNon p.s.d. covariance matrix...adjusting")
		Rbar = .makeposdef(Rbar)
	}
	
	res = .Call("copulaNormalC1", Rbar = Rbar, Z = Z, PACKAGE = "rmgarch")
	
	llh = res[[2]]
	
	if(is.finite(llh) | !is.na(llh) | !is.nan(llh)){
		assign("rmgarch_llh", llh, envir = mgarchenv)
	} else {
		llh = (get("rmgarch_llh", mgarchenv) + 0.1*(abs(get("rmgarch_llh", mgarchenv))))
	}
	
	ans = switch(arglist$returnType,
			lik = res[[1]],
			llh = llh,
			ALL = list(lik = res[[1]], llh = res[[2]], Rt = Rbar))
	return(ans)
}

copula.normalLLH2 = function(pars, arglist)
{
	.eps = 	.Machine$double.eps
	model = arglist$model
	umodel = arglist$umodel
	modelinc = model$modelinc
	estidx = arglist$estidx
	eidx = arglist$eidx
	midx = arglist$midx
	mpars = arglist$mpars
	idx = model$pidx	
	ipars = arglist$ipars
	m = arglist$m
	T = arglist$T
	mpars[which(eidx==1, arr.ind = TRUE)] = pars
	# assign the pars to the ipars (used in the DCC diltering)
	ipars[estidx,1] = mpars[which(eidx[,m+1]==1),m+1]
	trace = arglist$trace
	data = arglist$data
	n.start = model$modeldata$n.start	
	garch.llhvec = resids = H = matrix(0, nrow = T, ncol = m)
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, umodel$vt)
	specx = vector(mode = "list", length = m)
	for(i in 1:m){
		specx[[i]] = mspec@spec[[i]]
		setfixed(specx[[i]]) = as.list(mpars[which(midx[,i]==1), i])
	}
	flt = multifilter(multifitORspec = multispec(specx), data = xts(data, arglist$index[1:nrow(data)]), 
			out.sample = n.start, realizedVol = arglist$realizedVol[1:nrow(data),])
	garch.llhvec = sapply(flt@filter, FUN = function(x) x@filter$log.likelihoods)
	H = sigma(flt)
	resids = residuals(flt)
	stdres = resids/H
	
	T = dim(stdres)[1]
	transformation = model$modeldesc$transformation
	spd.control = model$spd.control
	ans = switch(transformation,
			parametric = .pparametric.filter(flt, stdres),
			empirical = .pempirical(stdres),
			spd = .pspd(stdres, spd.control))
	if(transformation == "spd") {
		ures = ans$ures
		sfit = ans$sfit
	} else{
		ures = ans
		sfit = NULL
	}
	# make tail adjustments in order to avoid problem with the quantile functions in
	# optimization
	if(any(ures > 0.99999)){
		xn = which(ures > 0.99999)
		ures[xn] = 0.99999
	}
	if(any(ures < .eps)){
		xn = which(ures < (1.5*.eps))
		ures[xn] = .eps
	}
	
	Z = qnorm(ures)
	
	if(modelinc[3]>0){
		Rbar = .Pconstruct(ipars[idx["C", 1]:idx["C", 2],1])
	} else{
		Rbar = cor(Z)
	}
	diag(Rbar) = 1
	
	eigR = eigen(Rbar)$values
	#print(eigR)
	if(any(is.complex(eigR)) || min(eigR) < 0){
		warning("\nNon p.s.d. covariance matrix...adjusting")
		Rbar = .makeposdef(Rbar)
	}
	
	res = .Call("copulaNormalC1", Rbar = Rbar, Z = Z, PACKAGE = "rmgarch")
	
	copula.llhvec = res[[1]]
	full.llhvec = apply(cbind(garch.llhvec, -copula.llhvec), 1, "sum")
	sol = list()
	sol$llhvec = full.llhvec
	sol$llh = sum(full.llhvec)
	
	ans = switch(arglist$returnType,
			lik = sol$llhvec,
			llh = sum(full.llhvec),
			ALL = list(lik = sol$llhvec, llh = sum(full.llhvec), Rt = Rbar))
	return(ans)
}

copula.normalLLH3 = function(arglist)
{
	.eps = 	.Machine$double.eps
	model = arglist$model
	umodel = arglist$umodel
	modelinc = model$modelinc
	estidx = arglist$estidx
	eidx = arglist$eidx
	midx = arglist$midx
	mpars = arglist$mpars
	idx = model$pidx
	ipars = arglist$ipars
	m = arglist$m
	T = arglist$T
	n.old = arglist$n.old
	filter.control = arglist$filter.control
	spd.control = arglist$spd.control
	flt = arglist$filterlist
	trace = arglist$trace
	data = arglist$data
	n.start = model$modeldata$n.start	
	garch.llhvec = resids = H = matrix(0, nrow = T, ncol = m)
	garch.llhvec = sapply(flt@filter, FUN = function(x) x@filter$log.likelihoods)
	H = sigma(flt)
	resids = residuals(flt)
	stdres = resids/H
	
	if(is.null(filter.control$n.old)) dcc.old = dim(stdres)[1] else dcc.old = n.old
	
	T = dim(stdres)[1]
	transformation = model$modeldesc$transformation
	spd.control = model$spd.control
	ans = switch(transformation,
			parametric = .pparametric.filter(flt, stdres),
			empirical = .pempirical.filter(stdres, dcc.old),
			spd = .pspd.filter(stdres, spd.control, dcc.old))
	if(transformation == "spd") {
		ures = ans$ures
		sfit = ans$sfit
	} else{
		ures = ans
		sfit = NULL
	}
	# make tail adjustments in order to avoid problem with the quantile functions in
	# optimization
	if(any(ures > 0.99999)){
		xn = which(ures > 0.99999)
		ures[xn] = 0.99999
	}
	if(any(ures < .eps)){
		xn = which(ures < (1.5*.eps))
		ures[xn] = .eps
	}
	
	Z = qnorm(ures)
	
	if(modelinc[3]>0){
		Rbar = .Pconstruct(ipars[idx["C", 1]:idx["C", 2],1])
	} else{
		Rbar = cor(Z[1:dcc.old, , drop = FALSE])
	}
	diag(Rbar) = 1
	
	eigR = eigen(Rbar)$values
	#print(eigR)
	if(any(is.complex(eigR)) || min(eigR) < 0){
		warning("\nNon p.s.d. covariance matrix...adjusting")
		Rbar = .makeposdef(Rbar)
	}
	
	res = .Call("copulaNormalC1", Rbar = Rbar, Z = Z, PACKAGE = "rmgarch")
	
	copula.llhvec = res[[1]]
	full.llhvec = apply(cbind(garch.llhvec, -copula.llhvec), 1, "sum")
	sol = list()
	sol$llhvec = full.llhvec
	sol$llh = sum(full.llhvec)
	
	ans = switch(arglist$returnType,
			lik = sol$llhvec,
			llh = sum(full.llhvec),
			ALL = list(lik = sol$llhvec, llh = sum(full.llhvec), Rt = Rbar))
	return(ans)
}

copula.tvstudentLLH1 = function(pars, arglist)
{
	mgarchenv = arglist$mgarchenv
	model = arglist$model
	modelinc = model$modelinc
	estidx = arglist$estidx
	idx = model$pidx
	ipars = arglist$ipars
	ipars[estidx, 1] = pars
	trace = arglist$trace
	m = arglist$m
	mx = model$maxdccOrder
	fit.control = arglist$fit.control
	sumdcc = sum(ipars[idx["dcca",1]:idx["dcca",2],1])+sum(ipars[idx["dccb",1]:idx["dccb",2],1])
	sumdccg = sum(ipars[idx["dccg",1]:idx["dccg",2],1])
	udata = arglist$ures
	T = dim(udata)[1]
	Z = matrix(rugarch:::qstd(udata, shape = ipars[idx["mshape",1], 1]), ncol = m)
	Qbar = cov(Z)
	# Take care of the Asymmetry Matrices
	if(modelinc[6]>0){
		Ibar = .asymI(Z)
		aZ = Ibar*Z
		Nbar = cov(aZ)
	} else{
		Ibar = .asymI(Z)
		# zero what is not needed
		aZ = Ibar*Z*0
		Nbar = matrix(0, m, m)
	}
	dtZ = log(matrix(rugarch:::dstd(Z, shape = ipars[idx["mshape",1], 1]), ncol = m))
	dtZ = rbind( matrix(0, nrow = mx, ncol = m), dtZ )
	Z = rbind( matrix(0, nrow = mx, ncol = m), Z )
	aZ = rbind( matrix(0, nrow = mx, ncol = m), aZ )
	
	arglist$Nbar = Nbar
	arglist$Qbar = Qbar
	
	
	if( fit.control$stationarity ){
		if(modelinc[6]>0){
			persist = .adcccon(pars, arglist)
		} else{
			persist = .dcccon(pars, arglist)
		}
		if( !is.na(persist) && persist >= 1 ) return(llh = get("rmgarch_llh", mgarchenv) + 0.1*(abs(get("rmgarch_llh", mgarchenv))))
	}
	res = .Call("copulaStudentC2", model = as.integer(modelinc), 
			pars = as.numeric(ipars[,1]), idx = as.integer(idx[,1]-1), 
			Qbar = Qbar, Nbar = Nbar, Z = Z, N = aZ, dtZ = dtZ, 
			epars = c(sumdcc, sumdccg, mx), PACKAGE = "rmgarch")
	
	llh = res[[3]]
	
	if(is.finite(llh) | !is.na(llh) | !is.nan(llh)){
		assign("rmgarch_llh", llh, envir = mgarchenv)
	} else {
		llh = (get("rmgarch_llh", mgarchenv) + 0.1*(abs(get("rmgarch_llh", mgarchenv))))
	}
	
	ans = switch(arglist$returnType,
			lik = res[[2]][(mx+1):(T+mx)],
			llh = llh,
			ALL = list(lik = res[[2]][(mx+1):(T+mx)], llh = res[[3]], Rt = res[[4]][(mx+1):(T+mx)], Qt = res[[1]][(mx+1):(T+mx)]))
	return(ans)
}

copula.tvstudentLLH2 = function(pars, arglist)
{
	.eps = 	.Machine$double.eps
	model = arglist$model
	umodel = arglist$umodel
	modelinc = model$modelinc
	estidx = arglist$estidx
	eidx = arglist$eidx
	midx = arglist$midx
	mpars = arglist$mpars
	idx = model$pidx	
	ipars = arglist$ipars
	m = arglist$m
	T = arglist$T
	mpars[which(eidx==1, arr.ind = TRUE)] = pars
	# assign the pars to the ipars (used in the DCC diltering)
	ipars[estidx,1] = mpars[which(eidx[,m+1]==1),m+1]
	trace = arglist$trace	
	data = arglist$data
	n.start = model$modeldata$n.start
	
	mx = model$maxdccOrder
	sumdcc = sum(ipars[idx["dcca",1]:idx["dcca",2],1])+sum(ipars[idx["dccb",1]:idx["dccb",2],1])
	sumdccg = sum(ipars[idx["dccg",1]:idx["dccg",2],1])
	garch.llhvec = resids = H = matrix(0, nrow = T, ncol = m)
	
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, umodel$vt)
	specx = vector(mode = "list", length = m)
	for(i in 1:m){
		specx[[i]] = mspec@spec[[i]]
		setfixed(specx[[i]]) = as.list(mpars[which(midx[,i]==1), i])
	}
	flt = multifilter(multifitORspec = multispec(specx), data = xts(data, arglist$index[1:nrow(data)]), 
			out.sample = n.start, realizedVol = arglist$realizedVol[1:nrow(data),])
	garch.llhvec = sapply(flt@filter, FUN = function(x) x@filter$log.likelihoods)
	H = sigma(flt)
	resids = residuals(flt)
	stdres = resids/H
	T = dim(stdres)[1]
	transformation = model$modeldesc$transformation
	spd.control = model$spd.control
	ans = switch(transformation,
			parametric = .pparametric.filter(flt, stdres),
			empirical = .pempirical(stdres),
			spd = .pspd(stdres, spd.control))
	if(transformation == "spd") {
		ures = ans$ures
		sfit = ans$sfit
	} else{
		ures = ans
		sfit = NULL
	}
	# make tail adjustments in order to avoid problem with the quantile functions in
	# optimization
	if(any(ures > 0.99999)){
		xn = which(ures > 0.99999)
		ures[xn] = 0.99999
	}
	if(any(ures < .eps)){
		xn = which(ures < (1.5*.eps))
		ures[xn] = .eps
	}	
	Z = matrix(rugarch:::qstd(ures, shape = ipars[idx["mshape",1], 1]), ncol = m)
	Qbar = cov(Z)
	# Take care of the Asymmetry Matrices
	if(modelinc[6]>0){
		Ibar = .asymI(Z)
		aZ = Ibar*Z
		Nbar = cov(aZ)
	} else{
		Ibar = .asymI(Z)
		# zero what is not needed
		aZ = Ibar*Z*0
		Nbar = matrix(0, m, m)
	}
	dtZ = log(matrix(rugarch:::dstd(Z, shape = ipars[idx["mshape",1], 1]), ncol = m))
	dtZ = rbind( matrix(0, nrow = mx, ncol = m), dtZ )
	Z = rbind( matrix(0, nrow = mx, ncol = m), Z )
	aZ = rbind( matrix(0, nrow = mx, ncol = m), aZ )
	
	sumdcc = sum(ipars[idx["dcca",1]:idx["dcca",2],1])+sum(ipars[idx["dccb",1]:idx["dccb",2],1])
	sumdccg = sum(ipars[idx["dccg",1]:idx["dccg",2],1])	

	
	res = .Call("copulaStudentC2", model = as.integer(modelinc), 
			pars = as.numeric(ipars[,1]), idx = as.integer(idx[,1]-1), 
			Qbar = Qbar, Nbar = Nbar, Z = Z, N = aZ, dtZ = dtZ, 
			epars = c(sumdcc, sumdccg, mx), PACKAGE = "rmgarch")

	copula.llhvec = res[[2]][(mx+1):(T+mx)]
	full.llhvec = apply(cbind(garch.llhvec, -copula.llhvec), 1, "sum")
	sol = list()
	sol$llhvec = full.llhvec
	sol$llh = sum(full.llhvec)
	
	ans = switch(arglist$returnType,
			lik = sol$llhvec,
			llh = sum(full.llhvec),
			ALL = list(lik = sol$llhvec, llh = sum(full.llhvec), Rt = res[[4]][(mx+1):(T+mx)], Qt = res[[1]][(mx+1):(T+mx)],
					Z = Z, Nbar = Nbar, Qbar = Qbar))
	return(ans)
}

copula.tvstudentLLH3 = function(arglist)
{
	.eps = 	.Machine$double.eps
	model = arglist$model
	umodel = arglist$umodel
	modelinc = model$modelinc
	estidx = arglist$estidx
	eidx = arglist$eidx
	midx = arglist$midx
	mpars = arglist$mpars
	idx = model$pidx
	ipars = arglist$ipars
	m = arglist$m
	T = arglist$T
	filter.control = arglist$filter.control
	spd.control = arglist$spd.control
	n.old = arglist$n.old
	flt = arglist$filterlist
	trace = arglist$trace
	data = arglist$data
	n.start = model$modeldata$n.start
	mx = model$maxdccOrder
	garch.llhvec = resids = H = matrix(0, nrow = T, ncol = m)
	garch.llhvec = sapply(flt@filter, FUN = function(x) x@filter$log.likelihoods)
	H = sigma(flt)
	resids = residuals(flt)
	stdres = resids/H

	if(is.null(filter.control$n.old)) dcc.old = dim(stdres)[1] else dcc.old = n.old
	
	T = dim(stdres)[1]
	transformation = model$modeldesc$transformation
	ans = switch(transformation,
			parametric = .pparametric.filter(flt, stdres),
			empirical = .pempirical.filter(stdres, dcc.old),
			spd = .pspd.filter(stdres, spd.control, dcc.old))
	if(transformation == "spd") {
		ures = ans$ures
		sfit = ans$sfit
	} else{
		ures = ans
		sfit = NULL
	}
	# make tail adjustments in order to avoid problem with the quantile functions in
	# optimization
	if(any(ures > 0.99999)){
		xn = which(ures > 0.99999)
		ures[xn] = 0.99999
	}
	if(any(ures < .eps)){
		xn = which(ures < (1.5*.eps))
		ures[xn] = .eps
	}
	Z = matrix(rugarch:::qstd(ures, shape = ipars[idx["mshape",1], 1]), ncol = m)
	dtZ = log(matrix(rugarch:::dstd(Z, shape = ipars[idx["mshape",1], 1]), ncol = m))	
	Qbar = cov(Z[1:dcc.old, , drop = FALSE])
	# Take care of the Asymmetry Matrices
	if(modelinc[6]>0){
		Ibar = .asymI(Z)
		aZ = Ibar*Z
		Nbar = cov(aZ[1:dcc.old, , drop = FALSE])
	} else{
		Ibar = .asymI(Z)
		# zero what is not needed
		aZ = Ibar*Z*0
		Nbar = matrix(0, m, m)
	}
	Z = rbind( matrix(0, nrow = mx, ncol = m), Z )
	aZ = rbind( matrix(0, nrow = mx, ncol = m), aZ )
	dtZ = rbind( matrix(0, nrow = mx, ncol = m), dtZ )
	
	sumdcc = sum(ipars[idx["dcca",1]:idx["dcca",2],1])+sum(ipars[idx["dccb",1]:idx["dccb",2],1])
	sumdccg = sum(ipars[idx["dccg",1]:idx["dccg",2],1])
	
	
	res = .Call("copulaStudentC2", model = as.integer(modelinc), 
			pars = as.numeric(ipars[,1]), idx = as.integer(idx[,1]-1), 
			Qbar = Qbar, Nbar = Nbar, Z = Z, N = aZ, dtZ = dtZ, 
			epars = c(sumdcc, sumdccg, mx), PACKAGE = "rmgarch")
	
	copula.llhvec = res[[2]][(mx+1):(T+mx)]
	full.llhvec = apply(cbind(garch.llhvec, -copula.llhvec), 1, "sum")
	sol = list()
	sol$llhvec = full.llhvec
	sol$llh = sum(full.llhvec)
	
	ans = switch(arglist$returnType,
			lik = sol$llhvec,
			llh = sum(full.llhvec),
			ALL = list(lik = sol$llhvec, llh = sum(full.llhvec), Rt = res[[4]][(mx+1):(T+mx)], Qt = res[[1]][(mx+1):(T+mx)],
					Z = Z, Nbar = Nbar, Qbar = Qbar))
	return(ans)
}

copula.studentLLH1 = function(pars, arglist)
{
	mgarchenv = arglist$mgarchenv
	model = arglist$model
	modelinc = model$modelinc
	estidx = arglist$estidx
	idx = model$pidx
	ipars = arglist$ipars
	ipars[estidx, 1] = pars
	trace = arglist$trace
	m = arglist$m
	fit.control = arglist$fit.control
	udata = arglist$ures
	T = dim(udata)[1]
	Z = matrix(rugarch:::qstd(udata, shape = ipars[idx["mshape",1], 1]), ncol = m)
	# Z = qt(udata, df = ipars[idx["mshape",1], 1])
	# dtZ = dt(Z, df = ipars[idx["mshape",1], 1], log = TRUE)
	dtZ = log(matrix(rugarch:::dstd(Z, shape = ipars[idx["mshape",1], 1]), ncol = m))
	if(modelinc[3]>0){
		Rbar = .Pconstruct(ipars[idx["C", 1]:idx["C", 2],1])
	} else{
		Rtau = cor.fk(Z)
		Rbar = sin(pi * Rtau/2)
	}
	diag(Rbar) = 1
	
	eigR = eigen(Rbar)$values
	#print(eigR)
	if(any(is.complex(eigR)) || min(eigR) < 0){
		warning("\nNon p.s.d. covariance matrix...adjusting")
		Rbar = .makeposdef(Rbar)
	}
	
	res = .Call("copulaStudentC1", pars = as.numeric(ipars[,1]), 
			idx = as.integer(idx[,1]-1), Rbar = Rbar, Z = Z, dtZ = dtZ, 
			PACKAGE = "rmgarch")
	
	llh = res[[2]]
	
	if(is.finite(llh) | !is.na(llh) | !is.nan(llh)){
		assign("rmgarch_llh", llh, envir = mgarchenv)
	} else {
		llh = (get("rmgarch_llh", mgarchenv) + 0.1*(abs(get("rmgarch_llh", mgarchenv))))
	}
	
	ans = switch(arglist$returnType,
			lik = res[[1]],
			llh = llh,
			ALL = list(lik = res[[1]], llh = res[[2]], Rt = Rbar))
	return(ans)
}

copula.studentLLH2 = function(pars, arglist)
{
	.eps = 	.Machine$double.eps
	model = arglist$model
	umodel = arglist$umodel
	modelinc = model$modelinc
	estidx = arglist$estidx
	eidx = arglist$eidx
	midx = arglist$midx
	mpars = arglist$mpars
	idx = model$pidx	
	ipars = arglist$ipars
	m = arglist$m
	T = arglist$T
	mpars[which(eidx==1, arr.ind = TRUE)] = pars
	# assign the pars to the ipars (used in the DCC diltering)
	ipars[estidx,1] = mpars[which(eidx[,m+1]==1),m+1]
	trace = arglist$trace
	data = arglist$data
	n.start = model$modeldata$n.start
	garch.llhvec = resids = H = matrix(0, nrow = T, ncol = m)
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, umodel$vt)
	specx = vector(mode = "list", length = m)
	for(i in 1:m){
		specx[[i]] = mspec@spec[[i]]
		setfixed(specx[[i]]) = as.list(mpars[which(midx[,i]==1), i])
	}
	flt = multifilter(multifitORspec = multispec(specx), data = xts(data, arglist$index[1:nrow(data)]), 
			out.sample = n.start, realizedVol = arglist$realizedVol[1:nrow(data),])
	garch.llhvec = sapply(flt@filter, FUN = function(x) x@filter$log.likelihoods)
	H = sigma(flt)
	resids = residuals(flt)
	stdres = resids/H
	T = dim(stdres)[1]
	transformation = model$modeldesc$transformation
	spd.control = model$spd.control
	ans = switch(transformation,
			parametric = .pparametric.filter(flt, stdres),
			empirical = .pempirical(stdres),
			spd = .pspd(stdres, spd.control))
	if(transformation == "spd") {
		ures = ans$ures
		sfit = ans$sfit
	} else{
		ures = ans
		sfit = NULL
	}
	# make tail adjustments in order to avoid problem with the quantile functions in
	# optimization
	if(any(ures > 0.99999)){
		xn = which(ures > 0.99999)
		ures[xn] = 0.99999
	}
	if(any(ures < .eps)){
		xn = which(ures < (1.5*.eps))
		ures[xn] = .eps
	}
	
	# Z = qt(ures, df = ipars[idx["mshape",1], 1])
	# dtZ = dt(Z, df = ipars[idx["mshape",1], 1], log = TRUE)
	Z = matrix(rugarch:::qstd(ures, shape = ipars[idx["mshape",1], 1]), ncol = m)
	dtZ = log(matrix(rugarch:::dstd(Z, shape = ipars[idx["mshape",1], 1]), ncol = m))
	
	
	if(modelinc[3]>0){
		Rbar = .Pconstruct(ipars[idx["C", 1]:idx["C", 2],1])
	} else{
		Rtau = cor.fk(Z)
		Rbar = sin(pi * Rtau/2)
	}
	diag(Rbar) = 1
	
	eigR = eigen(Rbar)$values
	#print(eigR)
	if(any(is.complex(eigR)) || min(eigR) < 0){
		warning("\nNon p.s.d. covariance matrix...adjusting")
		Rbar = .makeposdef(Rbar)
	}
	
	res = .Call("copulaStudentC1", pars = as.numeric(ipars[,1]), 
			idx = as.integer(idx[,1]-1), Rbar = Rbar, Z = Z, dtZ = dtZ, 
			PACKAGE = "rmgarch")
	
	copula.llhvec = res[[1]]
	full.llhvec = apply(cbind(garch.llhvec, -copula.llhvec), 1, "sum")
	sol = list()
	sol$llhvec = full.llhvec
	sol$llh = sum(full.llhvec)
	
	ans = switch(arglist$returnType,
			lik = sol$llhvec,
			llh = sum(full.llhvec),
			ALL = list(lik = sol$llhvec, llh = sum(full.llhvec), Rt = Rbar))
	return(ans)
}

copula.studentLLH3 = function(arglist)
{
	.eps = 	.Machine$double.eps
	model = arglist$model
	umodel = arglist$umodel
	modelinc = model$modelinc
	estidx = arglist$estidx
	eidx = arglist$eidx
	midx = arglist$midx
	mpars = arglist$mpars
	idx = model$pidx
	ipars = arglist$ipars
	m = arglist$m
	T = arglist$T
	n.old = arglist$n.old
	filter.control = arglist$filter.control
	spd.control = arglist$spd.control
	flt = arglist$filterlist
	trace = arglist$trace
	data = arglist$data
	n.start = model$modeldata$n.start
	garch.llhvec = resids = H = matrix(0, nrow = T, ncol = m)
	garch.llhvec = sapply(flt@filter, FUN = function(x) x@filter$log.likelihoods)
	H = sigma(flt)
	resids = residuals(flt)
	stdres = resids/H
	if(is.null(filter.control$n.old)) dcc.old = dim(stdres)[1] else dcc.old = n.old
	
	T = dim(stdres)[1]
	transformation = model$modeldesc$transformation
	spd.control = model$spd.control
	ans = switch(transformation,
			parametric = .pparametric.filter(flt, stdres),
			empirical = .pempirical.filter(stdres, dcc.old),
			spd = .pspd.filter(stdres, spd.control, dcc.old))
	if(transformation == "spd") {
		ures = ans$ures
		sfit = ans$sfit
	} else{
		ures = ans
		sfit = NULL
	}
	# make tail adjustments in order to avoid problem with the quantile functions in
	# optimization
	if(any(ures > 0.99999)){
		xn = which(ures > 0.99999)
		ures[xn] = 0.99999
	}
	if(any(ures < .eps)){
		xn = which(ures < (1.5*.eps))
		ures[xn] = .eps
	}
	
	#Z = qt(ures, df = ipars[idx["mshape",1], 1])
	#dtZ = dt(Z, df = ipars[idx["mshape",1], 1], log = TRUE)
	Z = matrix(rugarch:::qstd(ures, shape = ipars[idx["mshape",1], 1]), ncol = m)
	dtZ = log(matrix(rugarch:::dstd(Z, shape = ipars[idx["mshape",1], 1]), ncol = m))	
	
	if(modelinc[3]>0){
		Rbar = .Pconstruct(ipars[idx["C", 1]:idx["C", 2],1])
	} else{
		Rtau = cor.fk(Z[1:dcc.old, , drop = FALSE])
		Rbar = sin(pi * Rtau/2)
	}
	diag(Rbar) = 1
	
	eigR = eigen(Rbar)$values
	#print(eigR)
	if(any(is.complex(eigR)) || min(eigR) < 0){
		warning("\nNon p.s.d. covariance matrix...adjusting")
		Rbar = .makeposdef(Rbar)
	}
	
	res = .Call("copulaStudentC1", pars = as.numeric(ipars[,1]), 
			idx = as.integer(idx[,1]-1), Rbar = Rbar, Z = Z, dtZ = dtZ, 
			PACKAGE = "rmgarch")
	
	copula.llhvec = res[[1]]
	full.llhvec = apply(cbind(garch.llhvec, -copula.llhvec), 1, "sum")
	sol = list()
	sol$llhvec = full.llhvec
	sol$llh = sum(full.llhvec)
	
	ans = switch(arglist$returnType,
			lik = sol$llhvec,
			llh = sum(full.llhvec),
			ALL = list(lik = sol$llhvec, llh = sum(full.llhvec), Rt = Rbar))
	return(ans)
}
