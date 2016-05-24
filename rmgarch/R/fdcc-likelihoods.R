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
normal.fdccLLH1 = function(pars, arglist)
{
	# prepare inputs
	# rejoin fixed and pars
	mgarchenv = arglist$mgarchenv
	model = arglist$model
	modelinc = model$modelinc
	estidx = arglist$estidx
	idx = model$pidx
	ipars = arglist$ipars
	ipars[estidx, 1] = pars
	trace = arglist$trace
	m = arglist$m
	sidx = arglist$sidx
	Qbar = arglist$Qbar
	if(modelinc[4]>0){
		A = matrix(ipars[model$pidx["fdcca",1]:model$pidx["fdcca",2],1], ncol = modelinc[4], byrow=FALSE)
		A = apply(A, 2, function(x) sidx %*% x)
	} else{
		A = matrix(0, 1, 1)
	}
	if(modelinc[5]>0){
		B = matrix(ipars[model$pidx["fdccb",1]:model$pidx["fdccb",2],1], ncol = modelinc[5], byrow=FALSE)
		B = apply(B, 2, function(x) sidx %*% x)
		
	} else{
		B = matrix(0, 1, 1)
	}
	BB = (B %*% t(B))
	AA = (A %*% t(A))
	C = matrix(1, NROW(BB), NCOL(BB)) - BB - AA
	
	stdres = arglist$stdresid
	stdres = rbind( matrix(0, nrow = max(modelinc[4:5]), ncol = m), stdres )
	fit.control = arglist$fit.control
	H = arglist$H
	
	if( fit.control$stationarity ){
		persist = .fdcccon(pars, arglist)
		if( !is.na(persist) && any(persist) >= 1 ) 
			return(llh = get("rmgarch_llh", mgarchenv) + 0.1*(abs(get("rmgarch_llh", mgarchenv))))
	}
	
	res = .Call("fdccnormC1", model = as.integer(modelinc), A = A, B = B, C = C, Qbar = Qbar, 
			Z = stdres, PACKAGE = "rmgarch")
	
	llh = res[[3]]
	
	if(is.finite(llh) | !is.na(llh) | !is.nan(llh)){
		assign("rmgarch_llh", llh, envir = mgarchenv)
	} else {
		llh = (get("rmgarch_llh", mgarchenv) + 0.1*(abs(get("rmgarch_llh", mgarchenv))))
	}
	return(llh)
}

normal.fdccLLH2 = function(pars, arglist)
{
	# prepare inputs
	# rejoin fixed and pars
	mgarchenv = arglist$mgarchenv
	model = arglist$model
	umodel = arglist$umodel
	modelinc = model$modelinc
	idx = model$pidx
	mpars = arglist$mpars
	ipars = arglist$ipars
	m = arglist$m
	estidx = arglist$estidx
	eidx = arglist$eidx
	midx = arglist$midx
	sidx = arglist$sidx
	data = arglist$data
	n.start = model$modeldata$n.start
	# assign the pars to the matrix pars (used in the GARCH filtering)
	mpars[which(eidx==1, arr.ind = TRUE)] = pars
	# assign the pars to the ipars (used in the DCC diltering)
	ipars[estidx,1] = mpars[which(eidx[,m+1]==1),m+1]
	trace = arglist$trace
	T = arglist$T
	if(modelinc[4]>0){
		A = matrix(ipars[model$pidx["fdcca",1]:model$pidx["fdcca",2],1], ncol = modelinc[4], byrow=FALSE)
		A = apply(A, 2, function(x) sidx %*% x)
	} else{
		A = matrix(0, 1, 1)
	}
	if(modelinc[5]>0){
		B = matrix(ipars[model$pidx["fdccb",1]:model$pidx["fdccb",2],1], ncol = modelinc[5], byrow=FALSE)
		B = apply(B, 2, function(x) sidx %*% x)
	} else{
		B = matrix(0, 1, 1)
	}
	BB = (B %*% t(B))
	AA = (A %*% t(A))
	C = matrix(1, NROW(BB), NCOL(BB)) - BB - AA
	
	Qbar = arglist$Qbar
	# the dcc order
	mx = max(modelinc[4:5])
	m = dim(data)[2]
	H = matrix(0, nrow = T, ncol = m)
	resids = matrix(0, nrow = T, ncol = m)
	# simulate new H with which to standardize the dataset
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, umodel$vt)
	
	for(i in 1:m){
		specx = mspec@spec[[i]]
		setfixed(specx) = as.list(mpars[which(eidx[,i]==1), i])
		flt = ugarchfilter(spec = specx, data = xts(data[,i], arglist$index[1:nrow(data)]), out.sample = n.start, realizedVol = arglist$realizedVol[1:nrow(data),i])
		H[, i] = sigma(flt)
		resids[,i] = residuals(flt)
	}
	stdres = resids/H
	T = dim(stdres)[1]
	Qbar = cov(stdres)
	stdres = rbind( matrix(1, nrow = mx, ncol = m),  stdres )
	H = rbind( matrix(0, nrow = mx, ncol = m), H )
	
	res = .Call("fdccnormC2", model = as.integer(modelinc), A = A, B = B, C = C, 
			Qbar = Qbar,  H = H, Z = stdres, PACKAGE = "rmgarch")
	
	Qtout = res[[1]]
	likelihoods = res[[2]]
	llh = res[[3]]
	Rtout = res[[4]]
	Qtout = Qtout[(mx+1):(T+mx)]
	Rtout = Rtout[(mx+1):(T+mx)]
	likelihoods = likelihoods[(mx+1):(T+mx)]
	
	if(is.finite(llh) | !is.na(llh) | !is.nan(llh)){
		assign("rmgarch_llh", llh, envir = mgarchenv)
	} else {
		llh = (get("rmgarch_llh", mgarchenv) + 0.1*(abs(get("rmgarch_llh", mgarchenv))))
	}
	
	ans = switch(arglist$returnType,
			lik = likelihoods,
			llh = llh,
			all = list(lik = likelihoods, llh = llh, Rt = Rtout, Qt = Qtout))
	return(ans)
}


normalfilter.fdccLLH2 = function(pars, arglist)
{
	# prepare inputs
	# rejoin fixed and pars
	mgarchenv = arglist$mgarchenv
	model = arglist$model
	umodel = arglist$umodel
	modelinc = model$modelinc
	idx = model$pidx
	mpars = arglist$mpars
	ipars = arglist$ipars
	m = arglist$m
	estidx = arglist$estidx
	eidx = arglist$eidx
	sidx = arglist$sidx
	midx = arglist$midx
	data = arglist$data
	n.old = arglist$n.old
	dcc.old = arglist$dcc.old
	n.start = model$modeldata$n.start
	# assign the pars to the matrix pars (used in the GARCH filtering)
	mpars[which(eidx==1, arr.ind = TRUE)] = pars
	# assign the pars to the ipars (used in the DCC diltering)
	ipars[estidx,1] = mpars[which(eidx[,m+1]==1),m+1]
	trace = arglist$trace
	T = arglist$T
	Qbar = arglist$Qbar
	if(modelinc[4]>0){
		A = matrix(ipars[model$pidx["fdcca",1]:model$pidx["fdcca",2],1], ncol = modelinc[4], byrow=FALSE)
		A = apply(A, 2, function(x) sidx %*% x)
	} else{
		A = matrix(0, 1, 1)
	}
	if(modelinc[5]>0){
		B = matrix(ipars[model$pidx["fdccb",1]:model$pidx["fdccb",2],1], ncol = modelinc[5], byrow=FALSE)
		B = apply(B, 2, function(x) sidx %*% x)
		
	} else{
		B = matrix(0, 1, 1)
	}
	BB = (B %*% t(B))
	AA = (A %*% t(A))
	C = matrix(1, NROW(BB), NCOL(BB)) - BB - AA
	
	# the dcc order
	mx = max(modelinc[4:5])
	m = dim(data)[2]
	H = matrix(0, nrow = T, ncol = m)
	resids = matrix(0, nrow = T, ncol = m)
	# simulate new H with which to standardize the dataset
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, NULL)
	
	for(i in 1:m){
		specx = mspec@spec[[i]]
		setfixed(specx) = as.list(mpars[which(eidx[,i]==1), i])
		flt = ugarchfilter(spec = specx, data = xts(data[,i], arglist$index[1:nrow(data)]), out.sample = n.start, realizedVol = arglist$realizedVol[1:nrow(data),i])
		H[, i] = sigma(flt)
		resids[,i] = residuals(flt)
	}
	stdres = resids/H
	T = dim(stdres)[1]
	Qbar = cov(stdres[1:dcc.old, ])
	stdres = rbind( matrix(1, nrow = mx, ncol = m),  stdres )
	H = rbind( matrix(0, nrow = mx, ncol = m), H )
	
	res = .Call("fdccnormC2", model = as.integer(modelinc), A = A, B = B, C = C, 
			Qbar = Qbar,  H = H, Z = stdres, PACKAGE = "rmgarch")
	
	Qtout = res[[1]]
	likelihoods = res[[2]]
	llh = res[[3]]
	Rtout = res[[4]]
	Qtout = Qtout[(mx+1):(T+mx)]
	Rtout = Rtout[(mx+1):(T+mx)]
	likelihoods = likelihoods[(mx+1):(T+mx)]
	
	if(is.finite(llh) | !is.na(llh) | !is.nan(llh)){
		assign("rmgarch_llh", llh, envir = mgarchenv)
	} else {
		llh = (get("rmgarch_llh", mgarchenv) + 0.1*(abs(get("rmgarch_llh", mgarchenv))))
	}
	
	ans = switch(arglist$returnType,
			lik = likelihoods,
			llh = llh,
			all = list(lik = likelihoods, llh = llh, Rt = Rtout, Qt = Qtout))
	return(ans)
}