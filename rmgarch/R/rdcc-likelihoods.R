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


normal.dccLLH1 = function(pars, arglist)
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
	mx = model$maxdccOrder
	Qbar = arglist$Qbar
	Nbar = arglist$Nbar
	stdres = arglist$stdresid
	astdres = arglist$astdresid
	stdres = rbind( matrix(0, nrow = mx, ncol = m), stdres )
	astdres = rbind( matrix(0, nrow = mx, ncol = m), astdres )
	fit.control = arglist$fit.control
	H = arglist$H
	sumdcc = sum(ipars[idx["dcca",1]:idx["dcca",2],1])+sum(ipars[idx["dccb",1]:idx["dccb",2],1])
	sumdccg = sum(ipars[idx["dccg",1]:idx["dccg",2],1])
	
	if( fit.control$stationarity ){
		if(modelinc[5]>0){
			persist = .adcccon(pars, arglist)
		} else{
			persist = .dcccon(pars, arglist)
		}
		if( !is.na(persist) && persist >= 1 ) 
			return(llh = get("rmgarch_llh", mgarchenv) + 0.1*(abs(get("rmgarch_llh", mgarchenv))))
	}
	
	res = .Call("dccnormC1", model = as.integer(modelinc), pars = as.numeric(ipars[,1]), 
			idx = as.integer(idx[,1]-1), Qbar = Qbar, Nbar = Nbar, Z = stdres, 
			N = astdres, epars = c(sumdcc, sumdccg, mx), PACKAGE = "rmgarch")
	
	llh = res[[3]]
	
	if(is.finite(llh) | !is.na(llh) | !is.nan(llh)){
		assign("rmgarch_llh", llh, envir = mgarchenv)
	} else {
		llh = (get("rmgarch_llh", mgarchenv) + 0.1*(abs(get("rmgarch_llh", mgarchenv))))
	}
	return(llh)
}


normal.dccLLH2 = function(pars, arglist)
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
	data = arglist$data
	
	n.start = model$modeldata$n.start
	# assign the pars to the matrix pars (used in the GARCH filtering)
	mpars[which(eidx==1, arr.ind = TRUE)] = pars
	# assign the pars to the ipars (used in the DCC diltering)
	ipars[estidx,1] = mpars[which(eidx[,m+1]==1),m+1]
	trace = arglist$trace
	T = arglist$T
	Qbar = arglist$Qbar
	Nbar = arglist$Nbar
	# the dcc order
	mx = model$maxdccOrder
	m = dim(data)[2]
	H = matrix(0, nrow = T, ncol = m)
	resids = matrix(0, nrow = T, ncol = m)
	# simulate new H with which to standardize the dataset
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, NULL)
	
	for(i in 1:m){
		specx = mspec@spec[[i]]
		setfixed(specx) = as.list(mpars[which(midx[,i]==1), i])
		flt = ugarchfilter(spec = specx, data = xts(data[,i], arglist$index[1:nrow(data)]), out.sample = n.start, realizedVol = arglist$realizedVol[1:nrow(data),i])
		H[, i] = sigma(flt)
		resids[,i] = residuals(flt)
	}
	stdres = resids/H
	T = dim(stdres)[1]
	Qbar = cov(stdres)
	if(modelinc[5]>0){
		Ibar = .asymI(stdres)
		astdres = Ibar*stdres
		Nbar = cov(astdres)
	} else{
		Ibar = .asymI(stdres)
		# zero what is not needed
		astdres = Ibar*stdres*0
		Nbar = matrix(0, m, m)
	}
	
	#if(modelinc[5]>0){
	#	ipars[idx["dcca",1]:idx["dccg",2],1] = ipars[idx["dcca",1]:idx["dccg",2],1]^2
	#}
	sumdcc = sum(ipars[idx["dcca",1]:idx["dcca",2],1])+sum(ipars[idx["dccb",1]:idx["dccb",2],1])
	sumdccg = sum(ipars[idx["dccg",1]:idx["dccg",2],1])	
	stdres = rbind( matrix(1, nrow = mx, ncol = m),  stdres )
	astdres = rbind( matrix(1, nrow = mx, ncol = m), astdres )
	H = rbind( matrix(0, nrow = mx, ncol = m), H )
	
	res = .Call("dccnormC2", model = as.integer(modelinc), pars = as.numeric(ipars[,1]), 
			idx = as.integer(idx[,1]-1), Qbar = Qbar, Nbar = Nbar, H = H, Z = stdres, 
			N = astdres, epars = c(sumdcc, sumdccg, mx),
			PACKAGE = "rmgarch")
	
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


normalfilter.dccLLH2 = function(pars, arglist)
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
	Nbar = arglist$Nbar
	# the dcc order
	mx = model$maxdccOrder
	m = dim(data)[2]
	H = matrix(0, nrow = T, ncol = m)
	resids = matrix(0, nrow = T, ncol = m)
	# simulate new H with which to standardize the dataset
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, NULL)
	
	for(i in 1:m){
		specx = mspec@spec[[i]]
		setfixed(specx) = as.list(mpars[which(midx[,i]==1), i])
		flt = ugarchfilter(spec = specx, data =xts(data[,i], arglist$index[1:nrow(data)]), out.sample = n.start, n.old = n.old, realizedVol = arglist$realizedVol[1:nrow(data),i])
		H[, i] = sigma(flt)
		resids[,i] = residuals(flt)
	}
	
	stdres = resids/H

	T = dim(stdres)[1]
	Qbar = cov(stdres[1:dcc.old, ])
	if(modelinc[5]>0){
		Ibar = .asymI(stdres)
		astdres = Ibar*stdres
		Nbar = cov(astdres[1:dcc.old, ])
	} else{
		Ibar = .asymI(stdres)
		# zero what is not needed
		astdres = Ibar*stdres*0
		Nbar = matrix(0, m, m)
	}
	
	#if(modelinc[5]>0){
	#	ipars[idx["dcca",1]:idx["dccg",2],1] = ipars[idx["dcca",1]:idx["dccg",2],1]^2
	#}
	sumdcc = sum(ipars[idx["dcca",1]:idx["dcca",2],1])+sum(ipars[idx["dccb",1]:idx["dccb",2],1])
	sumdccg = sum(ipars[idx["dccg",1]:idx["dccg",2],1])	
	stdres = rbind( matrix(1, nrow = mx, ncol = m),  stdres )
	astdres = rbind( matrix(1, nrow = mx, ncol = m), astdres )
	H = rbind( matrix(0, nrow = mx, ncol = m), H )
	
	res = .Call("dccnormC2", model = as.integer(modelinc), pars = as.numeric(ipars[,1]), 
			idx = as.integer(idx[,1]-1), Qbar = Qbar, Nbar = Nbar, H = H, Z = stdres, 
			N = astdres, epars = c(sumdcc, sumdccg, mx),
			PACKAGE = "rmgarch")
	
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

student.dccLLH1 = function(pars, arglist)
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
	T = arglist$T
	m = arglist$m
	mx = model$maxdccOrder
	Qbar = arglist$Qbar
	Nbar = arglist$Nbar
	stdres = arglist$stdresid
	astdres = arglist$astdresid
	fit.control = arglist$fit.control
	H = arglist$H
	sumdcc = sum(ipars[idx["dcca",1]:idx["dcca",2],1])+sum(ipars[idx["dccb",1]:idx["dccb",2],1])
	sumdccg = sum(ipars[idx["dccg",1]:idx["dccg",2],1])
	stdres = rbind( matrix(0, nrow = mx, ncol = m), stdres )
	astdres = rbind( matrix(0, nrow = mx, ncol = m), astdres )
	
	if( fit.control$stationarity ){
		if(modelinc[5]>0){
			persist = .adcccon(pars, arglist)
		} else{
			persist = .dcccon(pars, arglist)
		}
		if( !is.na(persist) && persist >= 1 ) 
			return(llh = get("rmgarch_llh", mgarchenv) + 0.1*(abs(get("rmgarch_llh", mgarchenv))))
	}
	
	res = .Call("dccstudentC1", model = as.integer(modelinc), pars = as.numeric(ipars[,1]), 
			idx = as.integer(idx[,1]-1), Qbar = Qbar, Nbar = Nbar, Z = stdres, 
			N = astdres, epars = c(sumdcc, sumdccg, mx), PACKAGE = "rmgarch")
	
	llh = res[[3]]
	
	if(is.finite(llh) | !is.na(llh) | !is.nan(llh)){
		assign("rmgarch_llh", llh, envir = mgarchenv)
	} else {
		llh = (get("rmgarch_llh", mgarchenv) + 0.1*(abs(get("rmgarch_llh", mgarchenv))))
	}
	
	return(llh)
}

student.dccLLH2 = function(pars, arglist)
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
	data = arglist$data
	n.start = model$modeldata$n.start
	# assign the pars to the matrix pars (used in the GARCH filtering)
	mpars[which(eidx==1, arr.ind = TRUE)] = pars
	# assign the pars to the ipars (used in the DCC diltering)
	ipars[estidx,1] = mpars[which(eidx[,m+1]==1),m+1]
	trace = arglist$trace
	T = arglist$T
	Qbar = arglist$Qbar
	Nbar = arglist$Nbar
	mx = model$maxdccOrder
	m = dim(data)[2]
	H = matrix(0, nrow = T, ncol = m)
	resids = matrix(0, nrow = T, ncol = m)
	# simulate new H with which to standardize the dataset
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, NULL)
	
	for(i in 1:m){
		specx = mspec@spec[[i]]
		setfixed(specx) = as.list(mpars[which(midx[,i]==1), i])
		flt = ugarchfilter(spec = specx, data = xts(data[,i], arglist$index[1:nrow(data)]), out.sample = n.start, realizedVol = arglist$realizedVol[1:nrow(data),i])
		H[, i] = sigma(flt)
		resids[,i] = residuals(flt)
	}
	
	stdres = resids/H
	T = dim(stdres)[1]
	
	Qbar = cov(stdres)
	if(modelinc[5]>0){
		Ibar = .asymI(stdres)
		astdres = Ibar*stdres
		Nbar = cov(astdres)
	} else{
		Ibar = .asymI(stdres)
		# zero what is not needed
		astdres = Ibar*stdres*0
		Nbar = matrix(0, m, m)
	}
	sumdcc = sum(ipars[idx["dcca",1]:idx["dcca",2],1])+sum(ipars[idx["dccb",1]:idx["dccb",2],1])
	sumdccg = sum(ipars[idx["dccg",1]:idx["dccg",2],1])	
	stdres = rbind( matrix(1, nrow = mx, ncol = m),  stdres )
	astdres = rbind( matrix(1, nrow = mx, ncol = m), astdres )
	H = rbind( matrix(0, nrow = mx, ncol = m), H )
	
	res = .Call( "dccstudentC2", model = as.integer(modelinc), pars = as.numeric(ipars[,1]), 
			idx = as.integer(idx[,1]-1), Qbar = Qbar, Nbar = Nbar, H = H, 
			Z = stdres, N = astdres, epars = c(sumdcc, sumdccg, mx), PACKAGE = "rmgarch")
	
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

studentfilter.dccLLH2 = function(pars, arglist)
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
	Nbar = arglist$Nbar
	mx = model$maxdccOrder
	m = dim(data)[2]
	H = matrix(0, nrow = T, ncol = m)
	resids = matrix(0, nrow = T, ncol = m)
	# simulate new H with which to standardize the dataset
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, NULL)
	
	for(i in 1:m){
		specx = mspec@spec[[i]]
		setfixed(specx) = as.list(mpars[which(midx[,i]==1), i])
		flt = ugarchfilter(spec = specx, data = xts(data[,i], arglist$index[1:nrow(data)]), out.sample = n.start, n.old = n.old, realizedVol = arglist$realizedVol[1:nrow(data),i])
		H[, i] = sigma(flt)
		resids[,i] = residuals(flt)
	}
	
	stdres = resids/H
	T = dim(stdres)[1]
	
	Qbar = cov(stdres[1:dcc.old, ])
	if(modelinc[5]>0){
		Ibar = .asymI(stdres)
		astdres = Ibar*stdres
		Nbar = cov(astdres[1:dcc.old, ])
	} else{
		Ibar = .asymI(stdres)
		# zero what is not needed
		astdres = Ibar*stdres*0
		Nbar = matrix(0, m, m)
	}
	sumdcc = sum(ipars[idx["dcca",1]:idx["dcca",2],1])+sum(ipars[idx["dccb",1]:idx["dccb",2],1])
	sumdccg = sum(ipars[idx["dccg",1]:idx["dccg",2],1])	
	stdres = rbind( matrix(1, nrow = mx, ncol = m),  stdres )
	astdres = rbind( matrix(1, nrow = mx, ncol = m), astdres )
	H = rbind( matrix(0, nrow = mx, ncol = m), H )
	
	res = .Call( "dccstudentC2", model = as.integer(modelinc), pars = as.numeric(ipars[,1]), 
			idx = as.integer(idx[,1]-1), Qbar = Qbar, Nbar = Nbar, H = H, Z = stdres, 
			N = astdres, epars = c(sumdcc, sumdccg, mx), PACKAGE = "rmgarch")
	
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

laplace.dccLLH1 = function(pars, arglist)
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
	T = arglist$T
	m = arglist$m
	mx = model$maxdccOrder
	Qbar = arglist$Qbar
	Nbar = arglist$Nbar
	stdres = arglist$stdresid
	astdres = arglist$astdresid
	fit.control = arglist$fit.control
	H = arglist$H
	N = c(mx, T)
	
	#if(modelinc[5]>0){
	#	ipars[idx["dcca",1]:idx["dccg",2],1] = ipars[idx["dcca",1]:idx["dccg",2],1]^2
	#}
	sumdcc = sum(ipars[idx["dcca",1]:idx["dcca",2],1])+sum(ipars[idx["dccb",1]:idx["dccb",2],1])
	sumdccg = sum(ipars[idx["dccg",1]:idx["dccg",2],1])
	stdres = rbind( matrix(0, nrow = mx, ncol = m), stdres )
	astdres = rbind( matrix(0, nrow = mx, ncol = m), astdres )
	
	if( fit.control$stationarity ){
		if(modelinc[5]>0){
			persist = .adcccon(pars, arglist)
		} else{
			persist = .dcccon(pars, arglist)
		}
		if( !is.na(persist) && persist >= 1 ) 
			return(llh = get("rmgarch_llh", mgarchenv) + 0.1*(abs(get("rmgarch_llh", mgarchenv))))
	}
	
	res = .Call( "dcclaplaceC1", model = as.integer(modelinc), pars = as.numeric(ipars[,1]), 
			idx = as.integer(idx[,1]-1), Qbar = Qbar, Nbar = Nbar, Z = stdres, 
			N = astdres, epars = c(sumdcc, sumdccg, mx), PACKAGE = "rmgarch")
	
	llh = res[[3]]
	
	if(is.finite(llh) | !is.na(llh) | !is.nan(llh)){
		assign("rmgarch_llh", llh, envir = mgarchenv)
	} else {
		llh = (get("rmgarch_llh", mgarchenv) + 0.1*(abs(get("rmgarch_llh", mgarchenv))))
	}
	
	return(llh)
}

laplace.dccLLH2 = function(pars, arglist)
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
	data = arglist$data
	n.start = model$modeldata$n.start
	# assign the pars to the matrix pars (used in the GARCH filtering)
	mpars[which(eidx==1, arr.ind = TRUE)] = pars
	# assign the pars to the ipars (used in the DCC diltering)
	ipars[estidx,1] = mpars[which(eidx[,m+1]==1),m+1]
	trace = arglist$trace
	T = arglist$T
	Qbar = arglist$Qbar
	Nbar = arglist$Nbar
	mx = model$maxdccOrder
	m = dim(data)[2]
	H = matrix(0, nrow = T, ncol = m)
	resids = matrix(0, nrow = T, ncol = m)
	# simulate new H with which to standardize the dataset
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, NULL)
	
	for(i in 1:m){
		specx = mspec@spec[[i]]
		setfixed(specx) = as.list(mpars[which(midx[,i]==1), i])
		flt = ugarchfilter(spec = specx, data = xts(data[,i], arglist$index[1:nrow(data)]), out.sample = n.start, realizedVol = arglist$realizedVol[1:nrow(data),i])
		H[, i] = sigma(flt)
		resids[,i] = residuals(flt)
	}
	
	stdres = resids/H
	T = dim(stdres)[1]
	Qbar = cov(stdres)
	if(modelinc[5]>0){
		Ibar = .asymI(stdres)
		astdres = Ibar*stdres
		Nbar = cov(astdres)
	} else{
		Ibar = .asymI(stdres)
		# zero what is not needed
		astdres = Ibar*stdres*0
		Nbar = matrix(0, m, m)
	}
	sumdcc = sum(ipars[idx["dcca",1]:idx["dcca",2],1])+sum(ipars[idx["dccb",1]:idx["dccb",2],1])
	sumdccg = sum(ipars[idx["dccg",1]:idx["dccg",2],1])	
	stdres = rbind( matrix(1, nrow = mx, ncol = m),  stdres )
	astdres = rbind( matrix(1, nrow = mx, ncol = m), astdres )
	H = rbind( matrix(0, nrow = mx, ncol = m), H )
	
	res = .Call( "dcclaplaceC2", model = as.integer(modelinc), 
			pars = as.numeric(ipars[,1]), idx = as.integer(idx[,1]-1), 
			Qbar = Qbar, Nbar = Nbar, H = H, Z = stdres, N = astdres, 
			epars = c(sumdcc, sumdccg, mx), PACKAGE = "rmgarch")
	
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


laplacefilter.dccLLH2 = function(pars, arglist)
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
	Nbar = arglist$Nbar
	mx = model$maxdccOrder
	m = dim(data)[2]
	H = matrix(0, nrow = T, ncol = m)
	resids = matrix(0, nrow = T, ncol = m)
	# simulate new H with which to standardize the dataset
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, NULL)
	
	for(i in 1:m){
		specx = mspec@spec[[i]]
		setfixed(specx) = as.list(mpars[which(midx[,i]==1), i])
		flt = ugarchfilter(spec = specx, data = xts(data[,i], arglist$index[1:nrow(data)]), out.sample = n.start, realizedVol = arglist$realizedVol[1:nrow(data),i])
		H[, i] = sigma(flt)
		resids[,i] = residuals(flt)
	}
	
	stdres = resids/H
	T = dim(stdres)[1]
	Qbar = cov(stdres[1:dcc.old, ])
	if(modelinc[5]>0){
		Ibar = .asymI(stdres)
		astdres = Ibar*stdres
		Nbar = cov(astdres[1:dcc.old, ])
	} else{
		Ibar = .asymI(stdres)
		# zero what is not needed
		astdres = Ibar*stdres*0
		Nbar = matrix(0, m, m)
	}
	sumdcc = sum(ipars[idx["dcca",1]:idx["dcca",2],1])+sum(ipars[idx["dccb",1]:idx["dccb",2],1])
	sumdccg = sum(ipars[idx["dccg",1]:idx["dccg",2],1])	
	stdres = rbind( matrix(1, nrow = mx, ncol = m),  stdres )
	astdres = rbind( matrix(1, nrow = mx, ncol = m), astdres )
	H = rbind( matrix(0, nrow = mx, ncol = m), H )
	
	res = .Call( "dcclaplaceC2", model = as.integer(modelinc), 
			pars = as.numeric(ipars[,1]), idx = as.integer(idx[,1]-1), 
			Qbar = Qbar, Nbar = Nbar, H = H, Z = stdres, N = astdres, 
			epars = c(sumdcc, sumdccg, mx), PACKAGE = "rmgarch")
	
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