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

# A set of functions to compute VARX models for internal use fit/filter, forecast and simulate.
varxfit = function(X, p, constant = TRUE, exogen = NULL, robust = FALSE, gamma = 0.25, delta = 0.01, nc = 10, ns = 500, 
		postpad = c("none", "constant", "zero", "NA"), cluster = NULL)
{
	X = as.matrix(X)
	if (any(is.na(X))) stop("\nvarxfit:-->error: NAs in X.\n")
	if (ncol(X) < 2) stop("\nvarxfit:-->error: The matrix 'X' should contain at least two variables.\n")
	if (is.null(colnames(X))) colnames(X) = paste("X", 1:ncol(X), sep = "")
	colnames(X) = make.names(colnames(X))
	if(constant) ic = 1 else ic = 0
	obs = dim(X)[1]
	K = dim(X)[2]
	xsample = obs - p
	Xlags = embed(X, dimension = p + 1)[, -(1:K)]
	temp1 = NULL
	for (i in 1:p) {
		temp = paste(colnames(X), ".l", i, sep = "")
		temp1 = c(temp1, temp)
	}
	colnames(Xlags) = temp1
	Xend = X[-c(1:p), ]
	if(constant){
		rhs = cbind( Xlags, rep(1, xsample))
		colnames(rhs) <- c(colnames(Xlags), "const")
	} else{
		rhs = Xlags
		colnames(rhs) <- colnames(Xlags)		
	}
	if( !(is.null(exogen)) ) {
		exogen = as.matrix(exogen)
		if (!identical(nrow(exogen), nrow(X))) {
			stop("\nvarxfit:-->error: Different row size of X and exogen.\n")
		}
		XK = dim(exogen)[2]
		if (is.null(colnames(exogen))) colnames(exogen) = paste("exo", 1:ncol(exogen), sep = "")
		colnames(exogen) = make.names(colnames(exogen))
		tmp  = colnames(rhs)
		rhs =  cbind(rhs, exogen[-c(1:p), ])
		colnames(rhs) = c(tmp, colnames(exogen))
	} else{
		XK = 0
	}
	datamat = as.matrix(rhs)
	colnames(datamat) = colnames(rhs)
	postpad = tolower(postpad[1])
	if( robust ){
		sol = robustvar(data = X, exogen = exogen, constant = constant, lags = p, alpha = gamma, ns = ns, nc = nc, delta = delta, cluster = cluster)
		Bcoef = t(sol$betaR)
		Bcoef = cbind(Bcoef[,(ic+1):(p*K+ic)], if(constant) Bcoef[,1], if(!is.null(exogen)) Bcoef[,(ic+1+p*K):(ic+p*K+XK)] else NULL)
		colnames(Bcoef) = colnames(rhs)
		rownames(Bcoef) = colnames(X)
		xfitted = t( Bcoef %*% t( datamat ) )
		xresiduals = tail(X, obs - p) - xfitted
		sigma2 = diag( sol$sigmaR )
		if(postpad!="none"){
			if(postpad == "constant"){
				# pre-pad values with the constant
				xfitted = t( Bcoef %*% t( rbind(matrix(c(rep(0, p*K), if(constant) 1 else NULL, if(XK>0) rep(0, XK) else NULL), nrow = p, ncol=dim(Bcoef)[2], byrow = TRUE), datamat ) ) )
				xresiduals = X - xfitted
			} else if(postpad == "zero"){
				xfitted = t( Bcoef %*% t( rbind(matrix(rep(0, dim(Bcoef)[2]), nrow = p, ncol=dim(Bcoef)[2], byrow = TRUE), datamat ) ) )
				xresiduals = X - xfitted
			} else if(postpad == "NA"){
				xfitted = t( Bcoef %*% t( rbind(matrix(rep(NA, dim(Bcoef)[2]), nrow = p, ncol=dim(Bcoef)[2], byrow = TRUE), datamat ) ) )
				xresiduals = X - xfitted
			} else{
				# do nothing
				xfitted = t( Bcoef %*% t( datamat ) )
				xresiduals = tail(X, obs - p) - xfitted
			}
		}
		# I think we need to correct this to get the corrected values from the mlts function
		xp = solve( t(datamat)%*%datamat, diag(p*K + XK + ic)  )
		Bcov = lapply( as.list(sigma2), FUN = function(x) x * xp )
		# calculate standard errors and t-stats
		se = sapply(Bcov, FUN = function(x) sqrt(diag(x)))
		rownames(se) = colnames(rhs)
		tstat = t(Bcoef) / se
		pstat = 2*(1-pnorm(abs(tstat)))
	} else{
		Bcoef = matrix(NA, ncol = dim(datamat)[2], nrow = K)
		Bcoef =  t( apply(Xend, 2, FUN = function(x) solve(t(datamat) %*% datamat) %*% t(datamat) %*% x ) )
		colnames(Bcoef) = colnames(rhs)
		rownames(Bcoef) = colnames(X)
		xfitted = t( Bcoef %*% t( datamat ) )
		xresiduals = tail(X, obs - p) - xfitted
		# apply sigma on the unpadded values
		sigma2 = apply( xresiduals, 2, FUN = function(x) (t(x) %*% x)/(obs - p*K - XK - ic) )
		if(postpad!="none"){
			if(postpad == "constant"){
				# pre-pad values with the constant
				xfitted = t( Bcoef %*% t( rbind(matrix(c(rep(0, p*K), if(constant) 1 else NULL, if(XK>0) rep(0, XK) else NULL), nrow = p, ncol=dim(Bcoef)[2], byrow = TRUE), datamat ) ) )
				xresiduals = X - xfitted
			} else if(postpad == "zero"){
				xfitted = t( Bcoef %*% t( rbind(matrix(rep(0, dim(Bcoef)[2]), nrow = p, ncol=dim(Bcoef)[2], byrow = TRUE), datamat ) ) )
				xresiduals = X - xfitted
			} else if(postpad == "NA"){
				xfitted = t( Bcoef %*% t( rbind(matrix(rep(NA, dim(Bcoef)[2]), nrow = p, ncol=dim(Bcoef)[2], byrow = TRUE), datamat ) ) )
				xresiduals = X - xfitted
			} else{
				# do nothing
				xfitted = t( Bcoef %*% t( datamat ) )
				xresiduals = tail(X, obs - p) - xfitted
			}
		}
		xp = solve( t(datamat)%*%datamat, diag(p*K + XK + ic)  )
		Bcov = lapply( as.list(sigma2), FUN = function(x) x * xp )
		# calculate standard errors and t-stats
		se = sapply(Bcov, FUN = function(x) sqrt(diag(x)))
		rownames(se) = colnames(rhs)
		tstat = t(Bcoef) / se
		pstat = 2*(1-pnorm(abs(tstat)))
	}
	
	# We pad the p lags at the start with the constant mean for the fitted for consistency with ugarch.
	ans = list( Bcoef = Bcoef, xfitted = xfitted, xresiduals = xresiduals, Bcov = Bcov, se = se, tstat = tstat, pstat = pstat, lag = p,
			mxn = XK, constant = constant)
	return( ans )
}


varxfilter = function(X, p, Bcoef, exogen = NULL, postpad = c("none", "constant", "zero", "NA"))
{
	X = as.matrix(X)
	if(any(is.na(X))) stop("\nvarxfilter:-->error: NAs in X.\n")
	if(ncol(X) < 2) stop("\nvarxfilter:-->error: The matrix 'X' should contain at least two variables.\n")
	if(is.null(colnames(X))) colnames(X) = paste("X", 1:ncol(X), sep = "")
	colnames(X) = make.names(colnames(X))
	postpad = tolower(postpad[1])
	if(any(colnames(Bcoef)=="const")){
		constant = TRUE
		ic = 1
	} else{
		constant = FALSE
		ic = 0
	}
	obs = dim(X)[1]
	K = dim(X)[2]
	xsample = obs - p
	Xlags = embed(X, dimension = p + 1)[, -(1:K)]
	temp1 = NULL
	for (i in 1:p) {
		temp = paste(colnames(X), ".l", i, sep = "")
		temp1 = c(temp1, temp)
	}
	colnames(Xlags) = temp1
	Xend = X[-c(1:p), ]
	if(constant){
		rhs = cbind( Xlags, rep(1, xsample))
		colnames(rhs) <- c(colnames(Xlags), "const")
	} else{
		rhs = Xlags
		colnames(rhs) <- colnames(Xlags)
	}
	if( !(is.null(exogen)) ) {
		exogen = as.matrix(exogen)
		if (!identical(nrow(exogen), nrow(X))) {
			stop("\nvarxfit:-->error: Different row size of X and exogen.\n")
		}
		XK = dim(exogen)[2]
		if (is.null(colnames(exogen))) colnames(exogen) = paste("exo", 1:ncol(exogen), sep = "")
		colnames(exogen) = make.names(colnames(exogen))
		tmp  = colnames(rhs)
		rhs =  cbind(rhs, exogen[-c(1:p), ])
		colnames(rhs) = c(tmp, colnames(exogen))
	} else{
		XK = 0
	}
	datamat = as.matrix(rhs)
	colnames(datamat) = colnames(rhs)
	xfitted = t( Bcoef %*% t( datamat ) )
	xresiduals = tail(X, obs - p) - xfitted
	if(postpad!="none"){
		if(postpad == "constant"){
			# pre-pad values with the constant
			xfitted = t( Bcoef %*% t( rbind(matrix(c(rep(0, p*K), if(constant) 1 else NULL, if(XK>0) rep(0, XK) else NULL), nrow = p, ncol=dim(Bcoef)[2], byrow = TRUE), datamat ) ) )
			xresiduals = X - xfitted
		} else if(postpad == "zero"){
			xfitted = t( Bcoef %*% t( rbind(matrix(rep(0, dim(Bcoef)[2]), nrow = p, ncol=dim(Bcoef)[2], byrow = TRUE), datamat ) ) )
			xresiduals = X - xfitted
		} else if(postpad == "NA"){
			xfitted = t( Bcoef %*% t( rbind(matrix(rep(NA, dim(Bcoef)[2]), nrow = p, ncol=dim(Bcoef)[2], byrow = TRUE), datamat ) ) )
			xresiduals = X - xfitted
		} else{
			# do nothing
			xfitted = t( Bcoef %*% t( datamat ) )
			xresiduals = tail(X, obs - p) - xfitted
		}
	}
	ans = list( Bcoef = Bcoef, xfitted = xfitted, xresiduals = xresiduals, lag = p, constant = constant)
	return( ans )
}

# must return an array (n.ahead, n.assets, n.roll)
varxforecast = function(X, Bcoef, p, out.sample, n.ahead, n.roll, mregfor)
{
	X = as.matrix(X)
	X.orig = X
	n.roll = n.roll + 1
	m = ncol(X)
	n = nrow(X)
	X = X.orig[1:(n - out.sample), ,drop = F]
	n = nrow(X)
	if(any(colnames(Bcoef)=="const")){
		constant = TRUE
		ic = 1
	} else{
		constant = FALSE
		ic = 0
	}
	if(!is.null(mregfor)){
		mregfor = as.matrix(mregfor)
		mxn = ncol(mregfor)
	} else{
		mxn = 0
	}
	meanforc = array(NA, dim = c(n.ahead, m, n.roll))
	if(n.roll == 1){
		dmat1 = NULL		
		for(i in 0:(p-1))
		{
			dmat1 = cbind(dmat1, .lagx(X, n.lag = i, pad = 0) )
		}
		# for n.ahead > 1 and n.roll == 1
		dmat1 = as.numeric( tail(dmat1, 1) )
		for(i in 1:n.ahead){
			Z = cbind(matrix(dmat1[1 : (m * p)], nrow = 1), if(constant) matrix(1, ncol = 1, nrow = 1) else NULL, if(!is.null(mregfor)) matrix(as.numeric(mregfor[i, ]), nrow = 1 ) else NULL)
			meanforc[i, ,1] = Bcoef %*% t(Z)
			dmat1 = c(meanforc[i, ,1], dmat1)
		}
	} else{
		dmat2 = NULL
		N = dim(X.orig)[1] - out.sample
		for(i in 0:(p-1))
		{
			dmat2 = cbind(dmat2, .lagx(X.orig, n.lag = i, pad = 0) )
		}
		# n.roll
		dmat1 = dmat = as.numeric( tail(dmat2[1:N, ], 1) )
		for(i in 1:n.roll){
			for(j in 1:n.ahead){
				Z = cbind(matrix(dmat1[1 : (m * p)], nrow = 1), if(constant) matrix(1, ncol = 1, nrow = 1) else NULL, if(!is.null(mregfor)) matrix(as.numeric(mregfor[j+(i-1), ]), nrow = 1 ) else NULL)
				meanforc[j, ,i] = Bcoef %*% t(Z)
				dmat1 = c(meanforc[j, ,i], dmat1)
			}
			# move it one space along
			if( i < n.roll){
				dmat1 = c(dmat2[N+i, ], dmat)
				dmat = as.numeric( dmat2[N+i, ] )
			}
		}
	}
	return( meanforc )
}


varxsim = function(X, Bcoef, p, n.sim, n.start, prereturns, resids, mexsimdata)
{
	X = as.matrix(X)
	n = dim(X)[1]
	m = dim(X)[2]
	ns = n.start + n.sim
	if(any(colnames(Bcoef)=="const")){
		constant = TRUE
		ic = 1
	} else{
		constant = FALSE
		ic = 0
	}
	if(!is.null(mexsimdata)){
		mexsimdata = as.matrix(mexsimdata)
		mxn = ncol(mexsimdata)
	} else{
		mxn = 0
		Bcoef = Bcoef[,1:(m*p+1)]
	}
	dmat2 = NULL
	if( is.null(prereturns) ){
		for(i in 0:(p-1))
		{
			dmat2 = cbind(dmat2, .lagx(X, n.lag = i, pad = 0) )
		}
	} else{
		for(i in 0:(p-1))
		{
			dmat2 = cbind(dmat2, .lagx(prereturns, n.lag = i, pad = 0) )
		}
	}
	
	if( is.null(resids) ){
		resids = matrix(0, ncol = m, nrow = ns)
	}
	dmat = tail(dmat2, 1)
	meansim = matrix(NA, ncol = m, nrow = ns)
	for(i in 1:ns){
		Z = cbind(matrix(dmat[1 : (m * p)], nrow = 1), if(constant) matrix(1, ncol = 1, nrow = 1) else NULL, if(!is.null(mexsimdata)) matrix(as.numeric(mexsimdata[i, ]), nrow = 1 ) else NULL)
		meansim[i, ] = Bcoef %*% t(Z) + resids[i, ]
		if( i < ns ) dmat = c(meansim[i, ], dmat)
	}
	meansim = tail(meansim, n.sim)
	rownames(meansim) = 1:n.sim
	return( meansim )
}

varxsimXX = function(X, Bcoef, p, m.sim, prereturns, resids, mexsimdata)
{
	X = as.matrix(X)
	n = dim(X)[1]
	m = dim(X)[2]
	if(any(colnames(Bcoef)=="const")){
		constant = TRUE
		ic = 1
	} else{
		constant = FALSE
		ic = 0
	}
	if(!is.null(mexsimdata)){
		mexsimdata = as.matrix(mexsimdata)
		mxn = ncol(mexsimdata)
		# quick check to see if mxn = 1 (in which case it is inverted)
		if( mxn == m.sim){
			mexsimdata = t(mexsimdata)
			mxn = ncol(mexsimdata)
		}
		mxBcoef = matrix(Bcoef[,(m*p+ic+1):(m*p+ic+mxn)], ncol = mxn)
		Bcoef = Bcoef[,1:(m*p+ic)]
	} else{
		mxn = 0
		Bcoef = Bcoef[,1:(m*p+ic)]
	}
	dmat2 = NULL
	if( is.null(prereturns) ){
		for(i in 0:(p-1))
		{
			dmat2 = cbind(dmat2, .lagx(X, n.lag = i, pad = 0) )
		}
	} else{
		for(i in 0:(p-1))
		{
			dmat2 = cbind(dmat2, .lagx(prereturns, n.lag = i, pad = 0) )
		}
	}
	
	if( is.null(resids) ){
		resids = matrix(0, ncol = m, nrow = m.sim)
	}
	
	dmat = tail(dmat2, 1)
	#meansim = matrix(NA, ncol = m, nrow = m.sim)
	
	Z = cbind(matrix(dmat[1 : (m * p)], nrow = 1), if(constant) matrix(1, ncol = 1, nrow = 1) else NULL )
	meansim = matrix(Bcoef %*% t(Z), ncol = m, nrow = m.sim, byrow = TRUE)
	if(!is.null(mexsimdata)){
		for(i in 1:m.sim){
			meansim[i,] = meansim[i,] + matrix(mxBcoef %*% mexsimdata[i, ], ncol = m, byrow = TRUE)
		}
	}
	meansim = meansim + resids	
	return( meansim )
}
## VAR lag selection copied and amended from the 'vars' package
.varxselect = function (y, lag.max = 10, exogen = NULL) 
{
	y = as.matrix(y)
	if (any(is.na(y))) stop("\nNAs in Data.\n")
	colnames(y) = make.names(colnames(y))
	K = ncol(y)
	lag.max = abs(as.integer(lag.max))
	lag = abs(as.integer(lag.max + 1))
	ylagged = embed(y, lag)[, -c(1:K)]
	yendog = y[-c(1:lag.max), ]
	sample = nrow(ylagged)
	rhs = rep(1, sample)
	if (!(is.null(exogen))) {
		exogen = as.matrix(exogen)
		if (!identical(nrow(exogen), nrow(y))) { stop("\nDifferent row size of Data and exogenous regressors.\n") }
		if (is.null(colnames(exogen))) {
			colnames(exogen) = paste("exo", 1:ncol(exogen), 
					sep = "")
		}
		colnames(exogen) = make.names(colnames(exogen))
		rhs = cbind(rhs, exogen[-c(1:lag.max), ])
	}
	idx = seq(K, K * lag.max, K)
	detint = ncol(as.matrix(rhs))
	criteria = matrix(NA, nrow = 4, ncol = lag.max)
	rownames(criteria) = c("AIC", "HQ", "SC", "FPE")
	colnames(criteria) = paste(seq(1:lag.max))
	for (i in 1:lag.max) {
		ys.lagged = cbind(ylagged[, c(1:idx[i])], rhs)
		sampletot = nrow(y)
		nstar = ncol(ys.lagged)
		resids = resid(lm(yendog ~ -1 + ys.lagged))
		sigma.det = det(crossprod(resids)/sample)
		criteria[1, i] = log(sigma.det) + (2/sample) * (i * 
					K^2 + K * detint)
		criteria[2, i] = log(sigma.det) + (2 * log(log(sample))/sample) * 
				(i * K^2 + K * detint)
		criteria[3, i] = log(sigma.det) + (log(sample)/sample) * 
				(i * K^2 + K * detint)
		criteria[4, i] = ((sample + nstar)/(sample - nstar))^K * 
				sigma.det
	}
	order = apply(criteria, 1, which.min)
	return(list(selection = order, criteria = criteria))
}

.varcovres = function(X, mxn, p, resids, ic=1){
	n = dim(X)[1]
	m = dim(X)[2]
	# series x lags + constant  + exogenous
	nk = m * p + ic + mxn
	cov(resids) * (n - 1) / (n - nk)
}



################################################################################
# Multivariate Regression
# Input:
#   x: data-matrix (n,p)
#   y: data-matrix (n,q)
#   gamma: proportion of trimming 
#   arguments: 
#     - ns : contains number of subsets; default=5000
#     - nc : number of C-steps; default=10
#     - delta : critical value for Reweighted estimator, deafult=0.01
# Output:
#     beta : matrix (p,q) of MLTS-regression coefficients
#     sigma: matrix (q,q) containing MLTS-residual covariance
#     dres : residual distances (n,1) w.r.t. initial fit
#     betaR : matrix (p,q) of RMLTS-regression coefficients
#     sigmaR: matrix (q,q) containing RMLTS-residual covariance
#     dresR : residual distances (n,1) w.r.t. RMLTS
# Remark: if intercept needed, add a column of ones to the x-matrix
#
# Ref: Agullo,J., Croux, C., and Van Aelst, S. (2008) 
#      The Multivariate Least Trimmed Squares Estimator, 
#      Journal of multivariate analysis
#
# Author: Kristel Joossens
# 
#
###########
# EXAMPLE #
###########
# n=1000;
# p=5;
# q=2;
# x = cbind(rep(1,n),array(rnorm(n*(p-1)),dim=c(n,p-1)))
# beta = cbind(c(2,1,4,2,1),c(5,2,1,.5,2))
# y = x %*% beta +array(rnorm(n*q),dim=c(n,q))
# gamma = 0.25
# out = mlts(x,y,gamma)
################################################################################
mlts = function(x, y, gamma, ns = 500, nc = 10, delta = 0.01)
{ 
	d = dim(x)
	n = d[1]
	p = d[2]
	q = ncol(y) 
	h = floor(n*(1-gamma))+1
	obj0 = 1e10 
	for (i in 1:ns)
	{ 
		sorted = sort(runif(n), na.last = NA, index.return = TRUE)
		istart = sorted$ix[1:(p+q)]
		xstart = x[istart,]
		ystart = y[istart,]
		bstart = solve(t(xstart)%*%xstart,t(xstart)%*%ystart) 
		sigmastart = (t(ystart-xstart%*%bstart))%*%(ystart-xstart%*%bstart)/q
		for (j in 1:nc)
		{
			res  =  y - x %*% bstart
			tres = t(res)
			dist2 = colMeans(solve(sigmastart,tres)*tres)
			sdist2 = sort(dist2,na.last = NA,index.return = TRUE)
			idist2 = sdist2$ix[1:h]
			xstart = x[idist2,]
			ystart = y[idist2,]
			bstart = solve(t(xstart)%*%xstart,t(xstart)%*%ystart)
			sigmastart = (t(ystart-xstart%*%bstart))%*%(ystart-xstart%*%bstart)/(h-p)
		}
		obj = det(sigmastart)
		if(obj < obj0)
		{ 
			result.beta = bstart
			result.sigma = sigmastart
			obj0 = obj
		}
	}
	cgamma = (1-gamma)/pchisq(qchisq(1-gamma,q),q+2)
	result.sigma = cgamma * result.sigma
	res = y - x %*% result.beta
	tres = t(res)
	result.dres = colSums(solve(result.sigma,tres)*tres)
	result.dres = sqrt(result.dres)
	
	qdelta = sqrt(qchisq(1-delta,q))
	good  = (result.dres <= qdelta)
	xgood = x[good,]
	ygood = y[good,]
	result.betaR = solve(t(xgood)%*%xgood,t(xgood)%*%ygood)
	result.sigmaR = (t(ygood-xgood%*%result.betaR)) %*% (ygood-xgood%*%result.betaR)/(sum(good)-p)
	cdelta = (1-delta)/pchisq(qdelta^2,q+2)
	result.sigmaR = cdelta*result.sigmaR
	resR = y - x%*%result.betaR
	tresR = t(resR)
	result.dresR = colSums(solve(result.sigmaR,tresR)*tresR)
	result.dresR = sqrt(result.dresR)
	return( list( beta = result.beta, sigma = result.sigma, dres = result.dres,
			betaR = result.betaR, sigmaR = result.sigmaR, dresR = result.dresR ) )
}

mlts.parallel = function(x, y, gamma, ns = 500, nc = 10, delta = 0.01, cluster = NULL)
{
	d = dim(x)
	n = d[1]
	p = d[2]
	q = ncol(y) 
	h = floor(n*(1-gamma))+1
	obj0 = 1e10
	rseed = as.integer(runif(ns, 1, ns*ns*100))	
	.fn = function(x, y, nc, gamma, delta, rseed){
		result = list()
		n = NROW(x)
		p = NCOL(x)
		q = NCOL(y)
		set.seed(rseed)
		h = floor(n*(1-gamma))+1
		sorted = sort(runif(n), na.last = NA, index.return = TRUE)
		istart = sorted$ix[1:(p+q)]
		xstart = x[istart,]
		ystart = y[istart,]
		bstart = solve(t(xstart)%*%xstart,t(xstart)%*%ystart) 
		sigmastart = (t(ystart-xstart%*%bstart))%*%(ystart-xstart%*%bstart)/q
		for (j in 1:nc)
		{
			res  =  y - x %*% bstart
			tres = t(res)
			dist2 = colMeans(solve(sigmastart,tres)*tres)
			sdist2 = sort(dist2,na.last = NA,index.return = TRUE)
			idist2 = sdist2$ix[1:h]
			xstart = x[idist2,]
			ystart = y[idist2,]
			bstart = solve(t(xstart)%*%xstart,t(xstart)%*%ystart)
			sigmastart = (t(ystart-xstart%*%bstart))%*%(ystart-xstart%*%bstart)/(h-p)
		}
		result$obj = det(sigmastart)
		result$beta = bstart
		result$sigma = sigmastart
		return(result)
	}
	clusterExport(cluster, c(".fn", "y", "x", "rseed", "gamma", "delta", 
					"nc"), envir = environment())
	sol = parLapply(cluster, as.list(1:ns), fun = function(i){
				.fn(x, y, nc, gamma, delta, rseed[i])
			})
	
	d = sapply(sol, FUN = function(x) x$obj)
	idx = which(d==min(d))[1]
	result = sol[[idx]]
	cgamma = (1-gamma)/pchisq(qchisq(1-gamma,q),q+2)
	result$sigma = cgamma * result$sigma
	res = y - x %*% result$beta
	tres = t(res)
	result$dres = colSums(solve(result$sigma,tres)*tres)
	result$dres = sqrt(result$dres)
	qdelta = sqrt(qchisq(1-delta,q))
	good  = (result$dres <= qdelta)
	xgood = x[good,]
	ygood = y[good,]
	result$betaR = solve(t(xgood)%*%xgood,t(xgood)%*%ygood)
	result$sigmaR = (t(ygood-xgood%*%result$betaR)) %*% (ygood-xgood%*%result$betaR)/(sum(good)-p)
	cdelta = (1-delta)/pchisq(qdelta^2,q+2)
	result$sigmaR = cdelta*result$sigmaR
	resR = y - x%*%result$betaR
	tresR = t(resR)
	result$dresR = colSums(solve(result$sigmaR,tresR)*tresR)
	result$dresR = sqrt(result$dresR)
	return( list( beta = result$beta, sigma = result$sigma, dres = result$dres,
					betaR = result$betaR, sigmaR = result$sigmaR, dresR = result$dresR ) )
}

################################################################################
# Robust estimation of a VAR model using the robust MLST estimator.
# Robust lag length criteria for the give lag are reported.
# The user can as for the robust impulse response functions and their 
# analytic confidence bounds.
#
# Input arguments:
#   data  : multivariate time series (T x nvar) (nvar=number of variables)
#   lags  : order of the VAR model
# The three possible extra arguments can be
#   alpha : trimming portion to obtain the MLTS estimator 
#           the default value is 0.25
#   delta : trimming portion to obtain the RMLTS estimator (based on MLTS)
#           the default value is 0.01. If you want to redefine delta, you
#           have to redefine alpha as previous argument
#   nstep : number of steps: If are only interested in the effect of a unit
#           shock at an innovation to observations not far/far in the
#           future. You have to take nstep small(nstep=20)/large(nsetp=200)
# The combinations for the last 3 input arguments might be as follows
# nothing - alpha - alpha,delta - alpha,delta,nstep - 
#           nstep - nstep,alpha - nstep,alpha,delta
# eg. robustvar(data,lags,nstep,alpha)

# Output arguments:
#   beta : the estimated coefficients
#   AIC  : The Akaike criterion
#   HQ   : The Hannan Quinn criterion
#   SC   : The Schwarz criterion
# If nstep is not given, no IRF will be created.
# If nstep is given: figure 1 represent the  IRFs and their Analytic CB for
# VAR(lags)-model using OLS.
#   Z    : the estimated IRF
#   UB   : The 95# upper bound for the IRF
#   LB   : The 95# lower bound for the IRF
# robdis : contains the robust Mahalanobis distances

robustvar = function(data, exogen = NULL, constant = TRUE, lags = 2, alpha = 0.01, ns = 500, nc = 10, delta = 0.01, cluster)
{
	T = dim(data)[1]
	nvar  = dim(data)[2]
	ydata = data[(lags+1):T,]
	if(constant){
		xdata = matrix(1, nrow = T - lags, ncol = 1)
		for(i in 1:lags){ xdata = cbind(xdata, data[((1+lags)-i):(T-i), ])}
	} else{
		xdata = NULL
		for(i in 1:lags){ xdata = cbind(xdata, data[((1+lags)-i):(T-i), ])}
	}
	if(!is.null(exogen)) xdata = cbind(xdata, tail(exogen, T - lags))
	if(!is.null(cluster)){
		datamlts = mlts.parallel(x = as.matrix(xdata), y = as.matrix(ydata), 
				gamma = alpha, ns = ns, nc = nc, delta = delta, cluster = cluster)
	} else{
		datamlts = mlts(as.matrix(xdata), as.matrix(ydata), gamma = alpha, ns = ns, 
				nc = nc, delta = delta)  
	}
	U = datamlts$sigmaR
	logdetSigma = log(det(U))
	cdelta = (1-delta)/pchisq(qchisq(1 - delta, nvar),nvar+2)
	m = sum(datamlts$nooutlier)
	phi = nvar*(lags*nvar+1)
	ll  = -(T-lags)/2 *(log(2*pi)*nvar+log(det(U)))- (m-nvar)*nvar/(2*cdelta)
	res  = datamlts
	res$AIC =-2/(T-lags)*(ll-phi)
	res$HQ =-2/(T-lags)*(ll-log(log(T-lags))*phi)
	res$SC =-2/(T-lags)*(ll-log(T-lags)*phi/2)
	return( res )
}

# VAR unconditional mean
.umeanvar = function(Bcoef, p){
	#Bcoef[,11]%*%t(solve(diag(10) - Bcoef[,-11]))
	m = dim(Bcoef)[1]
	Id = diag(m)
	C = Bcoef[, p * m]
	n = dim(Bcoef)[2]
	if(n > (m+1)){
		muEX = Bcoef[, ( p * m + 2 ):n, drop = FALSE] 
		muC  = Bcoef[, p * m + 1] + apply(muEX, 1, "sum")
	} else{
		muEX = NULL
		muC  = Bcoef[, p * m + 1]
	}
	idx = t(rugarch:::.embed(1:(m*p), k=m, by=m, ascending = TRUE))
	for(i in 1:p) Id = Id - C[,idx[,i]]
	muC %*% t(solve(Id - C))
}

##############################################################################################