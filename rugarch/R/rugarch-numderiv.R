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

# Fitting procedure required functions:
.hessian2sided = function(f, x, ...)
{
	n = length(x)
	fx = f(x, ...)
	eps = .Machine$double.eps
	
	# Compute the stepsize (h)
	# h = eps^(1/3)*apply(as.data.frame(x), 1,FUN = function(z) max(abs(z), 1e-4))
	h = apply(as.data.frame(x), 1,FUN = function(z) max(abs(z*eps^(1/3)), 1e-9))
	xh = x+h
	h = xh-x
	if(length(h) == 1) ee = matrix(h, ncol = 1, nrow = 1) else ee = as.matrix(diag(h))
	
	# Compute forward and backward steps
	gp = vector(mode = "numeric", length = n)
	gp = apply(ee, 2, FUN = function(z) f(x+z, ...))
	gm = vector(mode="numeric",length=n)
	gm = apply(ee, 2, FUN = function(z) f(x-z, ...))
	H = h%*%t(h)
	Hm = H
	Hp = H
	# Compute "double" forward and backward steps
	for(i in 1:n){
		for(j in  i:n){
			Hp[i,j] = f(x+ee[,i]+ee[,j], ...)
			Hp[j,i] = Hp[i,j]
			Hm[i,j] = f(x-ee[,i]-ee[,j], ...)
			Hm[j,i] = Hm[i,j]
		}
	}
	#Compute the hessian
	for(i in 1:n){
		for(j in  i:n){
			H[i,j] = ( (Hp[i,j]-gp[i]-gp[j]+fx+fx-gm[i]-gm[j]+Hm[i,j]) /H[i,j] )/2
			H[j,i] = H[i,j]
		}
	}
	return(H)
}

# keep a copy with same name as it's called by other packages
.hessian2sidedcpp = function(f, x, ...)
{
	n = length(x)
	fx = f(x, ...)
	eps = .Machine$double.eps
	
	# Compute the stepsize (h)
	# h = eps^(1/3)*apply(as.data.frame(x), 1,FUN = function(z) max(abs(z), 1e-4))
	h = apply(as.data.frame(x), 1,FUN = function(z) max(abs(z*eps^(1/3)), 1e-9))
	xh = x+h
	h = xh-x
	if(length(h) == 1) ee = matrix(h, ncol = 1, nrow = 1) else ee = as.matrix(diag(h))
	
	# Compute forward and backward steps
	gp = vector(mode = "numeric", length = n)
	gp = apply(ee, 2, FUN = function(z) f(x+z, ...))
	gm = vector(mode="numeric",length=n)
	gm = apply(ee, 2, FUN = function(z) f(x-z, ...))
	H = h%*%t(h)
	Hm = H
	Hp = H
	# Compute "double" forward and backward steps
	for(i in 1:n){
		for(j in  i:n){
			Hp[i,j] = f(x+ee[,i]+ee[,j], ...)
			Hp[j,i] = Hp[i,j]
			Hm[i,j] = f(x-ee[,i]-ee[,j], ...)
			Hm[j,i] = Hm[i,j]
		}
	}
	#Compute the hessian
	for(i in 1:n){
		for(j in  i:n){
			H[i,j] = ( (Hp[i,j]-gp[i]-gp[j]+fx+fx-gm[i]-gm[j]+Hm[i,j]) /H[i,j] )/2
			H[j,i] = H[i,j]
		}
	}
	return(H)
}

# makes no difference in timings
.hessian2sidedcppOLD = function(f, x, ...)
{
	n = length(x)
	arglist = list(...)$arglist
	fx = function(y){
		ans = f(y, arglist)
		ifelse(!is.numeric(ans) | !is.finite(ans) | length(ans)==0, 1e10, ans)
	}
	eps = .Machine$double.eps
	# Compute the stepsize (h)
	#h = eps^(1/3)*apply(as.data.frame(x), 1,FUN = function(z) max(abs(z), 1e-2))
	# h = eps^(1/3) * sapply(x, FUN = function(z) max(abs(z), 1e-4))
	h = apply(as.data.frame(x), 1,FUN = function(z) max(abs(z*eps^(1/3)), 1e-9))
	xh = x+h
	h = xh-x
	if(length(h) == 1) ee = matrix(h, ncol = 1, nrow = 1) else ee = as.matrix(diag(h))
	# Compute forward and backward steps
	gp = vector(mode = "numeric", length = n)
	gp = apply(ee, 2, FUN = function(z) f(x+z, ...))
	gm = vector(mode="numeric",length=n)
	gm = apply(ee, 2, FUN = function(z) f(x-z, ...))
	H = h%*%t(h)
	# Hm = H
	# Hp = H
	sol = .Call("hessian2sided", fun = as.function(fx), x = as.numeric(x), H = as.matrix(H), 
			deps = as.matrix(ee), gminus = as.numeric(gm), gplus = as.numeric(gp), PACKAGE = "rugarch")
	return(sol)
}

# The following functions are based on on Kevin Sheppard's MFE toolbox
#---------------------------------------------------------------------
neweywestcv = function(data, nlag = NULL, center = TRUE)
{
	# Long-run covariance estimation using Newey-West (Bartlett) weights
	#  if nlag empty=NULL NLAG=min(floor(1.2*T^(1/3)),T)
	N = dim(as.matrix(data))[1]
	if(is.null(nlag)) nlag=min(floor(1.2*N^(1/3)),N)
	if(center) data = apply(data, 2, FUN = function(x) scale(x, center = TRUE, scale = FALSE))
	# weights
	bw = (nlag+1-(seq(0,nlag,by=1)))/(nlag+1)
	cv = 1/N * t(data)%*%data
	for(i in 1:nlag){
		gmi = 1/N * (t(data[(i+1):N,])%*%data[1:(N-i),])
		gpp = gmi + t(gmi)
		cv = cv + bw[i+1]*gpp
	}
	return(cv)
}

robustvcv = function(fun, pars, nlag = 0, hess, n, ...)
{
	k = length(pars)
	#h = apply(as.data.frame(pars), 1, FUN = function(x) max(abs(x*eps^(1/3)), 1e-7))
	h = apply(as.data.frame(pars), 1,FUN = function(z) max(abs(z*eps^(1/3)), 1e-9))
	hplus =  pars+h
	hminus = pars-h
	hparsminus = hparsplus = matrix(pars, ncol = k, nrow = k, byrow = T)
	diag(hparsplus) = hplus
	diag(hparsminus) = hminus
	
	likelihoodsminus = likelihoodsplus = zeros(n, k)
	likelihoodsplus =  apply(hparsplus, 1,  FUN = function(x) fun(x, ...))
	likelihoodsminus = apply(hparsminus, 1, FUN = function(x) fun(x, ...))
	
	scores = zeros(n,k)
	likpm = likelihoodsplus-likelihoodsminus
	scores = likpm/(2*repmat(t(h), n, 1))
	A = hess/n
	hess = A
	Ainv = try( solve(A), silent = TRUE )
	if( inherits(Ainv, "try-error")){
		info = 1
		vcv = NA
	} else{
		if(nlag>0){
			B = neweywestcv(scores, nlag)
			vcv = (Ainv%*%B%*%Ainv)/n
		} else{
			B = cov(scores)
			vcv = (Ainv%*%B%*%Ainv)/n
		}
		info = 0
	}
	return(list(vcv = vcv, scores = scores, info = info))
}
