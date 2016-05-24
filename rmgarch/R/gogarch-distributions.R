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
# Convolution Algorithms
#------------------------------------------------------------------------------
.intpolold = function(x,y,z)
{
	# linear interp of z to f(x,y), needs min(x)< z < max(x)
	n = length(z)
	j = 1
	i = 1
	r = NULL
	while(i <= n)
	{
		flag = i
		while(flag == i)
		{
			if(z[i] >= x[j] && z[i] < x[j+1])
			{
				r = c( r, y[j+1] - (x[j+1] - z[i]) * (y[j+1] - y[j]) / (x[j+1] - x[j]) )
				i = i + 1
			} else {
				j = j + 1
			}
		}
	}
	return(r)
}

.intpol = function(x, y, z)
{
	n = length(z)
	ans = approx(x, y = y, xout = z, method = "linear", n = length(z),
			yleft = min(z), yright = max(z), rule = 1, f = 0, ties = mean)$y
	return(ans)
	
}

.interval = function(x, y)
{
	# find the x: min(x)<min(y) and max(x)> max(y)
	i1 = which(x < min(y))
	i2 = which(x > max(y))
	i = seq(i1[dim(as.matrix(i1))[2]], i2[1])
	j = i[1]
	while(x[j] > y[1]){ j = j - 1 }
	k = max(i)
	while(x[k] < max(y)){ k = k + 1 }	
	i = c(if(i[1] - 1 != 0) j:(i[1] - 1) else NULL, i, if((max(i) + 1) < k) (max(i) + 1):k else NULL)
	i = sort(unique(i))
	return(i)
}

.wintpl = function(x, y, z, w = NULL)
{
	# interpolate with a window of size w on z.
	if(is.null(w)) w = length(z)
	m = length(z)
	n = floor(m/w)
	r = NULL
	for(i in 1:n)
	{
		dz = z[((i - 1) * w + 1):(i * w)]
		k = .interval(x, dz)
		zt = approx(x = x[k], y = y[k], xout = dz, method = "linear", n = length(dz),
				yleft = min(dz), yright = max(dz), rule = 1, f = 0, ties = mean)$y
		r = c(r, zt)
	}
	if(m%%w != 0)
	{
		e = m %% w 
		g = (n * w + 1):(n * w + e)
		dz = z[g]
		k = .interval(x, dz)
		zt = approx(x = x[k], y = y[k], xout = dz, method = "linear", n = length(dz),
				yleft = min(dz), yright = max(dz), rule = 1, f = 0, ties = mean)$y
		#.intpol(x[k], y[k], dz)
		r = c(r, zt)
	}
	return(r)
}

cfinv = function(z, f, step = 0.01, ...)
{
	pmax = 18
	p = 14 
	maxz = round(max(abs(z))) + 5
	while((maxz/step+1) > 2^(p-1)) { p = p + 1 }
	if(p > pmax) p = pmax
	if((maxz/step+1) > 2^(p - 1)) { step = (maxz + 1) * (1 + step/10)/(2^(p - 1)) }
	zs = sort(z)
	# run the fft
	n = 2^p
	x = seq(0, n - 1, by = 1) * step - (n * step/2)
	s = 1/(step * n)
	tt = 2 * pi * s * (seq(0, n - 1, by = 1) - n/2)
	sgn = rep(1, n)
	ds = seq(2, n, by = 2)
	sgn[ds] = -1 * rep(1, n/2)
	cf = f(tt, ...)
	phi = sgn * cf
	phi[n/2 + 1] = sgn[n/2 + 1]
	#zplan = planFFT(length(phi))
	#p = s * abs(FFT(phi, plan = zplan))
	# 9-June-2014 remove FFTW since maintainers unable to fix mavericks problem
	p = s * abs(fft(phi))
	pdf = .wintpl(x, p, zs, w = NULL)
	return(pdf)
}

# characteristic function of nig with independent margins
nigmvcf = function(z, alpha, beta, delta, mu)
{
	N = length(z)
	m = length(mu)
	x1 = 1i * z * sum(mu)
	zx = matrix(0, ncol = m, nrow = N)
	zx = apply(cbind(delta, alpha, beta), 1, FUN = function(x) x[1]*(sqrt(x[2]^2 - x[3]^2) - sqrt(x[2]^2 - (x[3] + 1i * z)^2)))
	x2 = apply(t(zx), 2, FUN = function(x) sum(x))
	ans = x1 + x2
	return(exp(ans))
}

# we use the function BesselK of the Bessel package of Martin Maechler for evaluating
# complex inputs to the bessel functions.
ghypmvcf = function(z, lambda, alpha, beta, delta, mu)
{
	N = length(z)
	m = length(mu)
	x0 = 1i * z * sum(mu)
	zx = matrix(0, ncol = m, nrow = N)
	zx = apply(cbind(lambda, alpha, beta, delta), 1, FUN = function(x) .ghypfn(x[1], x[2], x[3], x[4], z))
	x = apply(t(zx), 2, FUN = function(x) sum(x))
	ans = exp(x0 + x)
	return(ans)
}

.ghypfn = function(lambda, alpha, beta, delta, z)
{
	x1 = (lambda/2)*(log(alpha^2 - beta^2) - log(alpha^2 - (beta+1i*z)^2) )
	x2 = log(BesselK(delta * sqrt(alpha^2 - (beta + 1i*z)^2), abs(lambda))) - log(BesselK(delta * sqrt(alpha^2 - beta^2), abs(lambda)))
	x1 + x2
}


##############################################################################

#------------------------------------------------------------------------------
# .makeDNew, .makePNew, .makeQNew from the distr package adapted here for
# the convolution algorithm
#------------------------------------------------------------------------------
.makeDNew <- function(x, dx, h = NULL, standM = "sum"){
	dx <- (dx >= .Machine$double.eps)*dx
	if( length(dx) < length(x) ) dx <- c(0,dx)
	if (is.null(h)) h <- 1
	dx1 <- dx / h
	mfun <- approxfun
	## density
	df1 <- mfun(x = x, y = dx1, yleft = 0, yright = 0)
	if (standM == "sum")
		stand <- sum(dx)
	else   {
		stand <- try(integrate(df1, -Inf, Inf)$value, TRUE)
		if (is(stand,"try-error")){
			warning("'integrate()' threw an error ---result may be inaccurate.")
			stand <- sum(df1(x))*h*(x[2]-x[1])
		}
	}
	dfun <- function(x, log = FALSE)
	{ if (log){
			d0<-log(df1(x))-log(stand)
		} else{
			d0 <- df1(x) / stand
		}
		return (d0)}
	
	rm(x,dx1,h)
	return(dfun)
}

.primefun <- function(f,x, nm = NULL){
	
	h <- diff(x)
	l <- length(x)
	
	xm <- (x[-l]+x[-1])/2
	
	fxm <- f(xm)
	fx <- f(x)
	
	
	fxs  <- 2 * cumsum(fx) - fx - fx[1]
	fxsm <- 4 * cumsum(fxm)
	
	fxx <- c(0, (fxs[-1]+fxsm)* h / 6 )
	
	if (is.null(nm)) nm <- fxx[l]
	
	fx1 <- approxfun(x, fxx, yright = nm, yleft = 0)
	
	ffx <- function(u){
		ffy <- fx1(u) 
		ffy[u > max(x)] <- nm 
		ffy[u < min(x)] <- 0
		return(ffy)
	}
	
	return(ffx)
}

.csimpsum <- function(fx){
	l <- length(fx)
	l2 <- l%/%2
	if (l%%2 == 0) {
		fx <- c(fx[1:l2],(fx[l2]+fx[l2+1])/2,fx[(l2+1):l])
		l <- l+1}
	f.even <- fx[seq(l) %% 2 == 0]
	f.odd  <- fx[seq(l) %% 2 == 1]
	fs    <- 2 * cumsum(f.odd) - f.odd - f.odd[1]
	fsm   <- 4 * cumsum(f.even)
	ff <- c(0,(fs[2:(l2+1)]+fsm)/3 )
	ff
}

.makePd <- function(x,y, yleft, yright){
	stepfun(x = x, y = c(yleft, y))
}

.makeQd <- function(x,y, yleft, yright){
	force(y)
	force(x)
	f <- function(u) {
		q0 <- sapply(u, 
				function(z) y[min(sum(x < z-.Machine$double.eps) + 1,
									length(y)) ] )
		q0[.isEqual(u,0)] <- yleft
		q0[.isEqual(u,1)] <- yright
		return(q0)}
	return(f)
}
.makePNew <- function(x, dx, h = NULL, notwithLLarg = FALSE,
		Cont = TRUE, myPf = NULL, pxl = NULL, pxu = NULL){
	
	if (is.null (h)) h <- 0
	
	x.u <- x.l <- x
	if (Cont){
		mfun <- if (is.null (myPf)) approxfun else myPf
		l <- length(x)
		if ((l%%2==0)&& is.null(myPf)){
			l2 <- l/2
			if (is.null(pxl))
				x.l <- c(x[1:l2],(x[l2]+x[l2+1])/2,x[(l2+1):l])
			if (is.null(pxu))
				x.u <- c(x[1:l2],(x[l2]+x[l2+1])/2,x[(l2+1):l])
			l <- l+1
		}
		cfun <- .csimpsum
		if (is.null(pxl)&& is.null(myPf))
			x.l <- x.l[seq(l)%%2==1]
		if (is.null(pxu)&& is.null(myPf))
			x.u <- x.u[seq(l)%%2==1]
	}else    {
		mfun <- .makePd
		cfun <- cumsum
	}       
	
	p.l <- if(!is.null(pxl)) pxl else cfun(dx)
	
	nm <- max(p.l)
	p1.l <- mfun(x = x.l, y = p.l, yleft = 0, yright = nm)
	nm <- p1.l(max(x))
	if(notwithLLarg){
		ifElsePS <- substitute(if (lower.tail) p1.l(q) else 1 - p1.l(q))
	}else{
		p.u <- if(!is.null(pxu)) pxu else rev(cfun(rev(dx)))
		## continuity correction by h/2
		if (!Cont) p.u <- c(p.u[-1],0)
		p1.u <- mfun(x = x.u, y = p.u, yright = 0, yleft = nm)
		rm(p.u)
		ifElsePS <- substitute(if (lower.tail) p1.l(q) else p1.u(q))
	}
	pfun <- function(q, lower.tail = TRUE, log.p = FALSE){}
	body(pfun) <- substitute(
			{ p0 <- ifElsePC
				p0 <- if (log.p) log(p0)-log(nm) else p0/nm
				return(p0)
			}, list(ifElsePC = ifElsePS))
	rm(dx, p.l, notwithLLarg)
	return(pfun)
}
.isEqual <- function(p0, p1, tol = min(1e-05/2, .Machine$double.eps^.7)) abs(p0-p1)< tol

.isIn <- function(p0, pmat, tol = min(1e-05/2,.Machine$double.eps^.7))
{list1 <- lapply(1:nrow(pmat), function(x){ 
				(p0+tol > pmat[x,1]) & (p0-tol < pmat[x,2]) })
	apply(matrix(unlist(list1), ncol = nrow(pmat)), 1, any)}

.isEqual01<- function(x) .isEqual(x,0)|.isEqual(x,1)

.setEqual <- function(x, y, tol = 1e-7){
	### all elements of x equal to some element of y up tol are set to exactly
	###     the respective element of y
	x1 <- round(2*x/tol,0)
	y1 <- round(2*y/tol,0)
	z  <- x
	m  <- match(x1,y1)
	n.ina.m <- !is.na(m)
	z[n.ina.m] <- y[m[n.ina.m]]
	z
}
.makeQc <- function(x,y, yleft, yright){
	yl <- if(is.finite(yleft)) yleft  else y[1]
	yr <- if(is.finite(yright)) yright else y[length(y)]
	f1 <- approxfun(x = x, y = y, yleft = yl, yright = yr)
	f <- function(x) 
	{y1 <- f1(x)
		y1[.isEqual(x,0)] <- yleft
		y1[.isEqual(x,1)] <- yright
		return(y1)
	}
	return(f)
}

.makeQNew <- function(x, dx, h)
{
	pdnew=.makePNew(x, dx, h)
	notwithLLarg = FALSE
	yL = -Inf
	yR =  Inf
	px.l <- pdnew(x + 0.5*h)
	px.u <- pdnew(x + 0.5*h, lower.tail = FALSE)
	
	Cont=TRUE
	mfun <- .makeQc
	ix <- .isEqual01(px.l)
	if(!is.finite(yR)||Cont)
	{xx <- px.l[!ix]; yy <- x[!ix]}
	else  
	{xx <- px.l; yy <- x}
	q.l <- mfun(x = xx, y = yy, yleft = yL, yright = yR)
	rm(xx,yy)
	if(notwithLLarg){
		ifElseQS <- quote(if (lower.tail) q.l(p01) else q.l(1-p01))
	}else{
		#         px.u <- rev(px.u);
		x <- rev(x)
		if (Cont) px.u <- rev(px.u)
		ix <- .isEqual01(px.u)
		xx <- px.u[!ix]
		yy <- if (Cont) x[!ix] else x[rev(!ix)]
		q.u <- mfun(x = xx, y = yy, yleft = yR, yright = yL)
		rm(xx,yy)
		ifElseQS <- quote(if (lower.tail) q.l(p01) else q.u(p01))
	}
	qfun <- function(p, lower.tail = TRUE, log.p = FALSE){}
	body(qfun) <- substitute({
				if (log.p) p <- exp(p)
				if (any((p < -.Machine$double.eps)|(p > 1+.Machine$double.eps)))
					warning(gettextf("q method of %s produced NaN's ", objN))
				i01 <- (-.Machine$double.eps<=p)&(p<=1+.Machine$double.eps)
				p01 <- p[i01] ## only values in [0,1] are used
				q0  <- p*0
				q0[!i01] <- NaN
				q0[ i01] <- ifElseQC
				return(as.numeric(q0))
			}, list(ifElseQC = ifElseQS, objN = quote(.getObjName(1))))
	return(qfun)
}


.coskew.ind = function(skewval)
{
	n = length(skewval)
	cs = matrix(0, ncol = n^2, nrow = n)
	# diagonal tensor type indices using column type count
	indx = ((1:n)-1) * n^2 + ((1:n)-1) * n + (1:n)
	cs[indx] = skewval
	return(cs)
}

# N (pairs) = 2*(n^2)+n*(n-3)
# unique 2-pair combinations : factorial(n)/(factorial(n-2)*factorial(2))
# for each pair:
# Y = combn(1:n, 2)
# Tabulation:
# No. of unique pairs whose difference is 1: n-1 [{2,1},{3,2},...,{n,n-1}]
# No. of unique pairs whose difference is 2: n-2 [{3,1},{4,2},...,{n,n-2}]
# No. of unique pairs whose difference is n-1: 1 [{n,1}]


.cokurt.pairs = function(n)
{
	# Author: Alexios Ghalanos 
	# Copyright (C) 2013 (part of the rmgarch package)
	# Finds the locations of each i=j, k=l pair in an N x N^3 matrix representing 
	# the collapsed co-kurtosis tensor
	# As an example, the N=4 matrix below illustrates the indexing and the position of the pairs:
	# The code which follows generates the 1-dimensional index (based on columnwise counting) of
	# each pair and also returns the associates pair (for knowing which sigma values to multiply
	# to obtain the fourth co-moment).
	##
	## j:    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
	## k:    1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4     
	## l:    1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4     
	## i: 1  . . . . . x . . . . x . . . . x . x . . x . . . . . . . . . . . . . x . . . . . x . . . . . . . . . . x . . . . . . . . x . . .  
	## i: 2  . x . . x . . . . . . . . . . . x . . . . . . . . . x . . . . x . . . . . . x . . x . . . . . . . . . . . . . x . . . . . x . .
	## 1: 3  . . x . . . . . x . . . . . . . . . . . . . x . . x . . . . . . x . . . . x . . . . . . . . . x . . . . . . . . . . . x . . x .
	## i: 4  . . . . . . . . . . . . x . . . . . . . . . . x . . . . . . . . . . . . . . . . . . . x . . x . x . . . . . . . . . x . . . . .
	
	unique_pairs = factorial(n)/(factorial(n-2)*factorial(2))
	idx = vector(mode="list", length = unique_pairs)
	prs = t(combn(n, 2))
	prs = cbind(prs[,2], prs[,1], prs[,2] - prs[,1])
	prs = prs[order(prs[,1]-prs[,2]), ]
	colnames(prs) = c("ij","kl","d")
	ix = n+2
	# first pair is always : {2,1}
	idx[[1]] = cumsum(c(ix, (n-1)*n, (n-1), (n+1)*(n-1)^2, (n-1), (n-1)*n))
	for(i in 2:unique_pairs){
		d = prs[i,3]
		if(d == prs[i-1,3]){
			i1 = idx[[i-1]][1]+n^3+n^2+n+1
		} else{
			i1 = ix+n+1
			ix = i1
		}
		idx[[i]] = cumsum(c(i1, d*(n-1)*n, d*(n-1), d*(n+1)*(n-1)^2, d*(n-1), d*(n-1)*n))
	}
	return(list(index = idx, pairs = prs))
}

.cokurt.ind = function(sigma2, kurtvals)
{
	n = length(sigma2)
	x = .cokurt.pairs(n)
	ix = ((1:n)-1) * n^3  + ((1:n)-1) * n^2 + ((1:n)-1) * n + (1:n)
	z = rep(0, n*n^3)
	z[ix] = kurtvals
	for(i in 1:nrow(x$pairs)) z[x$index[[i]]] = sigma2[x$pairs[i,1]]*sigma2[x$pairs[i,2]]
	cs = matrix(z, ncol = n^3, nrow = n)
	return(cs)
}



.coskew.sigma = function(sigmas)
{
	n = length(sigmas)
	idx1 = sort(rep(1:n, n^2))
	idx2 = rep(sort(rep(1:n, n)), n)
	idx3 = rep(1:n, n^2)
	idx = cbind(idx1, idx2, idx3)
	zs = as.numeric(.Call("gogarchcssigma", idx = as.matrix(idx-1), SS = as.numeric(sigmas),
					PACKAGE = "rmgarch"))
	# apply(idx, 1, FUN = function(x) prod(sigmas[x]))
	cs = matrix(zs, ncol = n^2, nrow = n, byrow = TRUE)
	return(cs)
}

.cokurt.sigma = function(sigmas)
{
	n = length(sigmas)
	idx1 = sort(rep(1:n, n^3))
	idx2 = rep(sort(rep(1:n, n^2)), n)
	idx3 = rep(sort(rep(1:n, n)), n^2)
	idx4 = rep(1:n, n^3)
	idx = cbind(idx1, idx2, idx3, idx4)
	zs = as.numeric(.Call("gogarchcksigma", idx = as.matrix(idx-1), SS = as.numeric(sigmas),
					PACKAGE = "rmgarch"))
	# zs = apply(idx, 1, FUN = function(x) prod(sigmas[x]))
	cs = matrix(zs, ncol = n^3, nrow = n, byrow = TRUE)
	return(cs)
}

.dfft = function(object, index = 1)
{
	# change adjust index for the rolling forecast
	if(object@model$modtype=="goGARCHforecast" | object@model$modtype=="goGARCHroll"){
		if(object@model$n.ahead==1) index = index+1
	}
	switch(object@dist$dist,
			mvnorm = .dfft.mn(object, index),
			manig = .dfft.num(object, index),
			magh = .dfft.num(object, index))
}

.dfft.num = function(object, index = 1)
{
	if( object@dist$support.method == "user" ){
		n = dim(object@dist$y)[2]
		if(index>n) stop("\nindex outside dimensions of object\n", call. = FALSE)
		x = seq(object@dist$fft.support[1], object@dist$fft.support[2], by = object@dist$fft.by)
		return(.makeDNew(x = x, dx = object@dist$y[1:length(x),index], h = object@dist$fft.by, standM = "sum"))
	} else{
		n = length(object@dist$y)
		if(index>n) stop("\nindex outside dimensions of object\n", call. = FALSE)
		x = seq(object@dist$support.user[index, 1], object@dist$support.user[index, 2], by = object@dist$fft.by)
		return(.makeDNew(x = x, dx = object@dist$y[[index]], h = object@dist$fft.by, standM = "sum"))
	}
}

.dfft.mn = function(object, index)
{
	return(function(x) dnorm(x, mean = object@dist$y[index, 1], sd = sqrt(object@dist$y[index, 2])))
	
}

.pfft = function(object, index = 1)
{
	# change adjust index for the rolling forecast
	if(object@model$modtype=="goGARCHforecast" | object@model$modtype=="goGARCHroll"){
		if(object@model$n.ahead==1) index = index+1
	}
	switch(object@dist$dist,
			mvnorm = .pfft.mn(object, index),
			manig = .pfft.num(object, index),
			magh = .pfft.num(object, index))
}

.pfft.num = function(object, index = 1)
{
	if( object@dist$support.method == "user" ){
		n = dim(object@dist$y)[2]
		if(index>n) stop("\nindex outside dimensions of object\n", call. = FALSE)
		x = seq(object@dist$fft.support[1], object@dist$fft.support[2], by = object@dist$fft.by)
		return(.makePNew(x = x, dx = object@dist$y[1:length(x),index], h = object@dist$fft.by))
	} else{
		n = length(object@dist$y)
		if(index>n) stop("\nindex outside dimensions of object\n", call. = FALSE)
		x = seq(object@dist$support.user[index, 1], object@dist$support.user[index, 2], by = object@dist$fft.by)
		return(.makePNew(x = x, dx = object@dist$y[[index]], h = object@dist$fft.by))
	}
}

.pfft.mn = function(object, index)
{
	return(function(q) pnorm(q, mean = object@dist$y[index, 1], sd = sqrt(object@dist$y[index, 2])))
	
}

.qfft = function(object, index = 1)
{
	# change adjust index for the rolling forecast
	if(object@model$modtype=="goGARCHforecast" | object@model$modtype=="goGARCHroll"){
		if(object@model$n.ahead==1) index = index+1
	}
	switch(object@dist$dist,
			mvnorm = .qfft.mn(object, index),
			manig = .qfft.num(object, index),
			magh = .qfft.num(object, index))
}

.qfft.num = function(object, index = 1)
{
	if( object@dist$support.method == "user" ){
		n = dim(object@dist$y)[2]
		if(index>n) stop("\nindex outside dimensions of object\n", call. = FALSE)
		x = seq(object@dist$fft.support[1], object@dist$fft.support[2], by = object@dist$fft.by)
		return(.makeQNew(x = x, dx = object@dist$y[1:length(x),index], h = object@dist$fft.by))
	} else{
		n = length(object@dist$y)
		if(index>n) stop("\nindex outside dimensions of object\n", call. = FALSE)
		x = seq(object@dist$support.user[index, 1], object@dist$support.user[index, 2], by = object@dist$fft.by)
		return(.makeQNew(x = x, dx = object@dist$y[[index]], h = object@dist$fft.by))
	}
}

.qfft.mn = function(object, index)
{
	return(function(p) qnorm(p, mean = object@dist$y[index, 1], sd = sqrt(object@dist$y[index, 2])))
	
}