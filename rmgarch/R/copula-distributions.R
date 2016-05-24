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

# Implements copula methods for p, d, q, r for multivariate Normal amd Student.
# Future expansions should likely include the multivariate skew Normal and skew Student
# of Azzalini.

# General Setup:

# dcopula:
#-----------------------------------------------------------------------------------------------------------
# 	 
#	 c(u_1, ... , u_n) = multivariate_pdf( margin_quantile_1(u_1), ... , margin_quantile_n(u_n) )
# 	----------------------------------------------------------------------------------------------
# 					SumProduct_{1..n}[ margin_pdf_i( margin_quantile_i(u_i) ) ]
#
# where u_i is the uniform margin obtained by applying the univariate_cdf to the fitted data (each i can have
# a seperate distribution fitted)
#-----------------------------------------------------------------------------------------------------------


# pcopula:
#-----------------------------------------------------------------------------------------------------------
# C(u_1, ..., u_n) = multivariate_cdf( margin_quantile_1(u_1), ... , margin_quantile_n(u_n) )
#
# where u_i is the uniform margin obtained by applying the univariate_cdf to the fitted data (each i can have
# a seperate distribution fitted)
#-----------------------------------------------------------------------------------------------------------


# conversely, a multivariate density function multivariate_pdf() can be decomposed as:

# dinvcopula:
#-----------------------------------------------------------------------------------------------------------
# multivariate_pdf(x_1, ..., x_n) = c(margin_cdf_1(x_1), ... , margin_cdf_n(x_n)) * SumProduct_{1..n}[ univariate_pdf(x_i) ]
#
#-----------------------------------------------------------------------------------------------------------

# pinvcopula:
#-----------------------------------------------------------------------------------------------------------
# multivariate_cdf(x_1, ..., x_n) = c(margin_cdf_1(x_1), ... , margin_cdf_n(x_n))
#
#-----------------------------------------------------------------------------------------------------------

# Density of full Copula GARCH Model
dcgarch = function(mfit, cfit)
{
	# The full density of the model is given by:
	# log(copula density) + log(marginal GARCH densities)
	# where copula density = log( multivariate_f(q_1*(u_1),...,q_n(u_n)) ) ...
	# - (log(univariate_f(q_1*(u_1)) + log(univariate_f(q_n*(u_n)))
	garchllh = sum(sapply(mfit@fit, FUN = function(x) likelihood(x)))
	copllh = cfit$LLH
	llh = garchllh + copllh
	return(llh)
}


dcopula.gauss = function(U, Sigma, logvalue = FALSE){
	m = dim(Sigma)[2]
	# U = uniform(0,1)
	Z = apply(U, 2, FUN = function(x) qnorm(x))
	# Z = transormed to standard normal variates
	ans = .dmvnorm(Z, mean = rep(0, m), sigma = Sigma, log = TRUE) - apply(dnorm(Z, log = TRUE), 1,  "sum")
	if( !logvalue ) ans = exp( ans )
	return( ans )
}

#pcopula.gauss = function(U, Sigma, ...)
#{
#	require(mvtnorm)
#	m = dim(Sigma)[2]
#	U = matrix(U, ncol = m)
#	Z = apply(U, 2, FUN = function(x) qnorm(x))
#	mu = rep(0, m)
#	ans = apply(Z, 2, FUN = function(x) mvtnorm::pmvnorm(lower = rep(-Inf, m),  upper = x, mean = mu, sigma = Sigma, ...))
#	return( ans )
#}

rcopula.gauss = function(n, U, Sigma, ...)
{
	m = dim(Sigma)[2]
	mu = rep(0, m)
	ans = pnorm(.rmvnorm(n, mean = mu, sigma = Sigma))
	return( ans )
}

dcopula.student = function(U, Corr, df, logvalue = FALSE)
{
	m = dim(Corr)[2]
	Z = apply(U, 2, FUN = function(x) rugarch:::qstd(p = x , shape = df))
	mu = rep(0, m)
	ans = .dmvt(Z, delta = mu, sigma = Corr, df = df, log = TRUE) - apply(Z, 1, FUN = function(x) sum(rugarch:::dstd(x, shape = df, log = TRUE)))
	# .Call("dcopulaStudent", Z = Z, m = matrix(0, nrow = 1, ncol = m), sigma = Corr, df = df, dtZ = dt(Z, df = df, log = TRUE), 
	# PACKAGE = "rmgarch")
	if( !logvalue ) ans = exp(ans)
	return( ans )
}

#pcopula.student = function(U, Corr, df, ...)
#{
#	require(mvtnorm)
#	m = dim(Corr)[2]
#	U = matrix(U, ncol = m)
#	Z = apply(U, 2, FUN = function(x) qt(x, df = df))
#	mu = rep(0, m)
#	ans = apply(Z, 1, FUN = function(x) mvtnorm::pmvt(lower = rep(-Inf, m), upper = x, delta = mu, corr = Corr, df = df, ...))
#	return( ans )
#}

rcopula.student = function(n, U, Corr, df)
{
	m = dim(Corr)[2]
	mu = rep(0, m)
	ans = rugarch:::pstd(.rmvt(n, mu = mu, R = Corr, df = df), mu=0, sigma=1, shape = df)
	return ( ans )
}


#####################################################################################
# Copula Simulation
#------------------------------------------------------------------------------------
.sample.copula = function(model, Qbar, Rbar, Nbar, preQ = NULL, preZ = NULL, n.sim, n.start = 0, m.sim, rseed, cluster = NULL)
{
	timecopula = model$modeldesc$timecopula
	set.seed(rseed)
	if(timecopula){
		ans = .rtvcopula(model, Qbar, Rbar, Nbar, preQ = preQ, preZ = preZ, n.sim, n.start, m.sim, rseed, cluster = cluster)
	} else{
		ans = .rcopula(model, Rbar, n.sim, n.start, m.sim, rseed)
	}
	return(ans)
}


.rtvcopula = function(model, Qbar, Rbar, Nbar, preQ, preZ, n.sim, n.start = 0, m.sim, rseed, cluster = NULL)
{
	dccOrder = model$modelinc[4:5]
	cf = model$ipars[,1]
	idx = model$pidx
	mo = max( dccOrder )
	m = dim(Rbar)[1]
	z = array(NA,  dim = c(n.sim + n.start + mo, m, m.sim))
	if(model$modeldesc$distribution == "mvt"){
		for(i in 1:m.sim){
			set.seed(rseed[i])
			z[,,i] = rbind(preZ, .rmvt(n = (n.sim + n.start), mu = rep(0, m), 
					R = diag(m), df = as.numeric(cf["mshape"])))
		}
	} else{
		for(i in 1:m.sim){
			set.seed(rseed[i])
			z[,,i] = rbind(preZ, .rmvnorm(n = (n.sim + n.start), mean = rep(0, m), sigma = diag(m)))
		}
	}
	xseed = rseed+1
	simR = vector(mode = "list", length = m.sim)
	mtmp = vector(mode="list", length=m.sim)
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, require(rmgarch))
		clusterExport(cluster, c("model", "z", "preQ", "Rbar", "Nbar", 
						"mo", "n.sim", "n.start", "m", "xseed"), envir = environment())
		clusterExport(cluster, ".copuladccsimf", envir = environment())
		mtmp = parLapply(cluster, as.list(1:m.sim), fun = function(j){
					.copuladccsimf(model, Z = z[,,j], Qbar = Qbar, 
							preQ = preQ, Nbar = Nbar, Rbar = Rbar, mo = mo, 
							n.sim, n.start, m, rseed[j])
				})
	} else{
		for(i in 1:m.sim){
			mtmp[[i]] = .copuladccsimf(model, Z = z[,,i], Qbar = Qbar, preQ = preQ, 
					Nbar = Nbar, Rbar = Rbar, mo = mo, n.sim, n.start, m, rseed[i])
		}
	}
	simR = lapply(mtmp, FUN = function(x) if(is.matrix(x$R)) array(x$R, dim = c(m, m, n.sim)) else last(x$R, n.sim))
	Ures = vector(mode = "list", length = m)
	Usim = array(NA, dim = c(n.sim+n.start, m, m.sim))
	if(model$modeldesc$distribution  == "mvt"){
		for(i in 1:m) Ures[[i]] = matrix(rugarch:::pstd(sapply(mtmp, FUN = function(x) x$Z[,i]), mu=0, sigma=1, shape = cf["mshape"]), ncol = m.sim)
		for(i in 1:m.sim) Usim[,,i] = matrix(sapply(Ures, FUN = function(x) x[-(1:mo),i]), ncol = m)
	} else{
		for(i in 1:m) Ures[[i]] = pnorm(matrix(sapply(mtmp, FUN = function(x) x$Z[,i]), ncol = m.sim))
		for(i in 1:m.sim) Usim[,,i] = matrix(sapply(Ures, FUN = function(x) x[-(1:mo),i]), ncol = m)
	}
	rm(Ures)
	rm(mtmp)
	gc(verbose=FALSE)
	return(list(Usim = Usim, simR = simR))
}

.rcopula = function(model, Rbar, n.sim, n.start = 0, m.sim, rseed)
{
	nsim = n.sim + n.start
	m = dim(Rbar)[1]
	sim = array(data = NA, dim = c(nsim , m , m.sim))
	if(model$modeldesc$distribution == "mvt"){
		cf = model$ipars[,1]
		shape = cf["mshape"]
		for(i in 1:m.sim){
			set.seed(rseed[i])
			tmp = .rmvt(n = nsim, R = Rbar, df = shape, mu = rep(0, m))
			sim[,,i] = matrix(rugarch:::pstd(tmp, mu=0, sigma=1, shape = shape), nrow = nsim, ncol = m)
		}
	} else{
		for(i in 1:m.sim){
			set.seed(rseed[i])
			tmp = .rmvnorm(n = nsim, mean = rep(0, m), sigma = Rbar) 
			sim[,,i] = matrix(pnorm(tmp), nrow = nsim, ncol = m)
		}
	}
	return( sim )
}

.copuladccsimf = function(model, Z, Qbar, preQ, Rbar, Nbar, mo, n.sim, n.start, m, rseed){
	modelinc = model$modelinc
	ipars = model$pars
	idx = model$pidx
	n = n.sim + n.start + mo
	set.seed(rseed[1]+1)
	stdresid = matrix(rnorm(m * (n.sim+n.start+mo)), nrow = n.sim + n.start + mo, ncol = m)
	sumdcca = sum(ipars[idx["dcca",1]:idx["dcca",2],1])
	sumdccb = sum(ipars[idx["dccb",1]:idx["dccb",2],1])
	sumdcc = sumdcca + sumdccb
	sumdccg = sum(ipars[idx["dccg",1]:idx["dccg",2],1])
	res = switch(model$modeldesc$distribution,
	mvnorm = .Call("copuladccsimmvn", model = as.integer(modelinc), pars = as.numeric(ipars[,1]), 
			idx = as.integer(idx[,1]-1), Qbar = as.matrix(Qbar), preQ = as.matrix(preQ), 
			Rbar = as.matrix(Rbar), Nbar = as.matrix(Nbar), Z = as.matrix(Z),  Res = as.matrix(stdresid), 
			epars = c(sumdcc, sumdccg, mo), PACKAGE = "rmgarch"),
	mvt = .Call("copuladccsimmvt", model = as.integer(modelinc), pars = as.numeric(ipars[,1]), 
			idx = as.integer(idx[,1]-1), Qbar = as.matrix(Qbar), preQ = as.matrix(preQ), 
			Rbar = as.matrix(Rbar), Nbar = as.matrix(Nbar), Z = as.matrix(Z),  NZ = as.matrix(stdresid), 
			epars = c(sumdcc, sumdccg, mo), PACKAGE = "rmgarch"))
	Q = array(NA, dim = c(m, m, n.sim + n.start + mo))
	R = array(NA, dim = c(m, m, n.sim + n.start + mo))
	for(i in 1:(n.sim + n.start + mo)){
		R[,,i] = res[[2]][[i]]
		Q[,,i] = res[[3]][[i]]
	}
	Z = res[[3]]
	ans = list( Q = Q, R = R, Z = Z)
	rm(res)
	return( ans )
}
#------------------------------------------------------------------------------------

################################################################################
# some functions from mvtnorm package implemented locally here
.rmvnorm = function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean))){
	ev <- eigen(sigma, symmetric = TRUE)
	if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
		warning("sigma is numerically not positive definite")
	}
	retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% solve(ev$vectors)
	retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
	retval <- sweep(retval, 2, mean, "+")
	colnames(retval) <- names(mean)
	return( retval )
}


.dmvnorm = function (x, mean, sigma, log = FALSE) 
{
	if (is.vector(x)) {
		x <- matrix(x, ncol = length(x))
	}
	if (missing(mean)) {
		mean <- rep(0, length = ncol(x))
	}
	if (missing(sigma)) {
		sigma <- diag(ncol(x))
	}
	if (NCOL(x) != NCOL(sigma)) {
		stop("x and sigma have non-conforming size")
	}
	if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
			check.attributes = FALSE)) {
		stop("sigma must be a symmetric matrix")
	}
	if (length(mean) != NROW(sigma)) {
		stop("mean and sigma have non-conforming size")
	}
	distval <- mahalanobis(x, center = mean, cov = sigma)
	logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
	logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
	if(log) retval = logretval else retval = exp(logretval)
	return(retval)
}


.rmvt = function (n, R = diag(2), df = 5, mu = rep(0, nrow(R)))
{
	if (length(mu) != nrow(R)) stop("mu and sigma have non-conforming size")
	if (df == 0) return(.rmvnorm(n, mean = mu, sigma = R))
	rc = rchisq(n, df)
	RR = ((df-2)/df) * R
	v = sqrt(df/rc)
	retval = repmat(v,1,ncol(RR)) * .rmvnorm(n, mean = mu, sigma = RR)
	return(retval)
}


.dmvt = function (x, delta, sigma, df = 1, log = TRUE, type = "shifted") 
{
	if (df == 0) return(.dmvnorm(x, mean = delta, sigma = sigma, log = log))
	if (is.vector(x)) {
		x <- matrix(x, ncol = length(x))
	}
	if (missing(delta)) {
		delta <- rep(0, length = ncol(x))
	}
	if (missing(sigma)) {
		sigma <- diag(ncol(x))
	}
	if (NCOL(x) != NCOL(sigma)) {
		stop("x and sigma have non-conforming size")
	}
	if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
			check.attributes = FALSE)) {
		stop("sigma must be a symmetric matrix")
	}
	if (length(delta) != NROW(sigma)) {
		stop("mean and sigma have non-conforming size")
	}
	m <- NCOL(sigma)
	distval <- mahalanobis(x, center = delta, cov = sigma)
	logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
	logretval <- lgamma((m + df)/2) - (lgamma(df/2) + 0.5 * (logdet + 
					m * logb(pi * df))) - 0.5 * (df + m) * logb(1 + distval/df)
	if(log) retval = logretval else retval = exp(logretval)
	return(retval)
}
