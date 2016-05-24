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
#---------------------------------------------------------------------------------
# Locally imported copy of "klin.eval" from the klin package
incseq = function (a, b) 
{
	seq(a, b, length = max(0, b - a + 1))
}

.klin.eval = function (A, x, transpose = FALSE) 
{
	m <- sapply(A, nrow)
	n <- sapply(A, ncol)
	K <- length(A)
	if (K == 0) return(x)
	if(transpose) {
		tmp <- m
		m <- n
		n <- tmp
	}
	X <- x
	for(k in K:1) {
		X <- Matrix(as.vector(X), n[k], prod(c(n[incseq(1, k - 1)], m[incseq(k + 1, K)])))
		if(transpose) Y <- Matrix::crossprod(A[[k]], X) else Y <- A[[k]] %*% X
		X <- Matrix::t(Y)
	}
	return(as.vector(X))
}
#---------------------------------------------------------------------------------
fast_kron_M = function(rhs, lhs, n, p=3)
{
	Y = lapply(1:n, function(i){.klin.eval(rhs, lhs[i,], transpose = TRUE)})
	M = matrix(unlist(Y), nrow = n, ncol = n^p, byrow=TRUE)	
	return(M)
}

.kappaGH <-function(x, lambda = 1)
{    
	stopifnot(x >= 0)
	stopifnot(length(lambda) == 1)
	if (lambda == -0.5) {
		# NIG:
		kappa = 1/x
	} else {
		# GH:
		kappa = (
					besselK(x, lambda+1, expon.scaled = TRUE) /
					besselK(x, lambda, expon.scaled = TRUE) ) / x
	}
	return(kappa)
}
# ------------------------------------------------------------------------------
.deltaKappaGH<-function(x, lambda = 1)
{
	if (lambda == -0.5) {
		deltaKappa = .kappaGH(x, lambda+1) - .kappaGH(x, lambda)
	} else {
		deltaKappa = .kappaGH(x, lambda+1) - .kappaGH(x, lambda)
	}
	return(deltaKappa)
}
# ------------------------------------------------------------------------------
.paramGH = function(rho = 0 , zeta = 1, lambda = 1)
{
	# A function implemented by Diethelm Wuertz
	#   Change parameterizations to alpha(zeta, rho, lambda)
	Rho2 = 1 - rho^2
	alpha = zeta^2 * .kappaGH(zeta, lambda) / Rho2 
	alpha = alpha * ( 1 + rho^2 * zeta^2 * .deltaKappaGH(zeta, lambda) / Rho2)
	alpha = sqrt(alpha)  
	beta = alpha * rho
	delta = zeta / ( alpha * sqrt(Rho2) )
	mu = -beta * delta^2 * .kappaGH(zeta, lambda)	
	c(alpha = alpha, beta = beta, delta = delta, mu = mu)  
}
.nigtransform = function(zeta, rho)
{
	nigpars = t(apply(cbind(rho, zeta), 1, FUN = function(x) .paramGH(rho = x[1], zeta = x[2], lambda = -0.5)))
	colnames(nigpars) = c("alpha", "beta", "delta", "mu")
	return(nigpars)
}

.ghyptransform = function(zeta, rho, lambda)
{
	n = length(zeta)
	ghyppars = t(apply(cbind(rho, zeta), 1, FUN = function(x) .paramGH(rho = x[1], zeta = x[2], lambda = lambda)))
	ghyppars = cbind(ghyppars, rep(lambda, n))
	colnames(ghyppars) = c("alpha", "beta", "delta", "mu", "lambda")
	return(ghyppars)
}
.nigscale2 = function(mu, sigma, nigalpha, nigbeta, nigdelta, nigmu)
{
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	# alpha, beta, delta, mu [alpha = shape, beta = skew]
	xdensity[,1] = nigalpha/sigma
	xdensity[,2] = nigbeta/sigma
	xdensity[,3] = nigdelta*sigma
	xdensity[,4] = nigmu*sigma + mu
	colnames(xdensity) = c("alpha", "beta", "delta", "mu")
	return(xdensity)
}
.nigscale = function(mu, sigma, skew, shape)
{
	nigpars = t(apply(cbind(skew, shape), 1, FUN=function(x) .paramGH(rho = x[1], zeta = x[2], lambda = -0.5)))
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	# alpha, beta, delta, mu
	xdensity[,4] = nigpars[,1]/sigma
	xdensity[,3] = nigpars[,2]/sigma
	xdensity[,2] = nigpars[,3]*sigma
	xdensity[,1] = nigpars[,4]*sigma + mu
	colnames(xdensity) = c("mu", "sigma", "skew", "shape")
	return(xdensity)
}

.ghypscale = function(mu, sigma, skew, shape, lambda)
{
	ghpars = t(apply(cbind(skew, shape), 1, FUN=function(x) .paramGH(rho = x[1], zeta = x[2], lambda = lambda)))
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	# alpha, beta, delta, mu
	xdensity[,4] = ghpars[,1]/sigma
	xdensity[,3] = ghpars[,2]/sigma
	xdensity[,2] = ghpars[,3]*sigma
	xdensity[,1] = ghpars[,4]*sigma + mu
	colnames(xdensity) =c("mu", "sigma", "skew", "shape")
	return(xdensity)
}

.ghypscale2 = function(mu, sigma, ghalpha, ghbeta, ghdelta, ghmu)
{
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	# alpha, beta, delta, mu
	xdensity[,1] = ghalpha/sigma
	xdensity[,2] = ghbeta/sigma
	xdensity[,3] = ghdelta*sigma
	xdensity[,4] = ghmu*sigma + mu
	colnames(xdensity) = c("alpha", "beta", "delta", "mu")
	return(xdensity)
}


.normscale = function(mu, sigma)
{
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	xdensity[,1] = mu
	xdensity[,2] = sigma
	xdensity[,3] = 0
	xdensity[,4] = 0
	
	colnames(xdensity) = c("mu", "sigma", "skew", "shape")
	return(xdensity)
}
