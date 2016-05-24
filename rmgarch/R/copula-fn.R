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

##########################################################################################
# Some functions imported directly from QRMlib and mvtnorm and adjusted for use in this package
##########################################################################################

#------------------------------------------------------------------------------------
# Copula densities take the form: 
# dMultivariate/Product(dUnivariate)
# log(dMultivariate) - Sum(log(dUnivariate))
.Spearman = function(data) 
{
	Rho = cor(apply(data, 2, rank))
	return( Rho )
}

.Pconstruct <- function(theta){
	n <- length(theta)
	d <- (1 + sqrt(1+8*n))/2
	A <- matrix(0,nrow=d,ncol=d)
	A[lower.tri(A)] <- theta
	diag(A) <- rep(1,d)
	Q <- A %*% t(A)
	P <- cov2cor(Q)
	P
}

.Pdeconstruct <- function(P){
	A <- t(chol(P))
	Adiag <- diag(diag(A))
	Astar <- solve(Adiag) %*% A
	Astar[lower.tri(Astar)]
}

.vechR = function(Rho) 
{
	Rho[lower.tri(Rho)]
}

.ivechR = function(theta) 
{
	n = length(theta)
	n = (1 + sqrt(1+8*n))/2		
	Rho = matrix(NA, n, n)
	Rho[lower.tri(Rho)] = theta
	Rho[upper.tri(Rho)] = theta
	diag(Rho) = 1
	Rho
}


.vech2 = function(Rho) 
{
	A = t(chol(Rho))
	Adiag = diag(diag(A))
	Astar = solve(Adiag) %*% A
	return( Astar[lower.tri(Astar)] )
}

.ivech2 = function(theta)
{
	n = length(theta)
	d = (1 + sqrt(1+8*n))/2
	A = matrix(0, d, d)
	A[lower.tri(A)] = theta
	diag(A) = rep(1,d)
	Q = A %*% t(A)
	Rho = cov2cor(Q)
	return( Rho )
}

# about 3 times faster than "cor" function
# NO LONGER USED...replaced by cor.fk of pcaPP package
#.Kendall = function(data) 
#{
#	n <- dim(data)[1]
#	d <- dim(data)[2]
#	Rho <- matrix(0, nrow = d, ncol = d)
#	gr = as.matrix(combn(1:d, 2))
#	z = apply(gr, 2, FUN = function(x) Kendall(data[,x[1]], data[, x[2]])$tau)
#	Rho[lower.tri(Rho)] = z
#	Rho = Rho + t(Rho)
#	diag(Rho) = 1
#	nms <- dimnames(data)[[2]]
#	dimnames(Rho) <- list(nms, nms)
#	return( Rho )
#}

#.fit.kendall = function(Udata)
#{
#	Rho = .Kendall(Udata)
#	Rbar = sin(pi * Rho/2)
#	n = dim(Udata)[2]
#	dimn = n * (n - 1) / 2
#	vRbar = Rbar[lower.tri(Rbar)]
#	X = diag(dimn)
#	estimate = lm(vRbar ~ X - 1)$coef
#	A = matrix(0, ncol = n, nrow = n)
#	A[lower.tri(A)] = as.numeric(estimate)
#	A = A + t(A)
#	diag(A) = 1
#	return( A )
#}

########################################################################
# From corpcor package:
# Method by Higham 1988
.makeposdef = function(m)
{	
	d = dim(m)[1]
	es = eigen(m)
	esv = es$values
	delta =  2*d*max(abs(esv))*.Machine$double.eps 	
	tau = pmax(0, delta - esv)
	dm = es$vectors %*% diag(tau, d) %*% t(es$vectors)    	
	return( m +  dm )
}
#------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------

#######################################################################################
# copula likelihood functions (cpp coded)
#------------------------------------------------------------------------------------
.tvtcopulafn = function(pars, data, type = "LLH")
{
	nu = pars[1]
	alpha = pars[2]
	beta = pars[3]
	# this is the slowest part of the routine (student quantile)
	# qt is vectorized but no discernable difference in using apply instead
	Qdata = apply(data, 2, FUN = function(x) rugarch:::qstd(p = x, shape = nu))
	n = dim(data)[1]
	m = dim(data)[2]
	Qbar = cov(Qdata)
	Rbar = cov2cor(Qbar)
	ans = try(.Call("dccCopulaStudent", Qbar = Qbar, U = Qdata, Rbar = Rbar, dcca = alpha, 
					dccb = beta, tnu = nu, dccorder = c(1,1,1), dccsum = alpha + beta, 
					PACKAGE="rmgarch"), silent = TRUE)
	if(inherits(ans, "try-error")){
		ret = Inf
	} else {
		ret = switch(type, LLH = -ans[[3]], ALL = ans)
	}
	return(ret)
}

.tcopulafn = function(pars, data, Rbar, type = "LLH")
{
	nu = pars[1]
	Qdata = apply(data, 2, FUN = function(x) rugarch:::qstd(p = x ,shape = nu))
	#ans = try(.Call("staticCopulaStudent", U = Qdata, Rbar = Rbar, tnu = nu), silent = TRUE)
	ans = dcopula.student(U = Qdata, Corr = Rbar, df = nu, logvalue = TRUE)
		
	ret = switch(type, 
				LLH = -ans, 
				ALL = list(LLH = ans, Rbar = Rbar))
	return(ret)
}

.tvgausscopulafn = function(pars, data, type = "LLH")
{
	alpha = pars[2]
	beta = pars[3]
	Qdata = apply(data, 2, FUN = function(x) qnorm(p = x))
	n = dim(data)[1]
	m = dim(data)[2]
	Qbar = cov(Qdata)
	Rbar = cov2cor(Qbar)
	ans = try(.Call("dccCopulaNormal", Qbar = Qbar, U = Qdata, Rbar = Rbar, dcca = alpha, 
					dccb = beta, dccorder = c(1,1,1), dccsum = alpha + beta, 
					PACKAGE="rmgarch"), silent = TRUE)
	if(inherits(ans, "try-error")){
		ret = Inf
	} else {
		ret = switch(type, LLH = -ans[[3]], ALL = ans)
	}
	return(ret)
}

.gausscopulafn = function(pars, data, Rbar, type = "LLH")
{
	Qdata = apply(data, 2, FUN = function(x) qnorm(p = x))
	ans = dcopula.gauss(U = Qdata, Sigma = Rbar, logvalue = TRUE)
	ret = switch(type, 
				LLH = -ans, 
				ALL = list(LLH = ans, Rbar = Rbar))
	return(ret)
}
#------------------------------------------------------------------------------------



#####################################################################################
# Auxiliary Functions
#------------------------------------------------------------------------------------

.cor2cov = function(corr, sigmas)
{
	V = (sigmas)%o%(sigmas)*corr
	V
}