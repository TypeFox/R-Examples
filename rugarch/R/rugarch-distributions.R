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

# This file contains the location-scale invariant parameterization of the
# distributions used in the rugarch package. Since version 1.2-5, the 
# dpqr functions are based on the C coded functions and not R. Vectorization
# is handled internally without resorting to the Vectorize function which 
# was found to have a large overhead.

# Some of the original standardized distributions, implemented in R and found in 
# a number of the Rmetrics packages, were based on the work of Diethelm Wuertz 
# whose contribution is greatly appreciated. Some of the nig and gig C code is
# based on the work of others and this is carefully aknowledged in the source
# files.

# ------------------------------------------------------------------------------
# Skew Generalized Hyberolic Student's T 
# alpha = abs(beta)+1e-12, lambda = -nu/2
# ------------------------------------------------------------------------------
# Location-Scale Invariant Parametrization
.paramGHST = function(betabar, nu){
	# Alexios Ghalanos 2012
	# betabar is the skew parameter = beta*delta (parametrization 4 in Prause)
	# nu is the shape parameter
	# details can be found in the vignette
	delta = ( ((2 * betabar^2)/((nu-2)*(nu-2)*(nu-4))) + (1/(nu-2)) )^(-0.5)
	beta = betabar/delta
	mu = -( (beta * (delta^2))/(nu-2))
	return( c(mu, delta, beta, nu) )
}

dsghst = function(x, mu=0, sigma=1, skew=1, shape=8, log = FALSE)
{
	if(any(abs(skew)<1e-12)) skew[which(abs(skew)<1e-12)] = 1e-12
	n = c(length(x), length(mu), length(sigma), length(skew), length(shape))
	maxn = max(n)
	if(n[1]!=maxn) x = rep(x[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) skew  = rep(skew[1], maxn)
	if(n[5]!=maxn) shape = rep(shape[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_dghst", x = as.double(x), mu = as.double(mu), 
					sigma = as.double(sigma), skew = as.double(skew), 
					shape = as.double(shape), ans = ans, n = as.integer(maxn), 
					logr = as.integer(log), PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

rsghst = function(n, mu=0, sigma=1, skew=1, shape=8)
{
	if(any(abs(skew)<1e-12)) skew[which(abs(skew)<1e-12)] = 1e-12
	nn = c(length(mu), length(sigma), length(skew), length(shape))
	if(nn[1]!=n) mu    = rep(mu[1], n)
	if(nn[2]!=n) sigma = rep(sigma[1], n)
	if(nn[3]!=n) skew  = rep(skew[1], n)
	if(nn[4]!=n) shape = rep(shape[1], n)
	ans = double(n)
	sol = try(.C("c_rghst", n = as.integer(n), mu = as.double(mu), 
					sigma = as.double(sigma), skew = as.double(skew), 
					shape = as.double(shape), ans = ans, 
					PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

psghst = function(q, mu=0, sigma=1, skew=1, shape=8, lower.tail = TRUE, log.p = FALSE, ...){
	
	if(any(abs(skew)<1e-12)) skew[which(abs(skew)<1e-12)] = 1e-12
	n = c(length(q), length(mu), length(sigma), length(skew), length(shape))
	maxn = max(n)
	if(n[1]!=maxn) q = rep(q[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) skew  = rep(skew[1], maxn)
	if(n[5]!=maxn) shape = rep(shape[1], maxn)
	ans = double(maxn)
	for(i in 1:maxn){
		ans[i] = .psghst(q[i], mu=mu[i], sigma=sigma[i], skew=skew[i], 
				shape=shape[i], lower.tail = lower.tail, log.p = log.p)
	}
	return( ans )
}

.psghst = function(q, mu=0, sigma=1, skew=1, shape=8, lower.tail = TRUE, log.p = FALSE, ...){
	param = .paramGHST(skew, shape)
	# scale the parameters
	mux = param[1]*sigma + mu
	delta = param[2]*sigma
	beta = param[3]/sigma
	nu = param[4]
	ans = SkewHyperbolic::pskewhyp(q, mu = mux, delta = delta, beta = beta, nu = nu,
            log.p = log.p, lower.tail = lower.tail)
	return(ans)
}

qsghst = function(p, mu=0, sigma=1, skew=1, shape=8){
	if(any(abs(skew)<1e-12)) skew[which(abs(skew)<1e-12)] = 1e-12
	n = c(length(p), length(mu), length(sigma), length(skew), length(shape))
	maxn = max(n)
	if(n[1]!=maxn) p = rep(p[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) skew  = rep(skew[1], maxn)
	if(n[5]!=maxn) shape = rep(shape[1], maxn)
	ans = double(maxn)
	for(i in 1:maxn){
		ans[i] = .qsghst(p[i], mu=mu[i], sigma=sigma[i], skew=skew[i], 
				shape=shape[i])
	}
	return( ans )
}

.qsghst = function(p, mu=0, sigma=1, skew=1, shape=8, lower.tail = TRUE, ...){
	if (!lower.tail) {
		p <- 1 - p
		lower.tail == TRUE
	}
	param = .paramGHST(skew, shape)
	# scale the parameters
	mux = param[1]*sigma + mu
	delta = param[2]*sigma
	beta = param[3]/sigma
	nu = param[4]
	ans = SkewHyperbolic::qskewhyp(p, mu = mux, delta = delta, beta = beta, 
			nu = nu, lower.tail = lower.tail, method = c("spline","integrate")[1])
	return(ans)
}

ghstFit = function(x, control){
	f = function(pars, x){
		return(sum(-dsghst(x, mu=pars[1], sigma=pars[2], skew=pars[3], shape=pars[4], log = TRUE)))
	}
	x = as.numeric(x)
	x0 = c(mean(x), sd(x), 0.2, 8)
	# ARE THE shape LOWER BOUNDS TOO HIGH?
	# nu>6 (skewness exist) and nu>8 (kurtosis exist)
	fit = try(solnp(x0, fun = f, LB = c(-5, 1e-10, -10, 4.001), UB = c(5, 10, 10, 25), control = control, x = x), 
			silent = TRUE)
	# Add Names to $par
	names(fit$par) = c("mu", "sigma", "skew", "shape")
	# Return Value:
	return( fit )

}
# equivalence: 
# 1. dskewhyp(x, mu = (sigma*params[1]+0.5), delta = params[2]*sigma, beta = params[3]/sigma, nu = shape)
# 2. dskewhyp(z, mu = params[1], delta = params[2], beta = params[3], nu = shape)/sigma

# ------------------------------------------------------------------------------
# Skew Normal Distribution (Fernandez and Steel)
# ------------------------------------------------------------------------------
dsnorm = function(x, mu = 0, sigma = 1, skew = 1.5, log = FALSE)
{
	n = c(length(x), length(mu), length(sigma), length(skew))
	maxn = max(n)
	if(n[1]!=maxn) x = rep(x[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) skew  = rep(skew[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_dsnorm", x = as.double(x), mu = as.double(mu), 
					sigma = as.double(sigma), skew = as.double(skew), ans = ans, 
					n = as.integer(maxn), logr = as.integer(log), 
					PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

psnorm = function(q, mu = 0, sigma = 1, skew = 1.5)
{
	n = c(length(q), length(mu), length(sigma), length(skew))
	maxn = max(n)
	if(n[1]!=maxn) q = rep(q[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) skew  = rep(skew[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_psnorm", q = as.double(q), mu = as.double(mu), 
					sigma = as.double(sigma), skew = as.double(skew), ans = ans, 
					n = as.integer(maxn), PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

qsnorm = function(p, mu = 0, sigma = 1, skew = 1.5)
{
	n = c(length(p), length(mu), length(sigma), length(skew))
	maxn = max(n)
	if(n[1]!=maxn) p = rep(p[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) skew  = rep(skew[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_qsnorm", p = as.double(p), mu = as.double(mu), 
					sigma = as.double(sigma), skew = as.double(skew), ans = ans, 
					n = as.integer(maxn), PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

rsnorm = function(n, mu=0, sigma=1, skew=1.5)
{
	nn = c(length(mu), length(sigma), length(skew))
	if(nn[1]!=n) mu    = rep(mu[1], n)
	if(nn[2]!=n) sigma = rep(sigma[1], n)
	if(nn[3]!=n) skew  = rep(skew[1], n)
	ans = double(n)
	sol = try(.C("c_rsnorm", n = as.integer(n), mu = as.double(mu), 
					sigma = as.double(sigma), skew = as.double(skew), ans = ans, 
					PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

snormFit = function(x, control = list())
{   
	ctrl = .solnpctrl(control)
	start = c(mu = mean(x), sigma = sqrt(var(x)), xi = 1)	
	loglik = function(pars, x){ 
		f = -sum(log(dsnorm(x, pars[1], pars[2], pars[3])))
		f }
	fit = solnp(pars = start, fun = loglik, 
			LB = c(-Inf, 0, 0), UB = c(Inf, Inf, Inf), control = ctrl, x = x)
	names(fit$pars) = c("mu", "sigma", "skew")	
	return( fit )
}
# ------------------------------------------------------------------------------
# Normal Distribution (estimation function)
# ------------------------------------------------------------------------------
normFit <-function(x, control = list())
{   
	ctrl = .solnpctrl(control)
	start = c(mu = mean(x), sigma = sqrt(var(x)))
	loglik = function(pars, x){ 
		f = -sum(log(dnorm(x, pars[1], pars[2])))
		f }
	fit = solnp(pars = start, fun = loglik, LB = c(-Inf, 0), UB = c(Inf, Inf), 
			control = ctrl, x = x)
	names(fit$pars) = c("mu", "sigma")
	return( fit )
}   

# ------------------------------------------------------------------------------
# Generalized Error Distribution
# ------------------------------------------------------------------------------

dged = function(x, mu = 0, sigma = 1, shape = 2, log = FALSE)
{
	n = c(length(x), length(mu), length(sigma), length(shape))
	maxn = max(n)
	if(n[1]!=maxn) x = rep(x[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) shape = rep(shape[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_dged", x = as.double(x), mu = as.double(mu), 
					sigma = as.double(sigma), shape = as.double(shape), 
					ans = ans, n = as.integer(maxn), logr = as.integer(log), 
					PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

pged = function(q, mu = 0, sigma = 1, shape = 2)
{
	n = c(length(q), length(mu), length(sigma), length(shape))
	maxn = max(n)
	if(n[1]!=maxn) q = rep(q[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) shape = rep(shape[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_pged", q = as.double(q), mu = as.double(mu), 
					sigma = as.double(sigma), shape = as.double(shape), 
					ans = ans, n = as.integer(maxn),  
					PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

qged = function(p, mu = 0, sigma = 1, shape = 2)
{
	n = c(length(p), length(mu), length(sigma), length(shape))
	maxn = max(n)
	if(n[1]!=maxn) p = rep(p[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) shape = rep(shape[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_qged", p = as.double(p), mu = as.double(mu), 
					sigma = as.double(sigma), shape = as.double(shape), 
					ans = ans, n = as.integer(maxn),  
					PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

rged = function(n, mu = 0, sigma = 1, shape = 2)
{
	nn = c(length(mu), length(sigma), length(shape))
	if(nn[1]!=n) mu    = rep(mu[1], n)
	if(nn[2]!=n) sigma = rep(sigma[1], n)
	if(nn[3]!=n) shape = rep(shape[1], n)
	ans = double(n)
	sol = try(.C("c_rged", n = as.integer(n), mu = as.double(mu), 
					sigma = as.double(sigma), shape = as.double(shape), 
					ans = ans, PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

gedFit = function(x, control = list())
{   
	ctrl = .solnpctrl(control)
	start = c(mu = mean(x), sigma = sqrt(var(x)), shape = 2)
	loglik = function(pars, x){ 
		f = -sum(log(dged(x, pars[1], pars[2], pars[3])))
		f }
	fit = solnp(pars = start, fun = loglik, 
			LB = c(-Inf, 0, 0), UB = c(Inf, Inf, Inf), , control  = ctrl, x = x)	
	names(fit$pars) = c("mu", "sigma", "shape")	
	return( fit )
}
# ------------------------------------------------------------------------------
# Skewed Generalized Error Distribution (Fernandez and Steel)
# ------------------------------------------------------------------------------
dsged = function(x, mu = 0, sigma = 1, skew = 1.5, shape = 2, log = FALSE)
{
	n = c(length(x), length(mu), length(sigma), length(skew), length(shape))
	maxn = max(n)
	if(n[1]!=maxn) x = rep(x[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) skew  = rep(skew[1], maxn)
	if(n[5]!=maxn) shape = rep(shape[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_dsged", x = as.double(x), mu = as.double(mu), 
					sigma = as.double(sigma), skew = as.double(skew), 
					shape = as.double(shape), ans = ans, n = as.integer(maxn), 
					logr = as.integer(log), PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

psged = function(q, mu = 0, sigma = 1, skew = 1.5, shape = 2)
{
	n = c(length(q), length(mu), length(sigma), length(skew), length(shape))
	maxn = max(n)
	if(n[1]!=maxn) q = rep(q[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) skew  = rep(skew[1], maxn)
	if(n[5]!=maxn) shape = rep(shape[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_psged", q = as.double(q), mu = as.double(mu), 
					sigma = as.double(sigma), skew = as.double(skew), 
					shape = as.double(shape), ans = ans, 
					n = as.integer(maxn), PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

qsged = function(p, mu = 0, sigma = 1, skew = 1.5, shape = 2)
{
	n = c(length(p), length(mu), length(sigma), length(skew), length(shape))
	maxn = max(n)
	if(n[1]!=maxn) p = rep(p[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) skew  = rep(skew[1], maxn)
	if(n[5]!=maxn) shape = rep(shape[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_qsged", p = as.double(p), mu = as.double(mu), 
					sigma = as.double(sigma), skew = as.double(skew), 
					shape = as.double(shape), ans = ans, n = as.integer(maxn), 
					PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}


rsged = function(n, mu=0, sigma = 1, skew = 1.5, shape = 2)
{
	nn = c(length(mu), length(sigma), length(skew), length(shape))
	if(nn[1]!=n) mu    = rep(mu[1], n)
	if(nn[2]!=n) sigma = rep(sigma[1], n)
	if(nn[3]!=n) skew  = rep(skew[1], n)
	if(nn[4]!=n) shape = rep(shape[1], n)
	ans = double(n)
	sol = try(.C("c_rsged", n = as.integer(n), mu = as.double(mu), 
					sigma = as.double(sigma), skew = as.double(skew), 
					shape = as.double(shape), 
					ans = ans, PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

sgedFit = function(x, control = list())
{   
	ctrl = .solnpctrl(control)
	start = c(mu = mean(x), sigma = sqrt(var(x)), skew = 1, shape = 2)
	loglik = function(pars, x){ 
		f = -sum(log(dsged(x, pars[1], pars[2], pars[3], pars[4])))
		f }
	fit = solnp(pars = start, fun = loglik, 
			LB = c(-Inf, 0, 0, 0), UB = c(Inf, Inf, Inf, Inf), , control  = ctrl, x = x)
	names(fit$pars) = c("mu", "sigma", "skew", "shape")
	return( fit )
}
################################################################################
Heaviside<-function(x, a = 0) 
{   
	result = (sign(x-a) + 1)/2
	return( result )
}
# ------------------------------------------------------------------------------
# Student Distribution
# ------------------------------------------------------------------------------
dstd = function(x, mu = 0, sigma = 1, shape = 5, log = FALSE)
{
	n = c(length(x), length(mu), length(sigma), length(shape))
	maxn = max(n)
	if(n[1]!=maxn) x = rep(x[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) shape = rep(shape[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_dstd", x = as.double(x), mu = as.double(mu), 
					sigma = as.double(sigma), shape = as.double(shape), 
					ans = ans, n = as.integer(maxn), logr = as.integer(log), 
					PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}
	
pstd = function(q, mu = 0, sigma = 1, shape = 5)
{
	n = c(length(q), length(mu), length(sigma), length(shape))
	maxn = max(n)
	if(n[1]!=maxn) q = rep(q[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) shape = rep(shape[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_pstd", q = as.double(q), mu = as.double(mu), 
					sigma = as.double(sigma), shape = as.double(shape), 
					ans = ans, n = as.integer(maxn), 
					PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

qstd = function(p, mu = 0, sigma = 1, shape = 5)
{
	n = c(length(p), length(mu), length(sigma), length(shape))
	maxn = max(n)
	if(n[1]!=maxn) p = rep(p[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) shape = rep(shape[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_qstd", p = as.double(p), mu = as.double(mu), 
					sigma = as.double(sigma), shape = as.double(shape), 
					ans = ans, n = as.integer(maxn), 
					PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

rstd = function(n, mu = 0, sigma = 1, shape = 5)
{
	nn = c(length(mu), length(sigma), length(shape))
	if(nn[1]!=n) mu    = rep(mu[1], n)
	if(nn[2]!=n) sigma = rep(sigma[1], n)
	if(nn[3]!=n) shape = rep(shape[1], n)
	ans = double(n)
	sol = try(.C("c_rstd", n = as.integer(n), mu = as.double(mu), 
					sigma = as.double(sigma), shape = as.double(shape), 
					ans = ans, PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

stdFit = function(x, control = list())
{
	ctrl = .solnpctrl(control)
	start = c(mu = mean(x), sigma = sqrt(var(x)), shape = 4)
	loglik = function(pars, x){
		f = -sum(log(dstd(x, pars[1], pars[2], pars[3])))
		f }	
	fit = solnp(pars = start, fun = loglik,
			LB = c(-Inf, 0, 2.01), UB = c(Inf, Inf, Inf), control = ctrl, x = x)	
	names(fit$pars) = c("mu", "sigma", "shape")
	return( fit )
}

# ------------------------------------------------------------------------------
# Skewed Student Distribution (Fernandez and Steel)
# ------------------------------------------------------------------------------
dsstd = function(x, mu = 0, sigma = 1, skew = 1.5, shape = 5, log = FALSE)
{
	n = c(length(x), length(mu), length(sigma), length(skew), length(shape))
	maxn = max(n)
	if(n[1]!=maxn) x = rep(x[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) skew  = rep(skew[1], maxn)
	if(n[5]!=maxn) shape = rep(shape[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_dsstd", x = as.double(x), mu = as.double(mu), 
					sigma = as.double(sigma), skew = as.double(skew), 
					shape = as.double(shape), ans = ans, n = as.integer(maxn), 
					logr = as.integer(log), PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

psstd = function(q, mu = 0, sigma = 1, skew = 1.5, shape = 5)
{
	n = c(length(q), length(mu), length(sigma), length(skew), length(shape))
	maxn = max(n)
	if(n[1]!=maxn) q = rep(q[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) skew  = rep(skew[1], maxn)
	if(n[5]!=maxn) shape = rep(shape[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_psstd", q = as.double(q), mu = as.double(mu), 
					sigma = as.double(sigma), skew = as.double(skew), 
					shape = as.double(shape), ans = ans, n = as.integer(maxn),
					PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

qsstd = function(p, mu = 0, sigma = 1, skew = 1.5, shape = 5)
{
	n = c(length(p), length(mu), length(sigma), length(skew), length(shape))
	maxn = max(n)
	if(n[1]!=maxn) p = rep(p[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) skew  = rep(skew[1], maxn)
	if(n[5]!=maxn) shape = rep(shape[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_qsstd", p = as.double(p), mu = as.double(mu), 
					sigma = as.double(sigma), skew = as.double(skew), 
					shape = as.double(shape), ans = ans, n = as.integer(maxn),
					PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

rsstd = function(n, mu = 0, sigma = 1, skew = 1.5, shape = 5)
{
	nn = c(length(mu), length(sigma), length(skew), length(shape))
	if(nn[1]!=n) mu    = rep(mu[1], n)
	if(nn[2]!=n) sigma = rep(sigma[1], n)
	if(nn[3]!=n) skew  = rep(skew[1], n)
	if(nn[4]!=n) shape = rep(shape[1], n)
	ans = double(n)
	sol = try(.C("c_rsstd", n = as.integer(n), mu = as.double(mu), 
					sigma = as.double(sigma), skew = as.double(skew), 
					shape = as.double(shape), 
					ans = ans, PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

sstdFit = function(x, control = list())
{
	ctrl = .solnpctrl(control)
	start = c(mu = mean(x), sigma = sqrt(var(x)), skew = 1.5, shape = 4)	
	loglik = function(pars, x){
		f = -sum(log(dsstd(x, pars[1], pars[2], pars[3], pars[4])))
		f }	
	fit = solnp(pars = start, fun = loglik, LB = c(-Inf, 0, 0.1, 2.01), 
			UB = c(Inf, Inf, Inf, Inf), control = ctrl, x = x)
	names(fit$pars) = c("mu", "sigma", "skew", "shape")
	return( fit )
}
# ------------------------------------------------------------------------------
# Generalized Hyperbolic Distribution (standard representation)
# ------------------------------------------------------------------------------
.unitroot = function(f, interval, lower = min(interval), upper = max(interval), 
		tol = .Machine$double.eps^0.25, ...)
{   
	if (is.null(args(f))) {
		if (f(lower) * f(upper) >=0) return(NA)
	} else {
		if (f(lower, ...) * f(upper, ...) >= 0) return(NA)
	}
	ans = uniroot(f = f, interval = interval, lower = lower,
			upper = upper, tol = tol, ...)
	return( ans$root )
}
.kappaGH <-function(x, lambda = 1)
{    
	# A function implemented by Diethelm Wuertz
	#   Returns modified Bessel function ratio
	stopifnot(x >= 0)
	stopifnot(length(lambda) == 1)
	if (lambda == -0.5) {
		# NIG:
		kappa = 1/x
	} else {
		# GH:
		kappa = (besselK(x, lambda+1, expon.scaled = TRUE)/besselK(x, lambda, expon.scaled = TRUE) ) / x
	}
	return( kappa )
}
.deltaKappaGH<-function(x, lambda = 1)
{
	return( .kappaGH(x, lambda+1) - .kappaGH(x, lambda) )
}
.paramGH = function(rho = 0, zeta = 1, lambda = 1)
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
	return( c(alpha = alpha, beta = beta, delta = delta, mu = mu) )
}
.rgigjd = function(n, theta)
{	
	# A function implemented by Diethelm Wuertz	
	#	Original Version by David Scott
	lambda = theta[1]
	chi = theta[2]
	psi = theta[3]
	if (chi < 0) stop("chi can not be negative")
	if (psi < 0) stop("psi can not be negative")
	if ((lambda >= 0)&(psi==0)) stop("When lambda >= 0, psi must be > 0")
	if ((lambda <= 0)&(chi==0)) stop("When lambda <= 0, chi must be > 0")
	if (chi == 0) stop("chi = 0, use rgamma")
	if (psi == 0) stop("algorithm only valid for psi > 0")	
	alpha = sqrt(psi/chi)
	beta = sqrt(psi*chi)
	m = (lambda-1+sqrt((lambda-1)^2+beta^2))/beta
	g = function(y){
		0.5*beta*y^3 - y^2*(0.5*beta*m+lambda+1) + y*((lambda-1)*m-0.5*beta) + 0.5*beta*m
	}
	upper = m
	while (g(upper) <= 0) upper = 2*upper
	yM = uniroot(g, interval=c(0,m))$root
	yP = uniroot(g, interval=c(m,upper))$root
	a = (yP-m)*(yP/m)^(0.5*(lambda-1))*exp(-0.25*beta*(yP+1/yP-m-1/m))
	b = (yM-m)*(yM/m)^(0.5*(lambda-1))*exp(-0.25*beta*(yM+1/yM-m-1/m))
	c = -0.25*beta*(m+1/m) + 0.5*(lambda-1)*log(m)	
	output = numeric(n)
	for(i in 1:n){
		need.value = TRUE
		while(need.value==TRUE){
			R1 = runif (1)
			R2 = runif (1)
			Y = m + a*R2/R1 + b*(1-R2)/R1
			if (Y>0){
				if (-log(R1)>=-0.5*(lambda-1)*log(Y)+0.25*beta*(Y+1/Y)+c){
					need.value = FALSE
				}
			}
		}
		output[i] = Y
	}
	return( output/alpha )
}
.rgigjd1 = function(n, theta)
{	
	# A function implemented by Diethelm Wuertz
	# Description:
	# 	Modified version of rgigjd to generate random observations
	# 	from a generalised inverse Gaussian distribution in the
	# 	special case where lambda = 1.	
	#	Original Version by David Scott
	if (length(theta) == 2) theta = c(1, theta)
	lambda = 1
	chi = theta[2]
	psi = theta[3]
	if (chi < 0) stop("chi can not be negative")
	if (psi < 0) stop("psi can not be negative")	
	if (chi == 0) stop("chi = 0, use rgamma")
	if (psi == 0) stop("When lambda >= 0, psi must be > 0")	
	alpha = sqrt(psi/chi)
	beta = sqrt(psi*chi)
	m = abs(beta)/beta
	g = function(y){
		0.5*beta*y^3 - y^2*(0.5*beta*m+lambda+1) +
				y*(-0.5*beta) + 0.5*beta*m
	}
	upper = m
	while (g(upper)<=0) upper = 2*upper
	yM = uniroot(g,interval=c(0,m))$root
	yP = uniroot(g,interval=c(m,upper))$root
	a = (yP-m)*exp(-0.25*beta*(yP+1/yP-m-1/m))
	b = (yM-m)*exp(-0.25*beta*(yM+1/yM-m-1/m))
	c = -0.25*beta*(m+1/m)
	output = numeric(n)
	for(i in 1:n){
		need.value = TRUE
		while(need.value==TRUE){
			R1 = runif (1)
			R2 = runif (1)
			Y = m + a*R2/R1 + b*(1-R2)/R1
			if (Y>0){
				if (-log(R1)>=0.25*beta*(Y+1/Y)+c){
					need.value = FALSE
				}
			}
		}
		output[i] = Y
	}	
	return( output/alpha )
}

dgh = function(x, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1, log = FALSE)
{
	n = c(length(x), length(alpha), length(beta), length(delta), length(mu), length(lambda))
	maxn = max(n)
	if(n[1]!=maxn) x = rep(x[1], maxn)
	if(n[2]!=maxn) alpha = rep(alpha[1], maxn)
	if(n[3]!=maxn) beta  = rep(beta[1], maxn)
	if(n[4]!=maxn) delta = rep(delta[1], maxn)
	if(n[5]!=maxn) mu    = rep(mu[1], maxn)
	if(n[6]!=maxn) lambda= rep(lambda[1], maxn)
	if(any(alpha <= 0)) stop("alpha must be greater than zero")
	if(any(delta <= 0)) stop("delta must be greater than zero")
	if(any(abs(beta) >= alpha)) stop("abs value of beta must be less than alpha")
	ans = double(maxn)
	sol = try(.C("c_dgh", x = as.double(x), alpha = as.double(alpha), 
					beta = as.double(beta), delta = as.double(delta), 
					mu = as.double(mu), lambda = as.double(lambda), 
					ans = ans, n = as.integer(maxn), logr = as.integer(log), 
					PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

dgh_R = function(x, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1, log = FALSE)
{
	# A function implemented by Diethelm Wuertz
	# modified by Alexios Ghalanos
	n = c(length(x), length(alpha), length(beta), length(delta), length(mu), length(lambda))
	maxn = max(n)
	if(n[1]!=maxn) x = rep(x[1], maxn)
	if(n[2]!=maxn) alpha = rep(alpha[1], maxn)
	if(n[3]!=maxn) beta  = rep(beta[1], maxn)
	if(n[4]!=maxn) delta = rep(delta[1], maxn)
	if(n[5]!=maxn) mu    = rep(mu[1], maxn)
	if(n[6]!=maxn) lambda= rep(lambda[1], maxn)
	if(any(alpha <= 0)) stop("alpha must be greater than zero")
	if(any(delta <= 0)) stop("delta must be greater than zero")
	if(any(abs(beta) >= alpha)) stop("abs value of beta must be less than alpha")
	arg = delta*sqrt(alpha^2-beta^2)
	a = (lambda/2)*log(alpha^2-beta^2) - (
				log(sqrt(2*pi)) + (lambda-0.5)*log(alpha) + lambda*log(delta) +
				log(besselK(arg, lambda, expon.scaled = TRUE)) - arg )
	f = ((lambda-0.5)/2)*log(delta^2+(x - mu)^2)
	# Use exponential scaled form to prevent from overflows:
	arg = alpha * sqrt(delta^2+(x-mu)^2)
	k = log(besselK(arg, lambda-0.5, expon.scaled = TRUE)) - arg
	e = beta*(x-mu)
	ans = a + f + k + e
	if(!log) ans = exp(ans)
	return( ans )
}

pgh = function(q, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1)
{
	# A function implemented by Diethelm Wuertz
	# modified by Alexios Ghalanos
	n = c(length(q), length(alpha), length(beta), length(delta), length(mu), length(lambda))
	maxn = max(n)
	if(n[1]!=maxn) q = rep(q[1], maxn)
	if(n[2]!=maxn) alpha = rep(alpha[1], maxn)
	if(n[3]!=maxn) beta  = rep(beta[1], maxn)
	if(n[4]!=maxn) delta = rep(delta[1], maxn)
	if(n[5]!=maxn) mu    = rep(mu[1], maxn)
	if(n[6]!=maxn) lambda= rep(lambda[1], maxn)
	if(any(alpha <= 0)) stop("alpha must be greater than zero")
	if(any(delta <= 0)) stop("delta must be greater than zero")
	if(any(abs(beta) >= alpha)) stop("abs value of beta must be less than alpha")
	ans = rep(NA, maxn)
	for(i in 1:maxn){
		Integral = integrate(dgh, -Inf, q[i], stop.on.error = FALSE,
				alpha = alpha[i], beta = beta[i], delta = delta[i], mu = mu[i],
				lambda = lambda[i])
		ans[i] = as.numeric(unlist(Integral)[1])
	}
	return( ans )
}

qgh = function(p, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1)
{
	# A function implemented by Diethelm Wuertz
	# modified by Alexios Ghalanos
	n = c(length(p), length(alpha), length(beta), length(delta), length(mu), length(lambda))
	maxn = max(n)
	if(n[1]!=maxn) p = rep(p[1], maxn)
	if(n[2]!=maxn) alpha = rep(alpha[1], maxn)
	if(n[3]!=maxn) beta  = rep(beta[1], maxn)
	if(n[4]!=maxn) delta = rep(delta[1], maxn)
	if(n[5]!=maxn) mu    = rep(mu[1], maxn)
	if(n[6]!=maxn) lambda= rep(lambda[1], maxn)
	if(any(alpha <= 0)) stop("alpha must be greater than zero")
	if(any(delta <= 0)) stop("delta must be greater than zero")
	if(any(abs(beta) >= alpha)) stop("abs value of beta must be less than alpha")
	# Internal Function:
	.froot <- function(x, alpha, beta, delta, mu, lambda, p)
	{
		pgh(q = x, alpha = alpha, beta = beta, delta = delta,
				mu = mu, lambda = lambda) - p
	}
	# Quantiles:
	ans = rep(NA, maxn)
	for(i in 1:maxn){
		lower = -1
		upper = +1
		counter = 0
		iteration = NA
		while(is.na(iteration)){
			iteration = .unitroot(f = .froot, interval = c(lower,
							upper), alpha = alpha[i], beta = beta[i], delta = delta[i],
					mu = mu[i], lambda = lambda[i], p = p[i])
			counter = counter + 1
			lower = lower - 2^counter
			upper = upper + 2^counter
		}
		ans[i] = iteration
	}
	return( ans )
}

ghFit = function(x, control = list())
{
	ctrl = .solnpctrl(control)
	
	start = c(alpha = 1, beta = 0, delta = sqrt(var(x)), mu = mean(x), lambda = -2)
	eps = .Machine$double.eps^0.5
	BIG = 10000
	# Log-likelihood Function:
	loglik = function(pars, x){
		if(abs(pars[2])>pars[1]) return(100000)
		f =  -sum(log(dgh(x, pars[1], pars[2], pars[3], pars[4], pars[5], log = FALSE)))
		f 
	}
	con = function(pars, x){
		ans = abs(pars[2]) - pars[1]
		ans
	}
	fit = solnp(pars = start, fun = loglik, LB = c(eps, -BIG, eps, -BIG, -6), 
			UB = c(BIG,  BIG, BIG,  BIG,  6), ineqfun = con, ineqLB = -BIG, 
			ineqUB = -0.1, control = ctrl, x = x)
	names(fit$pars) = c("alpha", "beta", "delta", "mu", "lambda")
	return( fit )
	
}

rgh = function(n, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1)
{	# A function implemented by Diethelm Wuertz
	#	Original Version by David Scott
	# modified: Alexios Ghalanos
	nn = c(length(alpha), length(beta), length(delta), length(mu), length(lambda))
	maxn = n
	if(nn[1]!=maxn) alpha = rep(alpha[1], maxn)
	if(nn[2]!=maxn) beta  = rep(beta[1], maxn)
	if(nn[3]!=maxn) delta = rep(delta[1], maxn)
	if(nn[4]!=maxn) mu    = rep(mu[1], maxn)
	if(nn[5]!=maxn) lambda= rep(lambda[1], maxn)
	chi = delta^2
	psi = alpha^2 - beta^2
	if(any(alpha <= 0)) stop("alpha must be greater than zero")
	if(any(delta <= 0)) stop("delta must be greater than zero")
	if(any(abs(beta) >= alpha)) stop("abs value of beta must be less than alpha")
	V = cbind(lambda, chi, psi)
	X = apply(V, 1, function(x){
				ifelse(x[1] == 1, 
						.rgigjd1(1, c(x[1], x[2], x[3])), 
						.rgigjd(1, c(x[1], x[2], x[3])))
			})
	sigma = sqrt(as.numeric(X))
	Z = rnorm(n)
	Y = mu + beta*sigma^2 + sigma*Z
	return( Y )
}

# ------------------------------------------------------------------------------
# Standardized Generalized Hyperbolic Distribution
# ------------------------------------------------------------------------------
dsgh = function(x, mu = 0, sigma = 1, skew = 0, shape = 1, lambda = 1, log = FALSE)
{
	n = c(length(x), length(mu), length(sigma), length(skew), length(shape), length(lambda))
	maxn = max(n)
	if(n[1]!=maxn) x = rep(x[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) skew  = rep(skew[1], maxn)
	if(n[5]!=maxn) shape = rep(shape[1], maxn)
	if(n[6]!=maxn) lambda= rep(lambda[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_dghyp", x = as.double(x), mu = as.double(mu), 
					sigma = as.double(sigma), skew = as.double(skew), 
					shape = as.double(shape), lambda = as.double(lambda),
					ans = ans, n = as.integer(maxn), logr = as.integer(log), 
					PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

#dsgh_R = function(x, mu = 0, sigma = 1, zeta = 1, rho = 0, lambda = 1, log = FALSE)
#{
#	fun = function(x, mu = 0, sigma = 1, zeta = 1, rho = 0, lambda = 1, log = FALSE){
#		param = .paramGH(rho, zeta, lambda)
#		return( .dgh(x, param[1]/sigma, param[2]/sigma, param[3]*sigma, param[4]*sigma + mu, lambda, log) )
#	}
#	f = Vectorize( fun )
#	ans = f(x, mu, sigma, zeta, rho, lambda, log)
#	if(NCOL(ans)==1) ans = as.numeric(ans)
#	return( ans )	
#}

#system.time(dsgh_R(x, sigma=1.5, mu=2, lambda=-0.5))
#user  system elapsed 
#2.12    0.00    2.12 
# system.time(dsgh(x, sigma=1.5, mu=2, lambda=-0.5))
#user  system elapsed 
#0.03    0.00    0.03

psgh = function(q, mu = 0, sigma = 1, skew = 0, shape = 1, lambda = 1) 
{
	n = c(length(q), length(mu), length(sigma), length(skew), length(shape), length(lambda))
	maxn = max(n)
	if(n[1]!=maxn) q = rep(q[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) skew  = rep(skew[1], maxn)
	if(n[5]!=maxn) shape = rep(shape[1], maxn)
	if(n[6]!=maxn) lambda= rep(lambda[1], maxn)
	tmp = t(apply(cbind(skew, shape, lambda), 1, function(x) .paramGH(x[1], x[2], x[3])))
	ans = pgh((q-mu)/sigma, tmp[,1], tmp[,2], tmp[,3],tmp[,4], lambda)
	# equivalent: pgh(q, alpha/sigma, beta/sigma, delta*sigma, mu*sigma+mu, lambda)
	return( as.numeric( ans ) )
}

qsgh = function(p, mu = 0, sigma = 1, skew = 0, shape = 1, lambda = 1) 
{
	n = c(length(p), length(mu), length(sigma), length(skew), length(shape), length(lambda))
	maxn = max(n)
	if(n[1]!=maxn) p = rep(p[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) skew  = rep(skew[1], maxn)
	if(n[5]!=maxn) shape = rep(shape[1], maxn)
	if(n[6]!=maxn) lambda= rep(lambda[1], maxn)
	tmp = t(apply(cbind(skew, shape, lambda), 1, function(x) .paramGH(x[1], x[2], x[3])))
	ans = mu + sigma*as.numeric( qgh(p, tmp[,1], tmp[,2], tmp[,3], tmp[,4], lambda) )
	return( ans )
}

rsgh = function(n, mu = 0, sigma = 1, skew = 0, shape = 1, lambda = 1)
{
	nn = c(length(mu), length(sigma), length(skew), length(shape), length(lambda))
	if(nn[1]!=n) mu    = rep(mu[1], n)
	if(nn[2]!=n) sigma = rep(sigma[1], n)
	if(nn[3]!=n) skew  = rep(skew[1], n)
	if(nn[4]!=n) shape = rep(shape[1], n)
	if(nn[5]!=n) lambda= rep(lambda[1], n)
	ans = double(n)
	sol = try(.C("c_rghyp", n = as.integer(n), mu = as.double(mu), 
					sigma = as.double(sigma), skew = as.double(skew), 
					shape = as.double(shape), lambda = as.double(lambda),
					ans = ans, PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

rsgh_R = function(n, mu = 0, sigma = 1, skew = 0, shape = 1, lambda = 1) 
{
	nn = c(length(mu), length(sigma), length(skew), length(shape), length(lambda))
	if(nn[1]!=n) mu    = rep(mu[1], n)
	if(nn[2]!=n) sigma = rep(sigma[1], n)
	if(nn[3]!=n) skew  = rep(skew[1], n)
	if(nn[4]!=n) shape = rep(shape[1], n)
	if(nn[5]!=n) lambda= rep(lambda[1], n)
	tmp = t(apply(cbind(skew, shape, lambda), 1, function(x) .paramGH(x[1], x[2], x[3])))
	ans = mu + sigma*as.numeric( rgh(n, tmp[,1], tmp[,2], tmp[,3], tmp[,4], lambda) )
	return(ans)
}

# system.time(rsgh(1000, mu=0.1, sigma = 0.2))
# user  system elapsed 
# 0.02    0.00    0.02 
# system.time(rsgh_R(1000, mu=0.1, sigma = 0.2))
# user  system elapsed 
# 0.69    0.00    0.69 

sghFit = function(x, control=list(), ...) 
{
	x = as.vector(x)
	ctrl = .solnpctrl(control)
	start = c(mu = mean(x), sigma = sd(x), skew = 0, shape = 0.5, lambda = -2)
	loglik = function(pars, x){
		f =  sum(-log(dsgh(x, pars[1], pars[2], pars[3], pars[4], pars[5], log = FALSE)))
		f
	}
	fit = solnp(pars = start, fun = loglik, LB = c(-200, 1e-12, -0.999, 0.1, -4), 
			UB = c(200, 400, 0.999, 25, 4), control = ctrl, ..., x = x)
	names(fit$pars) <- c("mu", "sigma", "skew","shape","lambda")
	return(fit)
}
# ------------------------------------------------------------------------------
# Normal Inverse Gaussian (NIG) Distribution
# ------------------------------------------------------------------------------
.qnigC = function(p, alpha = 1, beta = 0, delta = 1, mu = 0)
{   
	if(alpha <= 0) stop("Invalid parameters: alpha <= 0.\n")
	if(alpha^2 <= beta^2) stop("Invalid parameters: alpha^2 <= beta^2.\n")
	if(delta <= 0) stop("Invalid parameters: delta <= 0.\n")
	if((sum(is.na(p)) > 0)) 
		stop("Invalid probabilities:\n",p,"\n")
	else 
	if(sum(p < 0)+sum(p > 1) > 0) stop("Invalid probabilities:\n",p,"\n")
	
	n <- length(p)
	q <- rep(0, n)
	retValues <- .C("qNIG",
			p = as.double(.CArrange(p,1,1,n)),
			i_mu = as.double(mu),
			i_delta = as.double(delta),
			i_alpha = as.double(alpha),
			i_beta = as.double(beta),
			i_n = as.integer(n),
			q = as.double(.CArrange(q, 1, 1, n)), PACKAGE="rugarch")
	quantiles <- retValues[[7]]
	quantiles[quantiles <= -1.78e+308] <- -Inf
	quantiles[quantiles >= 1.78e+308] <- Inf
	return( quantiles )
}
.CArrange <-function(obj, i, j, n)
{
	# Description:
	#   Arrange input matrices and vectors in a suitable way for the C program
	#   Matrices are transposed because the C program stores matrices by row 
	#   while R stores matrices by column
	
	# Arguments:
	#   i - length of first dimension
	#   j - length of second dimension
	#   n - length of third dimension
	
	# Value:
	#   out - transformed data set
	
	# Author: 
	#   Daniel Berg <daniel at nr.no> (Kjersti Aas <Kjersti.Aas at nr.no>)
	#   Date: 12 May 2005
	#   Version: 1.0.2	
	if(is.null(obj)) stop("Missing data")
	
	if(is.vector(obj)) {
		if(i==1 & j==1 & length(obj)==n) out <- as.double(obj)
		else stop("Unexpected length of vector")
	} else if(is.matrix(obj)) {
		if(nrow(obj) == i && ncol(obj) == j) out <- as.double(rep(t(obj), n))
		else stop("Unexpected dimensions of matrix")
	} else {
		stop("Unexpected object")
	}	
	return(out) 
}

dnig = function(x, alpha = 1, beta = 0, delta = 1, mu = 0, log = FALSE)
{   
	return( dgh(x = x, alpha = alpha, beta = beta, delta = delta, mu = mu, 
					lambda = -0.5, log = log) )
}
pnig = function(q, alpha = 1, beta = 0, delta = 1, mu = 0)
{   
	return( pgh(q = q, alpha = alpha, beta = beta, delta = delta, mu = mu, 
					lambda = -0.5) )
}
qnig = function(p, alpha = 1, beta = 0, delta = 1, mu = 0)
{
	n = c(length(p), length(alpha), length(beta), length(delta), length(mu))
	maxn = max(n)
	if(n[1]!=maxn) p = rep(p[1], maxn)
	if(n[2]!=maxn) alpha = rep(alpha[1], maxn)
	if(n[3]!=maxn) beta  = rep(beta[1], maxn)
	if(n[4]!=maxn) delta = rep(delta[1], maxn)
	if(n[5]!=maxn) mu    = rep(mu[1], maxn)
	ans = rep(NA, maxn)
	for(i in 1:maxn){
		ans[i] = .qnigC(p = p[i], alpha = alpha[i], beta = beta[i], delta = delta[i], mu = mu[i])
	}
	return( ans )
}
rnig = function(n, alpha = 1, beta = 0, delta = 1, mu = 0)
{   
	return( rgh(n, alpha = alpha, beta = beta, delta = delta, mu = mu, lambda=-0.5) )
}
nigFit = function(x, control = list())
{
	ctrl = .solnpctrl(control)
	start = c(alpha = 1, beta = 0, delta = sqrt(var(x)), mu = mean(x))
	loglik = function(pars, x){
		f = -sum(log(dnig(x, pars[1], pars[2], pars[3], pars[4])))
		f }
	con = function(pars, x){
		abs(pars[2])-pars[1]
	}
	fit = solnp(pars = start, fun = loglik, LB = c(0.1, -100, eps, -Inf), 
			UB = c(100, 100, Inf, Inf), ineqfun = con, ineqLB = -200, 
			ineqUB = 0,control = ctrl, x = x)
	names(fit$pars) = c("alpha", "beta", "delta", "mu")
	return( fit )
}

# ------------------------------------------------------------------------------
# Standardized Normal Inverse Gaussian (NIG) Distribution
# ------------------------------------------------------------------------------
.qsnigC = function(p, rho = 0, zeta = 1) 
{
	param = .paramGH(rho, zeta, lambda = -0.5)
	return( .qnigC(p, param[1], param[2], param[3], param[4]) )
}

dsnig = function(x, mu = 0, sigma = 1, skew = 0, shape = 1, log = FALSE) 
{
	n = c(length(x), length(mu), length(sigma), length(skew), length(shape))
	maxn = max(n)
	if(n[1]!=maxn) x = rep(x[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) skew  = rep(skew[1], maxn)
	if(n[5]!=maxn) shape = rep(shape[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_dsnig", x = as.double(x), mu = as.double(mu), 
					sigma = as.double(sigma), skew = as.double(skew), 
					shape = as.double(shape), ans = ans, n = as.integer(maxn), 
					logr = as.integer(log), PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}
psnig = function(q, mu = 0, sigma = 1, skew = 0, shape = 1) 
{
	return( psgh(q, mu, sigma, skew, shape, lambda = -0.5) )
}
qsnig = function(p, mu = 0, sigma = 1, skew = 0, shape = 1) 
{
	n = c(length(p), length(mu), length(sigma), length(skew), length(shape))
	maxn = max(n)
	if(n[1]!=maxn) p = rep(p[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) skew  = rep(skew[1], maxn)
	if(n[5]!=maxn) shape = rep(shape[1], maxn)
	ans = double(maxn)	
	for(i in 1:maxn){
		ans[i] = mu[i] + sigma[i]*.qsnigC(p, rho = skew[i], zeta = shape[i])
	}
	return( ans )
}
rsnig = function(n, mu = 0, sigma = 1, skew = 0, shape = 1)
{
	nn = c(length(mu), length(sigma), length(skew), length(shape))
	if(nn[1]!=n) mu    = rep(mu[1], n)
	if(nn[2]!=n) sigma = rep(sigma[1], n)
	if(nn[3]!=n) skew  = rep(skew[1], n)
	if(nn[4]!=n) shape = rep(shape[1], n)
	ans = double(n)
	sol = try(.C("c_rsnig", n = as.integer(n), mu = as.double(mu), 
					sigma = as.double(sigma), skew = as.double(skew), 
					shape = as.double(shape),
					ans = ans, PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

snigFit = function (x, control=list()) 
{   	
	x = as.vector(x)
	ctrl = .solnpctrl(control)
	start = c(mu = mean(x), sigma = sd(x), skew = 0, shape = 1)
	loglik = function(pars, x){
		f =  sum(-log(dsnig(x, pars[1], pars[2], pars[3], pars[4], log = FALSE)))
		f
	}
	fit = solnp(pars = start, fun = loglik, LB = c(-200, 1e-12, -0.999, 0.1), 
			UB = c(200, 400, 0.999, 25), control = ctrl, x = x)
	names(fit$pars) <- c("mu", "sigma", "skew","shape")
	return(fit)
}
# ------------------------------------------------------------------------------
# Johnson's SU Distribution (Rigby & Stasinopoulos parameterization)
# coded into C by Alexios Ghalanos
# ------------------------------------------------------------------------------
djsu = function(x, mu = 0, sigma = 1, skew = 1, shape = 0.5, log = FALSE)
{
	n = c(length(x), length(mu), length(sigma), length(skew), length(shape))
	maxn = max(n)
	if(n[1]!=maxn) x = rep(x[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) skew  = rep(skew[1], maxn)
	if(n[5]!=maxn) shape = rep(shape[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_djsu", x = as.double(x), mu = as.double(mu), 
					sigma = as.double(sigma), skew = as.double(skew), 
					shape = as.double(shape), ans = ans, n = as.integer(maxn), 
					logr = as.integer(log), PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

pjsu = function(q, mu = 0, sigma = 1, skew = 1, shape = 0.5)
{
	n = c(length(q), length(mu), length(sigma), length(skew), length(shape))
	maxn = max(n)
	if(n[1]!=maxn) q = rep(q[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) skew  = rep(skew[1], maxn)
	if(n[5]!=maxn) shape = rep(shape[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_pjsu", q = as.double(q), mu = as.double(mu), 
					sigma = as.double(sigma), skew = as.double(skew), 
					shape = as.double(shape), ans = ans, n = as.integer(maxn),
					PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

qjsu = function(p, mu = 0, sigma = 1, skew = 1, shape = 0.5)
{
	n = c(length(p), length(mu), length(sigma), length(skew), length(shape))
	maxn = max(n)
	if(n[1]!=maxn) p = rep(p[1], maxn)
	if(n[2]!=maxn) mu    = rep(mu[1], maxn)
	if(n[3]!=maxn) sigma = rep(sigma[1], maxn)
	if(n[4]!=maxn) skew  = rep(skew[1], maxn)
	if(n[5]!=maxn) shape = rep(shape[1], maxn)
	ans = double(maxn)
	sol = try(.C("c_qjsu", p = as.double(p), mu = as.double(mu), 
					sigma = as.double(sigma), skew = as.double(skew), 
					shape = as.double(shape), ans = ans, n = as.integer(maxn),  
					PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}
rjsu = function(n, mu = 0, sigma = 1, skew = 1, shape = 0.5)
{
	nn = c(length(mu), length(sigma), length(skew), length(shape))
	if(nn[1]!=n) mu    = rep(mu[1], n)
	if(nn[2]!=n) sigma = rep(sigma[1], n)
	if(nn[3]!=n) skew  = rep(skew[1], n)
	if(nn[4]!=n) shape = rep(shape[1], n)
	ans = double(n)
	sol = try(.C("c_rjsu", n = as.integer(n), mu = as.double(mu), 
					sigma = as.double(sigma), skew = as.double(skew), 
					shape = as.double(shape),
					ans = ans, PACKAGE="rugarch"), silent = TRUE)
	if(inherits(sol, 'try-error')){
		return(sol)
	} else{
		return(sol$ans)
	}
}

jsuFit = function(x, control = list())
{
	# a function implemented by Alexios Ghalanos
	ctrl = .solnpctrl(control)
	start = c(mean(x), sd(x), skew = 1, shape = 0.5)
	loglik = function(pars, x){
		-sum(log(djsu(x, pars[1], pars[2], pars[3], pars[4], log = FALSE)))
	}
	fit = solnp(pars = start, fun = loglik, LB = c(-Inf, 0, -20, 0.1),
			UB = c(Inf, Inf, 20, 100), control = ctrl, x = x)
	names(fit$pars) = c("mu", "sigma", "skew","shape")
	return( fit )
}
# ------------------------------------------------------------------------------
# Distribution Bounds
# ------------------------------------------------------------------------------
.DistributionBounds = function(distribution)
{
	ghlambda = 0
	ghlambda.LB = 0
	ghlambda.UB = 0
	if (distribution == "norm"){
		skew 	= 0
		skew.LB = 0
		skew.UB = 0
		shape 	= 0
		shape.LB = 0
		shape.UB = 0}
	if (distribution == "ged"){
		skew 	= 0
		skew.LB = 0
		skew.UB = 10
		shape 	= 2
		shape.LB = 0.1
		shape.UB = 50}
	if (distribution == "std"){
		skew 	= 0
		skew.LB = 0
		skew.UB = 0
		shape 	= 4
		shape.LB = 2.1
		shape.UB = 100}
	if (distribution == "snorm"){
		skew 	= 0.9
		skew.LB	= 0.1
		skew.UB	= 10
		shape 	= 0
		shape.LB = 0
		shape.UB = 0}
	if (distribution == "sged"){
		skew 	= 1
		skew.LB	= 0.01
		skew.UB	= 30
		shape 	= 2
		shape.LB = 0.1
		shape.UB = 60}
	if (distribution == "sstd"){
		skew 	= 1
		skew.LB = 0.01
		skew.UB = 30
		shape 	= 4
		shape.LB = 2.01
		shape.UB = 60}
	if (distribution == "nig"){
		skew 	= 0.2
		skew.LB = -0.99
		skew.UB	= 0.99
		shape 	= 0.4
		shape.LB = 0.01
		shape.UB = 25
		}
	if(distribution == "ghyp"){
		skew 	= 0.2
		skew.LB = -0.99
		skew.UB	= 0.99
		shape 	= 2
		shape.LB = 0.25
		shape.UB = 25
		ghlambda = -0.5
		ghlambda.LB = -6
		ghlambda.UB = 6
	}
	if(distribution == "jsu"){
		skew 	= 0
		skew.LB	= -20
		skew.UB	= 20
		shape 	= 1
		shape.LB = 0.1
		shape.UB = 10
	}
	if(distribution == "ghst"){
		skew 	= 0
		skew.LB	= -80
		skew.UB	= 80
		shape 	= 8.2
		shape.LB = 4.1
		shape.UB = 25
	}
	# johnson has 2 shape parameters. The second one we model with the "skew"
	# representation in rugarch
	skewed.dists = c("snorm", "sged", "sstd", "nig", "ghyp", "jsu", "ghst")
	shaped.dists = c("ged", "sged", "std", "sstd", "nig", "ghyp", "jsu", "ghst")
	skew0  = 0
	shape0 = 0
	if(any(skewed.dists == distribution)) include.skew=TRUE else include.skew=FALSE
	if(any(shaped.dists == distribution)) include.shape=TRUE else include.shape=FALSE
	if(distribution == "ghyp") include.ghlambda = TRUE else include.ghlambda = FALSE
	return( list(shape = shape, shape.LB = shape.LB, shape.UB = shape.UB, skew = skew,
					skew.LB = skew.LB, skew.UB = skew.UB, include.skew = include.skew, 
					include.shape = include.shape, skew0 = skew0, shape0 = shape0,
					include.ghlambda = include.ghlambda, ghlambda = ghlambda, 
					ghlambda.LB = ghlambda.LB, ghlambda.UB = ghlambda.UB) )
}

# ------------------------------------------------------------------------------
# Distribution Wrapper Functions
# ------------------------------------------------------------------------------
.makeSample = function(distribution, lambda = -0.5, skew, shape, n, seed)
{
	set.seed(seed)
	x = switch(distribution,
			norm = rnorm(n),
			snorm = rsnorm(n, skew = skew),
			std = rstd(n, shape = shape),
			sstd = rsstd(n, shape = shape, skew = skew),
			ged = rged(n, shape = shape),
			sged = rsged(n, shape = shape, skew = skew),
			nig = rsnig(n, shape = shape, skew = skew),
			ghyp = rsgh(n, shape = shape, skew = skew, lambda = lambda),
			ghst = rsghst(n, skew = skew, shape = shape),
			jsu = rjsu(n, skew = skew, shape = shape)
			)
	return(x)
}


# Scaling Transformation
.scaledist = function(dist, mu, sigma, lambda = -0.5, skew, shape)
{
	ans = switch(dist,
			norm  = .normscale(mu, sigma),
			snorm = .snormscale(mu, sigma, skew),
			std   = .stdscale(mu, sigma, shape),
			sstd  = .sstdscale(mu, sigma, skew, shape),
			nig   = .nigscale(mu, sigma, skew, shape),
			ghyp  = .ghypscale(mu, sigma, skew, shape, lambda),
			ged   = .gedscale(mu, sigma, shape),
			sged  = .sgedscale(mu, sigma, skew, shape),
			jsu   = .jsuscale(mu, sigma, skew, shape),
			ghst  = .ghstscale(mu, sigma, skew, shape)
			)
	return(ans)
}

# returns the time varying density parameters of the actual returns
# (rescaled from the (0,1) parametrization)
.nigtransform = function(rho, zeta)
{
	nigpars = t(apply(cbind(rho, zeta), 1, FUN = function(x) .paramGH(rho = x[1], zeta = x[2], lambda = -0.5)))
	colnames(nigpars) = c("alpha", "beta", "delta", "mu")
	return(nigpars)
}

.ghyptransform = function(rho, zeta, lambda)
{
	n = length(zeta)
	ghyppars = t(apply(cbind(rho, zeta), 1, FUN = function(x) .paramGH(rho = x[1], zeta = x[2], lambda = lambda)))
	ghyppars = cbind(ghyppars, rep(lambda, n))
	colnames(ghyppars) = c("alpha", "beta", "delta", "mu", "lambda")
	return(ghyppars)
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
	# technically: mu, delta, beta, alpha
	colnames(xdensity) = c("mu", "delta", "beta", "alpha")
	return(xdensity)
}

.ghstscale = function(mu, sigma, skew, shape)
{
	ghpars = t(apply(cbind(skew, shape), 1, FUN=function(x) .paramGHST(betabar = x[1], nu = x[2])))
	xdensity = matrix(0, ncol = 5, nrow = length(sigma))
	# alpha, beta, delta, mu
	xdensity[,5] = -shape/2
	xdensity[,4] = (abs(ghpars[,3])+1e-12)/sigma
	xdensity[,3] = ghpars[,3]/sigma
	xdensity[,2] = ghpars[,2]*sigma
	xdensity[,1] = ghpars[,1]*sigma + mu
	# technically: mu, delta, beta, alpha
	colnames(xdensity) = c("mu", "delta", "beta", "alpha", "lambda")
	return(xdensity)
}

.ghypscale = function(mu, sigma, skew, shape, lambda)
{
	if(length(lambda)>1){
		ghpars = t(apply(cbind(skew, shape, lambda), 1, FUN=function(x) .paramGH(rho = x[1], zeta = x[2], lambda = x[3])))
	} else{
		ghpars = t(apply(cbind(skew, shape), 1, FUN=function(x) .paramGH(rho = x[1], zeta = x[2], lambda = lambda)))
	}
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	# alpha, beta, delta, mu
	xdensity[,4] = ghpars[,1]/sigma
	xdensity[,3] = ghpars[,2]/sigma
	xdensity[,2] = ghpars[,3]*sigma
	xdensity[,1] = ghpars[,4]*sigma + mu
	# technically: mu, delta, beta, alpha
	colnames(xdensity) =c("mu", "delta", "beta", "alpha")
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

.snormscale = function(mu, sigma, skew)
{
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	xdensity[,1] = mu
	xdensity[,2] = sigma
	xdensity[,3] = skew
	xdensity[,4] = 0
	colnames(xdensity) = c("mu", "sigma", "skew", "shape")
	return(xdensity)
}

.stdscale = function(mu, sigma, shape)
{
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	xdensity[,1] = mu
	xdensity[,2] = sigma
	xdensity[,3] = 0
	xdensity[,4] = shape
	colnames(xdensity) = c("mu", "sigma", "skew", "shape")
	return(xdensity)
}

.sstdscale = function(mu, sigma, skew, shape)
{
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	xdensity[,1] = mu
	xdensity[,2] = sigma
	xdensity[,3] = skew
	xdensity[,4] = shape
	colnames(xdensity) = c("mu", "sigma", "skew", "shape")
	return(xdensity)
}


.gedscale = function(mu, sigma, shape)
{
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	xdensity[,1] = mu
	xdensity[,2] = sigma
	xdensity[,3] = 0
	xdensity[,4] = shape
	colnames(xdensity) = c("mu", "sigma", "skew", "shape")
	return(xdensity)
}

.sgedscale = function(mu, sigma, skew, shape)
{
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	xdensity[,1] = mu
	xdensity[,2] = sigma
	xdensity[,3] = skew
	xdensity[,4] = shape
	colnames(xdensity) = c("mu", "sigma", "skew", "shape")
	return(xdensity)
}


.jsuscale = function(mu, sigma, skew, shape)
{
	xdensity = matrix(0, ncol = 4, nrow = length(sigma))
	xdensity[,1] = mu
	xdensity[,2] = sigma
	xdensity[,3] = skew
	xdensity[,4] = shape
	colnames(xdensity) = c("mu", "sigma", "skew", "shape")
	return(xdensity)
}

#---------------------------------------------------------------------------------
# functions for export:
#---------------------------------------------------------------------------------
nigtransform = function(mu = 0, sigma = 1,  skew = 0, shape = 3)
{
	return(.scaledist(dist = "nig", mu = mu, sigma = sigma, skew = skew,
					shape = shape, lambda = -0.5))
}

ghyptransform = function(mu = 0, sigma = 1,  skew = 0, shape = 3, lambda = -0.5)
{
	return(.scaledist(dist = "ghyp", mu = mu, sigma = sigma, skew = skew,
					shape = shape, lambda = lambda))
}


fitdist = function(distribution = "norm", x, control=list()){
	valid.distributions = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig",
			"ghyp", "jsu", "ghst")
	if(!any(valid.distributions == distribution))
		stop("\nnot a valid distributions\n", call. = FALSE)
	ans = switch(distribution,
			norm = normFit(x, control),
			snorm = snormFit(x, control),
			std = stdFit(x, control),
			sstd = sstdFit(x, control),
			ged = gedFit(x, control),
			sged = sgedFit(x, control),
			nig = snigFit(x, control),
			ghyp = sghFit(x, control),
			jsu = jsuFit(x, control),
			ghst = ghstFit(x, control)
	)
	return(ans)
}
ddist = function(distribution = "norm", y, mu = 0, sigma = 1, lambda = -0.5, skew = 1, shape = 5)
{
	valid.distributions = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig",
			"ghyp", "jsu", "ghst")
	if(!any(valid.distributions == distribution))
		stop("\nnot a valid distributions\n", call. = FALSE)
	ans = switch(distribution,
		norm = dnorm(y, mean = mu, sd = sigma),
		snorm = dsnorm(y, mu = mu, sigma = sigma, skew = skew),
		std = dstd(y, mu = mu, sigma = sigma, shape = shape),
		sstd = dsstd(y, mu = mu, sigma = sigma, skew = skew, shape = shape),
		ged = dged(y, mu = mu, sigma = sigma, shape = shape),
		sged = dsged(y, mu = mu, sigma = sigma, skew = skew, shape = shape),
		nig = dsnig((y-mu)/sigma, skew = skew, shape = shape)/sigma,
		ghyp = dsgh((y-mu)/sigma, skew = skew, shape = shape, lambda = lambda)/sigma,
		jsu = djsu(y, mu = mu, sigma = sigma, skew = skew, shape = shape),
		ghst = dsghst(y, mu = mu, sigma = sigma, skew = skew, shape = shape)
		)
	return(ans)
}

pdist = function(distribution = "norm", q, mu = 0, sigma = 1, lambda = -0.5, skew = 1, shape = 5)
{
	valid.distributions = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig",
			"ghyp", "jsu","ghst")
	if(!any(valid.distributions == distribution))
		stop("\nnot a valid distributions\n", call. = FALSE)
	ans = switch(distribution,
			norm = pnorm(q, mean = mu, sd = sigma, log.p = FALSE),
			snorm = psnorm(q, mu = mu, sigma = sigma, skew = skew),
			std = pstd(q, mu = mu, sigma = sigma, shape = shape),
			sstd = psstd(q, mu = mu, sigma = sigma, skew = skew, shape = shape),
			ged = pged(q, mu = mu, sigma = sigma, shape = shape),
			sged = psged(q, mu = mu, sigma = sigma, skew = skew, shape = shape),
			nig = psnig((q-mu)/sigma, skew = skew, shape = shape),
			ghyp = psgh((q-mu)/sigma, skew = skew, shape = shape, lambda = lambda),
			jsu = pjsu(q, mu = mu, sigma = sigma, skew = skew, shape = shape),
			ghst = psghst(q, mu = mu, sigma = sigma, skew = skew, shape = shape)
	)
	return(ans)
}

qdist = function(distribution = "norm", p, mu = 0, sigma = 1, lambda = -0.5, skew = 1, shape = 5)
{
	valid.distributions = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig",
			"ghyp", "jsu","ghst")
	if(!any(valid.distributions == distribution))
		stop("\nnot a valid distributions\n", call. = FALSE)
	ans = switch(distribution,
			norm = qnorm(p, mean = mu, sd = sigma, log.p = FALSE),
			snorm = qsnorm(p, mu = mu, sigma = sigma, skew = skew),
			std = qstd(p, mu = mu, sigma = sigma, shape = shape),
			sstd = qsstd(p, mu = mu, sigma = sigma, skew = skew, shape = shape),
			ged = qged(p, mu = mu, sigma = sigma, shape = shape),
			sged = qsged(p, mu = mu, sigma = sigma, skew = skew, shape = shape),
			nig = qsnig(p, skew = skew, shape = shape)*sigma + mu,
			ghyp = qsgh(p, skew = skew, shape = shape, lambda = lambda)*sigma + mu,
			jsu = qjsu(p, mu = mu, sigma = sigma, skew = skew, shape = shape),
			ghst = qsghst(p, mu = mu, sigma = sigma, skew = skew, shape = shape),
	)
	return(ans)
}

# EQUAL:
# set.seed(10)
# rsnig(10, rho = skew, zeta = shape)*sigma + mu
# set.seed(10)
# rdist("nig", 10, mu = mu, sigma = sigma, skew = skew, shape = shape)
rdist = function(distribution = "norm", n, mu = 0, sigma = 1, lambda = -0.5, skew = 1, shape = 5)
{
	valid.distributions = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig",
			"ghyp", "jsu","ghst")
	if(!any(valid.distributions == distribution))
		stop("\nnot a valid distribution\n", call. = FALSE)
	ans = switch(distribution,
			norm = rnorm(n, mean = mu, sd = sigma),
			snorm = rsnorm(n, mu = mu, sigma = sigma, skew = skew),
			std = rstd(n, mu = mu, sigma = sigma, shape = shape),
			sstd = rsstd(n, mu = mu, sigma = sigma, skew = skew, shape = shape),
			ged = rged(n, mu = mu, sigma = sigma, shape = shape),
			sged = rsged(n, mu = mu, sigma = sigma, skew = skew, shape = shape),
			nig =  mu + sigma*rsnig(n, skew = skew, shape = shape),
			ghyp = mu + sigma*rsgh(n, skew = skew, shape = shape, lambda = lambda),
			jsu = rjsu(n, mu = mu, sigma = sigma, skew = skew, shape = shape),
			ghst = rsghst(n, mu = mu, sigma = sigma, skew = skew, shape = shape)
	)
	return(ans)
}


dskewness = function(distribution = "norm", skew = 1, shape = 5, lambda = -0.5)
{
	f = Vectorize(.dskewness)
	ans = f(distribution, skew, shape, lambda)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return(ans)
}
	
.dskewness = function(distribution = "norm", skew = 1, shape = 5, lambda = -0.5)
{
	valid.distributions = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig", "ghyp", "jsu","ghst")
	if(!any(valid.distributions == distribution))
		stop("\nnot a valid distribution\n", call. = FALSE)
	ans = switch(distribution,
			norm 	= 0,
			snorm 	= .snormskew(skew = skew),
			std 	= 0,
			sstd 	= .sstdskew(skew = skew, shape = shape),
			ged 	= 0,
			sged 	= .sgedskew(skew = skew, shape = shape),
			nig 	= .snigskew(skew = skew, shape = shape),
			ghyp 	= .sghypskew(skew = skew, shape = shape, lambda = lambda),
			jsu 	= .jsuskew(skew = skew, shape = shape),
			ghst	= .ghstskew(skew, shape)
	)
	return(as.numeric(ans))
}

dkurtosis = function(distribution = "norm", skew = 1, shape = 5, lambda = -0.5)
{
	f = Vectorize(.dkurtosis)
	ans = f(distribution, skew, shape, lambda)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return(ans)
}

.dkurtosis = function(distribution = "norm", skew = 1, shape = 5, lambda = -0.5)
{
	valid.distributions = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig", "ghyp", "jsu","ghst")
	if(!any(valid.distributions == distribution))
		stop("\nnot a valid distribution\n", call. = FALSE)
	ans = switch(distribution,
			norm 	= 0,
			snorm 	= 0,
			std 	= .stdexkurt(shape = shape),
			sstd 	= .sstdexkurt(skew = skew, shape = shape),
			ged 	= .gedexkurt(shape = shape),
			sged 	= .sgedexkurt(skew = skew, shape = shape),
			nig 	= .snigexkurt(skew = skew, shape = shape),
			ghyp 	= .sghypexkurt(skew = skew, shape = shape, lambda = lambda),
			jsu 	= .jsuexkurt(skew = skew, shape = shape),
			ghst	= .ghstexkurt(skew, shape)
	)
	return(as.numeric(ans))
}


# NIG Moments
.nigmu = function(alpha, beta, delta, mu){
	gm = sqrt(alpha^2 - beta^2)
	ans = mu + (delta*beta)/gm
	return(ans)
}
.nigsigma = function(alpha, beta, delta, mu){
	gm = sqrt(alpha^2 - beta^2)
	ans = sqrt((delta*alpha^2)/(gm^3))
	return(ans)
}

.snigskew = function(skew, shape){
	fun = function(skew, shape){
		pars = .paramGH(rho = skew, zeta = shape, lambda = -0.5)
		return(.nigskew(alpha = pars[1], beta = pars[2], delta = pars[3], mu = pars[4]))
	}
	f = Vectorize( fun )
	ans = f(skew, shape)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )
}

.nigskew = function(alpha, beta, delta, mu){
	gm = sqrt(alpha^2 - beta^2)
	ans = 3*beta/(alpha*sqrt(delta*gm))
	return(ans)
}

.snigexkurt = function(skew, shape){
	fun = function(skew, shape){
		pars = .paramGH(rho = skew, zeta = shape, lambda = -0.5)
		return(.nigexkurt(alpha = pars[1], beta = pars[2], delta = pars[3], mu = pars[4]))
	}
	f = Vectorize( fun )
	ans = f(skew, shape)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )
}

.nigexkurt = function(alpha, beta, delta, mu){
	gm = sqrt(alpha^2 - beta^2)
	ans = 3*(1+4*(beta^2)/(alpha^2))/(delta*gm)
	return(ans)
}

.ghypmu = function(lambda, alpha, beta, delta, mu){
	gm = sqrt(alpha^2 - beta^2)
	ans = mu + (delta * beta * besselK( delta * gm, lambda + 1) )/( gm * besselK( delta * gm, lambda) ) 
	return(ans)
}

.ghypsigma = function(lambda, alpha, beta, delta, mu){
	gm = sqrt(alpha^2 - beta^2)
	x1 = delta * besselK( delta * gm, lambda + 1) / ( gm * besselK( delta * gm, lambda) )
	x2 = ( (beta^2 * delta^2) / (gm^2) ) * ( ( besselK( delta * gm, lambda + 2) / besselK( delta * gm, lambda) ) 
				- ( besselK( delta * gm, lambda + 1)^2 / besselK( delta * gm, lambda)^2 ) )
	ans = sqrt(x1 + x2)
	return(ans)
}

.sghypskew = function(skew, shape, lambda){
	fun = function(skew, shape, lambda){
		pars = .paramGH(rho = skew, zeta = shape, lambda = lambda)
		return(.ghypskew(lambda = lambda, alpha = pars[1], beta = pars[2], delta = pars[3], mu = pars[4]))
	}
	f = Vectorize( fun )
	ans = f(skew, shape, lambda)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )
}

.ghypskew = function(lambda, alpha, beta, delta, mu){
	skew = ghypMom(3, lambda, alpha, beta, delta, mu, momType = "central")/(.ghypsigma(lambda, alpha, beta, delta, mu)^3)
	return(skew)
}

.sghypexkurt = function(skew, shape, lambda){
	fun = function(skew, shape, lambda){
		pars = .paramGH(rho = skew, zeta = shape, lambda = lambda)
		return(.ghypexkurt(lambda = lambda, alpha = pars[1], beta = pars[2], delta = pars[3], mu = pars[4]))
	}
	f = Vectorize( fun )
	ans = f(skew, shape, lambda)
	if(NCOL(ans)==1) ans = as.numeric(ans)
	return( ans )
}

.ghypexkurt = function(lambda, alpha, beta, delta, mu){
	kurt = ghypMom(4, lambda, alpha, beta, delta, mu, momType = "central")/(.ghypsigma(lambda, alpha, beta, delta, mu)^4) - 3
	return(kurt)
}

.norm2snorm1 = function(mu, sigma, skew)
{
	m1 = 2/sqrt(2 * pi)
	mu + m1 * (skew - 1/skew)*sigma
}

.norm2snorm2 = function(mu, sigma, skew)
{
	m1 = 2/sqrt(2 * pi)
	m2 = 1
	sigx = sqrt((1-m1^2)*(skew^2+1/skew^2) + 2*m1^2 - 1)
	
	sigx*sigma
}

.snormskew = function( skew )
{
	m1 = 2/sqrt(2 * pi)
	m2 = 1
	m3 = 4/sqrt(2 * pi)
	(skew - 1/skew) * ( ( m3 + 2 * m1^3 - 3 * m1 * m2 ) * ( skew^2 + (1/skew^2) ) + 3 * m1 * m2 - 4 * m1^3 )/
			( ( (m2 - m1^2) * (skew^2 + 1/skew^2) + 2 * m1^2 - m2) ^ (3/2) )
}

.snormexkurt = function( skew )
{
	0
}

.stdskew = function( shape )
{
	0
}

.stdexkurt = function( shape )
{
	ifelse(shape > 4, 6/(shape - 4), NA)
}

.sstdskew = function(skew, shape){
	# Theoretical moments based on bijection betweeen Fernandez and Steel verions
	# and Hansen's Generalized Skew-T (credit due to Michael Rockinger)
	if (shape > 2) {
		eta  = shape
		k2   = skew^2
		lda  = (k2-1)/(k2+1)
		ep1 = (eta+1)/2
		lnc = lgamma(ep1) - lgamma(eta/2) -0.5*log( pi*(eta-2))
		c   = exp(lnc)
		a   = 4*lda*c*(eta-2)/(eta-1)
		b   = sqrt(1+3*lda^2-a^2)
		my2 = 1+3*lda^2
		my3 = 16*c*lda*(1+lda^2)*((eta-2)^2)/((eta-1)*(eta-3))
		my4 = 3*(eta-2)*(1+10*lda^2+5*lda^4)/(eta-4)
		m3  = (my3-3*a*my2+2*a^3)/(b^3)
	} else{
		m3 = NA
	}
	return(m3)
}

.sstdexkurt = function( skew, shape )
{
	# Theoretical moments based on bijection betweeen Fernandez and Steel verions
	# and Hansen's Generalized Skew-T (credit due to Michael Rockinger)
	if(shape > 4 ){
		eta  = shape
		k2   = skew^2
		lda  = (k2-1)/(k2+1)
		ep1 = (eta+1)/2
		lnc = lgamma(ep1) - lgamma(eta/2) -0.5*log( pi*(eta-2))
		c   = exp(lnc)
		a   = 4*lda*c*(eta-2)/(eta-1)
		b   = sqrt(1+3*lda^2-a^2)
		my2 = 1+3*lda^2
		my3 = 16*c*lda*(1+lda^2)*((eta-2)^2)/((eta-1)*(eta-3))
		my4 = 3*(eta-2)*(1+10*lda^2+5*lda^4)/(eta-4)
		m3  = (my3-3*a*my2+2*a^3)/(b^3)
		m4  = -3 + (my4-4*a*my3+6*(a^2)*my2-3*a^4)/(b^4)
	} else{
		m4 = NA
	}
	return(m4)
}

.gedskew = function( shape )
{
	0
}

.gedexkurt = function( shape )
{
	( ( ( gamma(1/shape)/gamma(3/shape) )^2 ) * ( gamma(5/shape)/gamma(1/shape) ) ) - 3
}

.sgedskew = function( skew, shape )
{
	lambda = sqrt ( 2^(-2/shape) * gamma(1/shape) / gamma(3/shape) )
	m1 = ((2^(1/shape)*lambda)^1 * gamma(2/shape) / gamma(1/shape))
	m2 = 1
	m3 = ((2^(1/shape)*lambda)^3 * gamma(4/shape) / gamma(1/shape))
	(skew - 1/skew) * ( ( m3 + 2 * m1^3 - 3 * m1 * m2 ) * ( skew^2 + (1/skew^2) ) + 3 * m1 * m2 - 4 * m1^3 )/
			( ( (m2 - m1^2) * (skew^2 + 1/skew^2) + 2 * m1^2 - m2) ^ (3/2) )
	
}

.sgedexkurt= function( skew, shape )
{
	lambda = sqrt ( 2^(-2/shape) * gamma(1/shape) / gamma(3/shape) )
	m1 = ((2^(1/shape)*lambda)^1 * gamma(2/shape) / gamma(1/shape))
	m2 = 1
	m3 = ((2^(1/shape)*lambda)^3 * gamma(4/shape) / gamma(1/shape))
	m4 = ((2^(1/shape)*lambda)^4 * gamma(5/shape) / gamma(1/shape))
	cm4 = (-3 * m1^4 * (skew - 1/skew)^4) + 
			( 6 * m1^2 * (skew - 1/skew)^2 * m2*(skew^3 + 1/skew^3) )/(skew + 1/skew) - 
			( 4 * m1*(skew - 1/skew) * m3 * (skew^4 - 1/skew^4) )/(skew+1/skew) + 
			( m4 * (skew^5 + 1/skew^5) )/(skew + 1/skew)
	( cm4/( ( (m2 - m1^2) * (skew^2 + 1/skew^2) + 2 * m1^2 - m2) ^ 2 ) ) - 3
}

.jsuskew = function( mu = 0, sigma = 1, skew, shape )
{	
	Omega = -skew/shape
	w = exp(shape^-2)
	s3 = -0.25*sqrt(w)*( (w-1)^2 )*(w*(w+2)*sinh(3*Omega)+3*sinh(Omega))
	s3/(0.5*(w-1)*(w*cosh(2*Omega)+1))^(3/2)
}

.jsuexkurt = function( mu = 0, sigma = 1, skew, shape )
{
	Omega = -skew/shape
	w = exp(shape^-2)
	s4 = 0.125 * (w-1)^2*(w^2*(w^4+2*w^3+3*w^2-3)*cosh(4*Omega)+4*w^2*(w+2)*cosh(2*Omega)+3*(2*w+1))
	ans = s4/(0.5*(w-1)*(w*cosh(2*Omega)+1))^2
	return(ans - 3)
}

.ghstskew = function(skew, shape){
	if(shape<6){
		ans = NA
	} else{
		params = .paramGHST(nu = shape, betabar = skew)
		delta = params[2]
		beta = params[3]
		nu = params[4]
		beta2 = beta*beta
		delta2 = delta*delta
		ans = ( (2 * sqrt(nu - 4)*beta*delta)/( (2*beta2*delta2 + (nu-2)*(nu-4))^(3/2) ) ) * (3*(nu-2) + ((8*beta2*delta2)/(nu-6)))
	}
	return( ans )
}

.ghstexkurt = function(skew, shape){
	if(shape<8){
		ans = NA
	} else{
		params = .paramGHST(nu = shape, betabar = skew)
		delta = params[2]
		beta = params[3]
		nu = params[4]
		beta2 = beta*beta
		delta2 = delta*delta
		k1 = 6/( (2*beta2*delta2+(nu-2)*(nu-4))^2)
		k21 = (nu-2)*(nu-2)*(nu-4)
		k22 = (16*beta2*delta2*(nu-2)*(nu-4))/(nu-6)
		k23 = (8*(beta2^2)*(delta2^2)*(5*nu-22))/((nu-6)*(nu-8))
		ans = k1*(k21+k22+k23)
	}
	return( ans )
}
######################################################################################
# MGF
#mgf.sgh = function(u, rho, zeta, lambda)
#{
#	pars=.paramGH(rho, zeta, lambda)
#	alpha = unname(pars["alpha"])
#	beta = unname(pars["beta"])
#	delta = unname(pars["delta"])
#	mu = unname(pars["mu"])
#	d1 = alpha^2-beta^2
#	d2 = (alpha^2-(beta+u)^2)
#	mg = exp(u*mu)*( d1/d2 )^(lambda/2) * BesselK(lambda, delta*sqrt(d2))/BesselK(lambda, delta*sqrt(d1))
#	return(mg)
#}
#DD <- function(expr, name, order = 1) {
#	if(order < 1) stop("'order' must be >= 1")
#	if(order == 1) D(expr, name)
#	else DD(D(expr, name), name, order - 1)
#}
#mgf.sgh.expression = function(u, rho, zeta, lambda)
#{
#	pars=.paramGH(rho, zeta, lambda)
#	alpha = unname(pars["alpha"])
#	beta = unname(pars["beta"])
#	delta = unname(pars["delta"])
#	mu = unname(pars["mu"])
#	d1 = alpha^2-beta^2
#	d2 = (alpha^2-(beta+u)^2)
#	mg = paste("~exp(u*",mu,")*(", d1,"/(",alpha^2,"-(",beta,"+u)^2)^(",lambda/2,"))*besselK(",lambda,", ",delta,"*sqrt(",alpha^2,"-(",beta,"+u)^2))/besselK(",lambda,",",delta*sqrt(d1),")",sep="")
#
#	return(expression(mg))
#}