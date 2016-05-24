#################################################################################
##
##   R package Rsolnp by Alexios Ghalanos and Stefan Theussl Copyright (C) 2009-2013
##   This file is part of the R package Rsolnp.
##
##   The R package Rsolnp is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package Rsolnp is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

.eps = .Machine$double.eps

.subnpmsg = function(m){
	g1 = c("solnp-->")
	m1 = paste("\n",g1, "Redundant constraints were found. Poor\n",
			g1, "intermediate results may result.Suggest that you\n",
			g1, "remove redundant constraints and re-OPTIMIZE\n", sep = "")
	m2 = paste("\n",g1, "The linearized problem has no feasible\n",
			g1, "solution.  The problem may not be feasible.\n", sep = "")
	m3 = paste("\n",g1, "Minor optimization routine did not converge in the \n",
			g1, "specified number of minor iterations.  You may need\n",
			g1, "to increase the number of minor iterations.        \n", sep = "")
	ans = switch(m,
			m1 = m1,
			m2 = m2,
			m3 = m3)
	cat(ans)
}



.checkpars = function(pars, LB, UB, .env)
{
	if(is.null(pars))
		stop("\nsolnp-->error: must supply starting parameters\n", call. = FALSE)
	if(!is.null(LB)){
		if(length(pars) != length(LB))
			stop("\nsolnp-->error: LB length not equal to parameter length\n", call. = FALSE)
		if(is.null(UB)) UB = rep(.Machine$double.xmax/2, length(LB))
	} else{
		LB = NULL
	}
	if(!is.null(UB)){
		if(length(pars) != length(UB))
			stop("\nsolnp-->error: UB length not equal to parameter length\n", call. = FALSE)
		if(is.null(LB)) LB = rep(-.Machine$double.xmax/2, length(UB))
	} else{
		UB = NULL
	}
	if(!is.null(UB) && any(LB > UB))
		stop("\nsolnp-->error: UB must be greater than LB\n", call. = FALSE)
	
	if(!is.null(UB) && any(LB == UB))
		warning("\nsolnp-->warning: Equal Lower/Upper Bounds Found. Consider\n
						excluding fixed parameters.\n", call. = FALSE)
	# deal with infinite values as these are not accepted by solve.QP
	if(!is.null(LB) && !any(is.finite(LB))){
		idx = which(!is.finite(LB))
		LB[idx] = sign(LB[idx])*.Machine$double.xmax/2
	}
	if(!is.null(UB) && !any(is.finite(UB))){
		idx = which(!is.finite(UB))
		UB[idx] = sign(UB[idx])*.Machine$double.xmax/2
	}	
	assign(".LB", LB, envir = .env)
	assign(".UB", UB, envir = .env)
	return(1)
}

.checkfun = function(pars, fun, .env, ...)
{
	if(!is.function(fun)) stop("\nsolnp-->error: fun does not appear to be a function\n", call. = FALSE)
	val = fun(pars, ...)
	if(length(val) != 1)  stop("\nsolnp-->error: objective function returns value of length greater than 1!\n", call. = FALSE)
	
	assign(".solnp_fun", fun, envir = .env)
	ctmp = get(".solnp_nfn", envir =  .env)
	assign(".solnp_nfn", ctmp + 1, envir = .env)
	return(val)
}

# Might eventually use this, but really the user must take care of such problems
# in their own function/setup
.safefunx = function(pars, fun, .env, ...){
	xnames = get("xnames", envir = .env)
	names(pars) = xnames
	v  = fun(pars, ...)
	if(is.na(v) | !is.finite(v) | is.nan(v)) {
		warning(paste("\nsolnp-->warning: ", v , " detected in function call...check your function\n", sep = ""), immediate. = FALSE)
		v = 1e24
	}
	v
}

.checkgrad = function(pars, fun, .env, ...)
{
	n = length(pars)
	val = fun(pars, ...)
	if(length(val)!=n)
		stop("\nsolnp-->error: gradient vector length must be equal to length(pars)\n", call. = FALSE)
	assign(".solnp_gradfun", fun, envir = .env)
	return(val)
}

.checkhess = function(pars, fun, .env, ...)
{
	n = length(pars)
	val = fun(pars, ...)
	if(length(as.vector(val)) != (n*n))
		stop("\nsolnp-->error: hessian must be of length length(pars) x length(pars)\n", call. = FALSE)
	assign(".solnp_hessfun", fun, envir = .env)
	return(val)
}

.checkineq = function(pars, fun, ineqLB, ineqUB, .env, ...)
{
	xnames = get("xnames", envir = .env)
	val = fun(pars, ...)
	n = length(val)
	if(!is.null(ineqLB)){
		if(length(ineqLB) != n)
			stop("\nsolnp-->error: inequality function returns vector of different length to
							inequality lower bounds\n", call. = FALSE)
	} else{
		stop("\nsolnp-->error: inequality function given without lower bounds\n", call. = FALSE)
	}
	if(!is.null(ineqUB)){
		if(length(ineqUB) != n)
			stop("\nsolnp-->error: inequality function returns vector of different length to
							inequality upper bounds\n", call. = FALSE)
	} else{
		stop("\nsolnp-->error: inequality function given without upper bounds\n", call. = FALSE)
	}
	if(any(ineqLB > ineqUB))
		stop("\nsolnp-->error: ineqUB must be greater than ineqLB\n", call. = FALSE)

	assign(".ineqLB", ineqLB, envir = .env)
	assign(".ineqUB", ineqUB, envir = .env)
	.solnp_ineqfun = function(x, ...){
		names(x) = xnames
		fun(x, ...)
	}
	assign(".solnp_ineqfun", .solnp_ineqfun, envir = .env)
	return(val)
}


.checkeq = function(pars, fun, eqB, .env, ...)
{
	xnames = get("xnames", envir = .env)
	n = length(eqB)
	val = fun(pars, ...) - eqB
	if(length(val)!=n)
		stop("\nsolnp-->error: equality function returns vector of different length
						to equality value\n", call. = FALSE)
	.eqB = eqB
	assign(".eqB", .eqB, envir = .env)
	.solnp_eqfun = function(x, ...){
		names(x) = xnames
		fun(x, ...) - .eqB
	}
	assign(".solnp_eqfun", .solnp_eqfun, envir = .env)
	return(val)
}


# check the jacobian of inequality
.cheqjacineq = function(pars, fun, .env,  ...)
{
	# must be a matrix -> nrows = no.inequalities, ncol = length(pars)
	val = fun(pars, ...)
	.ineqLB = get(".ineqLB", envir = .env)
	.ineqUB = get(".ineqUB", envir = .env)
	if(!is.matrix(val))
		stop("\nsolnp-->error: Jacobian of Inequality must return a matrix type object\n", call. = FALSE)
	nd = dim(val)
	if(nd[2] != length(pars))
		stop("\nsolnp-->error: Jacobian of Inequality column dimension must be equal to length
						of parameters\n", call. = FALSE)
	if(nd[1] != length(.ineqUB))
		stop("\nsolnp-->error: Jacobian of Inequality row dimension must be equal to length
						of inequality bounds vector\n", call. = FALSE)
	# as in inequality function, transforms from a 2 sided inequality to a one sided inequality
	# (for the jacobian).
	.solnp_ineqjac = function(x, ...) { retval = fun(x, ...); rbind( retval ) }
	assign(".solnp_ineqjac", .solnp_ineqjac, envir = .env)
	return(val)
}

# check the jacobian of equality
.cheqjaceq = function(pars, fun, .env, ...)
{
	# must be a matrix -> nrows = no.equalities, ncol = length(pars)
	val = fun(pars, ...)
	.eqB = get(".eqB", envir = .env)
	if(!is.matrix(val))
		stop("\nsolnp-->error: Jacobian of Equality must return a matrix type object\n", call. = FALSE)
	nd = dim(val)
	if(nd[2] != length(pars))
		stop("\nsolnp-->error: Jacobian of Equality column dimension must be equal to length of parameters\n", call. = FALSE)
	if(nd[1] != length(.eqB))
		stop("\nsolnp-->error: Jacobian of Equality row dimension must be equal to length of equality bounds vector\n", call. = FALSE)
	assign(".solnp_eqjac", fun, envir = .env)
	return(val)
}

# reporting function
.report = function(iter, funv, pars)
{
	cat( paste( "\nIter: ", iter ," fn: ", format(funv, digits = 4, scientific = 5, nsmall = 4, zero.print = TRUE), "\t Pars: ", sep=""), 
			format(pars, digits = 4, scientific = 6, nsmall = 5, zero.print = TRUE) )
}

# finite difference gradient
.fdgrad = function(pars, fun, ...)
{
	if(!is.null(fun)){
		
		y0 = fun(pars, ...)
		nx = length(pars)
		grd = rep(0, nx)
		deltax = sqrt(.eps)
		for(i in 1:nx)
		{
			init = pars[i]
			pars[i]= pars[i] + deltax
			grd[i] = (fun(pars, ...) - y0) / deltax
			pars[i] = init
		}
	}
	else
	{
		grd = 0
	}
	return(grd)
}

# finite difference jacobian
.fdjac = function(pars, fun, ...)
{
	nx = length(pars)
	if(!is.null(fun))
	{
		y0 = fun(pars, ...)
		nf = length (y0)
		jac = matrix(0, nrow = nf, ncol= nx)
		deltax = sqrt (.eps)
		for(i  in 1:nx)
		{
			init = pars[i]
			pars[i]= pars[i] + deltax
			jac[,i] = (fun(pars, ...) - y0) / deltax
			pars[i] = init
		}
	} else{
		jac = rep(0, nx)
	}
	return(jac)
}

.emptygrad = function(pars, ...)
{
	matrix(0, nrow = 0, ncol = 1)
}

.emptyjac = function(pars, ...)
{
	#matrix(0, nrow = 0, ncol = length(pars))
	NULL
}

.emptyfun = function(pars, ...)
{
	NULL
}

.ineqlbfun = function(pars, .env, ...)
{
	LB = get(".solnp_LB", envir = .env)
	UB = get(".solnp_UB", envir = .env)
	.solnp_ineqfun = get(".solnp_ineqfun", envir = .env)
	res = c(pars - LB,  UB - pars)
	if(!is.null(.solnp_ineqfun)) res = c(.solnp_ineqfun(pars, ...), res)
	res
}

.ineqlbjac = function(pars, .env, ...)
{
	.solnp_ineqjac = get(".solnp_ineqjac", envir = .env)
	n = length(pars)
	res = rbind(diag(n), -diag(n))
	if(!is.null(.solnp_ineqjac)) res = rbind(.solnp_ineqjac(pars, ...), res)
	res
}

.solnpctrl = function(control){
	# parameters check is now case independent
	ans = list()
	params = unlist(control)
	if(is.null(params)) {
		ans$rho = 1
		ans$outer.iter = 400
		ans$inner.iter = 800
		ans$delta = 1.0e-7
		ans$tol = 1.0e-8
		ans$trace = 1
	} else{
		npar = tolower(names(unlist(control)))
		names(params) = npar
		if(any(substr(npar, 1, 3) == "rho")) ans$rho = as.numeric(params["rho"]) else ans$rho = 1
		if(any(substr(npar, 1, 10) == "outer.iter")) ans$outer.iter = as.numeric(params["outer.iter"]) else ans$outer.iter = 400
		if(any(substr(npar, 1, 10) == "inner.iter")) ans$inner.iter = as.numeric(params["inner.iter"]) else ans$inner.iter = 800
		if(any(substr(npar, 1, 5) == "delta")) ans$delta = as.numeric(params["delta"]) else ans$delta = 1.0e-7
		if(any(substr(npar, 1, 3) == "tol")) ans$tol = as.numeric(params["tol"]) else ans$tol = 1.0e-8
		if(any(substr(npar, 1, 5) == "trace")) ans$trace = as.numeric(params["trace"]) else ans$trace = 1
	}
	return(ans)
}

.zeros = function( n = 1, m = 1)
{
	if(missing(m)) m = n
	sol = matrix(0, nrow = n, ncol = m)
	return(sol)
}

.ones = function(n = 1, m = 1)
{
	if(missing(m)) m = n
	sol = matrix(1, nrow = n, ncol = m)
	return(sol)
}

.vnorm = function(x)
{
	sum((x)^2)^(1/2)
}

.solvecond = function(x)
{
	z = svd(x)$d
	if(any( z == 0 )) ret = Inf else ret = max( z ) / min( z )
	return(ret)
}
