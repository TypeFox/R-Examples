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
.dccsolver = function(solver, pars, fun, Ifn, ILB, IUB, gr, hessian, 
		control, LB, UB, ux = NULL, ci = NULL, mu = NULL, arglist)
{
	control = .dccgetcontrol(solver, control)
	retval = switch(solver,
			nlminb = .dccnlminbsolver(pars, fun, gr, hessian, control, LB, UB, arglist),
			solnp  = .dccsolnpsolver(pars, fun, Ifn, ILB, IUB, control, LB, UB, arglist),
			lbfgs  = .dcclbfgssolver(pars, fun, gr, control, LB, UB, arglist),
			gosolnp = .dccgosolnpsolver(pars, fun, Ifn, ILB, IUB, control, LB, UB, arglist))
	return(retval)
}

.dccsolnpsolver= function(pars, fun, Ifn, ILB, IUB, control, LB, UB, arglist)
{
	op <- options()
	options(warn = 0)
	ans = try(solnp(pars = pars, fun = fun, eqfun = NULL, 
					eqB = NULL, ineqfun = Ifn, ineqLB = ILB, 
					ineqUB = IUB, LB = LB, UB = UB, control = control, 
					arglist = arglist),
			silent = TRUE)
	if(inherits(ans,"try-error")){
		sol = list()
		sol$convergence = 1
		sol$message = ans
		sol$pars = rep(NA, length(pars))
		names(sol$pars) = names(pars)
	} else{
		sol = ans
	}
	hess = NULL
	options(op)
	return(list(sol = sol, hess = hess))
}
.dccgosolnpsolver = function(pars, fun, Ifn, ILB, IUB, gocontrol, LB, UB, arglist){
	control = .dccsolnpctrl(gocontrol)
	gocontrol = .dccgosolnpctrl(gocontrol)
	n.restarts = gocontrol$n.restarts
	rseed = gocontrol$rseed
	n.sim = gocontrol$n.sim
	cluster = arglist$cluster
	op <- options()
	options(warn = 0)
	# use the truncated normal distribution
	distr.opt = vector(mode = "list", length = length(pars))
	for(i in 1:length(pars)){
		distr.opt[[i]]$mean = pars[i]
		distr.opt[[i]]$sd = sqrt(pars[i]^2)*2
	}
	ans = try(gosolnp(pars = pars, fixed = NULL, fun = fun, eqfun = NULL, 
					eqB = NULL, ineqfun = Ifn, ineqLB = ILB, ineqUB = IUB, LB = LB, 
					UB = UB, control = control, distr = rep(2, length(LB)), 
					distr.opt = distr.opt, n.restarts = n.restarts, n.sim = n.sim, 
					cluster = cluster, rseed = rseed, arglist),
			silent = TRUE)
	if(inherits(ans,"try-error")){
		sol = list()
		sol$convergence = 1
		sol$message = ans
		sol$pars = rep(NA, length(pars))
		names(sol$pars) = names(pars)
	} else{
		sol = ans
	}
	hess = NULL
	options(op)
	return(list(sol = sol, hess = hess))
}
.dccnlminbsolver = function(pars, fun, gr, hessian, control, LB, UB, arglist){
	parscale = rep(1, length(pars))
	ans = try(nlminb(start = pars, objective = fun, gradient = gr, hessian = hessian,
					arglist = arglist, scale = 1/parscale, control = control, 
					lower = LB, upper = UB), silent = TRUE)
	if(inherits(ans, "try-error")){
		sol = list()
		sol$convergence = 1
		sol$message = ans
		sol$pars = rep(NA, length(pars))
		names(sol$pars) = names(pars)
	} else{
		sol = ans
		sol$pars = ans$par
		sol$par = NULL
	}
	hess = NULL
	return(list(sol = sol,hess = hess))
}

.dcclbfgssolver = function(pars, fun, gr, control, LB, UB, arglist){
	control$parscale = rep(1, length(pars))
	ans = try(optim(par = pars, fn = fun, gr = gr, arglist = arglist,
					method = "L-BFGS-B", lower = LB, upper = UB, control = control, 
					hessian = TRUE),silent=TRUE)
	if(inherits(ans, "try-error")){
		sol = list()
		sol$convergence = 1
		sol$message = ans
		sol$pars = rep(NA, length(pars))
		names(sol$pars) = names(pars)
	} else{
		sol = ans
		sol$pars = ans$par
		sol$par = NULL
	}
	hess = sol$hessian
	return(list(sol = sol, hess = hess))
}

# default control for solvers:
.dccgetcontrol = function(solver, control)
{
	ans = switch(solver,
			nlminb = .dccnlminbctrl(control),
			solnp = .dccsolnpctrl(control),
			gosolnp = .dccgosolnpctrl(control),
			lbfgs = .dcclbfgsctrl(control))
	return(ans)
}

.dccnlminbctrl = function(control)
{
	if(is.null(control$eval.max)) control$eval.max = 2000
	if(is.null(control$iter.max)) control$iter.max = 1500
	if(is.null(control$abs.tol)) control$abs.tol = 1e-20
	if(is.null(control$rel.tol)) control$rel.tol = 1e-10
	if(is.null(control$x.tol)) control$x.tol = 1.5e-8
	if(is.null(control$step.min)) control$step.min = 2.2e-14
	return(control)
}

.dcclbfgsctrl = function(control)
{
	if(is.null(control$REPORT)) control$REPORT = 10
	if(is.null(control$lmm)) control$lmm = 15
	if(is.null(control$pgtol)) control$pgtol = 1e-8
	if(is.null(control$factr)) control$factr = 1e-8
	return(control)
}

.dccsolnpctrl = function(control){
	# parameters check is now case independent
	ans = list()
	params = unlist(control)
	if(is.null(params)) {
		ans$rho = 1
		ans$outer.iter = 50
		ans$inner.iter = 1800
		ans$delta = 1.0e-7
		ans$tol = 1.0e-8
		ans$trace = 0
	} else{
		npar = tolower(names(unlist(control)))
		names(params) = npar
		if(any(substr(npar, 1, 3) == "rho")) ans$rho = as.numeric(params["rho"]) else ans$rho = 1
		if(any(substr(npar, 1, 5) == "outer.iter")) ans$outer.iter = as.numeric(params["outer.iter"]) else ans$outer.iter = 50
		if(any(substr(npar, 1, 5) == "inner.iter")) ans$inner.iter = as.numeric(params["inner.iter"]) else ans$inner.iter = 1000
		if(any(substr(npar, 1, 5) == "delta")) ans$delta = as.numeric(params["delta"]) else ans$delta = 1.0e-7
		if(any(substr(npar, 1, 3) == "tol")) ans$tol = as.numeric(params["tol"]) else ans$tol = 1.0e-8
		if(any(substr(npar, 1, 5) == "trace")) ans$trace = as.numeric(params["trace"]) else ans$trace = 0
	}
	return(ans)
}

.dccgosolnpctrl = function(control){
	# parameters check is now case independent
	ans = list()
	params = unlist(control)
	if(is.null(params)) {
		ans$n.restarts = 1
		ans$rseed
		ans$n.sim = 500
	} else{
		npar = tolower(names(unlist(control)))
		names(params) = npar
		if(any(substr(npar, 1, 10) == "n.restarts")) ans$n.restarts = as.numeric(params["n.restarts"]) else ans$n.restarts = 1
		if(any(substr(npar, 1, 5) == "rseed")) ans$rseed = as.numeric(params["rseed"]) else ans$rseed = NULL
		if(any(substr(npar, 1, 5) == "n.sim")) ans$n.sim = as.numeric(params["n.sim"]) else ans$n.sim = 500	
	}
	return(ans)
}

.adcccon = function(pars, arglist){
	ipars = arglist$ipars
	estidx = arglist$estidx
	idx = arglist$model$pidx
	Nbar = arglist$Nbar
	Qbar = arglist$Qbar
	ipars[estidx, 1] = pars
	dcca = ipars[idx["dcca", 1]:idx["dcca", 2], 1]
	dccb = ipars[idx["dccb", 1]:idx["dccb", 2], 1]
	dccg = ipars[idx["dccg", 1]:idx["dccg", 2], 1]
	Qbar2 = solve( .sqrtsymmat(Qbar) )
	delta = max( eigen( Qbar2 %*% Nbar %*% Qbar2, symmetric = TRUE, only.values = TRUE )$values )
	#print(delta)
	#if(!is.finite(delta) | is.na(delta) |is.nan(delta) | is.complex(delta)) delta = 1
	return( sum(dcca) + sum(dccb) + delta*sum(dccg) )
}

.dcccon = function(pars, arglist){
	ipars = arglist$ipars
	estidx = arglist$estidx
	idx = arglist$model$pidx
	ipars[estidx, 1] = pars
	dcca = ipars[idx["dcca", 1]:idx["dcca", 2], 1]
	dccb = ipars[idx["dccb", 1]:idx["dccb", 2], 1]
	return( sum(dcca) + sum(dccb) )
}