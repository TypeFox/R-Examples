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
# implements nlminb, lbgs and solnp
# only solnp implements true constraint (stationarity) optimization
.garchsolver = function(solver, pars, fun, Ifn, ILB, IUB, gr, hessian, parscale, 
		control, LB, UB, ux=NULL, ci=NULL, mu=NULL, arglist)
{
	gocontrol = control
	control = .getcontrol(solver, control)
	retval = switch(solver,
			hybrid = .hybridsolver(pars, fun, Ifn, ILB, IUB, gr, hessian, parscale, gocontrol, LB, UB, arglist),
			nlminb = .nlminbsolver(pars, fun, gr, hessian, parscale, control, LB, UB, arglist),
			solnp = .solnpsolver(pars, fun, Ifn, ILB, IUB, control, LB, UB, arglist),
			gosolnp = .gosolnpsolver(pars, fun, Ifn, ILB, IUB, gocontrol, LB, UB, arglist),
			lbfgs = .lbfgssolver(pars, fun, gr, parscale, control, LB, UB, arglist),
			nloptr = .nloptrsolver(pars, fun, control, LB, UB, arglist))
	return(retval)
}

.nlminbsolver = function(pars, fun, gr, hessian, parscale, control, LB, UB, arglist){
	ans = try(nlminb(start = pars, objective = fun, gradient = gr, hessian = hessian,
					arglist = arglist, scale = 1/parscale, control = control, 
					lower = LB, upper = UB), silent = TRUE)
	if(inherits(ans, "try-error")){
		sol = list()
		sol$convergence = 1
		sol$message = ans
		sol$par = rep(NA, length(pars))
		names(sol$par) = names(pars)
	} else{
		sol = ans
	}
	hess = NULL
	return(list(sol = sol,hess = hess))
}

.solnpsolver = function(pars, fun, Ifn, ILB, IUB, control, LB, UB, arglist){
	ans = try(solnp(pars, fun = fun, eqfun = NULL, eqB = NULL, ineqfun = Ifn, ineqLB = ILB, 
					ineqUB = IUB, LB = LB, UB = UB, control = control, arglist), silent = TRUE)
	if(inherits(ans,"try-error")){
		sol = list()
		sol$convergence = 1
		sol$message = ans
		sol$par = rep(NA, length(pars))
		names(sol$par) = names(pars)
		warning("\nrugarch-->warning: no convergence...\n")
	} else{
		sol = ans
	}
	hess = NULL
	return(list(sol = sol, hess = hess))
}

.hybridsolver = function(pars, fun, Ifn, ILB, IUB, gr, hessian, parscale, control, LB, UB, arglist){
	xcontrol = .getcontrol("solnp", control)
	ans = .solnpsolver(pars, fun, Ifn, ILB, IUB, xcontrol, LB, UB, arglist)
	if(ans$sol$convergence >= 1)
	{
		xcontrol = .getcontrol("nlminb", control)
		if(xcontrol$trace) cat("\nTrying nlminb solver...\n")
		ans = .nlminbsolver(pars, fun, gr, hessian, parscale, control=xcontrol, LB, UB, arglist = arglist)
		if(ans$sol$convergence>=1){
			if(xcontrol$trace) cat("\nTrying gosolnp solver...\n")
			ans = .gosolnpsolver(pars, fun, Ifn, ILB, IUB, control, LB, UB, arglist)
			if(ans$sol$convergence>=1){
				xcontrol = .getcontrol("nloptr", control)
				if(xcontrol$print_level) cat("\nLast try...nloptr solver...\n")
				ans = .nloptrsolver(pars, fun, xcontrol, LB, UB, arglist = arglist)
			}
		}
	}
	return(ans)
}

.gosolnpsolver = function(pars, fun, Ifn, ILB, IUB, gocontrol, LB, UB, arglist){
	control = .solnpctrl(gocontrol)
	gocontrol = .gosolnpctrl(gocontrol)
	n.restarts = gocontrol$n.restarts
	doparallel = gocontrol$parallel
	pkg = gocontrol$pkg
	cores = gocontrol$cores
	rseed = gocontrol$rseed
	n.sim = gocontrol$n.sim
	op <- options()
	options(warn = 0)
	# use the truncated normal distribution
	distr.opt = vector(mode = "list", length = length(pars))
	for(i in 1:length(pars)){
		distr.opt[[i]]$mean = pars[i]
		distr.opt[[i]]$sd = ifelse(pars[i]==0, 0.1, sqrt(pars[i]^2)*2)
	}
	if(doparallel){
		cl = makePSOCKcluster(as.integer(cores))
		clusterEvalQ(cl, library(rugarch))
		# set parallel mode on (pmode=1) to tell the likelihood function not to use
		# environment assignment
		arglist$pmode = 1
	} else{
		cl = NULL
	}
	ans = try(gosolnp(pars = pars, fixed = NULL, fun = fun, eqfun = NULL, 
			eqB = NULL, ineqfun = Ifn, ineqLB = ILB, ineqUB = IUB, LB = LB, 
			UB = UB, control = control, distr = rep(2, length(LB)), 
			distr.opt = distr.opt, n.restarts = n.restarts, n.sim = n.sim, 
			cluster = cl, rseed = rseed, arglist),
	silent = TRUE)
	if(doparallel) stopCluster(cl)
	# return to normal mode
	arglist$pmode = 0
	if(inherits(ans,"try-error")){
		sol = list()
		sol$convergence = 1
		sol$message = ans
		sol$par = rep(NA, length(pars))
		names(sol$par) = names(pars)
	} else{
		sol = ans
	}
	hess = NULL
	options(op)
	return(list(sol = sol, hess = hess))
}

.lbfgssolver = function(pars, fun, gr, parscale, control, LB, UB, arglist){
	control$parscale = parscale
	ans = try(optim(par = pars, fn = fun, gr = gr, arglist = arglist,
			method = "L-BFGS-B", lower = LB, upper = UB, control = control, 
			hessian = TRUE),silent=TRUE)
	if(inherits(ans, "try-error")){
		sol = list()
		sol$convergence = 1
		sol$message = ans
		sol$par = rep(NA, length(pars))
		names(sol$par) = names(pars)
	} else{
		sol = ans
	}
	hess = sol$hessian
	return(list(sol = sol, hess = hess))
}

# nloptr solver requires named (...)!
.nloptrsolver = function(pars, fun, control, LB, UB, arglist){
	ans = try(nloptr(x0 = pars, eval_f = fun,  eval_grad_f = NULL, 
					eval_g_ineq = NULL, lb = LB, ub = UB, eval_jac_g_ineq = NULL, 
					eval_g_eq = NULL, eval_jac_g_eq = NULL, opts = control, arglist = arglist), 
			silent=TRUE)
	if(inherits(ans, "try-error")){
		sol = list()
		sol$convergence = 1
		sol$message = ans
		sol$par = rep(NA, length(pars))
		names(sol$par) = names(pars)
	} else{
		sol = list()
		sol$convergence = 0
		sol$message = ans$message
		sol$par = ans$solution
	}
	hess = NULL
	return(list(sol = sol, hess = hess))
}

# default control for solvers:
.getcontrol = function(solver, control)
{
	ans = switch(solver,
		nlminb = .nlminbctrl(control),
		solnp = .solnpctrl(control),
		gosolnp = .gosolnpctrl(control),
		lbfgs = .lbfgsctrl(control),
		nloptr = .nloptrctrl(control))
	return(ans)
}

.nlminbctrl = function(control)
{
	if(is.null(control$eval.max)) control$eval.max = 2000
	if(is.null(control$trace)) control$trace = 0
	if(is.null(control$iter.max)) control$iter.max = 1500
	if(is.null(control$abs.tol)) control$abs.tol = 1e-20
	if(is.null(control$rel.tol)) control$rel.tol = 1e-10
	if(is.null(control$x.tol)) control$x.tol = 1.5e-8
	if(is.null(control$step.min)) control$step.min = 2.2e-14
	return(control)
}

.lbfgsctrl = function(control)
{
	if(is.null(control$REPORT)) control$REPORT = 10
	if(is.null(control$lmm)) control$lmm = 15
	if(is.null(control$pgtol)) control$pgtol = 1e-8
	if(is.null(control$factr)) control$factr = 1e-8
	return(control)
}

.solnpctrl = function(control){
	# parameters check is now case independent
	ans = list()
	params = unlist(control)
	if(is.null(params)) {
		ans$rho = 1
		ans$outer.iter = 50
		ans$inner.iter = 1800
		ans$delta = 1.0e-8
		ans$tol = 1.0e-8
		ans$trace = 0
	} else{
		npar = tolower(names(unlist(control)))
		names(params) = npar
		if(any(substr(npar, 1, 3) == "rho")) ans$rho = as.numeric(params["rho"]) else ans$rho = 1
		if(any(substr(npar, 1, 10) == "outer.iter")) ans$outer.iter = as.numeric(params["outer.iter"]) else ans$outer.iter = 50
		if(any(substr(npar, 1, 10) == "inner.iter")) ans$inner.iter = as.numeric(params["inner.iter"]) else ans$inner.iter = 1000
		if(any(substr(npar, 1, 5) == "delta")) ans$delta = as.numeric(params["delta"]) else ans$delta = 1.0e-8
		if(any(substr(npar, 1, 3) == "tol")) ans$tol = as.numeric(params["tol"]) else ans$tol = 1.0e-8
		if(any(substr(npar, 1, 5) == "trace")) ans$trace = as.numeric(params["trace"]) else ans$trace = 0
	}
	return(ans)
}

.gosolnpctrl = function(control){
	# parameters check is now case independent
	ans = list()
	params = unlist(control)
	if(is.null(params)) {
		ans$parallel = FALSE
		ans$pkg = "snowfall"
		ans$cores= 2
		ans$n.restarts = 1
		ans$rseed
		ans$n.sim = 500
		ans$trace = 0
	} else{
		npar = tolower(names(unlist(control)))
		names(params) = npar
		ans$parallel.control = list()
		if(any(substr(npar, 1, 8) == "parallel")) ans$parallel = as.logical(params["parallel"]) else ans$parallel = FALSE
		if(any(substr(npar, 1, 3) == "pkg")) ans$pkg = unname(params["pkg"]) else ans$pkg = "snowfall"
		if(any(substr(npar, 1, 5) == "cores")) ans$cores = unname(params["cores"]) else ans$cores = 2
		if(any(substr(npar, 1, 10) == "n.restarts")) ans$n.restarts = as.numeric(params["n.restarts"]) else ans$n.restarts = 1
		if(any(substr(npar, 1, 5) == "rseed")) ans$rseed = as.numeric(params["rseed"]) else ans$rseed = NULL
		if(any(substr(npar, 1, 5) == "n.sim")) ans$n.sim = as.numeric(params["n.sim"]) else ans$n.sim = 500
		if(any(substr(npar, 1, 5) == "trace")) ans$trace = as.numeric(params["trace"]) else ans$trace = 0
	}
	return(ans)
}

.nloptrctrl = function(control){
	#solver = c("NLOPT_LN_COBYLA", "NLOPT_LN_BOBYQA", "NLOPT_LN_PRAXIS", "NLOPT_LN_NELDERMEAD", "NLOPT_LN_SBPLX",
	# and AUGLAG + solver
	# subsolvers:
	xsub = c("NLOPT_LN_COBYLA", "NLOPT_LN_BOBYQA", "NLOPT_LN_PRAXIS", "NLOPT_LN_NELDERMEAD", "NLOPT_LN_SBPLX")
	ans = list()
	params = unlist(control)
	if(is.null(params)) {
		ans$ftol_rel = 1e-8
		ans$xtol_rel = 1e-6
		ans$maxeval = 25000
		ans$print_level = 0
		mainsolver = "NLOPT_LN_SBPLX"
		subsolver = NULL
	} else{
		npar = tolower(names(unlist(control)))
		names(params) = npar
		if(any(substr(npar, 1, 8) == "ftol_rel")) ans$ftol_rel = as.numeric(params["ftol_rel"]) else ans$ftol_rel = 1e-8
		if(any(substr(npar, 1, 8) == "xtol_rel")) ans$xtol_rel = as.numeric(params["xtol_rel"]) else ans$xtol_rel = 1e-6
		if(any(substr(npar, 1, 7) == "maxeval")) ans$maxeval = as.numeric(params["maxeval"]) else ans$maxeval = 25000
		if(any(substr(npar, 1, 11) == "print_level")) ans$print_level = as.numeric(params["print_level"]) else ans$print_level = 0
		if(any(substr(npar, 1, 6) == "solver")){
			if(abs(as.integer(params["solver"]))>5){
				mainsolver = "NLOPT_LN_AUGLAG"
				subsolver = xsub[min(5,abs(as.integer(params["solver"]))-5 )]
			} else{
				mainsolver = xsub[min(5,abs(as.integer(params["solver"])))]
				subsolver = NULL
			}
		} else{
			mainsolver = "NLOPT_LN_SBPLX"
			subsolver = NULL
		}
	}
	if(is.null(subsolver)){
		ans$algorithm = mainsolver
	} else{
		ans$algorithm = mainsolver
		ans$local_opts$algorithm = subsolver
		ans$local_opts$ftol_abs = ans$ftol_rel
		ans$local_opts$xtol_rel = ans$xtol_rel
		ans$local_opts$maxeval = 2000
		ans$local_opts$print_level = 0
		
	}
	return(ans)
}
#----------------------------------------------------------------------------------
.garchconbounds = function(){
	return(list(LB = eps,UB = 0.999))
}

.sgarchcon = function(pars, arglist){
	ipars = arglist$ipars
	estidx = arglist$estidx
	idx = arglist$model$pidx
	ipars[estidx, 1] = pars
	distribution = arglist$model$modeldesc$distribution
	con = .persistsgarch1(pars = ipars[,1], idx, distribution = distribution)
	return(con)
}

.igarchcon = function(pars, arglist){
	# this is an equality constraint
	ipars = arglist$ipars
	estidx = arglist$estidx
	idx = arglist$model$pidx
	ipars[estidx, 1] = pars
	modelinc = arglist$model$modelinc
	con = ifelse(modelinc[9]>1, sum(ipars[idx["alpha", 1]:idx["alpha", 2], 1]) + sum(ipars[idx["beta", 1]:(idx["beta", 2]-1), 1]) -
					ipars[idx["beta", 2],1], sum(ipars[idx["alpha", 1]:idx["alpha", 2], 1]) - ipars[idx["beta", 2],1])
	return(con)
}

.aparchcon = function(pars, arglist){
	ipars = arglist$ipars
	estidx = arglist$estidx
	idx = arglist$model$pidx
	ipars[estidx, 1] = pars
	distribution = arglist$model$modeldesc$distribution	
	con = .persistaparch1(pars = ipars[,1], idx = idx, distribution = distribution)
	if(is.na(con)) con = 1
	return(con)
}

.fgarchcon = function(pars, arglist){
	ipars = arglist$ipars
	estidx = arglist$estidx
	idx = arglist$model$pidx
	ipars[estidx, 1] = pars
	distribution = arglist$model$modeldesc$distribution	
	vsubmodel = arglist$model$modeldesc$vsubmodel
	con = .persistfgarch1(pars = ipars[,1], idx = idx, distribution = distribution, 
			submodel = vsubmodel)
	if(is.na(con)) con = 1
	return(con)
}

.realgarchcon = function(pars, arglist)
{
	ipars = arglist$ipars
	estidx = arglist$estidx
	idx = arglist$model$pidx
	ipars[estidx, 1] = pars
	modelinc = arglist$model$modelinc
	con = sum(ipars[idx["beta", 1]:idx["beta", 2], 1]) + ipars[idx["delta", 2], 1]*sum(ipars[idx["alpha", 1]:idx["alpha", 2], 1])
	return(con)
}

.gjrgarchcon = function(pars, arglist){
	ipars = arglist$ipars
	estidx = arglist$estidx
	idx = arglist$model$pidx
	ipars[estidx, 1] = pars
	distribution = arglist$model$modeldesc$distribution
	con = .persistgjrgarch1(pars = ipars[,1], idx = idx, distribution = distribution)
	if(is.na(con)) con = 1
	return(con)
}

.egarchcon = function(pars, arglist){
	ipars = arglist$ipars
	idx = arglist$model$pos.matrix
	betas  = pars[idx["beta",1]:idx["beta",2]]
	con = min(Mod(polyroot(c(1, -betas))))
	return(con)
}


# Constraints from Section 3.1 of Engle and Lee Paper
.csgarchcon = function(pars, arglist){
	ipars = arglist$ipars
	estidx = arglist$estidx
	idx = arglist$model$pidx
	ipars[estidx, 1] = pars
	modelinc = arglist$model$modelinc
	eta1 = ipars["eta11",1]
	eta2 = ipars["eta21",1]
	if(modelinc[9]>0) beta = sum(ipars[idx["beta",1]:idx["beta",2],1]) else beta = 0
	if(modelinc[8]>0) alpha = sum(ipars[idx["alpha",1]:idx["alpha",2],1]) else alpha = 0
	# \rho - alpha - beta > 0
	c1 = eta1 - alpha - beta
	# beta - \phi > 0
	c2 = beta - eta2
	return(c(c1, c2))
	#return( 1-  ( (alpha + beta)*(1-eta1) + eta1 ) )
}
