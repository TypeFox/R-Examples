#################################################################################
##
##   R package parma
##   Alexios Ghalanos Copyright (C) 2012-2013 (<=Aug)
##   Alexios Ghalanos and Bernhard Pfaff Copyright (C) 2013- (>Aug)
##   This file is part of the R package parma.
##
##   The R package parma is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package parma is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################

#  c("datatype", "benchmark", "targettype", "risk", "risktype", "leverage", "aux1", "aux2")
# targettype: 1 = inequality, 2 = equality
# risk: 1 (MAD), 2 (MiniMax), 3 (CVaR), 4 (CDaR), 5(EV), 6 (LPM), 7 (LPMUPM)
# risktype: 1 = minrisk, 2 = optimal (fractional)



################################################################
# General Penalty Functions
################################################################
penalty.minfun = function(w, optvars, uservars){
	m = optvars$wm
	pen = optvars$penalty
	rfun = optvars$rfun
	if(!is.null(dim(w)[2])){
		f = apply(w, 2, FUN = function(y){
					f = rfun(y, optvars, uservars) + 
							ifelse(optvars$index[3]==1 || !is.null(optvars$ineqfun), 
									pen * sum(pmax(parma.ineq.minfun(y, optvars, uservars), 0)^2), 0) + 
							pen * sum(parma.eq.minfun(y, optvars, uservars)^2)
					if(!is.null(optvars$max.pos)) f = f + pen * pmax(0, (sum((abs(y[optvars$widx])>0.001)) - optvars$max.pos) )^2
					return(f)
				})
	} else{
		f = rfun(w, optvars, uservars) + 
				ifelse(optvars$index[3]==1 || !is.null(optvars$ineqfun), 
						pen * sum(pmax(parma.ineq.minfun(w, optvars, uservars), 0)^2), 0) + 
				pen * sum(parma.eq.minfun(w, optvars, uservars)^2)
		if(!is.null(optvars$max.pos)) f = f + pen * pmax(0, sum((abs(w[optvars$widx])>0.001)) - optvars$max.pos)^2
	}
	return( f )
}


penalty.optfun = function(w, optvars, uservars){
	m = optvars$m
	pen = optvars$penalty
	rfun = optvars$rfun
	if(!is.null(dim(w)[2])){
		f = apply(w, 2, FUN = function(y){ 
					f = rfun(y, optvars, uservars) + 
							pen * sum(pmax(parma.ineq.optfun(y, optvars, uservars), 0)^2) + 
							pen * sum(parma.eq.optfun(y, optvars, uservars)^2)
					if(!is.null(optvars$max.pos)) f = f + pen * pmax(0, sum((abs(y[optvars$widx]/y[optvars$midx])>0.001)) - optvars$max.pos)^2
					return(f)
					
				})
	} else{
		f = rfun(w, optvars, uservars) + 
				pen * sum(pmax(parma.ineq.optfun(w, optvars, uservars), 0)^2) + 
				pen * sum(parma.eq.optfun(w, optvars, uservars)^2)
		if(!is.null(optvars$max.pos)) f = f + pen * pmax(0, sum((abs(w[optvars$widx]/w[optvars$midx])>0.001)) - optvars$max.pos)^2
	}
	return( f )
}

# custom penalty function for lpmupm (because of inequality containing the 
# nonlinear reward function of the Upper Partial Moment).
func.lpmupm = function(w, optvars, uservars){
	Data = optvars$Data
	lmoment = optvars$lmoment
	lthreshold = optvars$lthreshold
	port = Data %*% as.numeric(w[optvars$widx])
	f =  (mean(func.max.smooth( (w[optvars$midx]*lthreshold) - port)^lmoment))^(1/lmoment)
	return( f )
}

func.ineq.lpmupm = function(w, optvars, uservars){
	Data = optvars$Data
	umoment = optvars$umoment
	uthreshold = optvars$uthreshold
	N = optvars$N
	LB = optvars$LB
	UB = optvars$UB
	port = Data %*% as.numeric(w[optvars$widx])
	A1 =  1 - (mean(func.max.smooth(port - uthreshold*w[optvars$midx])^umoment))^(1/umoment)
	A2 = ( c( w[optvars$midx]*LB - w[optvars$widx], w[optvars$widx] - w[optvars$midx]*UB) )
	cons = c(A1, A2)
	if(!is.null(optvars$ineqfun)){
		n = length(optvars$ineqfun)
		fnlist = list(w, optvars, uservars)
		names(fnlist) = c("w", "optvars", "uservars")
		for(i in 1:n){
			cons = c(cons, do.call(optvars$ineqfun[[i]], args = fnlist))
		}
	}
	return( cons )
}

parma.eq.lpmupm = function(w, optvars, uservars){
	if(optvars$index[6]==0) cons = eqfun.budget.opt(w, optvars, uservars) else cons = eqfun.leverage.opt(w, optvars, uservars)
	if(!is.null(optvars$eqfun)){
		n = length(optvars$eqfun)
		fnlist = list(w, optvars, uservars)
		names(fnlist) = c("w", "optvars", "uservars")
		for(i in 1:n){
			cons = c(cons, do.call(optvars$eqfun[[i]], args = fnlist))
		}
	}
	return(cons)
}
penalty.lpmupm = function(w, optvars, uservars){
	pen = optvars$penalty
	if(!is.null(dim(w)[2])){
		f = apply(w, 2, FUN = function(y){
					f = func.lpmupm(y, optvars, uservars) + 
							pen * sum(pmax(func.ineq.lpmupm(y, optvars, uservars), 0)^2) + 
							pen * sum(parma.eq.lpmupm(y, optvars, uservars)^2)
					if(!is.null(optvars$max.pos)) f = f + pen * pmax(0, sum((abs(y[optvars$widx]/y[optvars$midx])>0.001)) - optvars$max.pos)^2
					return(f)})
	} else{
		f = func.lpmupm(w, optvars, uservars) + 
				pen * sum(pmax(func.ineq.lpmupm(w, optvars, uservars), 0)^2) + 
				pen * sum(parma.eq.lpmupm(w, optvars, uservars)^2)
		if(!is.null(optvars$max.pos)) f = f + pen * pmax(0, sum((abs(w[optvars$widx]/w[optvars$midx])>0.001)) - optvars$max.pos)^2
	}
	return(f)
}

gnlpport = function(optvars, uservars, solver = "cmaes", control = list(), ...)
{
	ctrl = .gnlp.ctrl(solver, control)
	risk = c("mad", "minimax", "cvar", "cdar", "ev", "lpm", "lpmupm")[optvars$index[4]]
	if(optvars$index[5] == 1){
		func = penalty.minfun
		rfun = switch(risk, 
				"cvar" = func.mincvar,
				"mad"  = func.minmad, 
				"lpm"  = func.minlpm,
				"minimax" = func.minminmax,
				"ev" = func.minvar)
		optvars$rfun = rfun
		solution = switch(solver, 
				"cmaes" = gnlp.mincmaes(func, optvars, uservars, ctrl),
				"crs"   = gnlp.mincrs(func, optvars, uservars, ctrl))
	} else{
		if(risk == "lpmupm"){
			func = penalty.lpmupm
		} else{
			func = penalty.optfun
			rfun = switch(risk, 
					"cvar" = func.optcvar,
					"mad"  = func.optmad, 
					"lpm"  = func.optlpm,
					"minimax" = func.optminmax,
					"ev" = func.optvar)
			optvars$rfun = rfun			
		}
		solution = switch(solver, 
				"cmaes" = gnlp.fraccmaes(func, optvars, uservars, ctrl),
				"crs"   = gnlp.fraccrs(func, optvars, uservars, ctrl))
	}
	return( solution )
}

################################################################################
.gnlp.ctrl = function(solver, control){
	ans = switch(solver, 
			"cmaes" = .cmaes.ctrl(control),
			"crs"   = .crs.ctrl(control))
	return(ans)
}

.cmaes.ctrl = function(control){
	cmaes.control(options = control$options, CMA = control$CMA)
}

.crs.ctrl = function(control){
	ctrl2 = list()
	ctrl2$algorithm = "NLOPT_GN_CRS2_LM"
	if( is.null(control$minf_max) ) ctrl2$minf_max = 0 else ctrl2$minf_max = control$minf_max
	if( is.null(control$ftol_rel) ) ctrl2$ftol_rel = NULL else ctrl2$ftol_rel = control$ftol_rel
	if( is.null(control$ftol_abs) ) ctrl2$ftol_abs = 1e-12 else  ctrl2$ftol_abs = control$ftol_abs
	if( is.null(control$xtol_rel) ) ctrl2$xtol_rel = 1e-11 else  ctrl2$xtol_rel = control$xtol_rel
	if( is.null(control$maxeval) ) ctrl2$maxeval = 100000 else  ctrl2$maxeval = as.integer( control$maxeval )
	if( is.null(control$maxtime) ) ctrl2$maxtime = 200 else  ctrl2$maxtime = control$maxtime
	if( is.null(control$print_level) ) ctrl2$print_level = 0 else  ctrl2$print_level = as.integer( control$print_level )
	return(ctrl2)
}

gnlp.mincmaes = function(func, optvars, uservars, ctrl){
	# use cmaes
	# refine once with nlminb
	sol = list()
	ans = try(cmaes(pars = optvars$x0, fun = func, optvars = optvars, uservars = uservars, 
					lower = optvars$LB, upper = optvars$UB, insigma = 0.1, 
					ctrl = ctrl), silent = TRUE)
	if(inherits(ans, 'try-error')){
		sol$status = 1
		sol$weights = rep(NA, optvars$wm)
		sol$risk = NA
		sol$reward = NA
		if(optvars$index[4]==3) sol$VaR = NA
		if(optvars$index[4]==4) sol$DaR = NA
	} else {
		ansf = nlminb(start = ans$par, func, lower = optvars$LB, upper = optvars$UB,  
				optvars = optvars, uservars = uservars, 
				control = list(trace= ifelse(is.null(ctrl$trace), 0, ctrl$trace), 
						eval.max = 2000, iter.max = 1000))
		sol$status = ansf$convergence
		sol$weights = ansf$par[optvars$widx]
		sol$risk = ansf$objective
		sol$reward = sum(sol$weights * optvars$mu)
		if(optvars$index[4]==3) sol$VaR = ansf$par[optvars$vidx]
		if(optvars$index[4]==4) sol$DaR = ansf$par[optvars$vidx]
	}
	return(sol)
}

gnlp.mincrs = function(func, optvars, uservars, ctrl){
	# use cmaes
	# refine once with nlminb
	sol = list()
	
	ans = try(nloptr( optvars$x0, eval_f = func, lb = optvars$LB, ub = optvars$UB, 
					opts = ctrl, optvars = optvars, uservars = uservars), 
			silent = TRUE)
	if(inherits(ans, 'try-error')){
		sol$status = 1
		sol$weights = rep(NA, optvars$wm)
		sol$risk = NA
		sol$reward = NA
		if(optvars$index[4]==3) sol$VaR = NA
		if(optvars$index[4]==4) sol$DaR = NA
	} else {
		ansf = nlminb(start = ans$solution, func, lower = optvars$LB, 
				upper = optvars$UB, optvars = optvars, uservars = uservars, 
				control = list(trace= ifelse(is.null(ctrl$trace), 0, ctrl$trace), 
						eval.max = 2000, iter.max = 1000))
		sol$status = ansf$convergence
		sol$weights = ansf$par[optvars$widx]
		sol$risk = ansf$objective
		sol$reward = sum(sol$weights * optvars$mu)
		if(optvars$index[4]==3) sol$VaR = ansf$par[optvars$vidx]
		if(optvars$index[4]==4) sol$DaR = ansf$par[optvars$vidx]
	}
	return(sol)
}

gnlp.fraccmaes = function(func, optvars, uservars, ctrl){
	# use cmaes
	# refine once with nlminb
	sol = list()
	ans = try(cmaes(pars = optvars$x0, fun = func, optvars = optvars, 
					uservars = uservars, lower = optvars$fLB, 
					upper = optvars$fUB, insigma = 0.1, ctrl = ctrl), 
			silent = TRUE)
	
	
	if(inherits(ans, 'try-error')){
		sol$status = 1
		sol$weights = rep(NA, optvars$wm)
		sol$risk = NA
		sol$reward = NA
		sol$multiplier = NA
		if(optvars$index[4]==3) sol$VaR = NA
		if(optvars$index[4]==4) sol$DaR = NA
	} else {
		ansf = nlminb(start = ans$par, func, lower = optvars$LB, upper = optvars$UB,  
				optvars = optvars, uservars = uservars, 
				control = list(trace= ifelse(is.null(ctrl$trace), 0, ctrl$trace), 
						eval.max = 2000, iter.max = 1000))
		sol$status = ansf$convergence
		sol$multiplier = ansf$par[optvars$midx]
		sol$weights = ansf$par[optvars$widx]/ansf$par[optvars$midx]
		sol$risk = ansf$objective/ansf$par[optvars$midx]
		sol$reward = sum(sol$weights * optvars$mu)
		if(optvars$index[4]==3) sol$VaR = ansf$par[optvars$vidx]/ansf$par[optvars$midx]
		if(optvars$index[4]==4) sol$DaR = ansf$par[optvars$vidx]/ansf$par[optvars$midx]
	}
	return(sol)	
}

gnlp.fraccrs = function(func, optvars, uservars, ctrl){
	ctrl2 = list()
	ctrl2$algorithm = "NLOPT_GN_CRS2_LM"
	sol = list()
	if( is.null(ctrl$minf_max) ) ctrl2$minf_max = 0 else ctrl2$minf_max = ctrl$minf_max
	if( is.null(ctrl$ftol_rel) ) ctrl2$ftol_rel = NULL else ctrl2$ftol_rel = ctrl$ftol_rel
	if( is.null(ctrl$ftol_abs) ) ctrl2$ftol_abs = 1e-12 else  ctrl2$ftol_abs = ctrl$ftol_abs
	if( is.null(ctrl$xtol_rel) ) ctrl2$xtol_rel = 1e-11 else  ctrl2$xtol_rel = ctrl$xtol_rel
	if( is.null(ctrl$maxeval) ) ctrl2$maxeval = 100000 else  ctrl2$maxeval = as.integer( ctrl$maxeval )
	if( is.null(ctrl$maxtime) ) ctrl2$maxtime = 200 else  ctrl2$maxtime = ctrl$maxtime
	if( is.null(ctrl$print_level) ) ctrl2$print_level = 0 else  ctrl2$print_level = as.integer( ctrl$print_level )
	
	ans = try(nloptr( optvars$x0, eval_f = func, lb = optvars$fLB, ub = optvars$fUB, 
					opts = ctrl, optvars = optvars, uservars = uservars), 
			silent = TRUE)
	if(inherits(ans, 'try-error')){
		sol$status = 1
		sol$weights = rep(NA, optvars$wm)
		sol$risk = NA
		sol$reward = NA
		sol$multiplier = NA
		if(optvars$index[4]==3) sol$VaR = NA
		if(optvars$index[4]==4) sol$DaR = NA
	} else {
		ansf = nlminb(start = ans$solution, func, lower = optvars$LB, 
				upper = optvars$UB, optvars = optvars, uservars = uservars, 
				control = list(trace= ifelse(is.null(ctrl$trace), 0, ctrl$trace), 
						eval.max = 2000, iter.max = 1000))
		sol$status = ansf$convergence
		sol$multiplier = ansf$par[optvars$midx]
		sol$weights = ansf$par[optvars$widx]/ansf$par[optvars$midx]
		sol$risk = ansf$objective/ansf$par[optvars$midx]
		sol$reward = sum(sol$weights * optvars$mu)
		if(optvars$index[4]==3) sol$VaR = ansf$par[optvars$vidx]/ansf$par[optvars$midx]
		if(optvars$index[4]==4) sol$DaR = ansf$par[optvars$vidx]/ansf$par[optvars$midx]
	}
	return(sol)
}
