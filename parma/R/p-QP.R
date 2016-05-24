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


# min risk QP and opt risk QP

.qp.min = function(optvars){
	S = as.matrix(optvars$S)	
	m = dim(S)[2]
	colnames(S) = NULL
	rownames(S) = NULL
	index = optvars$index
	sbench = optvars$benchmarkS
	#mbench = optvars$benchmarkM
	target = optvars$mutarget
	eq.mat = optvars$eq.mat
	eqB = optvars$eqB
	ineq.mat = optvars$ineq.mat
	ineq.LB  = optvars$ineq.LB
	ineq.UB  = optvars$ineq.UB
	LB = optvars$LB
	UB = optvars$UB
	budget = optvars$budget
	forecast = optvars$mu
	if(index[2] == 1){
		# benchmark case
		Saux = matrix(0, ncol = m+1, nrow = m+1)
		Saux[2:(m+1),2:(m+1)] = S
		Saux[1,] = sbench
		Saux[2:(m+1),1] = sbench[-1]
		if(any(eigen(Saux, TRUE, only.values=TRUE)$values<1e-10)){
			warning("\nparma: benchmark relative covariance matrix not PD...adjusting to make PD")
			Saux <- make.positive.definite(Saux, 1e-10)
		}
		forecast = matrix(forecast, ncol = m, nrow = 1)
		if(index[3] == 1){
			# inequality target
			if(!is.null(eq.mat)) eq.mat = cbind(matrix(0, nrow = dim(eq.mat)), eq.mat)
			if(!is.null(ineq.mat)) ineq.mat = cbind(matrix(0, nrow = dim(ineq.mat)), ineq.mat)
			.ineqx = parma.con.qpmin1(eq.mat, eqB = eqB, reward = NULL, rewardB = NULL, 
					ineq = rbind(cbind(0, forecast), ineq.mat), 
					ineqLB = c(target, ineq.LB), ineqUB = c(10000, ineq.UB), 
					LB = c(-1, LB), UB = c(-1, UB), budget, m+1)
		} else{
			if(!is.null(eq.mat)) eq.mat = cbind(matrix(0, nrow = dim(eq.mat)), eq.mat)
			if(!is.null(ineq.mat)) ineq.mat = cbind(matrix(0, nrow = dim(ineq.mat)), ineq.mat)
			.ineqx = parma.con.qpmin1(eq.mat, eqB = eqB, reward = c(0, forecast), 
					rewardB = target, ineq = ineq.mat, ineqLB = ineq.LB, 
					ineqUB = ineq.UB, LB = c(-1, LB), UB = c(-1, UB), budget, m+1)
		}
		dvec = rep(0, m+1)
	} else{
		Saux = S
		forecast = matrix(forecast, ncol = m, nrow = 1)
		if(index[3]==1){
			.ineqx = parma.con.qpmin2(eq.mat, eqB = eqB, reward = NULL, rewardB = NULL, 
					ineq = rbind(forecast, ineq.mat), ineqLB = c(target, ineq.LB), 
					ineqUB = c(50000, ineq.UB), LB = LB, UB = UB, budget, m)
		} else{
			.ineqx = parma.con.qpmin2(eq = eq.mat, eqB = eqB, reward = forecast, 
					rewardB = target, ineq = ineq.mat, ineqLB = ineq.LB, 
					ineqUB = ineq.UB, LB = LB, UB = UB, budget, m)
		}
		dvec = rep(0, m)
	}
	Amat = .ineqx$Amat
	bvec = .ineqx$bvec
	meq = .ineqx$meq
	Dmat = Saux
	return( list(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = meq, 
					forecast = forecast, S = S, benchvar = sbench[1], m = m, 
					index = index))
}

.qp.opt = function(optvars){
	S = as.matrix(optvars$S)	
	m = dim(S)[2]
	colnames(S) = NULL
	rownames(S) = NULL
	index = optvars$index
	sbench = optvars$benchmarkS
	#mbench = optvars$benchmarkM
	target = optvars$mutarget
	eq.mat = optvars$eq.mat
	eqB = optvars$eqB
	ineq.mat = optvars$ineq.mat
	ineq.LB  = optvars$ineq.LB
	ineq.UB  = optvars$ineq.UB
	LB = optvars$LB
	UB = optvars$UB
	budget = optvars$budget
	forecast = optvars$mu
	
	Saux = matrix(0, ncol = m+1, nrow = m+1)
	Saux[2:(m+1),2:(m+1)] = S
	Saux[1,] = sbench
	Saux[2:(m+1),1] = sbench[-1]
	if(any(eigen(Saux, TRUE, only.values=TRUE)$values<1e-10)){
		warning("\nparma: benchmark relative covariance matrix not PD...adjusting to make PD")
		Saux <- make.positive.definite(Saux, 1e-10)
	}
	# this guarantees PD of matrix when no benchmark is used:
	if(optvars$index[2]==0) Saux[1,1] = 1e-12
	# special constraints for optimal QP:
	.ineqx = parma.con.qpopt(eq = eq.mat, eqB = eqB, reward = forecast, 
			ineq = ineq.mat, ineqLB = ineq.LB, ineqUB = ineq.UB, LB = LB, UB = UB, 
			budget = budget, m = m+1)
	Amat = .ineqx$Amat
	bvec = .ineqx$bvec
	meq = .ineqx$meq
	Dmat = Saux
	dvec = rep(0, m+1)
	return( list(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = meq, 
					forecast = forecast, S = S, benchvar = sbench[1], 
					m = m, index = index))
}

# 2 solvers: QP for standard optimization and SOCP for quadratic tail constraints
qpport = function(optvars, ...){
	if(optvars$index[5]==1){
		setup = .qp.min(optvars)
		ans = qpmin.solver(setup, optvars, ...)
	} else{
		setup = .qp.opt(optvars)
		ans = try(qpopt.solver(setup, optvars, ...), silent = TRUE)			
	}
	sol = list()
	sol$weights = ans$weights
	sol$reward  = ans$reward
	sol$risk 	= as.numeric(ans$risk)
	# if(!is.null(optvars$SS)) 
	#	sol$riskbudget = sapply(optvars$SS, FUN = function(x) ans$weights %*% x %*% ans$weights)
	return( sol )
}



qpmin.solver = function(setup, optvars, ...){
	sol = try(solve.QP(Dmat = setup$Dmat, dvec = setup$dvec, Amat = t(setup$Amat),
					bvec = setup$bvec, meq = setup$meq), silent = TRUE)
	if(inherits(sol, "try-error")){
		warning(sol)
		status = "non-convergence"
		weights = rep(NA, setup$m)
		risk = NA
		reward = NA
	} else{
		status = 0
		if(setup$index[2]==1){
			weights = sol$solution[-1]
		} else{
			weights = sol$solution
		}
		risk = fun.qvar(sol$solution, setup$Dmat)
		reward = sum(weights * optvars$mu)
	}
	return(list(weights = weights, risk = risk, reward = reward, measure = "var"))
}

qpopt.solver = function(setup, optvars, ...){
	sol = try(solve.QP(Dmat = setup$Dmat, dvec = setup$dvec, Amat = t(setup$Amat),
					bvec = setup$bvec, meq = setup$meq), silent = TRUE)
	if(inherits(sol, "try-error")){
		warning(sol)
		status = "non-convergence"
		weights = rep(NA, setup$m)
		risk = NA
		reward = NA
		multiplier = NA
	} else{
		status = 0
		weights = -1*(sol$solution[-1]/sol$solution[1])
		risk = fun.qvar(-sol$solution/sol$solution[1], setup$Dmat)
		multiplier = sol$solution[1]
		reward = sum(weights * optvars$mu)
	}
	return(list(weights = weights, risk = risk, reward = reward, 
					multiplier = multiplier, measure = "var"))
}