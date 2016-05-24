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


#--------------------------------------------------------------------------------
# Rachev Ratio: Upper to Lower CVaR
# Research Reference (Rachev Ratio) : Biglova, A., Ortobelli, S., Rachev, S., Stoyanov (2004)
# Research Reference (Mixed Integer LP representation) : Konno, Tanaka, Yamamoto (2008)
#--------------------------------------------------------------------------------
# This expresses the fractional function as a MILP problem by applying the Charnes-Cooper Transformation
# Alexios Ghalanos 2010

# Completely unrealistic to expect to solve for scenario size>200 and alpha>0.01
# (at least with GLPK...cplex has more power but still very difficult).

milpmodrachev = function(optvars)
{
	Data = optvars$Data
	n = NROW(Data)
	m = NCOL(Data)
	colnames(Data) = NULL
	rownames(Data) = NULL
	probability = optvars$probability
	if(is.null(probability)) probability = rep(1/n, n)
	Data = scale(Data, scale = FALSE)
	
	upa = max(1, as.integer(optvars$upper.alpha * n))
	lpa = max(1, as.integer(optvars$lower.alpha * n))
	
	# We need specialized inequalities and parameter bounds inequalities multiplied by y0
	
	# S, alpha, theta, u_t, v_t, phi_t, psi_t, z_t, slack1_t, slack2_t, constant, y0
	
	# Constraint1: theta + sum(psi_t)/(1 - upper.alpha)*n = 1
	Amat = cbind(
			matrix(0, ncol = m, nrow = 1), 
			matrix(0, ncol = 1, nrow = 1), 
			matrix(1, ncol = 1, nrow = 1),
			matrix(0, ncol = n, nrow = 1),
			matrix(0, ncol = n, nrow = 1),
			matrix(0, ncol = n, nrow = 1),
			matrix(probability/upa, ncol = n, nrow = 1),
			matrix(0, ncol = n, nrow = 1),
			matrix(0, ncol = n, nrow = 1),
			matrix(0, ncol = n, nrow = 1),
			matrix(0, ncol = 1, nrow = 1),
			matrix(0, ncol = 1, nrow = 1)
	)
	
	# Constraint2: u_t - v_t + sum(w * S_t) + a = 0
	Amat = rbind( Amat, cbind( 
					Data,
					matrix(1, ncol = 1, nrow = n),
					matrix(0, ncol = 1, nrow = n),
					diag(n),
					-1 * diag(n),
					matrix(0, ncol = n, nrow = n),
					matrix(0, ncol = n, nrow = n),
					matrix(0, ncol = n, nrow = n),
					matrix(0, ncol = n, nrow = n),
					matrix(0, ncol = n, nrow = n),
					matrix(0, ncol = 1, nrow = n),
					matrix(0, ncol = 1, nrow = n)) )
	
	# Constraint3: phi_t - psi_t - sum(w * S_t) + theta = 0
	Amat = rbind( Amat, cbind( 
					-1 * Data,
					matrix(0, ncol = 1, nrow = n),
					matrix(1, ncol = 1, nrow = n),
					matrix(0, ncol = n, nrow = n),
					matrix(0, ncol = n, nrow = n),
					diag(n),
					-1 * diag(n),
					matrix(0, ncol = n, nrow = n),
					matrix(0, ncol = n, nrow = n),
					matrix(0, ncol = n, nrow = n),
					matrix(0, ncol = 1, nrow = n),
					matrix(0, ncol = 1, nrow = n)) )
	
	# Constraint4: phi_t - 1000*z_t + slack1_t = 0
	Amat = rbind( Amat, cbind( 
					matrix(0, ncol = m, nrow = n),
					matrix(0, ncol = 1, nrow = n),
					matrix(0, ncol = 1, nrow = n),
					matrix(0, ncol = n, nrow = n),
					matrix(0, ncol = n, nrow = n),
					diag(n),
					matrix(0, ncol = n, nrow = n),
					-1000 * diag(n),
					diag(n),
					matrix(0, ncol = n, nrow = n),
					matrix(0, ncol = 1, nrow = n),
					matrix(0, ncol = 1, nrow = n)) )
	
	# Constraint5: psi_t - 1000 + 1000*z_t + slack3_t = 0
	Amat = rbind( Amat, cbind( 
					matrix(0, ncol = m, nrow = n),
					matrix(0, ncol = 1, nrow = n),
					matrix(0, ncol = 1, nrow = n),
					matrix(0, ncol = n, nrow = n),
					matrix(0, ncol = n, nrow = n),
					matrix(0, ncol = n, nrow = n),
					diag(n),
					1000 * diag(n),
					matrix(0, ncol = n, nrow = n),
					diag(n),
					matrix(-1000, ncol = 1, nrow = n),
					matrix(0, ncol = 1, nrow = n) ) )
	
	# Constraint6: sum(z_t) = (1 - b) * n
	Amat = rbind( Amat, cbind( 
					matrix(0, ncol = m, nrow = 1),
					matrix(0, ncol = 1, nrow = 1),
					matrix(0, ncol = 1, nrow = 1),
					matrix(0, ncol = n, nrow = 1),
					matrix(0, ncol = n, nrow = 1),
					matrix(0, ncol = n, nrow = 1),
					matrix(0, ncol = n, nrow = 1),
					matrix(1, ncol = n, nrow = 1),
					matrix(0, ncol = n, nrow = 1),
					matrix(0, ncol = n, nrow = 1),
					matrix(0, ncol = 1, nrow = 1),
					matrix(0, ncol = 1, nrow = 1) ) )
	
	# Constraint7: sum(x_t) = budget*y0
	Amat = rbind( Amat, cbind(
					matrix(1, ncol = m, nrow = 1),
					matrix(0, ncol = 1, nrow = 1),
					matrix(0, ncol = 1, nrow = 1),
					matrix(0, ncol = n, nrow = 1),
					matrix(0, ncol = n, nrow = 1),
					matrix(0, ncol = n, nrow = 1),
					matrix(0, ncol = n, nrow = 1),
					matrix(0, ncol = n, nrow = 1),
					matrix(0, ncol = n, nrow = 1),
					matrix(0, ncol = n, nrow = 1),
					matrix(0, ncol = 1, nrow = 1),
					matrix(-optvars$budget, ncol = 1, nrow = 1) ) )
	
	# Constraint8: x<=UB.y0 i.e. x - UB.y0 <=0
	Amat = rbind( Amat, cbind( 
					diag(m),
					matrix(0, ncol = 1, nrow = m),
					matrix(0, ncol = 1, nrow = m),
					matrix(0, ncol = n, nrow = m),
					matrix(0, ncol = n, nrow = m),
					matrix(0, ncol = n, nrow = m),
					matrix(0, ncol = n, nrow = m),
					matrix(0, ncol = n, nrow = m),
					matrix(0, ncol = n, nrow = m),
					matrix(0, ncol = n, nrow = m),
					matrix(0, ncol = 1, nrow = m),
					matrix(-optvars$UB, ncol = 1, nrow = m) ) )
	
	# Constraint9: x - LB.y0 >=0
	Amat = rbind( Amat, cbind(
					diag(m),
					matrix(0, ncol = 1, nrow = m),
					matrix(0, ncol = 1, nrow = m),
					matrix(0, ncol = n, nrow = m),
					matrix(0, ncol = n, nrow = m),
					matrix(0, ncol = n, nrow = m),
					matrix(0, ncol = n, nrow = m),
					matrix(0, ncol = n, nrow = m),
					matrix(0, ncol = n, nrow = m),
					matrix(0, ncol = n, nrow = m),
					matrix(0, ncol = 1, nrow = m),
					matrix(-optvars$LB, ncol = 1, nrow = m) ) )
	
	# Constraint10: any user equalities
	if(!is.null(optvars$eq.mat)){
		eqm = dim(optvars$eq.mat)[1]
		Amat = rbind( Amat, cbind(
						optvars$eq.mat,
						matrix(0, ncol = 1, nrow = eqm),
						matrix(0, ncol = 1, nrow = eqm),
						matrix(0, ncol = n, nrow = eqm),
						matrix(0, ncol = n, nrow = eqm),
						matrix(0, ncol = n, nrow = eqm),
						matrix(0, ncol = n, nrow = eqm),
						matrix(0, ncol = n, nrow = eqm),
						matrix(0, ncol = n, nrow = eqm),
						matrix(0, ncol = n, nrow = eqm),
						matrix(0, ncol = 1, nrow = eqm),
						matrix(-optvars$eqB, ncol = 1, nrow = eqm) ) )
	} else{
		eqm = 0
	}
	
	# Constraint11: any user inequalities
	if(!is.null(optvars$ineq.mat)){
		ineqm = dim(optvars$ineq.mat)[1]
		Amat = rbind( Amat, cbind(
						optvars$ineq.mat,
						matrix(0, ncol = 1, nrow = ineqm),
						matrix(0, ncol = 1, nrow = ineqm),
						matrix(0, ncol = n, nrow = ineqm),
						matrix(0, ncol = n, nrow = ineqm),
						matrix(0, ncol = n, nrow = ineqm),
						matrix(0, ncol = n, nrow = ineqm),
						matrix(0, ncol = n, nrow = ineqm),
						matrix(0, ncol = n, nrow = ineqm),
						matrix(0, ncol = n, nrow = ineqm),
						matrix(0, ncol = 1, nrow = ineqm),
						matrix(-optvars$ineq.UB, ncol = 1, nrow = ineqm) ) )
		
		Amat = rbind( Amat, cbind(
						optvars$ineq.mat,
						matrix(0, ncol = 1, nrow = ineqm),
						matrix(0, ncol = 1, nrow = ineqm),
						matrix(0, ncol = n, nrow = ineqm),
						matrix(0, ncol = n, nrow = ineqm),
						matrix(0, ncol = n, nrow = ineqm),
						matrix(0, ncol = n, nrow = ineqm),
						matrix(0, ncol = n, nrow = ineqm),
						matrix(0, ncol = n, nrow = ineqm),
						matrix(0, ncol = n, nrow = ineqm),
						matrix(0, ncol = 1, nrow = ineqm),
						matrix(-optvars$ineq.LB, ncol = 1, nrow = ineqm) ) )
	} else{
		ineqm = 0
	}
	
	Amat = as.simple_triplet_matrix(Amat)
	
	rhs = c( 1, rep(0, n), rep(0, n), rep(0, 2*n), upa , 0, rep(0, m), rep(0, m), rep(0, eqm), rep(0, 2*ineqm) )	
	
	dir = c( "==", rep("==", 4*n), "==", rep("<=", 2 * ineqm), rep("==", eqm), "==", 
			rep("<=", m), rep(">=", m), rep("==", eqm), rep("<=", ineqm), rep(">=", ineqm) )
	
	# S, alpha, theta, u_t, v_t, phi_t, psi_t, z_t, slack1_t, slack2_t,  constant
	objL = c( rep(0, m), 1, 0, rep( probability/lpa, n ), rep(0, 6*n), 0 , 0 )
	
	xLB = c( rep(-Inf, m), 0, 0, rep(0, n), rep(0, n), rep(0, n), rep(0, n), rep(0, n), rep(0, 2*n), 1, 0.05)
	xUB = c( rep(Inf, m), Inf, Inf, rep(Inf, n), rep(Inf, n), rep(1, n), rep(1, 2 * n), rep(Inf, 2*n), 1, Inf)
	types = c(rep("C", m), "C", "C", rep("C", 4 * n), rep("B", n), rep("C", 2 * n), "C", "C")
	bounds = list( lower = list( ind = 1L:(4 + m + 7 * n), val = xLB ),
			upper = list( ind = 1L:(4 + m + 7 * n), val = xUB ) )
	
	ans = list(objL = objL, Amat = Amat, dir = dir, rhs = rhs, types = types, 
			bounds = bounds, dataidx = c(n, m), widx = 1:m)
	rm(Amat)
	rm(Data)
	gc()
	return(ans)
}


milplpm1.minrisk = function(optvars)
{
	# -[ Nxm ] -[ diag N ] + [1] = m+N+1
	Data = optvars$Data
	n = NROW(Data)
	m = NCOL(Data)
	colnames(Data) = NULL
	rownames(Data) = NULL
	probability = optvars$probability
	if(is.null(probability)) probability = rep(1/n, n)
	# -[ Nxm ] -[ diag N ] + [1] = m+N+1
	cons = lpmin.setup(optvars, cdimn = n+1)
	ineqcon = cons$ineqcon
	ineqB = cons$ineqB
	ineqn = cons$ineqn
	ineqm = cons$ineqm
	eqcon = cons$eqcon
	eqB = cons$eqB
	eqn = cons$eqn
	eqm = cons$eqm
	LB = optvars$LB
	UB = optvars$UB
	Amat = cbind( -Data, -diag(n), matrix(optvars$threshold, nrow = n, ncol = 1))
	Amat = rbind( Amat, ineqcon, eqcon )
	Amat = cbind( Amat, matrix(0, ncol = m, nrow = dim(Amat)[1]))
	Amat = rbind( Amat, cbind(-diag(m), matrix(0, ncol = n, nrow = m), matrix(0, ncol = 1, nrow = m), 
					diag(LB)) )
	Amat = rbind( Amat, cbind( diag(m), matrix(0, ncol = n, nrow = m), matrix(0, ncol = 1, nrow = m), 
					diag(-UB)) )
	Amat = rbind( Amat, cbind(matrix(0, ncol = m + n + 1, nrow = 1), matrix(1, ncol = m, nrow = 1) ) )	
	Amat = as.simple_triplet_matrix(Amat)
	
	rhs = rep(0, n)
	rhs = c( rhs, as.numeric(ineqB), as.numeric(eqB) , rep(0, 2*m), optvars$max.pos)
	dir = c( rep("<=", n), rep("<=", ineqn), rep("==", eqn), rep("<=", 2*m), "==")
	objL = c( rep(0, m), probability, 1 , rep(0, m))
	types = c( rep("C", m), rep("C", n), "C", rep("B", m))
	bounds = list( lower = list( ind = 1L:(m + n + 1 + m), val = c(LB,  rep(0, n), 1 , rep(0, m)) ),
			upper = list( ind = 1L:(m + n + 1 + m), val = c( UB, rep(Inf, n), 1, rep(1, m) ) ) )	
	ans = list(objL = objL, Amat = Amat, dir = dir, rhs = rhs, bounds = bounds, 
			problem.max = FALSE, types = types, 
			reward = optvars$mu, dataidx = c(n, m), widx = 1:m)
	rm(Amat)
	rm(Data)
	gc()
	return( ans )
}

milpmad.minrisk = function(optvars)
{
	Data = optvars$Data
	n = NROW(Data)
	m = NCOL(Data)
	colnames(Data) = NULL
	rownames(Data) = NULL
	probability = optvars$probability
	if(is.null(probability)) probability = rep(1/n, n)
	# -[ Nxm ] -[ diag N ] + [1] = m+N+1
	cons = lpmin.setup(optvars, cdimn = 2*n+1)
	ineqcon = cons$ineqcon
	ineqB = cons$ineqB
	ineqn = cons$ineqn
	ineqm = cons$ineqm
	eqcon = cons$eqcon
	eqB = cons$eqB
	eqn = cons$eqn
	eqm = cons$eqm
	LB = optvars$LB
	UB = optvars$UB
	
	Amat = cbind( -Data, matrix(optvars$benchmark, nrow = n, ncol = 1), diag(n), -diag(n) )
	Amat = rbind( Amat, ineqcon, eqcon )
	Amat = cbind( Amat, matrix(0, ncol = m, nrow = dim(Amat)[1]))

	Amat = rbind( Amat, cbind(-diag(m), matrix(0, ncol = 1, nrow = m), matrix(0, ncol = 2*n, nrow = m),
					diag(LB)) )
	Amat = rbind( Amat, cbind( diag(m), matrix(0, ncol = 1, nrow = m), matrix(0, ncol = 2*n, nrow = m),
					diag(-UB)) )
	Amat = rbind( Amat, cbind(matrix(0, ncol = m + 2*n + 1, nrow = 1), matrix(1, ncol = m, nrow = 1) ) )	
	Amat = as.simple_triplet_matrix(Amat)
	
	rhs = rep(0, n)
	
	rhs = c( rhs, as.numeric(ineqB), as.numeric(eqB) , rep(0, 2*m), optvars$max.pos)
	dir = c( rep("==", n), rep("<=", ineqn), rep("==", eqn) , rep("<=", 2*m), "==")
	objL = c( rep(0, m), 0, probability, probability, rep(0, m))
	bounds = list( lower = list( ind = 1L:(m + 2 * n + 1 + m), val = c(LB,  1, rep(0, 2 * n), rep(0, m) ) ),
			upper = list( ind = 1L:(m + 2 * n + 1 + m), val = c( UB, 1, rep(Inf, 2 * n), rep(1, m)) ) )	
	types = c( rep("C", m), "C", rep("C", 2*n), rep("B", m))
	
	ans = list(objL = objL, Amat = Amat, dir = dir, rhs = rhs, bounds = bounds, 
			problem.max = FALSE, reward = optvars$mu, 
			types = types, dataidx = c(n, m), widx = 1:m)
	rm(Amat)
	rm(Data)
	gc()
	return( ans )
}




milpminmax.minrisk = function(optvars)
{
	Data = optvars$Data
	n = NROW(Data)
	m = NCOL(Data)
	colnames(Data) = NULL
	rownames(Data) = NULL
	probability = optvars$probability
	if(is.null(probability)) probability = rep(1/n, n)
	# -[ Nxm ] -[ diag N ] + [1] = m+N+1
	cons = lpmin.setup(optvars, cdimn = 2)
	ineqcon = cons$ineqcon
	ineqB = cons$ineqB
	ineqn = cons$ineqn
	ineqm = cons$ineqm
	eqcon = cons$eqcon
	eqB = cons$eqB
	eqn = cons$eqn
	eqm = cons$eqm
	LB = optvars$LB
	UB = optvars$UB
	
	Amat = cbind( -Data, matrix(optvars$benchmark, nrow = n, ncol = 1), matrix(1, nrow = n, ncol = 1) )
	Amat = rbind( Amat, ineqcon, eqcon )
	Amat = cbind(Amat, matrix(0, ncol = m, nrow = dim(Amat)[1]))
	Amat = rbind( Amat, cbind(-diag(m), matrix(0, ncol = 2, nrow = m), diag(LB)) )
	Amat = rbind( Amat, cbind( diag(m), matrix(0, ncol = 2, nrow = m), diag(-UB)) )
	Amat = rbind( Amat, cbind(matrix(0, ncol = m + 2, nrow = 1), matrix(1, ncol = m, nrow = 1) ) )	
	Amat = as.simple_triplet_matrix(Amat)
	
	rhs = rep(0, n)
	rhs = c( rhs, as.numeric(ineqB), as.numeric(eqB), rep(0, 2*m), optvars$max.pos )
	
	dir = c( rep("<=", n), rep("<=", ineqn), rep("==", eqn), rep("<=", 2*m), "==")
	
	objL = c( rep(0, m), 0, -1 , rep(0, m))
	
	bounds = list( lower = list( ind = 1L:(m + 2 + m), val = c(LB,  1, -Inf, rep(0, m)) ),
			upper = list( ind = 1L:(m + 2 + m), val = c( UB, 1, Inf, rep(1, m)) ) )
	
	types = c( rep("C", m), rep("C", 2), rep("B", m))
	
	ans = list(objL = objL, Amat = Amat, dir = dir, rhs = rhs, bounds = bounds, 
			problem.max = FALSE, reward = optvars$mu,  types = types, 
			dataidx = c(n, m),  widx = 1:m, vidx = m+2)	
	rm(Amat)
	rm(Data)
	gc()
	return( ans )
}

milpcvar.minrisk = function(optvars)
{
	Data = optvars$Data
	n = NROW(Data)
	m = NCOL(Data)
	colnames(Data) = NULL
	rownames(Data) = NULL
	probability = optvars$probability
	if(is.null(probability)) probability = rep(1/n, n)
	# -[ Nxm ] -[ diag N ] + [1] = m+N+1
	cons = lpmin.setup(optvars, cdimn = n + 2)
	ineqcon = cons$ineqcon
	ineqB = cons$ineqB
	ineqn = cons$ineqn
	ineqm = cons$ineqm
	eqcon = cons$eqcon
	eqB = cons$eqB
	eqn = cons$eqn
	eqm = cons$eqm
	LB = optvars$LB
	UB = optvars$UB
	
	Amat = cbind( -Data,
			matrix(-1, nrow = n, ncol = 1), 
			-diag(n),
			matrix(optvars$benchmark, nrow = n, ncol = 1) )
	
	Amat = rbind( Amat, ineqcon, eqcon )
	Amat = cbind(Amat, matrix(0, ncol = m, nrow = dim(Amat)[1]))
	
	Amat = rbind( Amat, cbind(-diag(m), matrix(0, ncol = n+2, nrow = m), diag(LB)) )
	Amat = rbind( Amat, cbind( diag(m), matrix(0, ncol = n+2, nrow = m), diag(-UB)) )
	Amat = rbind( Amat, cbind(matrix(0, ncol = m + n + 2, nrow = 1), matrix(1, ncol = m, nrow = 1) ) )	
	Amat = as.simple_triplet_matrix(Amat)
	
	
	rhs = rep(0, n)
	rhs = c( rhs, as.numeric(ineqB), as.numeric(eqB), rep(0, 2*m), optvars$max.pos )
	
	dir = c( rep("<=", n), rep("<=", ineqn), rep("==", eqn), rep("<=", 2*m), "==")
	
	objL = c( rep(0, m), 1, probability/optvars$alpha, 0, rep(0, m))
	
	bounds = list( lower = list( ind = 1L:(m + n + 2 + m), val = c(LB,  -1, rep(0, n), 1, rep(0, m)) ),
			upper = list( ind = 1L:(m + n + 2 + m), val = c( UB, 1, rep(Inf, n), 1 , rep(1, m)) ) )
	
	types = c( rep("C", m), "C", rep("C", n), "C", rep("B", m))
	
	ans = list(objL = objL, Amat = Amat, dir = dir, rhs = rhs, bounds = bounds, 
			problem.max = FALSE, reward = optvars$mu, types = types, 
			dataidx = c(n, m), widx = 1:m, vidx = m+1)
	rm(Data)
	rm(Amat)
	gc()
	return( ans )	
}




milpcdar.minrisk = function(optvars)
{
	Data = optvars$Data
	n = NROW(Data)
	m = NCOL(Data)
	colnames(Data) = NULL
	rownames(Data) = NULL
	probability = optvars$probability
	if(is.null(probability)) probability = rep(1/n, n)
	# -[ Nxm ] -[ diag N ] + [1] = m+N+1
	cons = lpmin.setup(optvars, cdimn = 2*n + 2)
	ineqcon = cons$ineqcon
	ineqB = cons$ineqB
	ineqn = cons$ineqn
	ineqm = cons$ineqm
	eqcon = cons$eqcon
	eqB = cons$eqB
	eqn = cons$eqn
	eqm = cons$eqm
	LB = optvars$LB
	UB = optvars$UB
	
	xm = as.matrix( -diag(n) )
	idx = which(xm == -1, arr.ind = TRUE)
	xrow = idx[,1]
	xcol = idx[,2]
	xcol[2:length(xcol)] = xcol[2:length(xcol)] - 1
	diag(xm[ xrow[2:length(xrow)], xcol[2:length(xcol)] ]) = 1
	
	Amat = cbind( 
			-Data, 
			matrix(0, nrow = n, ncol = 1),
			matrix(0, nrow = n, ncol = n),
			as.matrix(xm),
			matrix(optvars$benchmark, nrow = n, ncol = 1))
	
	Amat = rbind( Amat, cbind(
					matrix(0, nrow = n, ncol = m), 
					matrix(-1, nrow = n, ncol = 1), 
					-1 * diag(n),
					diag(n),
					matrix(0, nrow = n, ncol = 1)) )
	
	Amat = rbind( Amat, ineqcon, eqcon )
	Amat = cbind(Amat, matrix(0, ncol = m, nrow = dim(Amat)[1]))
	
	Amat = rbind( Amat, cbind(-diag(m), matrix(0, ncol = 2*n+2, nrow = m), diag(LB)) )
	Amat = rbind( Amat, cbind( diag(m), matrix(0, ncol = 2*n+2, nrow = m), diag(-UB)) )
	Amat = rbind( Amat, cbind(matrix(0, ncol = m + 2*n + 2, nrow = 1), matrix(1, ncol = m, nrow = 1) ) )	
	Amat = as.simple_triplet_matrix(Amat)
	
	
	rhs = rep(0, 2 * n)
	rhs = c( rhs, as.numeric(ineqB), as.numeric(eqB), rep(0, 2*m), optvars$max.pos)
	
	dir = c( rep("<=", 2 * n), rep("<=", ineqn), rep("==", eqn), rep("<=", 2*m), "==")
	
	objL = c( rep(0, m), 1, probability/optvars$alpha, rep(0, n),  0, rep(0, m))
	
	bounds = list( lower = list( ind = 1L:(m + 2 * n + 2 + m), val = c( LB, -1, rep(0, 2 * n), 1, rep(0, m)) ),
			upper = list( ind = 1L:(m + 2 * n + 2 + m), val = c( UB, 1, rep(Inf, 2 * n), 1, rep(1, m) ) ) )	
	
	types = c( rep("C", m), "C", rep("C", 2*n), "C", rep("B", m))
	
	ans = list(objL = objL, Amat = Amat, dir = dir, rhs = rhs, bounds = bounds, 
			problem.max = FALSE, reward = optvars$mu, types = types, 
			dataidx = c(n, m), widx = 1:m, vidx = m+1)
	rm(Amat)
	rm(Data)
	gc()
	return( ans )	
}


milpport = function(optvars, ...)
{
	risk = c("mad", "minimax", "cvar", "cdar", "ev", "lpm", "lpmupm")[optvars$index[4]]
	setup = switch(risk,
				"mad"     = milpmad.minrisk(optvars),
				"minimax" = milpminmax.minrisk(optvars),
				"cvar"    = milpcvar.minrisk(optvars),
				"cdar"    = milpcdar.minrisk(optvars),
				"lpm"     = milplpm1.minrisk(optvars))
	sol = milpminsolver(setup, optvars, ...)
	
	#sol$risk = fun.risk(sol$weights, scenario, options, risk, benchmark)
	#sol$reward = sum(forecast * sol$weights)
	rm(setup)
	gc()
	return( sol )
}

milpminsolver = function(setup, optvars, ...)
{
	if( is.null(list(...)$verbose) ) verbose = FALSE else verbose = as.logical( list(...)$verbose )
	m = setup$dataidx[2]
	sol = try(Rglpk_solve_LP(obj = setup$objL, mat = setup$Amat, dir = setup$dir, 
					types = setup$types, rhs = setup$rhs, bounds = setup$bounds, 
					max = FALSE, verbose = verbose), silent = TRUE)
	if( inherits(sol, "try-error") ){
		status = "non-convergence"
		weights = rep(NA, m)
		risk = NA
		reward = NA
		if(optvars$index[4]==3) VaR = NA
		if(optvars$index[4]==4) DaR = NA
	} else{
		status  = sol$status
		weights = sol$solution[setup$widx]
		risk    = sol$optimum
		if( !is.null(setup$reward) ){
			reward = sum(setup$reward * weights)
		} else {
			reward = NULL
		}
		if(optvars$index[4]==3) VaR = sol$solution[setup$vidx]
		if(optvars$index[4]==4) DaR = sol$solution[setup$vidx]
	}
	rm(setup)
	gc()
	return( list(status = status, weights = weights, risk = risk, reward = reward, 
					multiplier = 1) )
}

# This is for the Rachev Ratio model
milpfracsolver = function(setup, optvars, ...)
{
	if( is.null(list(...)$verbose) ) verbose = FALSE else verbose = as.logical( list(...)$verbose )
	
	m = setup$dataidx[2]
	
	sol = try(Rglpk_solve_LP(obj = setup$objL, mat = setup$Amat, dir = setup$dir, 
					types = setup$type, rhs = setup$rhs, bounds = setup$bounds, 
					max = FALSE, verbose = verbose), silent = TRUE)
	ans = list()
	if( inherits(sol, "try-error") ){
		status = "non-convergence"
		weights = rep(NA, m)
		risk = NA
		reward = NA
		optimum = NA
		multiplier = NA
		VaR = NA
		DaR = NA
	} else{
		status = sol$status
		multiplier = sol$solution[setup$midx]
		tmp = sol$solution/multiplier
		weights = tmp[setup$widx]
		risk = sum(setup$objL * sol$solution)/multiplier
		reward = sum(optvars$mu * weights)
		optimum  = sol$optimum
		if(optvars$index[4]==3) VaR = sol$solution[setup$vidx]/multiplier else VaR = NA
		if(optvars$index[4]==4) DaR = sol$solution[setup$vidx]/multiplier else DaR = NA
	}
	
	
	
	rm(setup)
	
	return( list(status = status, weights = weights,  optimum = optimum) )
}