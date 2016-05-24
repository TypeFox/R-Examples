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

################################################################################
lpmin.setup = function(optvars, cdimn){
	Data = optvars$Data
	n = NROW(Data)
	m = NCOL(Data)
	colnames(Data) = NULL
	rownames(Data) = NULL
	index = optvars$index
	# create linear inequality constraints
	if(index[3] == 1){
		.ineqx = parma.ineq.lpmin(ineq = rbind(optvars$mu, optvars$ineq.mat), 
				ineqLB = c(optvars$mutarget, optvars$ineq.LB), 
				ineqUB = c(Inf, optvars$ineqUB), cdim = cdimn)
	} else{
		.ineqx = parma.ineq.lpmin(ineq = optvars$ineq.mat, ineqLB = optvars$ineq.LB, 
				ineqUB = optvars$ineq.UB, cdim = cdimn)
	}
	ineqcon = .ineqx$ineqcon
	ineqB = .ineqx$ineqB
	ineqn = .ineqx$ineqn
	ineqm = .ineqx$ineqm
	# create linear equality constraints
	if(index[3] == 1){
		.eqx = parma.eq.lpmin(optvars$eq.mat, optvars$eqB, reward = NULL, 
				rewardB = NULL, cdim = cdimn, optvars$budget, m)
	} else{
		.eqx = parma.eq.lpmin(optvars$eq.mat, optvars$eqB, reward = optvars$mu, 
				rewardB = optvars$mutarget, cdim = cdimn, 
				optvars$budget, m)
	}
	eqcon = .eqx$eqcon
	eqB = .eqx$eqB
	eqn = .eqx$eqn
	eqm = .eqx$eqm
	return(list(ineqcon = ineqcon, ineqB = ineqB, ineqn = ineqn, ineqm = ineqm,
					eqcon = eqcon, eqB = eqB, eqn = eqn, eqm = eqm))
}

lpopt.setup = function(optvars, cdimn){	
	Data = optvars$Data
	n = NROW(Data)
	m = NCOL(Data)
	colnames(Data) = NULL
	rownames(Data) = NULL
	index = optvars$index
	# create linear inequality constraints
	.ineqx = parma.ineq.lpopt(ineq = optvars$ineq.mat, ineqLB = optvars$ineq.LB, 
			ineqUB = optvars$ineq.UB, LB = optvars$LB, UB = optvars$UB, 
			cdim = cdimn)
	ineqcon = .ineqx$ineqcon
	ineqB = .ineqx$ineqB
	ineqn = .ineqx$ineqn
	ineqm = .ineqx$ineqm
	
	.eqx = parma.eq.lpopt(optvars$eq.mat, optvars$eqB, reward = NULL, 
			rewardB = NULL, cdim = cdimn, optvars$budget, m)
	eqcon = .eqx$eqcon
	eqB = .eqx$eqB
	eqn = .eqx$eqn
	eqm = .eqx$eqm

	return(list(ineqcon = ineqcon, ineqB = ineqB, ineqn = ineqn, ineqm = ineqm,
					eqcon = eqcon, eqB = eqB, eqn = eqn, eqm = eqm))
}

################################################################################
# LPM1
################################################################################
lplpm1.minrisk = function(optvars)
{
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
	Amat = as.simple_triplet_matrix( rbind( Amat, ineqcon, eqcon ) )
	rhs = rep(0, n)
	rhs = c( rhs, as.numeric(ineqB), as.numeric(eqB) )
	dir = c( rep("<=", n), rep("<=", ineqn), rep("==", eqn) )
	objL = c( rep(0, m), probability, 1 )
	bounds = list( lower = list( ind = 1L:(m + n + 1), val = c(LB,  rep(0, n), 1 ) ),
			upper = list( ind = 1L:(m + n + 1), val = c( UB, rep(Inf, n), 1 ) ) )	
	ans = list(objL = objL, Amat = Amat, dir = dir, rhs = rhs, bounds = bounds, 
			problem.max = FALSE, reward = optvars$mu, 
			dataidx = c(n, m), widx = 1:m)
	rm(Amat)
	rm(Data)
	gc()
	return( ans )
}

# there is no benchmark in this problem since the presence of the threshold makes it impossible
# to have them both
lplpm1.optrisk = function(optvars)
{
	Data = optvars$Data
	n = NROW(Data)
	m = NCOL(Data)
	colnames(Data) = NULL
	rownames(Data) = NULL
	probability = optvars$probability
	if(is.null(probability)) probability = rep(1/n, n)
	# -[ Nxm ] -[ diag N ] + [1] = m+N+1
	cons = lpopt.setup(optvars, cdimn = n+1)
	ineqcon = cons$ineqcon
	ineqB = cons$ineqB
	ineqn = cons$ineqn
	ineqm = cons$ineqm
	eqcon = cons$eqcon
	eqB = cons$eqB
	eqn = cons$eqn
	eqm = cons$eqm
	
	Amat = cbind( -Data, -diag(n), matrix(optvars$threshold, ncol = 1, nrow = n) )
	Amat = rbind( Amat, 
			cbind( matrix(optvars$mu, ncol = m, nrow = 1), 
					matrix(0, ncol = n, nrow = 1),
					matrix(0, ncol = 1, nrow = 1) ) )
	Amat = rbind( Amat, ineqcon, eqcon )
	Amat = as.simple_triplet_matrix(Amat)
	rhs = rep(0, n)
	rhs = c( rhs, 1, as.numeric(ineqB), as.numeric(eqB) )
	dir = c( rep("<=", n), "==", rep("<=", ineqn), rep("==", eqn) )
	objL = c( rep(0, m), probability, 0 )
	bounds = list( lower = list( ind = 1L:(m + n + 1), val = c( rep(-10000, m),  rep(0, n), 0.05 ) ),
			upper = list( ind = 1L:(m + n + 1), val = c( rep(10000, m), rep(Inf, n), 10000*m ) ) )
	ans = list(objL = objL, Amat = Amat, dir = dir, rhs = rhs, bounds = bounds, 
			problem.max = FALSE, reward = optvars$mu, 
			dataidx = c(n, m), widx = 1:m, midx = length(objL))
	rm(Amat)
	rm(Data)
	gc()
	return( ans )
}

################################################################################
# MAD
################################################################################
lpmad.minrisk = function(optvars)
{
	Data = optvars$Data
	n = NROW(Data)
	m = NCOL(Data)
	colnames(Data) = NULL
	rownames(Data) = NULL
	probability = optvars$probability
	if(is.null(probability)) probability = rep(1/n, n)
	# -[ Nxm ] -[ diag N ] + [1] = m+N+1
	cons = lpmin.setup(optvars, cdimn = 2*n + 1)
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
	Amat = as.simple_triplet_matrix( rbind( Amat, ineqcon, eqcon ) )
	rhs = rep(0, n)
	rhs = c( rhs, as.numeric(ineqB), as.numeric(eqB) )
	dir = c( rep("==", n), rep("<=", ineqn), rep("==", eqn) )
	objL = c( rep(0, m), 0, probability, probability )
	bounds = list( lower = list( ind = 1L:(m + 2 * n + 1), val = c(LB,  1, rep(0, 2 * n) ) ),
			upper = list( ind = 1L:(m + 2 * n + 1), val = c( UB, 1, rep(Inf, 2 * n) ) ) )	
	ans = list(objL = objL, Amat = Amat, dir = dir, rhs = rhs, bounds = bounds, 
			problem.max = FALSE, reward = optvars$mu, 
			dataidx = c(n, m), widx = 1:m)
	rm(Amat)
	rm(Data)
	gc()
	return( ans )
}


lpmad.optrisk = function(optvars)
{
	Data = optvars$Data
	n = NROW(Data)
	m = NCOL(Data)
	colnames(Data) = NULL
	rownames(Data) = NULL
	probability = optvars$probability
	if(is.null(probability)) probability = rep(1/n, n)
	# -[ Nxm ] -[ diag N ] + [1] = m+N+1
	cons = lpopt.setup(optvars, cdimn = 2*n + 1)
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
	Amat = cbind( -Data, diag(n) , -diag(n), matrix(optvars$benchmark, ncol = 1, nrow = n) )
	Amat = rbind( Amat, 
			cbind( matrix(optvars$mu, ncol = m, nrow = 1), 
					matrix(0, ncol = n, nrow = 1),
					matrix(0, ncol = n, nrow = 1),
					matrix(0, ncol = 1, nrow = 1) ) )
	Amat = rbind( Amat, ineqcon, eqcon )
	Amat = as.simple_triplet_matrix(Amat)
	rhs = rep(0, n)
	rhs = c( rhs, 1, as.numeric(ineqB), as.numeric(eqB) )
	dir = c( rep("==", n), "==", rep("<=", ineqn), rep("==", eqn) )
	objL = c( rep(0, m), probability, probability, 0 )
	bounds = list( lower = list( ind = 1L:(m + 2 * n + 1), val = c( rep(-10000, m),  rep(0, 2 * n), 0.05 ) ),
			upper = list( ind = 1L:(m + 2 * n + 1), val = c( rep(10000, m), rep(Inf, 2 * n), 10000*m ) ) )
	ans = list(objL = objL, Amat = Amat, dir = dir, rhs = rhs, bounds = bounds, 
			problem.max = FALSE, reward = optvars$mu, 
			dataidx = c(n, m), widx = 1:m, midx = length(objL))
	rm(Amat)
	rm(Data)
	gc()
	return( ans )
}

################################################################################
# MinMax
################################################################################

lpminmax.minrisk = function(optvars)
{
	Data = optvars$Data
	n = NROW(Data)
	m = NCOL(Data)
	colnames(Data) = NULL
	rownames(Data) = NULL
	# Probability not used in MinMax problem
	# probability = optvars$probability
	# if(is.null(probability)) probability = rep(1/n, n)
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
	Amat = as.simple_triplet_matrix(Amat)
	rhs = rep(0, n)
	rhs = c( rhs, as.numeric(ineqB), as.numeric(eqB) )
	
	dir = c( rep("<=", n), rep("<=", ineqn), rep("==", eqn) )
	
	objL = c( rep(0, m), 0, -1 )
	
	bounds = list( lower = list( ind = 1L:(m + 2), val = c(LB,  1, -Inf ) ),
			upper = list( ind = 1L:(m + 2), val = c( UB, 1, Inf ) ) )
		
	ans = list(objL = objL, Amat = Amat, dir = dir, rhs = rhs, bounds = bounds, 
			problem.max = FALSE, reward = optvars$mu, 
			dataidx = c(n, m), widx = 1:m, vidx = length(objL))	
	rm(Amat)
	rm(Data)
	gc()
	return( ans )
}


lpminmax.optrisk = function(optvars)
{
	Data = optvars$Data
	n = NROW(Data)
	m = NCOL(Data)
	colnames(Data) = NULL
	rownames(Data) = NULL
	# probability = optvars$probability
	# if(is.null(probability)) probability = rep(1/n, n)
	# -[ Nxm ] -[ diag N ] + [1] = m+N+1
	cons = lpopt.setup(optvars, cdimn = 2)
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
	
	Amat = cbind( -Data, matrix(1, ncol = 1, nrow = n), matrix(optvars$benchmark, ncol = 1, nrow = n) )
	Amat = rbind( Amat, 
			cbind( matrix(optvars$mu, ncol = m, nrow = 1), 
					matrix(0, nrow = 1, ncol = 1),
					matrix(0, nrow = 1, ncol = 1) ) )
	Amat = rbind( Amat, ineqcon, eqcon )
	Amat = as.simple_triplet_matrix(Amat)
	rhs = rep(0, n)
	rhs = c( rhs, 1, as.numeric(ineqB), as.numeric(eqB) )
	
	dir = c( rep("<=", n), "==", rep("<=", ineqn), rep("==", eqn) )
	
	objL = c( rep(0, m), -1, 0 )
	
	bounds = list( lower = list( ind = 1L:(m + 2), val = c(rep(-10000, m),  -Inf, 0.05 ) ),
			upper = list( ind = 1L:(m + 2), val = c( rep(10000, m), Inf, 10000*m ) ) )
		
	ans = list(objL = objL, Amat = Amat, dir = dir, rhs = rhs, bounds = bounds, 
			problem.max = FALSE, reward = optvars$mu, 
			dataidx = c(n, m), widx = 1:m, midx = m+2, vidx = m+1)	
	rm(Amat)
	rm(Data)
	gc()
	return( ans )
}

################################################################################
# CVaR
################################################################################
lpcvar.minrisk = function(optvars)
{
	Data = optvars$Data
	n = NROW(Data)
	m = NCOL(Data)
	colnames(Data) = NULL
	rownames(Data) = NULL
	probability = optvars$probability
	if(is.null(probability)) probability = rep(1/n, n)
	# -[ Nxm ] -[ diag N ] + [1] = m+N+1
	cons = lpmin.setup(optvars, cdimn = n+2)
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
	Amat = as.simple_triplet_matrix( rbind( Amat, ineqcon, eqcon ) )
	rhs = rep(0, n)
	rhs = c( rhs, as.numeric(ineqB), as.numeric(eqB) )
	
	dir = c( rep("<=", n), rep("<=", ineqn), rep("==", eqn) )
	
	objL = c( rep(0, m), 1, probability/optvars$alpha, 0)
	
	bounds = list( lower = list( ind = 1L:(m + n + 2), val = c(LB,  -1, rep(0, n), 1 ) ),
			upper = list( ind = 1L:(m + n + 2), val = c( UB, 1, rep(Inf, n), 1 ) ) )
	
	ans = list(objL = objL, Amat = Amat, dir = dir, rhs = rhs, bounds = bounds, 
			problem.max = FALSE, reward = optvars$mu, 
			dataidx = c(n, m), widx = 1:m, vidx = m+1)
	rm(Data)
	rm(Amat)
	gc()
	return( ans )
}

lpcvar.optrisk = function(optvars)
{
	Data = optvars$Data
	n = NROW(Data)
	m = NCOL(Data)
	colnames(Data) = NULL
	rownames(Data) = NULL
	probability = optvars$probability
	if(is.null(probability)) probability = rep(1/n, n)
	# -[ Nxm ] -[ diag N ] + [1] = m+N+1
	cons = lpopt.setup(optvars, cdimn = n+2)
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
	
	Amat = cbind( -Data, matrix(-1, nrow = n, ncol = 1), -diag(n) , matrix(optvars$benchmark, ncol = 1, nrow = n) )
	Amat = rbind( Amat, 
			cbind( matrix(optvars$mu, ncol = m, nrow = 1), 
					matrix(0, ncol = 1, nrow = 1), 
					matrix(0, ncol = n, nrow = 1), 
					matrix(0, ncol = 1, nrow = 1) ) )
	Amat = rbind( Amat, ineqcon, eqcon )
	Amat = as.simple_triplet_matrix(Amat)
	rhs = rep(0, n)
	rhs = c( rhs, 1, as.numeric(ineqB), as.numeric(eqB) )
	
	dir = c( rep("<=", n), "==", rep("<=", ineqn), rep("==", eqn) )
	
	objL = c( rep(0, m), 1, probability/optvars$alpha, 0 )
	
	bounds = list( lower = list( ind = 1L:(m + n + 2), val = c( rep(-10000, m),  0, rep(0, n), 0.05 ) ),
			upper = list( ind = 1L:(m + n + 2), val = c( rep(10000, m), 100000, rep(Inf, n), 10000*m ) ) )
		
	ans = list(objL = objL, Amat = Amat, dir = dir, rhs = rhs, bounds = bounds, 
			problem.max = FALSE, reward = optvars$mu, 
			dataidx = c(n, m), widx = 1:m, vidx = m+1, midx = length(objL))
	rm(Data)
	rm(Amat)
	gc()
	return( ans )
}
################################################################################
# CDaR
################################################################################
lpcdar.minrisk = function(optvars)
{
	Data = optvars$Data
	n = NROW(Data)
	m = NCOL(Data)
	colnames(Data) = NULL
	rownames(Data) = NULL
	probability = optvars$probability
	if(is.null(probability)) probability = rep(1/n, n)
	# -[ Nxm ] -[ diag N ] + [1] = m+N+1
	cons = lpmin.setup(optvars, cdimn = 2*n+2)
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
	
	Amat = rbind( Amat, ineqcon, eqcon)
	Amat = as.simple_triplet_matrix(Amat)
	rhs = rep(0, 2 * n)
	rhs = c( rhs, as.numeric(ineqB), as.numeric(eqB) )
	
	dir = c( rep("<=", 2 * n), rep("<=", ineqn), rep("==", eqn) )
	
	objL = c( rep(0, m), 1, probability/optvars$alpha, rep(0, n),  0)
	
	bounds = list( lower = list( ind = 1L:(m + 2 * n + 2), val = c( LB, -1, rep(0, 2 * n), 1) ),
			upper = list( ind = 1L:(m + 2 * n + 2), val = c( UB, 1, rep(Inf, 2 * n), 1 ) ) )	
	
	ans = list(objL = objL, Amat = Amat, dir = dir, rhs = rhs, bounds = bounds, 
			problem.max = FALSE, reward = optvars$mu, 
			dataidx = c(n, m), widx = 1:m, vidx = m+1)
	rm(Data)
	rm(Amat)
	gc()
	return( ans )	
}

lpcdar.optrisk = function(optvars)
{
	Data = optvars$Data
	n = NROW(Data)
	m = NCOL(Data)
	colnames(Data) = NULL
	rownames(Data) = NULL
	probability = optvars$probability
	if(is.null(probability)) probability = rep(1/n, n)
	# -[ Nxm ] -[ diag N ] + [1] = m+N+1
	cons = lpopt.setup(optvars, cdimn = 2*n+2)
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
			matrix(optvars$benchmark, nrow = n, ncol = 1) )
	
	Amat = rbind( Amat, cbind(
					matrix(0, nrow = n, ncol = m),
					matrix(-1, nrow = n, ncol = 1), 
					-1 * diag(n),
					diag(n),
					matrix(0, nrow = n, ncol = 1) ) )
	
	Amat = rbind( Amat, cbind(
					matrix(optvars$mu, nrow = 1, ncol = m),
					matrix(0, nrow = 1, ncol = 1), 
					matrix(0, nrow = 1, ncol = n),
					matrix(0, nrow = 1, ncol = n),
					matrix(0, nrow = 1, ncol = 1) ) )
	
	
	Amat = rbind( Amat, ineqcon, eqcon)
	Amat = as.simple_triplet_matrix(Amat)
	rhs = c(rep(0, 2 * n), 1)
	rhs = c( rhs, as.numeric(ineqB), as.numeric(eqB) )
	
	dir = c( rep("<=", 2 * n), "==", rep("<=", ineqn), rep("==", eqn) )
	
	objL = c( rep(0, m) , 1, probability/optvars$alpha, rep(0, n), 0)
	
	bounds = list( lower = list( ind = 1L:(m + 2 * n + 2), val = c( rep(-10000, m), -100000, rep(0, 2 * n), 1e-12 ) ),
			upper = list( ind = 1L:(m + 2 * n + 2), val = c( rep( 10000, m), 100000, rep(Inf, 2 * n), 10000*m ) ) )	
	
	ans = list(objL = objL, Amat = Amat, dir = dir, rhs = rhs, bounds = bounds, 
			problem.max = FALSE, reward = optvars$mu, 
			dataidx = c(n, m), widx = 1:m, vidx = m+1, midx = length(objL))
	rm(Data)
	rm(Amat)
	gc()
	return( ans )
}
################################################################################

#risk: 1 (MAD), 2 (MiniMax), 3 (CVaR), 4 (CDaR), 5(EV), 6 (LPM), 7 (LPMUPM)
lpport = function(optvars, solver, ...)
{
	valid.solvers = c("glpk", "symphony")
	solver = match.arg(tolower(solver), valid.solvers)
	#if(solver=="symphony"){
	#	if(!requireNamespace(Rsymphony, quietly=TRUE)) stop("\nPackage 'Rsymphony' must be installed") 
	#}
	risk = c("mad", "minimax", "cvar", "cdar", "ev", "lpm", "lpmupm")[optvars$index[4]]
	if(optvars$index[5] == 1){
		setup = switch(risk,
				"mad"     = lpmad.minrisk(optvars),
				"minimax" = lpminmax.minrisk(optvars),
				"cvar"    = lpcvar.minrisk(optvars),
				"cdar"    = lpcdar.minrisk(optvars),
				"lpm"     = lplpm1.minrisk(optvars))
		sol = lpminsolver(setup, solver, optvars, ...)
	} else{
		setup = switch(risk,
				'mad'  = lpmad.optrisk(optvars),
				'minimax' = lpminmax.optrisk(optvars),
				'cvar' = lpcvar.optrisk(optvars),
				'cdar' = lpcdar.optrisk(optvars),
				'lpm' = lplpm1.optrisk(optvars))
		sol = lpoptsolver(setup, solver, optvars, ...)
	}
	#sol$risk = fun.risk(sol$weights, scenario, options, risk, benchmark)	
	rm(setup)
	gc()
	return( sol )
}

lpminsolver = function(setup, solver, optvars, ...)
{
	if( is.null(list(...)$verbose) ) verbose = FALSE else verbose = as.logical( list(...)$verbose )
	m = setup$dataidx[2]
	sol = switch(solver,
			glpk = try( Rglpk_solve_LP(obj = setup$objL, mat = setup$Amat, dir = setup$dir, 
							rhs = setup$rhs, bounds = setup$bounds, max = setup$problem.max, 
							verbose = verbose), silent = FALSE ),
			symphony =  try( Rsymphony::Rsymphony_solve_LP(obj = setup$objL, 
							mat = setup$Amat, dir = setup$dir, rhs = setup$rhs, 
							bounds = setup$bounds, max = setup$problem.max), 
					silent = FALSE ))
	if( inherits(sol, "try-error") ){
		status = "non-convergence"
		weights = rep(NA, m)
		risk = NA
		reward = NA
		VaR = NA
		DaR = NA
	} else{
		status  = sol$status
		weights = sol$solution[setup$widx]
		# GLPK has a documentation bug, statting that it returns objval when it
		# in fact returns "optimum".
		risk    = switch(solver, glpk = sol$optimum, symphony = sol$objval)
		if( !is.null(setup$reward) ){
			reward = sum(optvars$mu * weights)
		} else {
			reward = NULL
		}
		if(optvars$index[4]==6) risk = risk - 1
		if(optvars$index[4]==3) VaR = sol$solution[setup$vidx] else VaR = NA
		if(optvars$index[4]==4) DaR = sol$solution[setup$vidx] else DaR = NA
	}
	rm(setup)
	gc()
	return( list(status = status, weights = weights, risk = risk, reward = reward, 
					multiplier = 1, VaR = VaR, DaR = DaR) )
}

lpoptsolver = function(setup, solver, optvars, ...)
{
	if( is.null(list(...)$verbose) ) verbose = FALSE else verbose = as.logical( list(...)$verbose )
	
	m = setup$dataidx[2]
	
	sol = switch(solver,
			glpk = try( Rglpk_solve_LP(obj = setup$objL, mat = setup$Amat, dir = setup$dir, 
							rhs = setup$rhs, bounds = setup$bounds, max = setup$problem.max, 
							verbose = verbose), silent = FALSE ),
			symphony =  try( Rsymphony::Rsymphony_solve_LP(obj = setup$objL, 
							mat = setup$Amat, dir = setup$dir, rhs = setup$rhs, 
							bounds = setup$bounds, max = setup$problem.max), 
					silent = FALSE ))
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
		optimum  = switch(solver, glpk = sol$optimum, symphony = sol$objval)
		if(optvars$index[4]==3) VaR = sol$solution[setup$vidx]/multiplier else VaR = NA
		if(optvars$index[4]==4) DaR = sol$solution[setup$vidx]/multiplier else DaR = NA
	}
	rm(setup)
	
	return( list(status = status, weights = weights, risk = risk, reward = reward, 
					optimum = optimum, multiplier = multiplier, VaR = VaR, DaR = DaR) )
}