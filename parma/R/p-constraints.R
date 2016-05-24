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

# constraints for the different formulations
################################################################
# LP: Fractional Formulation Constraints
# Based on Charnes/Cooper Method
# (Optimal Risk/Reward Portfolios)
################################################################
parma.ineq.lpopt = function(ineq, ineqLB, ineqUB, LB, UB, cdim)
{
	m = length(LB)
	if( is.matrix(ineq) ) {
		ineqn = dim(ineq)[1]
		ineqcon = cbind(
				rbind(ineq, -ineq),
				matrix(0, nrow = 2 * ineqn, ncol = cdim) )
		ineqcon = rbind( ineqcon, 
				cbind(diag(m), matrix(0, nrow = m, ncol = cdim) ),
				cbind(-1 * diag(m), matrix(0, nrow = m, ncol = cdim) ) )
		ineqcon[1:ineqn, m + cdim] = -1 * ineqUB
		ineqcon[(ineqn + 1):(2 * ineqn), m + cdim] = 1 * ineqLB
		ineqcon[(2 * ineqn + 1):(2 * ineqn + m), m + cdim] =  -1 * UB
		ineqcon[(2 * ineqn + m + 1):(2 * ineqn + 2 * m), m + cdim] = 1 * LB
		ineqB = c( rep(0, 2 * ineqn), rep(0, 2 * m) )
		nn = 2 * ineqn + 2 * m
	} else{
		ineqn = 0
		ineqcon = rbind(
				cbind(diag(m), matrix(0, nrow = m, ncol = cdim) ),
				cbind(-1 * diag(m), matrix(0, nrow = m, ncol = cdim) ) )		
		ineqcon[1:m, m + cdim] =  -1 * UB
		ineqcon[(m + 1):(2 * m), m + cdim] = 1 * LB
		ineqB = c( rep(0, 2 * m) )
		nn = 2 * m
	}
	return( list( ineqcon = ineqcon, ineqB = ineqB, ineqn = nn, ineqm = m ) )
}

# rewardB here is NULL/we deal with fractional objectives in the main problem formulation
parma.eq.lpopt = function(eq, eqB, reward, rewardB, cdim, budget = 1, m = NULL)
{
	eqcon = NULL
	eqn = 0
	if( !is.null(m) ) eqm = m else eqm = length( reward )
	eqB = NULL
	if( is.matrix(eq) ){
		eqn = dim(eq)[1]
		eqm = dim(eq)[2]
		if( !is.null(budget) ){
			# if sum of weights is missing we ipomse summation constraint
			eqcon = rbind( eqcon, 
					cbind( matrix(1, nrow =  1, ncol = eqm), matrix(0, nrow =  1, ncol = cdim)),
					cbind( eq, matrix(0, nrow =  eqn, ncol = cdim)))
			eqcon[ , eqm + cdim] = c(-budget, -eqB)
			eqn = eqn + 1
			eqm = eqm
			eqB = c( 0, rep(0, eqn) )
		} else{
			eqcon = rbind( eqcon, 
					cbind( eq, matrix(0, nrow = eqn, ncol = cdim) ) )
			eqcon[ , eqm + cdim] = -eqB
			eqn = eqn
			eqB = rep(0, eqn)
		}
	} else{
		if( !is.null(budget) ){
			eqcon = rbind( eqcon, 
					cbind( matrix(1, nrow =  1, ncol = eqm), matrix(0, nrow =  1, ncol = cdim)))
			eqcon[ , eqm + cdim] = -budget
			eqn = 1
			eqm = eqm
			eqB = 0
		}
	}
	return( list( eqcon = eqcon, eqB = eqB, eqn = eqn, eqm = eqm ) )
}

################################################################
# LP: Standard Problem Formulation Constraints
# (Standard Constrained Portfolios)
################################################################
parma.ineq.lpmin = function(ineq, ineqLB, ineqUB, cdim)
{
	if( is.matrix(ineq) ) {
		ineqn = dim(ineq)[1]
		ineqm = dim(ineq)[2]
		ineqcon = cbind(
				rbind(ineq, -ineq), 
				matrix(0, nrow = 2 * ineqn, ncol = cdim) )
		ineqB = c( as.numeric(ineqUB), as.numeric(-1 * ineqLB) )
		nn = 2 * ineqn
	} else{
		ineqcon = NULL
		nn = 0
		ineqm = 0
		ineqB = NULL		
	}
	return( list( ineqcon = ineqcon, ineqB = ineqB, ineqn = nn, ineqm = ineqm ) )
}

parma.eq.lpmin = function(eq, eqB, reward, rewardB, cdim, budget = 1, m = NULL)
{
	mbench = 0
	if( !is.null(rewardB) ){
		eqcon = cbind( matrix(reward, nrow = 1), 
				matrix(mbench, nrow = 1, ncol = 1), 
				matrix(0, nrow = 1, ncol = cdim - 1) )
		eqn = 1
		if( !is.null(m) ) eqm = m else eqm = length( reward )
		eqB = rewardB
		if(is.matrix(eq)){
			eqn = dim(eq)[1]
			eqm = dim(eq)[2]
			if( budget ){
				# if sum of weights is missing we impose summation constraint
				eqcon = rbind(  cbind( matrix(1, nrow =  1, ncol = eqm), matrix(0, nrow =  1, ncol = cdim) ),
						cbind( eq, matrix(0, nrow =  eqn, ncol = cdim) ), 
						eqcon )
				eqn = eqn + 2
				eqm = eqm
				eqB = c( 1, eqB, rewardB )
			} else{
				eqcon = rbind( cbind( eq, matrix(0, nrow = eqn, ncol = cdim) ), eqcon )
				eqn = eqn + 1
				eqB = c( eqB, rewardB )
			}
		} else{
			if( !is.null(budget) ){
				eqcon = rbind( cbind( matrix(1, nrow =  1, ncol = eqm), matrix(0, nrow =  1, ncol = cdim) ),
						eqcon )
				eqn = eqn + 1
				eqm = eqm
				eqB = c( budget, rewardB )
			}
		}
	} else{
		eqcon = NULL
		eqn = 0
		if( !is.null(m) ) eqm = m else eqm = length( reward )
		if( is.matrix(eq) ){
			eqn = dim(eq)[1]
			eqm = dim(eq)[2]
			if( !is.null(budget) ){
				# if sum of weights is missing we impose summation constraint
				eqcon = rbind(  cbind( matrix(1, nrow =  1, ncol = eqm), matrix(0, nrow =  1, ncol = cdim) ),
						cbind( eq, matrix(0, nrow =  eqn, ncol = cdim) ) )
				eqn = eqn + 1
				eqm = eqm
				eqB = c(budget, eqB)
			} else{
				eqcon = rbind( cbind( eq, matrix(0, nrow = eqn, ncol = cdim) ),  eqcon )
				eqn = eqn
				eqB = eqB
			}
		} else{
			if( !is.null(budget) ){
				eqcon = rbind( cbind( matrix(1, nrow =  1, ncol = eqm), matrix(0, nrow =  1, ncol = cdim) ), eqcon )
				eqn = 1
				eqm = eqm
				eqB = budget
			}
		}
	}
	
	return( list( eqcon = eqcon, eqB = eqB, eqn = eqn, eqm = eqm ) )
}

################################################################
# QP constraints
################################################################
parma.con.qpmin1 = function(eq, eqB, reward, rewardB, ineq, ineqLB, ineqUB, LB, UB, budget, m)
{
	Amat = NULL
	bvec = NULL
	meq = 0
	dvec = rep(0, m)
	# equalities
	if( !is.null(rewardB) )
	{
		Amat = rbind( Amat, matrix( reward, ncol = m, nrow = 1) )
		meq = meq + 1
		bvec = c(bvec, rewardB)
	}
	
	if( !is.null(budget) )
	{
		Amat = rbind( Amat, matrix( c(0, rep(1,m-1)), ncol = m, nrow = 1))
		meq = meq + 1
		bvec = c(bvec, budget)
	}
	
	if( !is.null(eq) ){
		neq = length(eqB)
		Amat = rbind( Amat, eq)
		meq = meq + neq
		bvec = c(bvec, eqB)
	}
	
	# inequalities ( >= )
	if( !is.null(ineq) ){
		nineq = length(ineqLB)
		Amat = rbind(Amat, 
				ineq, -ineq)
		bvec = c(bvec, ineqLB, -ineqUB)
	}
	# Upper and Lower Bounds
	Amat = rbind(Amat, diag(m), -diag(m))
	bvec = c(bvec, if(!is.null(LB)) LB else rep(-1000, m), if(!is.null(UB)) -UB else rep(-1000, m))
	
	return( list(Amat = Amat, bvec = bvec, meq = meq) )
}

parma.con.qpmin2 = function(eq, eqB, reward, rewardB, ineq, ineqLB, ineqUB, LB, UB, budget, m)
{
	Amat = NULL
	bvec = NULL
	meq = 0
	dvec = rep(0, m)
	# equalities
	if( !is.null(rewardB) )
	{
		Amat = rbind( Amat, matrix( reward, ncol = m, nrow = 1) )
		meq = meq + 1
		bvec = c(bvec, rewardB)
	}
	
	if( !is.null(budget) )
	{
		Amat = rbind( Amat, matrix( 1, ncol = m, nrow = 1))
		meq = meq + 1
		bvec = c(bvec, budget)
	}
	
	if( !is.null(eq) ){
		neq = length(eqB)
		Amat = rbind( Amat, eq)
		meq = meq + neq
		bvec = c(bvec, eqB)
	}
	
	# inequalities ( >= )
	if( !is.null(ineq) ){
		nineq = length(ineqLB)
		Amat = rbind(Amat, 
				ineq, -ineq)
		bvec = c(bvec, ineqLB, -ineqUB)
	}
	# Upper and Lower Bounds
	Amat = rbind(Amat, diag(m), -diag(m))
	bvec = c(bvec, if(!is.null(LB)) LB else rep(-1000, m), if(!is.null(UB)) -UB else rep(-1000, m))
	
	return( list(Amat = Amat, bvec = bvec, meq = meq) )
}

parma.con.qpopt = function(eq, eqB, reward, ineq, ineqLB, ineqUB, LB, UB, budget, m)
{
	Amat = NULL
	bvec = NULL
	meq = 0
	dvec = rep(0, m)
	# equalities
	mbench = 0
	Amat = rbind( Amat, matrix( c(mbench, reward), ncol = m, nrow = 1) )
	meq = meq + 1
	bvec = c(bvec, 1)
	Amat = rbind( Amat, matrix( c(budget, rep(1, m-1)), ncol = m, nrow = 1))
	meq = meq + 1
	bvec = c(bvec, 0)
	if( !is.null(eq) ){
		neq = length(eqB)
		Amat = rbind( Amat, cbind(eqB, eq) )
		meq = meq + neq
		bvec = c(bvec, rep(0, neq))
	}
	# inequalities ( >= )
	Amat = rbind(Amat, matrix(c(-1, rep(0, m-1)), ncol = m, nrow = 1))
	Amat = rbind(Amat, cbind( LB, diag(m-1)))
	Amat = rbind(Amat, cbind(-UB,-diag(m-1)))
	
	bvec = c(bvec, rep(0, 2*(m-1) + 1))
	nineq = 2*(m)
	if( !is.null(ineq) ){
		nineq = nineq+length(ineqLB)
		Amat = rbind(Amat, 
				cbind(ineqLB, ineq), cbind(-ineqUB,-ineq))
		bvec = c(bvec, rep(0, 2*length(ineqLB)))
	}
	# Upper and Lower Bounds	
	return( list(Amat = Amat, bvec = bvec, meq = meq) )
}
###############################################################################
# General NLP Constraint Functions
###############################################################################
parma.ineq.minfun = function(w, optvars, uservars){
	cons = NULL
	if(!is.null(optvars$ineqfun)){
		n = length(optvars$ineqfun)
		if(optvars$index[3]==1) cons = ineqfun.target.min(w, optvars, uservars) else cons = NULL
		fnlist = list(w, optvars, uservars)
		names(fnlist) = c("w", "optvars", "uservars")
		for(i in 1:n){
			cons = c(cons, do.call(optvars$ineqfun[[i]], args = fnlist))
		}
	} else{
		if(optvars$index[3]==1) cons = ineqfun.target.min(w, optvars, uservars) else cons = NULL
	}
	if(optvars$index[4]==2) cons = c(cons, ineqfun.minminmax(w, optvars, uservars))	
	return(cons)
}

parma.ineq.mingrad = function(w, optvars, uservars){
	cons = NULL
	if(!is.null(optvars$ineqgrad)){
		n = length(optvars$ineqgrad)
		if(optvars$index[3]==1) cons = ineqjac.target.min(w, optvars, uservars) else cons = NULL
		fnlist = list(w, optvars, uservars)
		names(fnlist) = c("w", "optvars", "uservars")
		for(i in 1:n){
			cons = rbind(cons, do.call(optvars$ineqgrad[[i]], args = fnlist))
		}
	} else{
		if(optvars$index[3]==1) cons = ineqjac.target.min(w, optvars, uservars) else cons = NULL
	}
	if(optvars$index[4]==2){
		if(is.null(cons)) cons = ineqjac.minminmax(w, optvars, uservars) else cons = rbind(cons, ineqjac.minminmax(w, optvars, uservars))
	}
	return(cons)
}

parma.eq.minfun = function(w, optvars, uservars){
	if(!is.null(optvars$eqfun)){
		n = length(optvars$eqfun)
		if(optvars$index[6]==0) cons = eqfun.budget.min(w, optvars, uservars) else cons = eqfun.leverage.min(w, optvars, uservars)
		if(optvars$index[3]==2) cons = c(cons, eqfun.target.min(w, optvars, uservars))
		fnlist = list(w, optvars, uservars)
		names(fnlist) = c("w", "optvars", "uservars")
		for(i in 1:n){
			cons = c(cons, do.call(optvars$eqfun[[i]], args = fnlist))
		}
	} else{
		if(optvars$index[6]==0) cons = eqfun.budget.min(w, optvars, uservars) else cons = eqfun.leverage.min(w, optvars, uservars)
		if(optvars$index[3]==2) cons = c(cons, eqfun.target.min(w, optvars, uservars))
	}
	return(cons)
}

parma.eq.mingrad = function(w, optvars, uservars){
	if(!is.null(optvars$eqgrad)){
		n = length(optvars$eqgrad)
		if(optvars$index[6]==0) cons = eqjac.budget.min(w, optvars, uservars) else cons = eqjac.leverage.min(w, optvars, uservars)
		if(optvars$index[3]==2) cons = rbind(cons, eqjac.target.min(w, optvars, uservars))
		fnlist = list(w, optvars, uservars)
		names(fnlist) = c("w", "optvars", "uservars")
		for(i in 1:n){
			cons = rbind(cons, do.call(optvars$eqgrad[[i]], args = fnlist))
		}
	} else{
		if(optvars$index[6]==0) cons = eqjac.budget.min(w, optvars, uservars) else cons = eqjac.leverage.min(w, optvars, uservars)
		if(optvars$index[3]==2) cons = rbind(cons, eqjac.target.min(w, optvars, uservars))
	}
	return(cons)
}


# Special consideration for MiniMax optimization
parma.ineq.optfun = function(w, optvars, uservars){
	cons = NULL
	if(!is.null(optvars$ineqfun)){
		n = length(optvars$ineqfun)
		cons = c(cons, ineqfun.bounds.opt(w, optvars, uservars))
		fnlist = list(w, optvars, uservars)
		names(fnlist) = c("w", "optvars", "uservars")
		for(i in 1:n){
			cons = c(cons, do.call(optvars$ineqfun[[i]], args = fnlist))
		}
	} else{
		cons = ineqfun.bounds.opt(w, optvars, uservars)
	}
	if(optvars$index[4]==2) cons = c(cons, ineqfun.optminmax(w, optvars, uservars))
	
	return(cons)
}

parma.ineq.optgrad = function(w, optvars, uservars){
	cons = NULL
	if(!is.null(optvars$ineqgrad)){
		n = length(optvars$ineqgrad)
		cons = ineqjac.bounds.opt(w, optvars, uservars)
		fnlist = list(w, optvars, uservars)
		names(fnlist) = c("w", "optvars", "uservars")
		for(i in 1:n){
			cons = rbind(cons, do.call(optvars$ineqgrad[[i]], args = fnlist))
		}
	} else{
		cons = ineqjac.bounds.opt(w, optvars, uservars)
	}
	if(optvars$index[4]==2) cons = rbind(cons, ineqjac.optminmax(w, optvars, uservars))
	return(cons)
}


parma.eq.optfun = function(w, optvars, uservars){
	# equality target is a given in fractional programming
	cons = eqfun.target.opt(w, optvars, uservars)
	if(!is.null(optvars$eqfun)){
		n = length(optvars$eqfun)
		if(optvars$index[6]==0) cons = c(cons, eqfun.budget.opt(w, optvars, uservars)) else cons = c(cons, eqfun.leverage.opt(w, optvars, uservars))
		fnlist = list(w, optvars, uservars)
		names(fnlist) = c("w", "optvars", "uservars")
		for(i in 1:n){
			cons = c(cons, do.call(optvars$eqfun[[i]], args = fnlist))
		}
	} else{
		if(optvars$index[6]==0) cons = c(cons, eqfun.budget.opt(w, optvars, uservars)) else cons = c(cons, eqfun.leverage.opt(w, optvars, uservars))
	}
	return(cons)
}

parma.eq.optgrad = function(w, optvars, uservars){
	cons = eqjac.target.opt(w, optvars, uservars)
	
	if(!is.null(optvars$eqgrad)){
		n = length(optvars$eqgrad)
		if(optvars$index[6]==0) cons = rbind(cons, eqjac.budget.opt(w, optvars, uservars)) else cons = rbind(cons, eqjac.leverage.opt(w, optvars, uservars))
		fnlist = list(w, optvars, uservars)
		names(fnlist) = c("w", "optvars", "uservars")
		for(i in 1:n){
			cons = rbind(cons, do.call(optvars$eqgrad[[i]], args = fnlist))
		}
	} else{
		if(optvars$index[6]==0) cons = rbind(cons, eqjac.budget.opt(w, optvars, uservars)) else cons = rbind(cons, eqjac.leverage.opt(w, optvars, uservars))
	}
	return(cons)
}

parma.ineq.lpmupm = function(w, optvars, uservars){
	if(!is.null(optvars$ineqfun)){
		n = length(optvars$ineqfun)
		cons = func.ineq.lpmupm(w, optvars, uservars)
		fnlist = list(w, optvars, uservars)
		names(fnlist) = c("w", "optvars", "uservars")
		for(i in 1:n){
			cons = c(cons, do.call(optvars$ineqfun[[i]], args = fnlist))
		}
	} else{
		cons = func.ineq.lpmupm(w, optvars, uservars)
	}
	return(cons)
}