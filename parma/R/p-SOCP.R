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
# SOCP formulation:
# || A{i}*x + b{i} || <= C{i}*x + d{i}
# p-norm constraints (leverage)

# Types of problems:
# 1. Benchmark Relative
# 2. Leverage
# 3. QCQP
# 4. Minrisk (equality and inequality reward)
# 5. MaxReward (inequality risk)
# 6. Fractional

# Minimize Risk given reward equality

parma.eqcon.qpsocp = function(S, reward, rewardB, benchmarkS = NULL, eq = NULL, eqB = NULL, ineq = NULL, 
		ineqLB = NULL, ineqUB = NULL, budget = NULL, leverage = NULL, LB, UB, eqSlack = 1e-5)
{
	m = length(LB)	
	
	if(!is.null(leverage)){
		psize = m
		plev = rep(0, m)
		pLB =  rep(0, m)
		pUB = abs(UB)
		uselev = TRUE
	} else{
		psize = 0
		plev = NULL
		pUB = NULL
		pLB = NULL
		uselev = FALSE
	}
	if(!is.null(benchmarkS))
	{
		bsize = 1
		bvalue = 0
		useBench = TRUE
		bLB =   1.00000
		bUB =  -0.99999
		widx = 2:(m+1)
	} else{
		widx = 1:m
		bsize = 0
		bvalue = NULL
		useBench = FALSE
		bLB = NULL
		bUB = NULL
	}
	# setup:
	C = NULL
	d = NULL
	A = NULL
	b = NULL
	N = NULL
	# Objective Function
	# min t
	#  f = [0,...,0, 1]	
	f = c(rep(0, m+bsize), plev, 1)
	# x   : [w_1, ..., w_m, t]
	# w_i : weight of asset i
	# t   : auxiliiary variable used for risk minimization
	
	# Risk (LHS): ||A{1}*x|| <= t
	# A{1} : [ Sigma^(1/2) 0]
	## [ Sqrt.Sigma_i_i, ..., Sqrt.Sigma_i_m,  0] + [0]
	## [ .............., ..., ..............., 0] + [0]
	#  [ Sqrt.Sigma_m_i, ..., Sqrt.Sigma_m_m,  0] + [0]
	# Constraint (RHS): C{1}*x+d
	# C{1}: [0,...,0, 1]
	# C{1}*x : t
	if(useBench){
		sqrtS = cbind(matrix(benchmarkS, ncol=1), rbind(matrix(0, nrow=1, ncol = m), S))
		sqrtS[1,2:(m+1)] = benchmarkS[2:(m+1)]
		if(any(eigen(sqrtS, TRUE, only.values=TRUE)$values<1e-10)){
			warning("\nparma: benchmark relative covariance matrix not PD...adjusting to make PD")
			sqrtS <- make.positive.definite(sqrtS, 1e-10)
		}
		sqrtS = .sqrtm(sqrtS)
	} else{
		if(any(eigen(S, TRUE, only.values=TRUE)$values<1e-10)){
			warning("\nparma: benchmark relative covariance matrix not PD...adjusting to make PD")
			S <- make.positive.definite(S, 1e-10)
		}
		sqrtS = .sqrtm(S)
	}
	A = rbind( A, cbind(sqrtS, matrix(0, ncol = psize+1, nrow=m+bsize)) )
	b = c( b, rep(0, m+bsize) )
	C = rbind(C, matrix(c(rep(0, m+bsize), plev, 1), ncol = m+bsize+psize+1) )
	d = c( d, 0 )
	N = c( N, m+bsize)
	
	# Eq Constraints
	if( !is.null( eq ) ){
		eqn = dim( eq )[1]
		C = rbind( C, 
				cbind(if(useBench) matrix(1, ncol = bsize, nrow = eqn) else NULL, eq, 
						matrix(0, ncol = psize+1, nrow = eqn)), 
				cbind(if(useBench) matrix(1, ncol = bsize, nrow = eqn) else NULL, -eq, 
						matrix(0, ncol = psize+1, nrow = eqn)) )
		# might need to add some slack
		d = c( d, -eqB*(1-eqSlack), eqB*(1+eqSlack))
		b = c( b, rep( 0, 2 * eqn ) )
		A = rbind( A, matrix( 0, nrow = 2 * eqn, ncol = m+bsize+psize+1 ) )
		N = c( N, rep(1, 2 * eqn) )
	}
	# Reward Equality Constraint
	if(!is.null(reward)){
		C = rbind( C, 
				matrix( c(bvalue,  reward, plev, 0), ncol = m+bsize+psize+1 ), 
				matrix( c(bvalue, -reward, plev, 0), ncol = m+bsize+psize+1 ) )
		d = c( d, - rewardB*(1-eqSlack) , rewardB*(1+eqSlack))
		b = c(b, 0, 0 )
		A = rbind( A, matrix( 0, nrow = 2, ncol = m+bsize+psize+1) )
		N = c( N, 1, 1 )
	}
	# Budget Constraint
	if( !is.null(budget) ){
		C = rbind( C, 
				matrix( c(bvalue,rep( 1,m),plev,0), ncol = m+bsize+psize+1 ), 
				matrix( c(bvalue,rep(-1,m),plev,0), ncol = m+bsize+psize+1 ) )
		d = c( d, - budget*(1-eqSlack) , budget*(1+eqSlack))
		b = c(b, 0, 0 )
		A = rbind( A, matrix( 0, nrow = 2, ncol = m+bsize+psize+1) )
		N = c( N, 1, 1 )
	}
	if( !is.null(leverage) ){
		A = rbind(A, matrix(0, ncol = m+bsize+psize+1, nrow = 2*m))
		b = c(b, rep(0, 2*m))
		# x+s>=0
		C = rbind(C, cbind( if(useBench)  matrix(0, ncol=1, nrow = m) else NULL, 
						diag(m), diag(m), matrix(0, ncol=1, nrow = m)))
		# -x+s>=0
		C = rbind(C, cbind(if(useBench)   matrix(0, ncol=1, nrow = m) else NULL, 
					   -diag(m), diag(m), matrix(0, ncol=1, nrow = m)))
		# -s-x>=0
		d = c(d, rep(0, 2*m))
		A = rbind(A, matrix(0, ncol = m+bsize+psize+1, nrow = 2))
		b = c(b, rep(0, 2))
		C = rbind(C, cbind(matrix(0, ncol=m+bsize, nrow=1), matrix(-1, ncol=m, nrow=1), matrix(0, ncol=1, nrow=1)))
		C = rbind(C, cbind(matrix(0, ncol=m+bsize, nrow=1), matrix( 1, ncol=m, nrow=1), matrix(0, ncol=1, nrow=1)))
		d = c(d, (1+eqSlack)*leverage, -(1-eqSlack)*leverage)
		N = c( N, rep(1, 2*m+2) )
	}
	
	# Ineq Constraints
	if( !is.null(ineq) ){
		ineqn = dim(ineq)[1]
		C = rbind( C, 
				cbind(if(useBench)  matrix(0, ncol=1, nrow = ineqn) else NULL, ineq, matrix(0, ncol = psize+1, nrow = ineqn)), 
				cbind(if(useBench)  matrix(0, ncol=1, nrow = ineqn) else NULL,-ineq, matrix(0, ncol = psize+1, nrow = ineqn)) )		
		d = c(d, -ineqLB, ineqUB)
		b = c(b, rep( 0, 2 * ineqn ) )
		A = rbind( A, matrix(0, nrow = 2 * ineqn, ncol = bsize+m+psize+1) )
		N = c(N, rep(1, 2 * ineqn) )
	}
	# Parameter Bounds
	C = rbind( C, diag(m+bsize+psize+1), -diag(m+bsize+psize+1) )
	d = c(d, c(bLB, -LB, pLB, 0), c(bUB, UB, pUB, 1))
	b = c(b, rep( 0, 2 * (m+bsize+psize+1) ) )
	A = rbind( A, matrix(0, nrow = 2 * (m+bsize+psize+1), ncol = (m+bsize+psize+1)) )
	N = c(N, rep(1, 2 * (m+bsize+psize+1)) )
	return( list( f = f, A = A, C = C, d = d, b = b, A = A, N = N, widx = widx, m = m) )
}

# Minimize Risk given reward equality
parma.ineqcon.qpsocp = function(S, reward, rewardB, benchmarkS = NULL, eq = NULL, eqB = NULL, ineq = NULL, 
		ineqLB = NULL, ineqUB = NULL, budget = NULL, leverage = NULL, LB, UB, eqSlack = 1e-5)
{
	m = length(LB)	
	
	if(!is.null(leverage)){
		psize = m
		plev = rep(0, m)
		pLB =  rep(0, m)
		pUB = abs(UB)
		uselev = TRUE
	} else{
		psize = 0
		plev = NULL
		pUB = NULL
		pLB = NULL
		uselev = FALSE
	}
	if(!is.null(benchmarkS))
	{
		bsize = 1
		bvalue = 0
		useBench = TRUE
		bLB =   1.00000
		bUB =  -0.99999
		widx = 2:(m+1)
	} else{
		widx = 1:m
		bsize = 0
		bvalue = NULL
		useBench = FALSE
		bLB = NULL
		bUB = NULL
	}
	# setup:
	C = NULL
	d = NULL
	A = NULL
	b = NULL
	N = NULL
	# Objective Function
	# min t
	#  f = [0,...,0, 1]	
	f = c(rep(0, m+bsize), plev, 1)
	# x   : [w_1, ..., w_m, t]
	# w_i : weight of asset i
	# t   : auxiliiary variable used for risk minimization
	
	# Risk (LHS): ||A{1}*x|| <= t
	# A{1} : [ Sigma^(1/2) 0]
	## [ Sqrt.Sigma_i_i, ..., Sqrt.Sigma_i_m,  0] + [0]
	## [ .............., ..., ..............., 0] + [0]
	#  [ Sqrt.Sigma_m_i, ..., Sqrt.Sigma_m_m,  0] + [0]
	# Constraint (RHS): C{1}*x+d
	# C{1}: [0,...,0, 1]
	# C{1}*x : t
	if(useBench){
		sqrtS = cbind(matrix(benchmarkS, ncol=1), rbind(matrix(0, nrow=1, ncol = m), S))
		sqrtS[1,2:(m+1)] = benchmarkS[2:(m+1)]
		if(any(eigen(sqrtS, TRUE, only.values=TRUE)$values<1e-10)){
			warning("\nparma: benchmark relative covariance matrix not PD...adjusting to make PD")
			sqrtS <- make.positive.definite(sqrtS, 1e-10)
		}
		sqrtS = .sqrtm(sqrtS)
	} else{
		if(any(eigen(S, TRUE, only.values=TRUE)$values<1e-10)){
			warning("\nparma: benchmark relative covariance matrix not PD...adjusting to make PD")
			S <- make.positive.definite(S, 1e-10)
		}
		sqrtS = .sqrtm(S)
	}
	A = rbind( A, cbind(sqrtS, matrix(0, ncol = psize+1, nrow=m+bsize)) )
	b = c( b, rep(0, m+bsize) )
	C = rbind(C, matrix(c(rep(0, m+bsize), plev, 1), ncol = m+bsize+psize+1) )
	d = c( d, 0 )
	N = c( N, m+bsize)
	
	# Eq Constraints
	if( !is.null( eq ) ){
		eqn = dim( eq )[1]
		C = rbind( C, 
				cbind(if(useBench) matrix(1, ncol = bsize, nrow = eqn) else NULL, eq, 
						matrix(0, ncol = psize+1, nrow = eqn)), 
				cbind(if(useBench) matrix(1, ncol = bsize, nrow = eqn) else NULL, -eq, 
						matrix(0, ncol = psize+1, nrow = eqn)) )
		# might need to add some slack
		d = c( d, -eqB*(1-eqSlack), eqB*(1+eqSlack))
		b = c( b, rep( 0, 2 * eqn ) )
		A = rbind( A, matrix( 0, nrow = 2 * eqn, ncol = m+bsize+psize+1 ) )
		N = c( N, rep(1, 2 * eqn) )
	}
	# Reward Inequality Constraint
	if(!is.null(reward)){
		C = rbind( C, matrix( c(bvalue,  reward, plev, 0), ncol = m+bsize+psize+1 ) )		
		d = c( d, - rewardB)
		b = c(b, 0)
		A = rbind( A, matrix( 0, nrow = 1, ncol = m+bsize+psize+1) )
		N = c( N, 1)
	}
	# Budget Constraint
	if( !is.null(budget) ){
		C = rbind( C, 
				matrix( c(bvalue,rep( 1,m),plev,0), ncol = m+bsize+psize+1 ), 
				matrix( c(bvalue,rep(-1,m),plev,0), ncol = m+bsize+psize+1 ) )
		d = c( d, - budget*(1-eqSlack) , budget*(1+eqSlack))
		b = c(b, 0, 0 )
		A = rbind( A, matrix( 0, nrow = 2, ncol = m+bsize+psize+1) )
		N = c( N, 1, 1 )
	}
	if( !is.null(leverage) ){
		A = rbind(A, matrix(0, ncol = m+bsize+psize+1, nrow = 2*m))
		b = c(b, rep(0, 2*m))
		# x+s>=0
		C = rbind(C, cbind( if(useBench)  matrix(0, ncol=1, nrow = m) else NULL, 
						diag(m), diag(m), matrix(0, ncol=1, nrow = m)))
		# -x+s>=0
		C = rbind(C, cbind(if(useBench)   matrix(0, ncol=1, nrow = m) else NULL, 
						-diag(m), diag(m), matrix(0, ncol=1, nrow = m)))
		# -s-x>=0
		d = c(d, rep(0, 2*m))
		A = rbind(A, matrix(0, ncol = m+bsize+psize+1, nrow = 2))
		b = c(b, rep(0, 2))
		C = rbind(C, cbind(matrix(0, ncol=m+bsize, nrow=1), matrix(-1, ncol=m, nrow=1), matrix(0, ncol=1, nrow=1)))
		C = rbind(C, cbind(matrix(0, ncol=m+bsize, nrow=1), matrix( 1, ncol=m, nrow=1), matrix(0, ncol=1, nrow=1)))
		d = c(d, (1+eqSlack)*leverage, -(1-eqSlack)*leverage)
		N = c( N, rep(1, 2*m+2) )
	}
	
	# Ineq Constraints
	if( !is.null(ineq) ){
		ineqn = dim(ineq)[1]
		C = rbind( C, 
				cbind(if(useBench)  matrix(0, ncol=1, nrow = ineqn) else NULL, ineq, matrix(0, ncol = psize+1, nrow = ineqn)), 
				cbind(if(useBench)  matrix(0, ncol=1, nrow = ineqn) else NULL,-ineq, matrix(0, ncol = psize+1, nrow = ineqn)) )		
		d = c(d, -ineqLB, ineqUB)
		b = c(b, rep( 0, 2 * ineqn ) )
		A = rbind( A, matrix(0, nrow = 2 * ineqn, ncol = bsize+m+psize+1) )
		N = c(N, rep(1, 2 * ineqn) )
	}
	# Parameter Bounds
	C = rbind( C, diag(m+bsize+psize+1), -diag(m+bsize+psize+1) )
	d = c(d, c(bLB, -LB, pLB, 0), c(bUB, UB, pUB, 1))
	b = c(b, rep( 0, 2 * (m+bsize+psize+1) ) )
	A = rbind( A, matrix(0, nrow = 2 * (m+bsize+psize+1), ncol = (m+bsize+psize+1)) )
	N = c(N, rep(1, 2 * (m+bsize+psize+1)) )
	return( list( f = f, A = A, C = C, d = d, b = b, A = A, N = N, widx = widx, m = m) )
}

parma.eqcon.qcqpsocp = function(S, Q = NULL, qB = NULL, reward, rewardB, benchmarkS = NULL, 
		eq = NULL, eqB = NULL, ineq = NULL, ineqLB = NULL, ineqUB = NULL, budget = NULL, 
		leverage = NULL, LB, UB, eqSlack = 1e-5)
{
	# setup:
	m = length(LB)
	if(!is.null(leverage)){
		psize = m
		plev = rep(0, m)
		pLB =  rep(0, m)
		pUB = abs(UB)
		uselev = TRUE
	} else{
		psize = 0
		plev = NULL
		pUB = NULL
		pLB = NULL
		uselev = FALSE
	}
	if(!is.null(benchmarkS))
	{
		bsize = 1
		bvalue = 0
		useBench = TRUE
		bLB =   1.00000
		bUB =  -0.99999
		widx = 2:(m+1)
	} else{
		widx = 1:m
		bsize = 0
		bvalue = NULL
		useBench = FALSE
		bLB = NULL
		bUB = NULL
	}
	# setup:
	C = NULL
	d = NULL
	A = NULL
	b = NULL
	N = NULL
	# Objective Function
	# min t
	#  f = [0,...,0, 1]	
	f = c(rep(0, m+bsize), plev, 1)
	
	# x   : [w_1, ..., w_m, t]
	# w_i : weight of asset i
	# t   : auxiliiary variable used for risk minimization
	
	# Risk (LHS): ||A{1}*x|| <= t
	# A{1} : [ Sigma^(1/2) 0]
	## [ Sqrt.Sigma_i_i, ..., Sqrt.Sigma_i_m,  0] + [0]
	## [ .............., ..., ..............., 0] + [0]
	#  [ Sqrt.Sigma_m_i, ..., Sqrt.Sigma_m_m,  0] + [0]
	# Constraint (RHS): C{1}*x+d
	# C{1}: [0,...,0, 1]
	# C{1}*x : t
	if(useBench){
		sqrtS = cbind(matrix(benchmarkS, ncol=1), rbind(matrix(0, nrow=1, ncol = m), S))
		sqrtS[1,2:(m+1)] = benchmarkS[2:(m+1)]
		if(any(eigen(sqrtS, TRUE, only.values=TRUE)$values<1e-10)){
			warning("\nparma: benchmark relative covariance matrix not PD...adjusting to make PD")
			sqrtS <- make.positive.definite(sqrtS, 1e-10)
		}
		sqrtS = .sqrtm(sqrtS)
	} else{
		if(any(eigen(S, TRUE, only.values=TRUE)$values<1e-10)){
			warning("\nparma: benchmark relative covariance matrix not PD...adjusting to make PD")
			S <- make.positive.definite(S, 1e-10)
		}
		sqrtS = .sqrtm(S)
	}
	A = rbind( A, cbind(sqrtS, matrix(0, ncol = psize+1, nrow=m+bsize)) )
	b = c( b, rep(0, m+bsize) )
	C = rbind(C, matrix(c(rep(0, m+bsize), plev, 1), ncol = m+bsize+psize+1) )
	d = c( d, 0 )
	N = c( N, m+bsize)

	# Q is for the quadratic constraint matrix bounded by r
	if(!is.null(Q)){
		v = length(Q)
		for(i in 1:v){
			vm = nrow(Q[[i]])
			if(useBench){
				sqrtQ = .sqrtm( Q[[i]] )
				sqrtQ = cbind(matrix(0, ncol = 1, nrow = vm), sqrtQ)
				sqrtQ = rbind(matrix(0, nrow = 1, ncol = m + 1))
				A = rbind( A, cbind(sqrtQ, matrix(0, ncol = psize+1, nrow=vm)) )
			} else{
				A = rbind( A, cbind(.sqrtm( Q[[i]] ), matrix(0, ncol = psize+1, nrow=vm)) )
			}
			b = c( b, rep(0, vm+bsize) )
			C = rbind(C, matrix(rep(0, bsize+psize+m+1), ncol = bsize+psize+m+1) )
			d = c( d,  qB[i] )
			N = c( N, vm+bsize )
		}
	}
	# Eq Constraints
	if( !is.null( eq ) ){
		eqn = dim( eq )[1]
		C = rbind( C, 
				cbind(if(useBench) matrix(1, ncol = bsize, nrow = eqn) else NULL, eq, 
						matrix(0, ncol = psize+1, nrow = eqn)), 
				cbind(if(useBench) matrix(1, ncol = bsize, nrow = eqn) else NULL, -eq, 
						matrix(0, ncol = psize+1, nrow = eqn)) )
		# might need to add some slack
		d = c( d, -eqB*(1-eqSlack), eqB*(1+eqSlack))
		b = c( b, rep( 0, 2 * eqn ) )
		A = rbind( A, matrix( 0, nrow = 2 * eqn, ncol = m+bsize+psize+1 ) )
		N = c( N, rep(1, 2 * eqn) )
	}
	# Reward Equality Constraint
	if(!is.null(reward)){
		C = rbind( C, 
				matrix( c(bvalue,  reward, plev, 0), ncol = m+bsize+psize+1 ), 
				matrix( c(bvalue, -reward, plev, 0), ncol = m+bsize+psize+1 ) )
		d = c( d, - rewardB*(1-eqSlack) , rewardB*(1+eqSlack))
		b = c(b, 0, 0 )
		A = rbind( A, matrix( 0, nrow = 2, ncol = m+bsize+psize+1) )
		N = c( N, 1, 1 )
	}
	# Budget Constraint
	if( !is.null(budget) ){
		C = rbind( C, 
				matrix( c(bvalue,rep( 1,m),plev,0), ncol = m+bsize+psize+1 ), 
				matrix( c(bvalue,rep(-1,m),plev,0), ncol = m+bsize+psize+1 ) )
		d = c( d, - budget*(1-eqSlack) , budget*(1+eqSlack))
		b = c(b, 0, 0 )
		A = rbind( A, matrix( 0, nrow = 2, ncol = m+bsize+psize+1) )
		N = c( N, 1, 1 )
	}
	if( !is.null(leverage) ){
		A = rbind(A, matrix(0, ncol = m+bsize+psize+1, nrow = 2*m))
		b = c(b, rep(0, 2*m))
		# x+s>=0
		C = rbind(C, cbind( if(useBench)  matrix(0, ncol=1, nrow = m) else NULL, 
						diag(m), diag(m), matrix(0, ncol=1, nrow = m)))
		# -x+s>=0
		C = rbind(C, cbind(if(useBench)   matrix(0, ncol=1, nrow = m) else NULL, 
						-diag(m), diag(m), matrix(0, ncol=1, nrow = m)))
		# -s-x>=0
		d = c(d, rep(0, 2*m))
		A = rbind(A, matrix(0, ncol = m+bsize+psize+1, nrow = 2))
		b = c(b, rep(0, 2))
		C = rbind(C, cbind(matrix(0, ncol=m+bsize, nrow=1), matrix(-1, ncol=m, nrow=1), matrix(0, ncol=1, nrow=1)))
		C = rbind(C, cbind(matrix(0, ncol=m+bsize, nrow=1), matrix( 1, ncol=m, nrow=1), matrix(0, ncol=1, nrow=1)))
		d = c(d, (1+eqSlack)*leverage, -(1-eqSlack)*leverage)
		N = c( N, rep(1, 2*m+2) )
	}
	
	# Ineq Constraints
	if( !is.null(ineq) ){
		ineqn = dim(ineq)[1]
		C = rbind( C, 
				cbind(if(useBench)  matrix(0, ncol=1, nrow = ineqn) else NULL, ineq, matrix(0, ncol = psize+1, nrow = ineqn)), 
				cbind(if(useBench)  matrix(0, ncol=1, nrow = ineqn) else NULL,-ineq, matrix(0, ncol = psize+1, nrow = ineqn)) )		
		d = c(d, -ineqLB, ineqUB)
		b = c(b, rep( 0, 2 * ineqn ) )
		A = rbind( A, matrix(0, nrow = 2 * ineqn, ncol = bsize+m+psize+1) )
		N = c(N, rep(1, 2 * ineqn) )
	}
	# Parameter Bounds
	C = rbind( C, diag(m+bsize+psize+1), -diag(m+bsize+psize+1) )
	d = c(d, c(bLB, -LB, pLB, 0), c(bUB, UB, pUB, 1))
	b = c(b, rep( 0, 2 * (m+bsize+psize+1) ) )
	A = rbind( A, matrix(0, nrow = 2 * (m+bsize+psize+1), ncol = (m+bsize+psize+1)) )
	N = c(N, rep(1, 2 * (m+bsize+psize+1)) )
	return( list( f = f, A = A, C = C, d = d, b = b, A = A, N = N, widx = widx, m = m) )
}

parma.ineqcon.qcqpsocp = function(S, Q = NULL, qB = NULL, reward, rewardB, benchmarkS = NULL, 
		eq = NULL, eqB = NULL, ineq = NULL, ineqLB = NULL, ineqUB = NULL, budget = NULL, 
		leverage = NULL, LB, UB, eqSlack = 1e-5)
{
	# setup:
	m = length(LB)
	if(!is.null(leverage)){
		psize = m
		plev = rep(0, m)
		pLB =  rep(0, m)
		pUB = abs(UB)
		uselev = TRUE
	} else{
		psize = 0
		plev = NULL
		pUB = NULL
		pLB = NULL
		uselev = FALSE
	}
	if(!is.null(benchmarkS))
	{
		bsize = 1
		bvalue = 0
		useBench = TRUE
		bLB =   1.00000
		bUB =  -0.99999
		widx = 2:(m+1)
	} else{
		widx = 1:m
		bsize = 0
		bvalue = NULL
		useBench = FALSE
		bLB = NULL
		bUB = NULL
	}
	# setup:
	C = NULL
	d = NULL
	A = NULL
	b = NULL
	N = NULL
	# Objective Function
	# min t
	#  f = [0,...,0, 1]	
	f = c(rep(0, m+bsize), plev, 1)
	
	# x   : [w_1, ..., w_m, t]
	# w_i : weight of asset i
	# t   : auxiliiary variable used for risk minimization
	
	# Risk (LHS): ||A{1}*x|| <= t
	# A{1} : [ Sigma^(1/2) 0]
	## [ Sqrt.Sigma_i_i, ..., Sqrt.Sigma_i_m,  0] + [0]
	## [ .............., ..., ..............., 0] + [0]
	#  [ Sqrt.Sigma_m_i, ..., Sqrt.Sigma_m_m,  0] + [0]
	# Constraint (RHS): C{1}*x+d
	# C{1}: [0,...,0, 1]
	# C{1}*x : t
	if(useBench){
		sqrtS = cbind(matrix(benchmarkS, ncol=1), rbind(matrix(0, nrow=1, ncol = m), S))
		sqrtS[1,2:(m+1)] = benchmarkS[2:(m+1)]
		if(any(eigen(sqrtS, TRUE, only.values=TRUE)$values<1e-10)){
			warning("\nparma: benchmark relative covariance matrix not PD...adjusting to make PD")
			sqrtS <- make.positive.definite(sqrtS, 1e-10)
		}
		sqrtS = .sqrtm(sqrtS)
	} else{
		if(any(eigen(S, TRUE, only.values=TRUE)$values<1e-10)){
			warning("\nparma: benchmark relative covariance matrix not PD...adjusting to make PD")
			S <- make.positive.definite(S, 1e-10)
		}
		sqrtS = .sqrtm(S)
	}
	A = rbind( A, cbind(sqrtS, matrix(0, ncol = psize+1, nrow=m+bsize)) )
	b = c( b, rep(0, m+bsize) )
	C = rbind(C, matrix(c(rep(0, m+bsize), plev, 1), ncol = m+bsize+psize+1) )
	d = c( d, 0 )
	N = c( N, m+bsize)
	
	# Q is for the quadratic constraint matrix bounded by r
	if(!is.null(Q)){
		v = length(Q)
		for(i in 1:v){
			vm = nrow(Q[[i]])
			if(useBench){
				sqrtQ = .sqrtm( Q[[i]] )
				sqrtQ = cbind(matrix(0, ncol = 1, nrow = vm), sqrtQ)
				sqrtQ = rbind(matrix(0, nrow = 1, ncol = m + 1))
				A = rbind( A, cbind(sqrtQ, matrix(0, ncol = psize+1, nrow=vm)) )
			} else{
				A = rbind( A, cbind(.sqrtm( Q[[i]] ), matrix(0, ncol = psize+1, nrow=vm)) )
			}
			b = c( b, rep(0, vm+bsize) )
			C = rbind(C, matrix(rep(0, bsize+psize+m+1), ncol = bsize+psize+m+1) )
			d = c( d,  qB[i] )
			N = c( N, vm+bsize )
		}
	}
	# Eq Constraints
	if( !is.null( eq ) ){
		eqn = dim( eq )[1]
		C = rbind( C, 
				cbind(if(useBench) matrix(1, ncol = bsize, nrow = eqn) else NULL, eq, 
						matrix(0, ncol = psize+1, nrow = eqn)), 
				cbind(if(useBench) matrix(1, ncol = bsize, nrow = eqn) else NULL, -eq, 
						matrix(0, ncol = psize+1, nrow = eqn)) )
		# might need to add some slack
		d = c( d, -eqB*(1-eqSlack), eqB*(1+eqSlack))
		b = c( b, rep( 0, 2 * eqn ) )
		A = rbind( A, matrix( 0, nrow = 2 * eqn, ncol = m+bsize+psize+1 ) )
		N = c( N, rep(1, 2 * eqn) )
	}
	# Reward Inequality Constraint
	if(!is.null(reward)){
		C = rbind( C, matrix( c(bvalue,  reward, plev, 0), ncol = m+bsize+psize+1 ) )		
		d = c( d, - rewardB)
		b = c(b, 0)
		A = rbind( A, matrix( 0, nrow = 1, ncol = m+bsize+psize+1) )
		N = c( N, 1)
	}
	# Budget Constraint
	if( !is.null(budget) ){
		C = rbind( C, 
				matrix( c(bvalue,rep( 1,m),plev,0), ncol = m+bsize+psize+1 ), 
				matrix( c(bvalue,rep(-1,m),plev,0), ncol = m+bsize+psize+1 ) )
		d = c( d, - budget*(1-eqSlack) , budget*(1+eqSlack))
		b = c(b, 0, 0 )
		A = rbind( A, matrix( 0, nrow = 2, ncol = m+bsize+psize+1) )
		N = c( N, 1, 1 )
	}
	if( !is.null(leverage) ){
		A = rbind(A, matrix(0, ncol = m+bsize+psize+1, nrow = 2*m))
		b = c(b, rep(0, 2*m))
		# x+s>=0
		C = rbind(C, cbind( if(useBench)  matrix(0, ncol=1, nrow = m) else NULL, 
						diag(m), diag(m), matrix(0, ncol=1, nrow = m)))
		# -x+s>=0
		C = rbind(C, cbind(if(useBench)   matrix(0, ncol=1, nrow = m) else NULL, 
						-diag(m), diag(m), matrix(0, ncol=1, nrow = m)))
		# -s-x>=0
		d = c(d, rep(0, 2*m))
		A = rbind(A, matrix(0, ncol = m+bsize+psize+1, nrow = 2))
		b = c(b, rep(0, 2))
		C = rbind(C, cbind(matrix(0, ncol=m+bsize, nrow=1), matrix(-1, ncol=m, nrow=1), matrix(0, ncol=1, nrow=1)))
		C = rbind(C, cbind(matrix(0, ncol=m+bsize, nrow=1), matrix( 1, ncol=m, nrow=1), matrix(0, ncol=1, nrow=1)))
		d = c(d, (1+eqSlack)*leverage, -(1-eqSlack)*leverage)
		N = c( N, rep(1, 2*m+2) )
	}
	
	# Ineq Constraints
	if( !is.null(ineq) ){
		ineqn = dim(ineq)[1]
		C = rbind( C, 
				cbind(if(useBench)  matrix(0, ncol=1, nrow = ineqn) else NULL, ineq, matrix(0, ncol = psize+1, nrow = ineqn)), 
				cbind(if(useBench)  matrix(0, ncol=1, nrow = ineqn) else NULL,-ineq, matrix(0, ncol = psize+1, nrow = ineqn)) )		
		d = c(d, -ineqLB, ineqUB)
		b = c(b, rep( 0, 2 * ineqn ) )
		A = rbind( A, matrix(0, nrow = 2 * ineqn, ncol = bsize+m+psize+1) )
		N = c(N, rep(1, 2 * ineqn) )
	}
	# Parameter Bounds
	C = rbind( C, diag(m+bsize+psize+1), -diag(m+bsize+psize+1) )
	d = c(d, c(bLB, -LB, pLB, 0), c(bUB, UB, pUB, 1))
	b = c(b, rep( 0, 2 * (m+bsize+psize+1) ) )
	A = rbind( A, matrix(0, nrow = 2 * (m+bsize+psize+1), ncol = (m+bsize+psize+1)) )
	N = c(N, rep(1, 2 * (m+bsize+psize+1)) )
	return( list( f = f, A = A, C = C, d = d, b = b, A = A, N = N, widx = widx, m = m) )
}


parma.ineqcon.qpmaxrewardsocp = function(S, reward, riskB, benchmarkS = NULL, 
		eq = NULL, eqB = NULL, ineq = NULL, ineqLB = NULL, ineqUB = NULL, 
		budget = NULL, leverage = NULL, LB, UB, eqSlack = 1e-5)
{
	# setup:
	m = length(LB)
	
	if(!is.null(leverage)){
		psize = m
		plev = rep(0, m)
		pLB =  rep(0, m)
		pUB = abs(UB)
		uselev = TRUE
	} else{
		psize = 0
		plev = NULL
		pUB = NULL
		pLB = NULL
		uselev = FALSE
	}
	if(!is.null(benchmarkS))
	{
		bsize = 1
		bvalue = 0
		useBench = TRUE
		bLB =   1.00000
		bUB =  -0.99999
		widx = 2:(m+1)
	} else{
		widx = 1:m
		bsize = 0
		bvalue = NULL
		useBench = FALSE
		bLB = NULL
		bUB = NULL
	}
	# setup:
	C = NULL
	d = NULL
	A = NULL
	b = NULL
	N = NULL
	
	f = c(if(useBench) 0 else NULL, -reward, plev)
	if(useBench){
		sqrtS = cbind(matrix(benchmarkS, ncol=1), rbind(matrix(0, nrow=1, ncol = m), S))
		sqrtS[1,2:(m+1)] = benchmarkS[2:(m+1)]
		if(any(eigen(sqrtS, TRUE, only.values=TRUE)$values<1e-10)){
			warning("\nparma: benchmark relative covariance matrix not PD...adjusting to make PD")
			sqrtS <- make.positive.definite(sqrtS, 1e-10)
		}
		sqrtS = .sqrtm(sqrtS)
	} else{
		if(any(eigen(S, TRUE, only.values=TRUE)$values<1e-10)){
			warning("\nparma: benchmark relative covariance matrix not PD...adjusting to make PD")
			S <- make.positive.definite(S, 1e-10)
		}
		sqrtS = .sqrtm(S)
	}
	A = rbind( A, cbind(sqrtS, if(uselev) matrix(0, ncol = psize, nrow=m) else NULL) )
	b = c( b, rep(0, m+bsize) )
	C = rbind(C, matrix(c(rep(0, m+bsize), plev), ncol = m+bsize+psize) )
	d = c( d, riskB )
	N = c( N, m+bsize )	
	# Eq Constraints
	if( !is.null( eq ) ){
		eqn = dim( eq )[1]
		C = rbind( C, 
				cbind(if(useBench) matrix(1, ncol = bsize, nrow = eqn) else NULL, eq, if(uselev) matrix(0, ncol = psize, nrow = eqn) else NULL), 
				cbind(if(useBench) matrix(1, ncol = bsize, nrow = eqn) else NULL,-eq, if(uselev) matrix(0, ncol = psize, nrow = eqn) else NULL) )
		# might need to add some slack
		d = c( d, -eqB*(1-eqSlack), eqB*(1+eqSlack))
		b = c( b, rep( 0, 2 * eqn ) )
		A = rbind( A, matrix( 0, nrow = 2 * eqn, ncol = m+bsize+psize) )
		N = c( N, rep(1, 2 * eqn) )
	}
	# Budget Constraint
	if( !is.null(budget) ){
		C = rbind( C, 
				matrix( c(bvalue,rep( 1,m),plev), ncol = m+bsize+psize ), 
				matrix( c(bvalue,rep(-1,m),plev), ncol = m+bsize+psize) )
		d = c( d, - budget*(1-eqSlack) , budget*(1+eqSlack))
		b = c(b, 0, 0 )
		A = rbind( A, matrix( 0, nrow = 2, ncol = m+bsize+psize) )
		N = c( N, 1, 1 )
	}
	if( !is.null(leverage) ){
		A = rbind(A, matrix(0, ncol = m+bsize+psize, nrow = 2*m))
		b = c(b, rep(0, 2*m))
		# x+s>=0
		C = rbind(C, 
				cbind(if(useBench)  matrix(0, ncol=1, nrow = m) else NULL, 
						diag(m), diag(m)))
		# -x+s>=0
		C = rbind(C, 
				cbind(if(useBench)  matrix(0, ncol=1, nrow = m) else NULL, 
						-diag(m), diag(m)))
		# -s-x>=0
		d = c(d, rep(0, 2*m))
		A = rbind(A, matrix(0, ncol = m+bsize+psize, nrow = 2))
		b = c(b, rep(0, 2))
		C = rbind(C, 
				cbind(if(useBench)  matrix(0, ncol=1, nrow = m) else NULL, 
						matrix( 0, ncol=m, nrow=1), 
						matrix(-1, ncol=m, nrow=1)))
		C = rbind(C, 
				cbind(if(useBench)  matrix(0, ncol=1, nrow = m) else NULL, 
						matrix( 0, ncol=m, nrow=1), 
						matrix( 1, ncol=m, nrow=1)))
		d = c(d, (1+eqSlack)*leverage, -(1-eqSlack)*leverage)
		N = c( N, rep(1, 2*m+2) )
	}
	
	# Ineq Constraints
	if( !is.null(ineq) ){
		ineqn = dim(ineq)[1]
		C = rbind( C, 
				cbind(if(useBench)  matrix(0, ncol=1, nrow = ineqn) else NULL,  ineq, if(uselev) matrix(0, ncol = psize, nrow = ineqn) else NULL), 
				cbind(if(useBench)  matrix(0, ncol=1, nrow = ineqn) else NULL, -ineq, if(uselev) matrix(0, ncol = psize, nrow = ineqn) else NULL) )		
		d = c(d, -ineqLB, ineqUB)
		b = c(b, rep( 0, 2 * ineqn ) )
		A = rbind( A, matrix(0, nrow = 2 * ineqn, ncol = bsize+psize+m) )
		N = c(N, rep(1, 2 * ineqn) )
	}
	# Parameter Bounds
	C = rbind( C, diag(m+bsize+psize), -diag(m+bsize+psize) )
	d = c(d, c(bLB, -LB, pLB), c(bUB, UB, pUB))
	b = c(b, rep( 0, 2 * (m+bsize+psize) ) )
	A = rbind( A, matrix(0, nrow = 2 * (m+bsize+psize), ncol = (m+bsize+psize)) )
	N = c(N, rep(1, 2 * (m+bsize+psize)) )
	return( list( f = f, A = A, C = C, d = d, b = b, A = A, N = N, widx = widx , m = m) )
}

# maximize return over risk
parma.ineqcon.qpoptsocp = function(S, reward, eq = NULL, eqB = NULL, benchmarkS = NULL, 
		ineq = NULL, ineqLB = NULL, ineqUB = NULL, budget = NULL, leverage = NULL, LB, UB, 
		eqSlack = 1e-5)
{
	m = length(LB)
	
	if(!is.null(leverage)){
		psize = m
		plev = rep(0, m)
		pLB =  rep(0, m)
		pUB = rep(1e10, m)
		useLev = TRUE
	} else{
		psize = 0
		plev = NULL
		pUB = NULL
		pLB = NULL
		useLev = FALSE
	}
	if(!is.null(benchmarkS))
	{
		bsize = 1
		bvalue = 0
		useBench = TRUE
	} else{
		bsize = 0
		bvalue = NULL
		useBench = FALSE
	}
	midx = 1
	widx = 2:(m+1)
	# setup:
	C = NULL
	d = NULL
	A = NULL
	b = NULL
	N = NULL
	# Objective Function
	# min t
	#  f = [0,...,0, 1]	
	f = c(0, -reward, plev)
	# x   : [w_1, ..., w_m, t]
	# w_i : weight of asset i
	# t   : auxiliiary variable used for risk minimization
	
	# Risk (LHS): ||A{1}*x|| <= t
	# A{1} : [ Sigma^(1/2) 0]
	## [ Sqrt.Sigma_i_i, ..., Sqrt.Sigma_i_m,  0] + [0]
	## [ .............., ..., ..............., 0] + [0]
	#  [ Sqrt.Sigma_m_i, ..., Sqrt.Sigma_m_m,  0] + [0]
	# Constraint (RHS): C{1}*x+d
	# C{1}: [0,...,0, 1]
	# C{1}*x : t
	if(useBench){
		benchmarkS = as.numeric(benchmarkS)
		sqrtS = cbind(matrix(benchmarkS, ncol=1), rbind(matrix(0, nrow=1, ncol = m), S))
		sqrtS[1,2:(m+1)] = benchmarkS[2:(m+1)]
		if(any(eigen(sqrtS, TRUE, only.values=TRUE)$values<1e-12)){
			warning("\nparma: benchmark relative covariance matrix not PD...adjusting to make PD")
			sqrtS <- make.positive.definite(sqrtS, 1e-12)
		}
		sqrtS[1,] = -sqrtS[1,]
		sqrtS = .sqrtm(sqrtS)
	} else{
		if(any(eigen(S, TRUE, only.values=TRUE)$values<1e-12)){
			warning("\nparma: benchmark relative covariance matrix not PD...adjusting to make PD")
			S <- make.positive.definite(S, 1e-12)
		}
		sqrtS = .sqrtm(S)
	}
	A = rbind( A, cbind(if(!useBench) matrix(0, ncol=1, nrow=m) else NULL, sqrtS, if(useLev) matrix(0, ncol = psize, nrow=m+bsize) else NULL) )
	b = c( b, rep(0, m+bsize) )
	C = rbind(C, matrix(c(rep(0, m+1), plev), ncol = m+psize+1) )
	d = c( d, 1 )
	N = c( N, m+bsize )
	
	# Eq Constraints
	if( !is.null( eq ) ){
		eqn = dim( eq )[1]
		eqx1 = cbind(matrix(0, ncol = 1, nrow = eqn), eq, if(useLev) matrix(0, ncol = psize, nrow = eqn) else NULL)
		eqx1[1:eqn, 1] = -eqB*(1-eqSlack)
		eqx2 = cbind(matrix(0, ncol = 1, nrow = eqn),-eq, if(useLev) matrix(0, ncol = psize, nrow = eqn) else NULL)
		eqx1[1:eqn, 1] =  eqB*(1+eqSlack)
		C = rbind( C, eqx1,  eqx2)
		# might need to add some slack
		d = c( d, rep(0, 2*eqn) )
		b = c( b, rep(0, 2*eqn ) )
		A = rbind( A, matrix( 0, nrow = 2 * eqn, ncol = m+psize+1 ) )
		N = c( N, rep(1, 2*eqn) )
	}
	# Budget Constraint
	if( !is.null(budget) ){
		C = rbind( C, 
				matrix( c(-budget*(1-eqSlack), rep( 1,m),plev), ncol = m+psize+1 ), 
				matrix( c( budget*(1+eqSlack), rep(-1,m),plev), ncol = m+psize+1 ) )
		d = c(d, 0, 0)
		b = c(b, 0, 0)
		A = rbind( A, matrix( 0, nrow = 2, ncol = m+psize+1) )
		N = c(N, 1, 1)
	}
	if( !is.null(leverage) ){
		A = rbind(A, matrix(0, ncol = m+psize+1, nrow = 2*m))
		b = c(b, rep(0, 2*m))
		# x+s>=0
		C = rbind(C, cbind( matrix(0, ncol=1, nrow=m), diag(m), diag(m) ))
		# -x+s>=0
		C = rbind(C, cbind( matrix(0, ncol=1, nrow=m),-diag(m), diag(m) ))
		# -s-x>=0
		d = c(d, rep(0, 2*m))
		A = rbind(A, matrix(0, ncol = m+psize+1, nrow = 2))
		b = c(b, rep(0, 2))
		C = rbind(C, cbind(matrix( (1+eqSlack)*leverage, ncol=1, nrow=1), matrix(0, ncol=m, nrow=1), matrix(-1, ncol=m, nrow=1) ))
		C = rbind(C, cbind(matrix(-(1-eqSlack)*leverage, ncol=1, nrow=1), matrix(0, ncol=m, nrow=1), matrix( 1, ncol=m, nrow=1) ))
		d = c(d, 0, 0)
		N = c( N, rep(1, 2*m+2) )
	}
	
	# Ineq Constraints
	if( !is.null(ineq) ){
		ineqn = dim(ineq)[1]
		ineqx1 = cbind(matrix(0, ncol = 1, nrow = ineqn), ineq, if(useLev) matrix(0, ncol = psize, nrow = ineqn) else NULL)
		ineqx1[1:ineqn, 1] = -ineqLB
		ineqx2 = cbind(matrix(0, ncol = 1, nrow = ineqn),-ineq, if(useLev) matrix(0, ncol = psize, nrow = ineqn) else NULL)
		ineqx2[1:ineqn, 1] = ineqUB
		C = rbind( C, ineqx1, ineqx2)
		d = c(d, rep(0, 2 * ineqn))
		b = c(b, rep(0, 2 * ineqn))
		A = rbind( A, matrix(0, nrow = 2 * ineqn, ncol = psize+m+1) )
		N = c(N, rep(1, 2 * ineqn))
	}
	# Parameter Bounds
	lB = diag(m+psize+1)
	lB[2:(m+1),1] = -LB
	uB = -diag(m+psize+1)
	uB[2:(m+1),1] = UB
	C = rbind( C, lB, uB )
	d = c(d, c(0, rep(0, m), pLB), c(50000, rep(0, m), pUB))
	b = c(b, rep( 0, 2 * (m+psize+1) ) )
	A = rbind( A, matrix(0, nrow = 2 * (m+psize+1), ncol = (m+psize+1)) )
	N = c(N, rep(1, 2 * (m+psize+1)) )
	return( list( f = f, A = A, C = C, d = d, b = b, A = A, N = N, widx = widx, midx = midx, m = m) )
}

##################################################################################
.socp.min = function(optvars, eqSlack = 1e-5)
{
	# inequality vs equality
	# QCQP
	idx = optvars$index
	if(idx[3]==1){
		if(!is.null(optvars$Q)){
			setup = parma.ineqcon.qcqpsocp(S = optvars$S, Q = optvars$Q, qB = optvars$qB, 
					reward = optvars$mu, rewardB = optvars$mutarget, 
					benchmarkS = if(idx[2]>0) optvars$benchmarkS else NULL, eq = optvars$eq.mat, 
					eqB = optvars$eqB, ineq = optvars$ineq.mat, 
					ineqLB = optvars$ineq.LB, ineqUB = optvars$ineq.UB, 
					budget = optvars$budget, leverage = optvars$leverage, 
					LB = optvars$LB, UB = optvars$UB, eqSlack = eqSlack)
		} else{
			setup = parma.ineqcon.qpsocp(S = optvars$S, reward = optvars$mu, 
					rewardB = optvars$mutarget, benchmarkS = if(idx[2]>0) optvars$benchmarkS else NULL, 
					eq = optvars$eq.mat, eqB = optvars$eqB, ineq = optvars$ineq.mat, 
					ineqLB = optvars$ineq.LB, ineqUB = optvars$ineq.UB, 
					budget = optvars$budget, leverage = optvars$leverage, 
					LB = optvars$LB, UB = optvars$UB, eqSlack = eqSlack)
		}
	} else {
		if(!is.null(optvars$Q)){
			setup = parma.eqcon.qcqpsocp(S = optvars$S, Q = optvars$Q, qB = optvars$qB, 
					reward = optvars$mu, rewardB = optvars$mutarget, 
					benchmarkS = if(idx[2]>0) optvars$benchmarkS else NULL, eq = optvars$eq.mat, 
					eqB = optvars$eqB, ineq = optvars$ineq.mat, 
					ineqLB = optvars$ineq.LB, ineqUB = optvars$ineq.UB, 
					budget = optvars$budget, leverage = optvars$leverage, 
					LB = optvars$LB, UB = optvars$UB, eqSlack = eqSlack)
		} else{
			setup = parma.eqcon.qpsocp(S = optvars$S, reward = optvars$mu, 
					rewardB = optvars$mutarget, benchmarkS = if(idx[2]>0) optvars$benchmarkS else NULL, 
					eq = optvars$eq.mat, eqB = optvars$eqB, ineq = optvars$ineq.mat, 
					ineqLB = optvars$ineq.LB, ineqUB = optvars$ineq.UB, 
					budget = optvars$budget, leverage = optvars$leverage, 
					LB = optvars$LB, UB = optvars$UB, eqSlack = eqSlack)
		}
	}
	return(setup)
}

.socp.max = function(optvars, eqSlack = 1e-5)
{
	idx = optvars$index
	setup = parma.ineqcon.qpmaxrewardsocp(S = optvars$S, reward = optvars$mu, 
			riskB = optvars$riskB, benchmarkS = if(idx[2]>0) optvars$benchmarkS else NULL, 
			eq = optvars$eq.mat, eqB = optvars$eqB, ineq = optvars$ineq.mat, 
			ineqLB = optvars$ineq.LB, ineqUB = optvars$ineq.UB, 
			budget = optvars$budget, leverage = optvars$leverage, 
			LB = optvars$LB, UB = optvars$UB, eqSlack = eqSlack)
	return(setup)
}

.socp.opt = function(optvars, eqSlack = 1e-5)
{
	idx = optvars$index
	setup = parma.ineqcon.qpoptsocp(S = optvars$S, reward = optvars$mu, 
			benchmarkS = if(idx[2]>0) optvars$benchmarkS else NULL, eq = optvars$eq.mat, 
			eqB = optvars$eqB, 
			ineq = optvars$ineq.mat, ineqLB = optvars$ineq.LB, ineqUB = optvars$ineq.UB, 
			budget = optvars$budget, leverage = optvars$leverage, 
			LB = optvars$LB, UB = optvars$UB, eqSlack = eqSlack)
	return(setup)
}
socpport = function(optvars, control = list(abs.tol = 1e-8, rel.tol = 1e-8, 
				Nu=2, max.iter=5250, BigM.K = 4, BigM.iter = 15), eqSlack = 1e-5,...){

	if(optvars$index[5]==1){
		setup = .socp.min(optvars, eqSlack)
		ans = socpmin.solver(setup, optvars, control, ...)
	} else if(optvars$index[5]==2){	
		setup = .socp.opt(optvars, eqSlack)
		ans = try(socpopt.solver(setup, optvars, control, ...), silent = TRUE)			
	} else{
		setup = .socp.max(optvars, eqSlack)
		ans = socpmax.solver(setup, optvars, control, ...)
	}
	sol = list()
	sol$weights = ans$weights
	sol$reward  = ans$reward
	sol$risk 	= as.numeric(ans$risk)
	sol$status  = ans$status
	# if(!is.null(optvars$SS)) 
	#	sol$riskbudget = sapply(optvars$SS, FUN = function(x) ans$weights %*% x %*% ans$weights)
	return( sol )
}
socpmin.solver = function(setup, optvars, control = list(abs.tol = 1e-8, rel.tol = 1e-8, 
				Nu=2, max.iter=5250, BigM.K = 4, BigM.iter = 15), ...)
{
	sol = try(Socp(setup$f, setup$A, setup$b, setup$C, setup$d, setup$N, 
					x = NULL, z = NULL, w = NULL, control=control), silent = TRUE)
	if(inherits(sol, "try-error")){
		warning(sol)
		status = "non-convergence"
		weights = rep(NA, setup$m)
		risk = NA
		reward = NA
	} else{
		weights = sol$x[setup$widx]
		risk = norm(setup$A[1:setup$N[1],] %*% sol$x, "f")^2
		reward = sum(weights * optvars$mu)
		status = sol$message
	}
	return(list(weights = weights, risk = risk, reward = reward, measure = "var",
					status = status))
}

socpmax.solver = function(setup, optvars, control = list(abs.tol = 1e-8, rel.tol = 1e-8, 
				Nu=2, max.iter=5250, BigM.K = 4, BigM.iter = 15), ...)
{
	sol = try(Socp(setup$f, setup$A, setup$b, setup$C, setup$d, setup$N, 
					x = NULL, z = NULL, w = NULL, control=control), silent = TRUE)
	if(inherits(sol, "try-error")){
		warning(sol)
		status = "non-convergence"
		weights = rep(NA, setup$m)
		risk = NA
		reward = NA
	} else{
		weights = sol$x[setup$widx]
		risk = norm(setup$A[1:setup$N[1],] %*% sol$x, "f")^2
		reward = sum(weights * optvars$mu)
		status = sol$message
		
	}
	return(list(weights = weights, risk = risk, reward = reward, measure = "var",
					status = status))
}


socpopt.solver = function(setup, optvars, control = list(abs.tol = 1e-8, rel.tol = 1e-8, 
				Nu=2, max.iter=5250, BigM.K = 4, BigM.iter = 15), ...)
{
	sol = try(Socp(setup$f, setup$A, setup$b, setup$C, setup$d, setup$N, 
					x = NULL, z = NULL, w = NULL, control=control), silent = TRUE)
	if(inherits(sol, "try-error")){
		warning(sol)
		status = "non-convergence"
		weights = rep(NA, setup$m)
		risk = NA
		reward = NA
		multiplier = NA
	} else{
		weights = sol$x[setup$widx]/sol$x[setup$midx]
		risk = norm(setup$A[1:setup$N[1],] %*% sol$x/sol$x[setup$midx], "f")^2
		reward = sum(weights * optvars$mu)
		multiplier = sol$x[setup$midx]
		status = sol$message
	}
	return(list(weights = weights, risk = risk, reward = reward, multiplier = multiplier, 
					measure = "var", status = status))
}