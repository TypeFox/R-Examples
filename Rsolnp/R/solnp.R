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

# Based on the original solnp by Yinyu Ye
# http://www.stanford.edu/~yyye/Col.html


#----------------------------------------------------------------------------------
# The Function SOLNP solves nonlinear programs in standard form:
#
#        minimize              J(P)
#        subject to            EC(P)  =0
#                   IB(:,1)<=  IC(P)  <=IB(:,2)
#                   PB(:,1)<=    P    <=PB(:,2).
#where
#
#  J       : Cost objective scalar function
#  EC      : Equality constraint vector function
#  IC      : Inequality constraint vector function
#  P       : Decision parameter vector
#  IB, PB  : lower and upper bounds for IC and P.
#----------------------------------------------------------------------------------

# control list
#           RHO  : penalty parameter
#           MAJIT: maximum number of major iterations
#           MINIT: maximum number of minor iterations
#           DELTA: relative step size in forward difference evaluation
#           TOL  : tolerance on feasibility and optimality
# defaults RHO=1, MAJIT=10, MINIT=10, DELTA=1.0e-5, TOL=1.0e-4

solnp = function(pars, fun, eqfun = NULL, eqB = NULL, ineqfun = NULL, ineqLB = NULL, ineqUB = NULL, LB = NULL, UB = NULL, control = list(), ...)
{
	# start timer
	tic = Sys.time()
	xnames = names(pars)
	# get environment
	.solnpenv <- environment()
	assign("xnames", xnames, envir = .solnpenv)
	# initiate function count
	assign(".solnp_nfn", 0, envir = .solnpenv)
	assign(".solnp_errors", 0, envir = .solnpenv)
	
	# index of function indicators
	# [1] length of pars
	# [2] has function gradient?
	# [3] has hessian?
	# [4] has ineq?
	# [5] ineq length
	# [6] has jacobian (inequality)
	# [7] has eq?
	# [8] eq length
	# [9] has jacobian (equality)
	# [10] has upper / lower bounds
	# [11] has either lower/upper bounds or ineq
	
	
	ind = rep(0, 11)
	np = ind[1]  = length(pars)
	# lower parameter bounds - indicator
	# lpb[1]=1 means lower/upper bounds present
	# lpb[2]=1 means lower/upper bounds OR inequality bounds present
	
	# do parameter and LB/UB checks
	check1 = .checkpars(pars, LB, UB, .solnpenv)
	
	# .LB and .UB assigned
	
	.LB = get(".LB", envir = .solnpenv)
	.UB = get(".UB", envir = .solnpenv)
	
	
	if( !is.null(.LB) || !is.null(.UB) ) ind[10] = 1
	
	# do function checks and return starting value
	funv = .checkfun(pars, fun, .solnpenv, ...)
	#.solnp_fun assigned
	.solnp_fun = get(".solnp_fun", envir = .solnpenv)
		
	# Analytical Gradient Functionality not yet implemented in subnp function
	
	# gradient and hessian checks
	#if(!is.null(grad)){
	#	gradv = .checkgrad(pars, grad, .solnpenv, ...)
	#	ind[2] = 1
	#} else{
	#	.solnp_gradfun = function(pars, ...) .fdgrad(pars, fun = .solnp_fun, ...)
		ind[2] = 0
	#	gradv = .solnp_gradfun(pars, ...)
	#}
	# .solnp_gradfun(pars, ...) assigned

	.solnp_hessfun = NULL
	ind[3] = 0
	#hessv = NULL
	# .solnp_hessfun(pars, ...) assigned

	# do inequality checks and return starting values
	
	if(!is.null(ineqfun)){
		ineqv 	= .checkineq(pars, ineqfun, ineqLB, ineqUB, .solnpenv, ...)
		ind[4] 	= 1
		nineq 	= length(ineqLB)
		ind[5] 	= nineq
		
		# check for infinities/nans
		.ineqLBx = .ineqLB
		.ineqUBx = .ineqUB
		.ineqLBx[!is.finite(.ineqLB)] = -1e10
		.ineqUBx[!is.finite(.ineqUB)] =  1e10
		ineqx0 	= (.ineqLBx + .ineqUBx)/2
		#if(!is.null(ineqgrad)){
		#	ineqjacv = .cheqjacineq(pars, gradineq, .ineqUB, .ineqLB, .solnpenv, ...)
		#	ind[6] = 1
		#} else{
		# .solnp_ineqjac = function(pars, ...) .fdjac(pars, fun = .solnp_ineqfun, ...)
		ind[6] = 0
		#ineqjacv = .solnp_ineqjac(pars, ...)
		#}
	} else{
		.solnp_ineqfun = function(pars, ...) .emptyfun(pars, ...)
		# .solnp_ineqjac = function(pars, ...) .emptyjac(pars, ...)
		ineqv 	= NULL
		ind[4] 	= 0
		nineq 	= 0
		ind[5] 	= 0
		ind[6] 	= 0
		ineqx0 	= NULL
		.ineqLB = NULL
		.ineqUB = NULL
	}
	# .solnp_ineqfun and .solnp_ineqjac assigned
	# .ineqLB and .ineqUB assigned
	.solnp_ineqfun = get(".solnp_ineqfun", envir = .solnpenv)
	.ineqLB = get(".ineqLB", envir = .solnpenv)
	.ineqUB = get(".ineqUB", envir = .solnpenv)


	# equality checks
	if(!is.null(eqfun)){
		eqv 	= .checkeq(pars, eqfun, eqB, .solnpenv, ...)
		ind[7] 	= 1
		.eqB = get(".eqB", envir = .solnpenv)
		neq 	= length(.eqB)
		ind[8] 	= neq
		#if(!is.null(eqgrad)){
		#	eqjacv = .cheqjaceq(pars, gradeq, .solnpenv, ...)
		#	ind[9] = 1
		#} else{
		#	.solnp_eqjac = function(pars, ...) .fdjac(pars, fun = .solnp_eqfun, ...)
		#	eqjacv = .solnp_eqjac(pars, ...)
			ind[9] = 0
		#}
	} else {
		eqv = NULL
		#eqjacv = NULL
		.solnp_eqfun = function(pars, ...) .emptyfun(pars, ...)
		#.solnp_eqjac = function(pars, ...) .emptyjac(pars, ...)
		ind[7] 	= 0
		neq 	= 0
		ind[8] 	= 0
		ind[9] 	= 0
	}
	# .solnp_eqfun(pars, ...) and .solnp_eqjac(pars, ...) assigned
	# .solnp_eqB assigned
	.solnp_eqfun = get(".solnp_eqfun", envir = .solnpenv)

	if( ind[ 10 ] || ind [ 4 ]) ind[ 11 ] = 1
		
	# parameter bounds (pb)
	pb  = rbind( cbind(.ineqLB, .ineqUB), cbind(.LB, .UB) )
	
	# check control list
	ctrl  = .solnpctrl( control )
	rho   = ctrl[[ 1 ]]
	# maxit = outer iterations
	maxit = ctrl[[ 2 ]]
	# minit = inner iterations
	minit = ctrl[[ 3 ]]
	delta = ctrl[[ 4 ]]
	tol   = ctrl[[ 5 ]]
	trace = ctrl[[ 6 ]]
	
	# total constraints (tc) = no.inequality constraints + no.equality constraints
	tc = nineq + neq
	
	# initialize fn value and inequalities and set to NULL those not needed
	j  = jh = funv
	tt = 0 * .ones(3, 1)
	
	if( tc > 0 ) {
		# lagrange multipliers (lambda)
		lambda = 0 * .ones(tc, 1)
		# constraint vector = [1:neq 1:nineq]
		constraint = c(eqv, ineqv)
		if( ind[4] ) {
			tmpv = cbind(constraint[ (neq + 1):tc ] - .ineqLB, .ineqUB - constraint[ (neq + 1):tc ] )
			testmin = apply( tmpv, 1, FUN = function( x ) min(x[ 1 ], x[ 2 ]) )
			if( all(testmin > 0) ) ineqx0 = constraint[ (neq + 1):tc ]
			constraint[ (neq + 1):tc ] = constraint[ (neq + 1):tc ] - ineqx0
		}
		tt[ 2 ] = .vnorm(constraint)
		if( max(tt[ 2 ] - 10 * tol, nineq, na.rm = TRUE) <= 0 ) rho = 0
	} else{
		lambda = 0
	}
	# starting augmented parameter vector
	p  = c(ineqx0, pars)
	hessv  = diag(np + nineq)
	mu = np
	.solnp_iter = 0
	ob = c(funv, eqv, ineqv)
	
	while( .solnp_iter < maxit ){
		.solnp_iter = .solnp_iter + 1
		.subnp_ctrl = c(rho, minit, delta, tol, trace)
		
		# make the scale for the cost, the equality constraints, the inequality
		# constraints, and the parameters
		if( ind[7] ) {
			# [1 neq]
			vscale = c( ob[ 1 ], rep(1, neq) * max( abs(ob[ 2:(neq + 1) ]) ) )
		} else {
			vscale = 1
		}
		
		if( !ind[ 11 ] ) {
			vscale = c(vscale, p)
		} else {
			# [ 1 neq np]
			vscale = c(vscale, rep( 1, length.out = length(p) ) )
		}
		
		vscale = apply( matrix(vscale, ncol = 1), 1, FUN = function( x ) min( max( abs(x), tol ), 1/tol ) )
		
		res   = .subnp(pars = p, yy = lambda, ob = ob, hessv = hessv, lambda = mu, vscale = vscale, 
				ctrl = .subnp_ctrl, .env = .solnpenv, ...)
		if(get(".solnp_errors", envir =  .solnpenv) == 1){
			maxit = .solnp_iter
		}
		p  = res$p
		lambda  = res$y
		hessv  = res$hessv
		mu = res$lambda
		temp = p[ (nineq + 1):(nineq + np) ]
		funv = .safefunx(temp, .solnp_fun, .env = .solnpenv, ...)
		ctmp = get(".solnp_nfn", envir =  .solnpenv)
		assign(".solnp_nfn", ctmp + 1, envir = .solnpenv)
		
		tempdf = cbind(temp, funv)
		
		if( trace ){
			.report(.solnp_iter, funv, temp)
		}
		
		eqv = .solnp_eqfun(temp, ...)		
		ineqv = .solnp_ineqfun(temp, ...)
		
		ob = c(funv, eqv, ineqv)
		
		tt[ 1 ] = (j - ob[ 1 ]) / max(abs(ob[ 1 ]), 1)
		j = ob[ 1 ]
		
		if( tc > 0 ){
			constraint = ob[ 2:(tc + 1) ]
			
			if( ind[ 4 ] ){
				tempv = rbind( constraint[ (neq + 1):tc ] - pb[ 1:nineq, 1 ],
				              pb[ 1:nineq, 2 ] - constraint[ (neq + 1):tc ] )
				              
				if( min(tempv) > 0 ) {
					p[ 1:nineq ] = constraint[ (neq + 1):tc ]
				}
				
				constraint[ (neq + 1):tc ] = constraint[ (neq + 1):tc ] - p[ 1:nineq ]
			}
			
			tt[ 3 ] = .vnorm(constraint)
			
			if( tt[ 3 ] < 10 * tol ) { 
				rho = 0
				mu  = min(mu, tol)
			}
			
			if( tt[ 3 ] < 5 * tt[ 2 ]) {
				rho = rho/5
			}
			
			if( tt[ 3 ] > 10 * tt[ 2 ]) {
				rho = 5 * max( rho, sqrt(tol) )
			}
			
			if( max( c( tol + tt[ 1 ], tt[ 2 ] - tt[ 3 ] ) ) <= 0 ) { 
				lambda = 0
				hessv = diag( diag ( hessv ) )
			}

			tt[ 2 ] = tt[ 3 ]
		}
		
		if( .vnorm( c(tt[ 1 ], tt[ 2 ]) ) <= tol ) {
			maxit = .solnp_iter
		}
		
		jh = c(jh, j)
	}
	
	if( ind[ 4 ] ) {
		ineqx0 = p[ 1:nineq ]
	}
	
	p = p[ (nineq + 1):(nineq + np) ]
	
	if(get(".solnp_errors", envir =  .solnpenv) == 1){
		convergence = 2
		if( trace ) cat( paste( "\nsolnp--> Solution not reliable....Problem Inverting Hessian.\n", sep="" ) )
	} else{
		if( .vnorm( c(tt[ 1 ], tt[ 2 ]) ) <= tol ) {
			convergence = 0
			if( trace ) cat( paste( "\nsolnp--> Completed in ", .solnp_iter, " iterations\n", sep="" ) )
		} else{
			convergence = 1
			if( trace ) cat( paste( "\nsolnp--> Exiting after maximum number of iterations\n",
							"Tolerance not achieved\n", sep="" ) )
		}
	}
	# end timer
	ctmp = get(".solnp_nfn", envir =  .solnpenv)
	toc = Sys.time() - tic
	names(p) = xnames
	ans = list(pars = p, convergence = convergence, values = as.numeric(jh), lagrange = lambda, 
			hessian = hessv, ineqx0 = ineqx0, nfuneval = ctmp, outer.iter = .solnp_iter, 
			elapsed = toc, vscale = vscale)
	return( ans )
}
