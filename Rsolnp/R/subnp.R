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

# Based on the original subnp Yinyu Ye
# http://www.stanford.edu/~yyye/Col.html
.subnp = function(pars, yy, ob, hessv, lambda, vscale, ctrl, .env, ...)
{
	.solnp_fun = get(".solnp_fun", envir = .env)
	.solnp_eqfun = get(".solnp_eqfun", envir = .env)
	.solnp_ineqfun = get(".solnp_ineqfun", envir = .env)
	ineqLB = get(".ineqLB", envir = .env)
	ineqUB = get(".ineqUB", envir = .env)
	LB = get(".LB", envir = .env)
	UB = get(".UB", envir = .env)
	#.solnp_gradfun = get(".solnp_gradfun", envir = .env)
	#.solnp_eqjac = get(".solnp_eqjac", envir = .env)
	#.solnp_ineqjac = get(".solnp_ineqjac", envir = .env)
	ind = get("ind", envir = .env)
	
	# pars [nineq + np]	
	rho   = ctrl[ 1 ]
	maxit = ctrl[ 2 ]
	delta = ctrl[ 3 ]
	tol   = ctrl[ 4 ]
	trace = ctrl[ 5 ]
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
	
	
	neq   = ind[ 8 ]
	nineq = ind[ 5 ]
	np    = ind[ 1 ]
	ch    = 1
	alp   = c(0,0,0)
	nc    = neq + nineq
	npic  = np + nineq
	p0    = pars
	
	# pb [ 2 x (nineq + np) ]
	pb    = rbind( cbind(ineqLB, ineqUB), cbind(LB,UB) )
	sob   = numeric()
	ptt   = matrix()
	sc    = numeric()
	
	# scale the cost, the equality constraints, the inequality constraints, 
	# the parameters (inequality parameters AND actual parameters), 
	# and the parameter bounds if there are any
	# Also make sure the parameters are no larger than (1-tol) times their bounds
	# ob [ 1 neq nineq]
	
	ob = ob / vscale[ 1:(nc + 1) ]
	# p0 [np]
	p0 = p0 / vscale[ (neq + 2):(nc + np + 1) ]
	if( ind[ 11 ] ) {
		
		if( !ind[ 10 ] ) {
			mm = nineq
		} else {
			mm = npic
		}
		
		pb = pb / cbind(vscale[ (neq + 2):(neq + mm + 1) ], vscale[ (neq + 2):(neq + mm + 1) ])
	}

	# scale the lagrange multipliers and the Hessian
	
	if( nc > 0) {
		# yy [total constraints = nineq + neq]
		# scale here is [tc] and dot multiplied by yy
		yy = vscale[ 2:(nc + 1) ] * yy / vscale[ 1 ]
	}
	# yy = [zeros 3x1]
	
	# h is [ (np+nineq) x (np+nineq) ]
	#columnvector %*% row vector (size h) then dotproduct h then dotdivide scale[1]
	hessv = hessv * (vscale[ (neq + 2):(nc + np + 1) ] %*% t(vscale[ (neq + 2):(nc + np + 1)]) ) / vscale[ 1 ]
	# h[ 8x8 eye]
	j = ob[ 1 ]
	
	if( ind[4] ) {
		
		if( !ind[7] ) {
			# [nineq x (nineq+np) ]
			a = cbind( -diag(nineq), matrix(0, ncol = np, nrow = nineq) ) 
		} else {
			# [ (neq+nineq) x (nineq+np)]
			a = rbind( cbind( 0 * .ones(neq, nineq), matrix(0, ncol = np, nrow = neq) ), 
					cbind( -diag(nineq), matrix(0, ncol = np, nrow = nineq) ) )
		}
		
	}
	if( ind[7] && !ind[4] ) {
		a = .zeros(neq, np)
	}
	
	if( !ind[7] && !ind[4] ) {
		a = .zeros(1, np)
	}
	
	# gradient
	g= 0 * .ones(npic, 1)
	p = p0 [ 1:npic ]
	if( nc > 0 ) {
		# [ nc ]
		constraint = ob[ 2:(nc + 1) ]
		# constraint [5 0 11 3x1]
		# gradient routine
		for( i in 1:np ) {
			# scale the parameters (non ineq)
			p0[ nineq + i ] = p0[ nineq + i ] + delta
			tmpv = p0[ (nineq + 1):npic ] * vscale[ (nc + 2):(nc + np + 1) ]
			funv 	= .safefunx(tmpv, .solnp_fun, .env, ...)
			eqv 	= .solnp_eqfun(tmpv, ...)
			ineqv 	= .solnp_ineqfun(tmpv, ...)
			ctmp = get(".solnp_nfn", envir =  .env)
			assign(".solnp_nfn", ctmp + 1, envir = .env)

			ob = c(funv, eqv, ineqv) / vscale[ 1:(nc + 1) ]
			g[ nineq + i ]   = (ob[ 1 ] - j) / delta
			a[ , nineq + i ] = (ob[ 2:(nc + 1) ] - constraint) / delta
			p0[ nineq + i ]  = p0[ nineq + i ] - delta
		}
		
		if( ind[4] ) {
			constraint[ (neq + 1):(neq + nineq) ] = constraint[ (neq + 1):(neq + nineq) ] - p0[ 1:nineq ]
		}
		
		# solver messages
		if( .solvecond(a) > 1 / .eps ) { 
			if( trace ) .subnpmsg( "m1" )
		}

		# a(matrix) x columnvector - columnvector
		# b [nc,1]
		b  = a %*% p0 - constraint
		ch = -1
		alp[ 1 ] = tol - max( abs(constraint) )
		
		if( alp[ 1 ] <= 0 ) {
			ch = 1
			
			if( !ind[11] ) {
				# a %*% t(a) gives [nc x nc]
				# t(a) %*% above gives [(np+nc) x 1]
				p0 = p0 - t(a) %*% solve(a %*% t(a), constraint)
				alp[ 1 ] = 1
			}

		}
		
		if( alp[ 1 ] <= 0 ) {
			# this expands p0 to [nc+np+1]
			p0[ npic + 1 ] = 1
			a  = cbind(a, -constraint)
			# cx is rowvector
			cx = cbind(.zeros(1,npic), 1)
			dx = .ones(npic + 1, 1)
			go = 1 
			minit = 0
			
			while( go >= tol ) {
				minit = minit + 1
				# gap [(nc + np) x 2]
				gap = cbind(p0[ 1:mm ] - pb[ , 1 ], pb[ , 2 ] - p0[ 1:mm ] )
				# this sorts every row
				gap = t( apply( gap, 1, FUN=function( x ) sort(x) ) )
				dx[ 1:mm ] = gap[ , 1 ]
				# expand dx by 1
				dx[ npic + 1 ] = p0[ npic + 1 ]
				
				if( !ind[10] ) {
					dx[ (mm + 1):npic ] = max( c(dx[ 1:mm ], 100) ) * .ones(npic - mm, 1)
				}
				# t( a %*% diag( as.numeric(dx) ) ) gives [(np+nc + 1 (or more) x nc]
				# dx * t(cx) dot product of columnvectors
				# qr.solve returns [nc x 1]
				
				# TODO: Catch errors here
				y = try( qr.solve( t( a %*% diag( as.numeric(dx) , length(dx), length(dx) ) ), dx * t(cx) ), silent = TRUE)
				if(inherits(y, "try-error")){
					p = p0 * vscale[ (neq + 2):(nc + np + 1) ]  # unscale the parameter vector
					if( nc > 0 ) {
						y = 0 # unscale the lagrange multipliers
					}
					hessv = vscale[ 1 ] * hessv / (vscale[ (neq + 2):(nc + np + 1) ] %*% t(vscale[ (neq + 2):(nc + np + 1) ]) )
					ans = list(p = p, y = y, hessv = hessv, lambda = lambda)
					assign(".solnp_errors", 1, envir = .env)
					return(ans)
				}
				v = dx * ( dx *(t(cx) - t(a) %*% y) )
				
				if( v[ npic + 1 ] > 0 ) {
					z = p0[ npic + 1 ] / v[ npic + 1 ]
					
					for( i in 1:mm ) {
					
						if( v[ i ] < 0 ) {
							z = min(z, -(pb[ i, 2 ] - p0[ i ]) / v[ i ])
						} else if( v[ i ] > 0 ) { 
							z = min( z, (p0[ i ] - pb[ i , 1 ]) / v[ i ]) 
						}
					}
					
					if( z >= p0[ npic + 1 ] / v[ npic + 1 ] ) {
						p0 = p0 - z * v
					} else {
						p0 = p0 - 0.9 * z * v 
					}
					go = p0[ npic + 1 ]
					
					if( minit >= 10 ) {
						go = 0 
					}
					
				} else {
					go = 0
					minit = 10
				}
				
			}
			
			if( minit >= 10 ) {
				if( trace ) .subnpmsg( "m2" )
			}
			
			a = matrix(a[ , 1:npic ], ncol = npic)
			b = a %*% p0[ 1:npic ]
		}
		
	}
	
	p = p0 [ 1:npic ]
	y = 0
	
	if( ch > 0 ) {
		
		tmpv = p[ (nineq + 1):npic ] * vscale[ (nc + 2):(nc + np + 1) ]
		funv = .safefunx(tmpv, .solnp_fun, .env,...)
		eqv = .solnp_eqfun(tmpv, ...)
		ineqv = .solnp_ineqfun(tmpv, ...)
		ctmp = get(".solnp_nfn", envir =  .env)
		assign(".solnp_nfn", ctmp + 1, envir = .env)
		ob = c(funv, eqv, ineqv) / vscale[ 1:(nc + 1) ]
	}
	
	j = ob[ 1 ]
	
	if( ind[4] ) {
		ob[ (neq + 2):(nc + 1) ] = ob[ (neq + 2):(nc + 1) ] - p[ 1:nineq ]

	}
	
	if( nc > 0 ) {
		ob[ 2:(nc + 1) ] = ob[ 2:(nc + 1) ] - a %*% p + b
		j = ob[ 1 ] - t(yy) %*% matrix(ob[ 2:(nc + 1) ], ncol=1) + rho * .vnorm(ob[ 2:(nc + 1) ]) ^ 2
	}
	
	minit = 0
	while( minit < maxit ) {
		minit = minit + 1
		
		if( ch > 0 ) {
		
			for( i in 1:np ) {
				
				p[ nineq + i ] = p[ nineq + i ] + delta
				tmpv = p[ (nineq + 1):npic ] * vscale[ (nc + 2):(nc + np + 1) ]
				funv 	= .safefunx(tmpv, .solnp_fun, .env, ...)
				eqv 	= .solnp_eqfun(tmpv, ...)
				ineqv 	= .solnp_ineqfun(tmpv, ...)
				ctmp = get(".solnp_nfn", envir =  .env)
				assign(".solnp_nfn", ctmp + 1, envir = .env)
				obm = c(funv, eqv, ineqv) / vscale[ 1:(nc + 1) ]
				
				if( ind[4] ) {
					obm[ (neq + 2):(nc + 1)] = obm[ (neq + 2):(nc + 1) ] - p[ 1:nineq ]
				}
				
				if( nc > 0 ) {
					
					obm[ 2:(nc + 1) ] = obm[ 2:(nc + 1) ] - a %*% p + b
					obm = obm[ 1 ] - t(yy) %*% obm[ 2:(nc + 1) ] + rho * .vnorm(obm[ 2:(nc + 1 ) ]) ^ 2
				}
				
				g[ nineq + i ] = (obm - j) / delta
				p[ nineq + i ] = p[ nineq + i ] - delta
			}
			
			if( ind[4] ) {
				g[ 1:nineq ] = 0
			}
			
		}
		
		if( minit > 1 ) {
			yg = g - yg
			sx = p - sx
			sc[ 1 ] = t(sx) %*% hessv %*% sx
			sc[ 2 ] = t(sx) %*% yg
			
			if( (sc[ 1 ] * sc[ 2 ]) > 0 ) {
				sx = hessv %*% sx
				hessv  = hessv - ( sx %*% t(sx) ) / sc[ 1 ] + ( yg %*% t(yg) ) / sc[ 2 ]
			}
			
		}
		
		dx = 0.01 * .ones(npic, 1)
		if( ind[11] ) {
			
			gap = cbind(p[ 1:mm ] - pb[ , 1 ], pb[ , 2 ] - p[ 1:mm ])
			gap = t( apply( gap, 1, FUN = function( x ) sort(x) ) )
			gap = gap[ , 1 ] + sqrt(.eps) * .ones(mm, 1)
			dx[ 1:mm, 1 ] = .ones(mm, 1) / gap
			if( !ind[10] ){
				dx[ (mm + 1):npic, 1 ] = min (c( dx[ 1:mm, 1 ], 0.01) ) * .ones(npic - mm, 1)
			}
			
		}

		go = -1
		lambda = lambda / 10
		while( go <= 0 ) {
			cz = try(chol( hessv + lambda * diag( as.numeric(dx * dx), length(dx), length(dx) ) ),  silent = TRUE)
			if(inherits(cz, "try-error")){
				p = p * vscale[ (neq + 2):(nc + np + 1) ]  # unscale the parameter vector
				if( nc > 0 ) {
					y = 0 # unscale the lagrange multipliers
				}
				hessv = vscale[ 1 ] * hessv / (vscale[ (neq + 2):(nc + np + 1) ] %*% t(vscale[ (neq + 2):(nc + np + 1) ]) )
				ans = list(p = p, y = y, hessv = hessv, lambda = lambda)
				assign(".solnp_errors", 1, envir = .env)
				return(ans)
			}
			cz = try(solve(cz), silent = TRUE)
			if(inherits(cz, "try-error")){
				p = p * vscale[ (neq + 2):(nc + np + 1) ]  # unscale the parameter vector
				if( nc > 0 ) {
					y = 0 # unscale the lagrange multipliers
				}
				hessv = vscale[ 1 ] * hessv / (vscale[ (neq + 2):(nc + np + 1) ] %*% t(vscale[ (neq + 2):(nc + np + 1) ]) )
				ans = list(p = p, y = y, hessv = hessv, lambda = lambda)
				assign(".solnp_errors", 1, envir = .env)
				return(ans)
			}
			yg = t(cz) %*% g
			
			if( nc == 0 ) {
				u = -cz %*% yg
			} else{
				y = try( qr.solve(t(cz) %*% t(a), yg), silent = TRUE )
				if(inherits(y, "try-error")){
					p = p * vscale[ (neq + 2):(nc + np + 1) ]  # unscale the parameter vector
					if( nc > 0 ) {
						# y = vscale[ 1 ] * y / vscale[ 2:(nc + 1) ] # unscale the lagrange multipliers
						y = 0
					}
					hessv = vscale[ 1 ] * hessv / (vscale[ (neq + 2):(nc + np + 1) ] %*% t(vscale[ (neq + 2):(nc + np + 1) ]) )
					ans = list(p = p, y = y, hessv = hessv, lambda = lambda)
					assign(".solnp_errors", 1, envir = .env)
					return(ans)
				}
				
				u = -cz %*% (yg - ( t(cz) %*% t(a) ) %*% y)
			}
			
			p0 = u[ 1:npic ] + p
			if( !ind[ 11 ] ) {
				go = 1
			} else {
				go = min( c(p0[ 1:mm ] - pb[ , 1 ], pb[ , 2 ] - p0[ 1:mm ]) )
				lambda = 3 * lambda
			}
			
		}
		
		alp[ 1 ] = 0
		ob1 = ob
		ob2 = ob1
		sob[ 1 ] = j
		sob[ 2 ] = j
		ptt = cbind(p, p)
		alp[ 3 ] = 1.0
		ptt = cbind(ptt, p0)
		tmpv = ptt[ (nineq + 1):npic, 3 ] * vscale[ (nc + 2):(nc + np + 1) ]
		funv 	= .safefunx(tmpv, .solnp_fun, .env, ...)
		eqv 	= .solnp_eqfun(tmpv, ...)
		ineqv 	= .solnp_ineqfun(tmpv, ...)
		ctmp = get(".solnp_nfn", envir =  .env)
		assign(".solnp_nfn", ctmp + 1, envir = .env)
		
		ob3 = c(funv, eqv, ineqv) / vscale[ 1:(nc + 1) ]
		sob[ 3 ] = ob3[ 1 ]
		
		if( ind[4] ) {
			ob3[ (neq + 2):(nc + 1) ] = ob3[ (neq + 2):(nc + 1) ] - ptt[ 1:nineq, 3 ]
		}
		
		if( nc > 0 ) {
			ob3[ 2:(nc + 1) ] = ob3[ 2:(nc + 1) ] - a %*% ptt[ , 3 ] + b
			sob[ 3 ] = ob3[ 1 ] - t(yy) %*% ob3[ 2:(nc + 1) ] + rho * .vnorm(ob3[ 2:(nc + 1) ]) ^ 2
		}
		
		go = 1
		while( go > tol ) {
			alp[ 2 ] = (alp[ 1 ] + alp[ 3 ]) / 2
			ptt[ , 2 ] = (1 - alp[ 2 ]) * p + alp[ 2 ] * p0
			tmpv = ptt[ (nineq + 1):npic, 2 ] * vscale[ (nc + 2):(nc + np + 1) ]
			funv 	= .safefunx(tmpv, .solnp_fun, .env, ...)
			eqv 	= .solnp_eqfun(tmpv, ...)
			ineqv 	= .solnp_ineqfun(tmpv, ...)
			ctmp = get(".solnp_nfn", envir =  .env)
			assign(".solnp_nfn", ctmp + 1, envir = .env)
			
			ob2 = c(funv, eqv, ineqv) / vscale[ 1:(nc + 1) ]
			
			sob[ 2 ] = ob2[ 1 ]

			if( ind[4] ) {
				ob2[ (neq + 2):(nc + 1) ] = ob2[ (neq + 2):(nc + 1) ] - ptt[ 1:nineq , 2 ]
			}
			
			if( nc > 0 ) {
				ob2[ 2:(nc + 1) ] = ob2[ 2:(nc + 1) ] - a %*% ptt[ , 2 ] + b
				sob[ 2 ] = ob2[ 1 ] - t(yy) %*% ob2[ 2:(nc + 1) ] + rho * .vnorm(ob2[ 2:(nc + 1) ]) ^ 2
			}
			
			obm = max(sob)
			
			if( obm < j ) {
				obn = min(sob)
				go = tol * (obm - obn) / (j - obm)
			}

			condif1 = sob[ 2 ] >= sob[ 1 ]
			condif2 = sob[ 1 ] <= sob[ 3 ] && sob[ 2 ] < sob[ 1 ]
			condif3 = sob[ 2 ] <  sob[ 1 ] && sob[ 1 ] > sob[ 3 ]
			
			if( condif1 ) {
				sob[ 3 ] = sob[ 2 ]
				ob3 = ob2
				alp[ 3 ] = alp[ 2 ]
				ptt[ , 3 ] = ptt[ , 2 ]
			}
			
			if( condif2 ) {
				sob[ 3 ] = sob[ 2 ]
				ob3 = ob2
				alp[ 3 ] = alp[ 2 ]
				ptt[ , 3 ] = ptt[ , 2 ]
			}
			
			if( condif3 ) {
				sob[ 1 ] = sob[ 2 ]
				ob1 = ob2
				alp[ 1 ] = alp[ 2 ]
				ptt[ , 1 ] = ptt[ , 2 ]
			}
			
			if( go >= tol ) {
				go = alp[ 3 ] - alp[ 1 ]
			}
			
		}
		
		sx = p
		yg = g
		ch = 1
		obn = min(sob)
		if( j <= obn ) {
			maxit = minit
		}
		
		reduce = (j - obn) / ( 1 + abs(j) )
		
		if( reduce < tol ) {
			maxit = minit
		}
		
		condif1 = sob[ 1 ] <  sob[ 2 ]
		condif2 = sob[ 3 ] <  sob[ 2 ] && sob[ 1 ] >= sob[ 2 ]
		condif3 = sob[ 1 ] >= sob[ 2 ] && sob[ 3 ] >= sob[ 2 ]
		
		if( condif1 ) {
			j = sob[ 1 ]
			p = ptt[ , 1 ]
			ob = ob1
		}
		
		if( condif2 ) {
			j = sob [ 3 ]
			p = ptt[ , 3 ]
			ob = ob3
		}
		
		if( condif3 ) {
			j = sob[ 2 ]
			p = ptt[ , 2 ]
			ob = ob2
		}
		
	}
	
	p = p * vscale[ (neq + 2):(nc + np + 1) ]  # unscale the parameter vector
	
	if( nc > 0 ) {
		y = vscale[ 1 ] * y / vscale[ 2:(nc + 1) ] # unscale the lagrange multipliers
	}
	
	hessv = vscale[ 1 ] * hessv / (vscale[ (neq + 2):(nc + np + 1) ] %*% t(vscale[ (neq + 2):(nc + np + 1) ]) )	
	if( reduce > tol ) {
		if( trace ) .subnpmsg( "m3" )
	}

	ans = list(p = p, y = y, hessv = hessv, lambda = lambda)
	return( ans )
}


.fdgrad = function(i, p, delta, np, vscale, constraint, j, nineq, npic, nc, .solnp_fun, 
		.solnp_eqfun, .solnp_ineqfun, .env, ...)
{
	ans = list()
	px = p
	px[ nineq + i ] = px[ nineq + i ] + delta
	tmpv = px[ (nineq + 1):npic ] * vscale[ (nc + 2):(nc + np + 1) ]
	funv = .safefunx(tmpv, .solnp_fun, .env, ...)
	eqv 	= .solnp_eqfun(tmpv, ...)
	ineqv 	= .solnp_ineqfun(tmpv, ...)
	ctmp = get(".solnp_nfn", envir =  .env)
	assign(".solnp_nfn", ctmp + 1, envir = .env)
	ob = c(funv, eqv, ineqv) / vscale[ 1:(nc + 1) ]
	ans$p	= px
	ans$ob 	= ob
	ans$g   = (ob[ 1 ] - j) / delta
	ans$a   = (ob[ 2:(nc + 1) ] - constraint) / delta
	return( ans )
}


