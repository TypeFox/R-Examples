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
.parmautility = function(U = c("CARA", "Power"), method = c("moment", "scenario"), 
		scenario = NULL, M1 = NULL, M2 =  NULL, M3 = NULL, M4 = NULL, RA = 1, 
		budget = 1, LB = rep(0, length(M1)), UB = rep(1, length(M1)))
{
	if(tolower(method[1])=="scemario"){
		stop("\nparma: scenario based utility optimization not yet implemented")
	} else{
		M1 = as.numeric(M1)
		m = length(M1)
		
		if(tolower(U) == "power") stop("\nparma: Power Utility optimization not yet implemented")
		if(dim(M2)[1]!=m) stop("\nparma: M2 dimension 1 does not match M1 length")
		if(dim(M2)[2]!=m) stop("\nparma: M2 dimension 2 does not match M1 length")
		if(!is.null(M3) && !is.null(M4)){
			mx = "c4"
			mom = 4
		} else{
			mx = "c2"
			mom = 2
		}
		if(mx == "c4"){
			if(dim(M3)[1]!=m) stop("\nparma: M3 dimension 2 does not match M1 length")
			if(dim(M3)[2]!= m^2) stop("\nparma: M3 dimension 2 does not match M1^2 length")
			if(dim(M4)[1]!=m) stop("\nparma: M4 dimension 2 does not match M1 length")
			if(dim(M4)[2]!= m^3) stop("\nparma: M4 dimension 2 does not match M1^3 length")
		}
		if(is.null(LB)){
			LB = rep(0, m)
			warning("\nparma: no LB provided...setting Lower Bounds to zero.")
		}
		if(is.null(UB)){
			UB = rep(1, m)
			warning("\nparma: no UB provided...setting Upper Bounds to 1.")
		}
		if(any(UB<LB)) stop("\nparma: UB must be greater than LB.")
		if(length(LB)!=m) LB = rep(LB[1], m)
		if(length(UB)!=m) UB = rep(UB[1], m)
		sol = switch(mx,
				c2 = try(nloptr(x0 = rep(1/m, m), 
								eval_f = cara2.fn,
								eval_grad_f = cara2.gr, 
								lb = LB, 
								ub = UB, 
								eval_g_eq = cara2.eq, 
								eval_jac_g_eq = cara2.eqgr, 
								opts = .nloptr.ctrl(), 
								l=RA, 
								M1 = M1, 
								M2 = M2, budget = budget)),
				c4 = try(nloptr(x0 = rep(1/m, m), 
								eval_f = cara4.fn,
								eval_grad_f = cara4.gr, 
								lb = LB, 
								ub = UB, 
								eval_g_eq = cara4.eq, 
								eval_jac_g_eq = cara4.eqgr, 
								opts = .nloptr.ctrl(), 
								l=RA, 
								M1 = M1, 
								M2 = M2, 
								M3 = M3, 
								M4 = M4, budget = budget), silent = TRUE))
		if( inherits(sol, "try-error") ){
			status = "non-convergence"
			weights = rep(NA, m)
			utility = NA
			reward = NA
		} else{
			status = sol$status
			weights = sol$solution[1:m]
			utility = -sol$objective
			reward  = sum(sol$solution[1:m] * M1)
		}
		model = list(method = method, utility = U, moments = 4)
		solution = list(status = status, weights = weights, utility = utility, reward = reward, risk = NA)
	}
	ans = new("parmaPort",
			solution = solution,
			model = model)
	return(ans)
}
# grad(cara4.fn,  weights, method = "Richardson", method.args = list(), l=RA, M1 = M1, M2 = M2, M3 = M3, M4 = M4)
# cara4.gr(weights, l=RA, M1 = M1, M2 = M2, M3 = M3, M4 = M4)
# cara4.fn(rep(1/15,15), l=RA, M1 = M1, M2 = M2, M3 = M3, M4 = M4)
# cara4.fn(weights, l=RA, M1 = M1, M2 = M2, M3 = M3, M4 = M4)
parma4CARA = function(M1, M2, M3, M4, RA = 1, budget = 1, LB, UB){
	m = length(M1)
	sol = try(nloptr(x0 = rep(1/m, m), 
			eval_f = cara4.fn,
			eval_grad_f = cara4.gr, 
			lb = LB, 
			ub = UB, 
			eval_g_eq = cara4.eq, 
			eval_jac_g_eq = cara4.eqgr, 
			opts = .nloptr.ctrl(), 
			l=RA, 
			M1 = M1, 
			M2 = M2, 
			M3 = M3, 
			M4 = M4, budget = budget), silent = TRUE)
	return(sol)
}

parma2CARA = function(M1, M2, RA = 1, budget = 1, LB, UB){
	m = length(M1)
	sol = try(nloptr(x0 = rep(1/m, m), 
			eval_f = cara2.fn,
			eval_grad_f = cara2.gr, 
			lb = LB, 
			ub = UB, 
			eval_g_eq = cara2.eq, 
			eval_jac_g_eq = cara2.eqgr, 
			opts = .nloptr.ctrl(), 
			l=RA, 
			M1 = M1, 
			M2 = M2))
	return(sol)
}

cara2.eq = function(w, l, M1, M2, budget){
	sum(w)-budget
}

cara2.eqgr = function(w, l, M1, M2, budget){
	matrix(1, ncol = length(w), nrow = 1)
}

cara2.fn = function(w, l, M1, M2, budget){
	n = length(w)
	w = matrix(w, ncol = 1)
	pm1 = t(w) %*% M1
	pm2 = t(w) %*% M2 %*% w
	pm3 = 0
	pm4 = 0
	l2 = (l^2)/2
	l3 = 0
	l4 = 0
	f = -exp(-l*pm1) * (1 + l2*pm2 + l3*0 + l4*0)
	return( as.numeric( -f ) )
}


cara2.gr = function(w, l, M1, M2, budget){
	n = length(w)
	w = matrix(w, ncol = 1)
	pm1 = as.numeric(t(w) %*% M1)
	pm2 = t(w) %*% M2 %*% w
	pm3 = 0
	pm4 = 0
	l2 = (l^2)/2
	l3 = (l^3)/factorial(3)
	l4 = (l^4)/factorial(4)
	dl2 = (l^2)/2
	dl3 = (l^3)/factorial(3)
	dl4 = (l^4)/factorial(4)
	dm1 = M1
	dm2 = 2 * M2 %*% w
	dm3 = 0
	dm4 = 0
	g = (l*exp(-l*pm1))*(M1%*%(1 + l2*pm2 + l3*pm3 + l4*pm4))-exp(-l*pm1)*(dl2*dm2 + dl3*dm3 + dl4*dm4)
	return( as.numeric(-g))
}


cara4.eq = function(w, l, M1, M2, M3, M4, budget){
	sum(w)-budget
}

cara4.eqgr = function(w, l, M1, M2, M3, M4, budget){
	matrix(1, ncol = length(w), nrow = 1)
}



cara4.fn = function(w, l, M1, M2, M3, M4, budget){
	n = length(w)
	w = matrix(w, ncol = 1)
	pm1 = t(w) %*% M1
	pm2 = t(w) %*% M2 %*% w
	pm3 = t(w) %*% M3 %*% kronecker(w, w)
	pm4 = t(w) %*% M4 %*% kronecker(w, kronecker(w, w))
	l2 = (l^2)/2
	l3 = (l^3)/factorial(3)
	l4 = (l^4)/factorial(4)
	f = -exp(-l*pm1) * (1 + l2*pm2 + l3*pm3 + l4*pm4)
	return( as.numeric( -f ) )
}



cara4.gr = function(w, l, M1, M2, M3, M4, budget){
	n = length(w)
	w = matrix(w, ncol = 1)
	pm1 = as.numeric(t(w) %*% M1)
	pm2 = t(w) %*% M2 %*% w
	pm3 = t(w) %*% M3 %*% kronecker(w, w)
	pm4 = t(w) %*% M4 %*% kronecker(w, kronecker(w, w))
	l2 = (l^2)/2
	l3 = (l^3)/factorial(3)
	l4 = (l^4)/factorial(4)
	dl2 = (l^2)/2
	dl3 = (l^3)/factorial(3)
	dl4 = (l^4)/factorial(4)
	dm1 = M1
	dm2 = 2 * M2 %*% w
	dm3 = 3 * M3 %*% kronecker(w, w)
	dm4 = 4 * M4 %*% kronecker(w, kronecker(w, w))
	g = (l*exp(-l*pm1))*(M1%*%(1 + l2*pm2 + l3*pm3 + l4*pm4))-exp(-l*pm1)*(dl2*dm2 + dl3*dm3 + dl4*dm4)
	return( as.numeric(-g))
}

