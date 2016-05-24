
oed_for_slope_over_intercept_homo = function(n, xmin, xmax, theta, f_hetero = NULL, ...){
	rho_to_design(n, xmin, xmax, compute_rho_star_homo(xmin, xmax, theta))
}

rho_to_design = function(n, xmin, xmax, rho){
	num_xmin = round(n * rho)
	if (num_xmin == n){ #just in case...
		num_xmin = num_xmin - 1
	} else if (num_xmin == 0){
		num_xmin = 1
	}
	c(rep(xmin, num_xmin), rep(xmax, n - num_xmin))	
}

#plot_rho_star_by_theta = function(xmin, xmax, thetas = seq(0, 10, length.out = 100), ...){
#	rhostars = array(NA, length(thetas))		
#	for (i in 1 : length(thetas)){
#		rhostars[i] = compute_rho_star_homo(xmin, xmax, theta = thetas[i])
#	}		
#	plot(thetas, rhostars, type = "l", ...)
#}

compute_rho_star_homo = function(xmin = 1, xmax = 10, theta = 1){
	(1 + theta * xmax) / (2 + theta * (xmax + xmin))
}

oed_for_slope_over_intercept_hetero = function(n, xmin, xmax, theta, f_hetero, MaxIter = 6000, MaxFunEvals = 6000, TolFun = 1e-6, NUM_RAND_STARTS = 50){
	#define the objection function in local scope so I can use theta and f_hetero without fear of conflict elsewhere
	Q_prop = function(xs){
		xs = as.numeric(xs)
		#create design matrix
		X = cbind(rep(1, length(xs)), xs)
		#do some intermediate calculations
		XtXinv = solve(t(X) %*% X)
		XtSigmaX = t(X) %*% diag(f_hetero(xs)) %*% X
		var_B = XtXinv %*% XtSigmaX %*% XtXinv
		del_g_beta = t(as.matrix(c(-theta, 1)))
		#asymptotic variance
		as.numeric(del_g_beta %*% var_B %*% t(del_g_beta))
	}
#	f_hetero = function(x){2 - (x - xmin) / (xmax - xmin)}
#	f_hetero = function(x){1}
	options = optimset(MaxIter = MaxIter, MaxFunEvals = MaxFunEvals, TolFun = TolFun)
	#use the homoskedastic as a starting point
	
	x_starts = rbind(
			seq(xmin, xmax, length.out = n)
	)
	for (rho in seq(0, 1, length.out = 2 * n)){
		x_starts = rbind(x_starts, rho_to_design(n, xmin, xmax, rho))
	}
	for (i in 1 : NUM_RAND_STARTS){
		x_starts = rbind(x_starts, runif(n, xmin, xmax))
	}
	x_starts = unique(x_starts)
	
	sols = list()
	#run the Nelder-Mead search from different starting points
	for (i in 1 : nrow(x_starts)){
#		cat(i, "\n")
		sols[[i]] = fminbnd(Q_prop, as.numeric(x_starts[i, ]), rep(xmin - TolFun, n), rep(xmax + TolFun, n), options, verbose = F)
	}
	
	sols_vals = lapply(sols, function(sol){neldermead.get(this = sol, key = "fopt")})
	sol = sols[[which.min(sols_vals)]]
	
	sols_vecs = round(matrix(unlist(lapply(sols, function(sol){sort(neldermead.get(this = sol, key = "xopt")[,1])})), nrow = length(sols), byrow = TRUE), 2)
	sols_vecs
	#return the result as a sorted array for convenience
	sort(neldermead.get(this = sol, key = "xopt")[, 1])	
}