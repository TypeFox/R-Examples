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

#---------------------------------------------------------------------------------
# optimization by randomly restarted parameters using simulated parameter strategy
# Alexios Ghalanos 2010
#---------------------------------------------------------------------------------

# allowed distributions:
# 1: uniform (no confidence in the location of the parameter...somewhere in LB-UB space)
# 2: truncnorm (high confidence in the location of the parameter)
# 3: normal (Uncertainty in Lower and Upper bounds, but some idea about the dispersion about the location)
# ...

gosolnp = function(pars = NULL, fixed = NULL, fun, eqfun = NULL, eqB = NULL, ineqfun = NULL, ineqLB = NULL,
		ineqUB = NULL, LB = NULL, UB = NULL, control = list(), distr = rep(1, length(LB)), distr.opt = list(),
		n.restarts = 1, n.sim = 20000, cluster = NULL, rseed = NULL, ...)
{
	if( !is.null(pars) ) gosolnp_parnames = names(pars) else gosolnp_parnames = NULL
	if(is.null(control$trace)) trace = FALSE else trace = as.logical(control$trace)
	if(is.null(control$eval.type)) parmethod = 1 else parmethod = as.integer(min(abs(control$eval.type),2))
	if(parmethod == 0) parmethod = 1
	control$eval.type = NULL
	# use a seed to initialize random no. generation
	if(is.null(rseed)) rseed = as.numeric(Sys.time()) else rseed = as.integer(rseed)
	# function requires both upper and lower bounds
	if(is.null(LB))
		stop("\ngosolnp-->error: the function requires lower parameter bounds\n", call. = FALSE)
	if(is.null(UB))
		stop("\ngosolnp-->error: the function requires upper parameter bounds\n", call. = FALSE)
	# allow for fixed parameters (i.e. non randomly chosen), but require pars vector in that case
	if(!is.null(fixed) && is.null(pars))
		stop("\ngosolnp-->error: you need to provide a pars vector if using the fixed option\n", call. = FALSE)
	if(!is.null(pars)) n = length(pars) else n = length(LB)

	np = 1:n

	if(!is.null(fixed)){
		# make unique
		fixed = unique(fixed)
		# check for violations in indices
		if(any(is.na(match(fixed, np))))
			stop("\ngosolnp-->error: fixed indices out of bounds\n", call. = FALSE)
	}
	# check distribution options
	# truncated normal
	if(any(distr == 2)){
		d2 = which(distr == 2)
		for(i in 1:length(d2)) {
			if(is.null(distr.opt[[d2[i]]]$mean))
				stop(paste("\ngosolnp-->error: distr.opt[[,",d2[i],"]] missing mean\n", sep = ""), call. = FALSE)
			if(is.null(distr.opt[[d2[i]]]$sd))
				stop(paste("\ngosolnp-->error: distr.opt[[,",d2[i],"]] missing sd\n", sep = ""), call. = FALSE)
		}
	}
	#  normal
	if(any(distr == 3)){
		d3 = which(distr == 3)
		for(i in 1:length(d3)) {
			if(is.null(distr.opt[[d3[i]]]$mean))
				stop(paste("\ngosolnp-->error: distr.opt[[,",d3[i],"]] missing mean\n", sep = ""), call. = FALSE)
			if(is.null(distr.opt[[d3[i]]]$sd))
				stop(paste("\ngosolnp-->error: distr.opt[[,",d3[i],"]] missing sd\n", sep = ""), call. = FALSE)
		}
	}
	# setup cluster exports:
	if( !is.null(cluster) ){
		clusterExport(cluster, c("gosolnp_parnames", "fun", "eqfun",
						"eqB", "ineqfun", "ineqLB", "ineqUB", "LB", "UB"), envir = environment())
		if( !is.null(names(list(...))) ){
			# evaluate promises
			xl = names(list(...))
			for(i in 1:length(xl)){
				eval(parse(text=paste(xl[i],"=list(...)[[i]]",sep="")))
			}
			clusterExport(cluster, names(list(...)), envir = environment())
		}
		clusterEvalQ(cluster, require(Rsolnp))
	}
	# initiate random search
	gosolnp_rndpars = switch(parmethod,
			.randpars(pars = pars, fixed = fixed, fun = fun, eqfun = eqfun,
					eqB = eqB, ineqfun = ineqfun, ineqLB = ineqLB,
					ineqUB = ineqUB, LB = LB, UB = UB, distr = distr,
					distr.opt = distr.opt, n.restarts = n.restarts,
					n.sim = n.sim, trace = trace, rseed = rseed,
					gosolnp_parnames = gosolnp_parnames, cluster = cluster, ...),
			.randpars2(pars = pars, fixed = fixed, fun = fun, eqfun = eqfun,
					eqB = eqB, ineqfun = ineqfun, ineqLB = ineqLB,
					ineqUB = ineqUB, LB = LB, UB = UB, distr = distr,
					distr.opt = distr.opt, n.restarts = n.restarts,
					n.sim = n.sim, rseed = rseed, trace = trace,
					gosolnp_parnames = gosolnp_parnames, cluster = cluster, ...))

	gosolnp_rndpars = gosolnp_rndpars[,1:n, drop = FALSE]
	# initiate solver restarts
	if( trace ) cat("\ngosolnp-->Starting Solver\n")
	solution = vector(mode = "list", length = n.restarts)
	if( !is.null(cluster) )
	{
		clusterExport(cluster, c("gosolnp_rndpars"), envir = environment())
		solution = parLapply(cluster, as.list(1:n.restarts), fun = function(i) {
					xx = gosolnp_rndpars[i,]
					names(xx) = gosolnp_parnames
					ans = try(solnp(pars = xx, fun = fun, eqfun = eqfun,
									eqB = eqB, ineqfun = ineqfun,
									ineqLB = ineqLB, ineqUB = ineqUB,
									LB = LB, UB = UB,
									control = control, ...), silent = TRUE)
					if(inherits(ans, "try-error")){
						ans = list()
						ans$values = 1e10
						ans$convergence = 0
						ans$pars = rep(NA, length(xx))
					}
					return( ans )
				})
	} else {
		solution = lapply(as.list(1:n.restarts), FUN = function(i){
					xx = gosolnp_rndpars[i,]
					names(xx) = gosolnp_parnames
					ans = try(solnp(pars = xx, fun = fun, eqfun = eqfun,
									eqB = eqB, ineqfun = ineqfun,
									ineqLB = ineqLB, ineqUB = ineqUB,
									LB = LB, UB = UB,
									control = control, ...), silent = TRUE)
					if(inherits(ans, "try-error")){
						ans = list()
						ans$values = 1e10
						ans$convergence = 0
						ans$pars = rep(NA, length(xx))
					}
					return( ans )
				})
	}
	if(n.restarts>1){
		best = sapply(solution, FUN = function(x) if(x$convergence!=0) NA else x$values[length(x$values)])
		if(all(is.na(best)))
			stop("\ngosolnp-->Could not find a feasible starting point...exiting\n", call. = FALSE)
		nb = which(best == min(best, na.rm = TRUE))[1]
		solution = solution[[nb]]
		if( trace ) cat("\ngosolnp-->Done!\n")
		solution$start.pars = gosolnp_rndpars[nb,]
		names(solution$start.pars) = gosolnp_parnames
		solution$rseed = rseed
	} else{
		solution = solution[[1]]
		solution$start.pars = gosolnp_rndpars[1,]
		names(solution$start.pars) = gosolnp_parnames
		solution$rseed = rseed
	}
	return(solution)
}



startpars = function(pars = NULL, fixed = NULL, fun, eqfun = NULL, eqB = NULL,
		ineqfun = NULL, ineqLB = NULL, ineqUB = NULL, LB = NULL, UB = NULL,
		distr = rep(1, length(LB)), distr.opt = list(), n.sim = 20000, cluster = NULL,
		rseed = NULL, bestN = 15, eval.type = 1, trace = FALSE, ...)
{
	if( !is.null(pars) ) gosolnp_parnames = names(pars) else gosolnp_parnames = NULL
	if(is.null(eval.type)) parmethod = 1 else parmethod = as.integer(min(abs(eval.type),2))
	if(parmethod == 0) parmethod = 1
	eval.type = NULL
	#trace = FALSE
	# use a seed to initialize random no. generation
	if(is.null(rseed)) rseed = as.numeric(Sys.time()) else rseed = as.integer(rseed)
	# function requires both upper and lower bounds
	if(is.null(LB))
		stop("\nstartpars-->error: the function requires lower parameter bounds\n", call. = FALSE)
	if(is.null(UB))
		stop("\nstartpars-->error: the function requires upper parameter bounds\n", call. = FALSE)

	# allow for fixed parameters (i.e. non randomly chosen), but require pars vector in that case
	if(!is.null(fixed) && is.null(pars))
		stop("\nstartpars-->error: you need to provide a pars vector if using the fixed option\n", call. = FALSE)
	if(!is.null(pars)) n = length(pars) else n = length(LB)

	np = seq_len(n)

	if(!is.null(fixed)){
		# make unique
		fixed = unique(fixed)
		# check for violations in indices
		if(any(is.na(match(fixed, np))))
			stop("\nstartpars-->error: fixed indices out of bounds\n", call. = FALSE)
	}

	# check distribution options
	# truncated normal
	if(any(distr == 2)){
		d2 = which(distr == 2)
		for(i in 1:length(d2)) {
			if(is.null(distr.opt[[d2[i]]]$mean))
				stop(paste("\nstartpars-->error: distr.opt[[,",d2[i],"]] missing mean\n", sep = ""), call. = FALSE)
			if(is.null(distr.opt[[d2[i]]]$sd))
				stop(paste("\nstartpars-->error: distr.opt[[,",d2[i],"]] missing sd\n", sep = ""), call. = FALSE)
		}
	}
	#  normal
	if(any(distr == 3)){
		d3 = which(distr == 3)
		for(i in 1:length(d3)) {
			if(is.null(distr.opt[[d3[i]]]$mean))
				stop(paste("\nstartpars-->error: distr.opt[[,",d3[i],"]] missing mean\n", sep = ""), call. = FALSE)
			if(is.null(distr.opt[[d3[i]]]$sd))
				stop(paste("\nstartpars-->error: distr.opt[[,",d3[i],"]] missing sd\n", sep = ""), call. = FALSE)
		}
	}

	# setup cluster exports:
	if( !is.null(cluster) ){
		clusterExport(cluster, c("gosolnp_parnames", "fun", "eqfun",
						"eqB", "ineqfun", "ineqLB", "ineqUB", "LB", "UB"), envir = environment())
		if( !is.null(names(list(...))) ){
			# evaluate promises
			xl = names(list(...))
			for(i in 1:length(xl)){
			  eval(parse(text = paste(xl[i], "=list(...)", "[[" , i, "]]", sep = "")))
			}
			clusterExport(cluster, names(list(...)), envir = environment())
		}
		if( !is.null(names(list(...))) ) parallel::clusterExport(cluster, names(list(...)), envir = environment())
		clusterEvalQ(cluster, require(Rsolnp))
	}

	# initiate random search
	gosolnp_rndpars = switch(parmethod,
			.randpars(pars = pars, fixed = fixed, fun = fun, eqfun = eqfun,
					eqB = eqB, ineqfun = ineqfun, ineqLB = ineqLB, ineqUB = ineqUB,
					LB = LB, UB = UB, distr = distr, distr.opt = distr.opt,
					n.restarts = as.integer(bestN), n.sim = n.sim, trace = trace,
					rseed = rseed, gosolnp_parnames = gosolnp_parnames,
					cluster = cluster, ...),
			.randpars2(pars = pars, fixed = fixed, fun = fun, eqfun = eqfun,
					eqB = eqB, ineqfun = ineqfun, ineqLB = ineqLB, ineqUB = ineqUB,
					LB = LB, UB = UB, distr = distr, distr.opt = distr.opt,
					n.restarts = as.integer(bestN), n.sim = n.sim, trace = trace,
					rseed = rseed, gosolnp_parnames = gosolnp_parnames,
					cluster = cluster, ...))
	return(gosolnp_rndpars)
}


.randpars = function(pars, fixed, fun, eqfun, eqB,  ineqfun, ineqLB, ineqUB,
		LB, UB, distr, distr.opt, n.restarts, n.sim, trace = TRUE, rseed,
		gosolnp_parnames, cluster, ...)
{
	if( trace ) cat("\nCalculating Random Initialization Parameters...")
	N = length(LB)
	gosolnp_rndpars = matrix(NA, ncol = N, nrow = n.sim * n.restarts)
	if(!is.null(fixed)) for(i in 1:length(fixed)) gosolnp_rndpars[,fixed[i]] = pars[fixed[i]]
	nf = 1:N
	if(!is.null(fixed)) nf = nf[-c(fixed)]
	m = length(nf)
	set.seed(rseed)
	for(i in 1:m){
		j = nf[i]
		gosolnp_rndpars[,j] = switch(distr[j],
				.distr1(LB[j], UB[j], n.restarts*n.sim),
				.distr2(LB[j], UB[j], n.restarts*n.sim, mean = distr.opt[[j]]$mean, sd = distr.opt[[j]]$sd),
				.distr3(n.restarts*n.sim, mean = distr.opt[[j]]$mean, sd = distr.opt[[j]]$sd)
		)
	}

	if( trace ) cat("ok!\n")

	if(!is.null(ineqfun)){
		if( trace ) cat("\nExcluding Inequality Violations...\n")
		ineqv = matrix(NA, ncol = length(ineqLB), nrow = n.restarts*n.sim)
		# ineqv = t(apply(rndpars, 1, FUN = function(x) ineqfun(x)))
		if(length(ineqLB) == 1){
			ineqv = apply(gosolnp_rndpars, 1, FUN = function(x){
						names(x) = gosolnp_parnames
						ineqfun(x, ...)} )
			lbviol = sum(ineqv<ineqLB)
			ubviol = sum(ineqv>ineqUB)
			if( lbviol > 0 | ubviol > 0 ){
				vidx = c(which(ineqv<ineqLB), which(ineqv>ineqUB))
				vidx = unique(vidx)
				gosolnp_rndpars = gosolnp_rndpars[-c(vidx),,drop=FALSE]
				lvx = length(vidx)
			} else{
				vidx = 0
				lvx = 0
			}
		} else{
			ineqv = t(apply(gosolnp_rndpars, 1, FUN = function(x){
								names(x) = gosolnp_parnames
								ineqfun(x, ...)} ))

			# check lower and upper violations
			lbviol = apply(ineqv, 1, FUN = function(x) sum(any(x<ineqLB)))
			ubviol = apply(ineqv, 1, FUN = function(x) sum(any(x>ineqUB)))
			if( any(lbviol > 0) | any(ubviol > 0) ){
				vidx = c(which(lbviol>0), which(ubviol>0))
				vidx = unique(vidx)
				gosolnp_rndpars = gosolnp_rndpars[-c(vidx),,drop=FALSE]
				lvx = length(vidx)

			} else{
				vidx = 0
				lvx = 0
			}
		}
		if( trace ) cat(paste("\n...Excluded ", lvx, "/",n.restarts*n.sim, " Random Sequences\n", sep = ""))
	}
	# evaluate function value
	if( trace ) cat("\nEvaluating Objective Function with Random Sampled Parameters...")
	if( !is.null(cluster) ){
		nx = dim(gosolnp_rndpars)[1]
		clusterExport(cluster, c("gosolnp_rndpars", ".safefun"), envir = environment())
		evfun = parLapply(cluster, as.list(1:nx), fun = function(i){
					.safefun(gosolnp_rndpars[i, ], fun, gosolnp_parnames, ...)
				})
		evfun = as.numeric( unlist(evfun) )
	} else{
		evfun = apply(gosolnp_rndpars, 1, FUN = function(x) .safefun(x, fun, gosolnp_parnames, ...))
	}
	if( trace ) cat("ok!\n")
	if( trace ) cat("\nSorting and Choosing Best Candidates for starting Solver...")
	z = sort.int(evfun, index.return = T)
	ans = gosolnp_rndpars[z$ix[1:n.restarts],,drop = FALSE]
	prtable = cbind(ans, z$x[1:n.restarts])
	if( trace ) cat("ok!\n")
	colnames(prtable) = c(paste("par", 1:N, sep = ""), "objf")
	if( trace ){
		cat("\nStarting Parameters and Starting Objective Function:\n")
		if(n.restarts == 1) print(t(prtable), digits = 4) else print(prtable, digits = 4)
	}
	return(prtable)
}

# form a barrier function before passing the parameters
.randpars2 = function(pars, fixed, fun, eqfun, eqB,  ineqfun, ineqLB, ineqUB, LB,
		UB, distr, distr.opt, n.restarts, n.sim, rseed, trace = TRUE,
		gosolnp_parnames, cluster, ...)
{
	if( trace ) cat("\nCalculating Random Initialization Parameters...")
	N = length(LB)
	gosolnp_idx = "a"
	gosolnp_R = NULL
	if(!is.null(ineqfun) && is.null(eqfun) ){
		gosolnp_idx = "b"
		gosolnp_R = 100
	}
	if( is.null(ineqfun) && !is.null(eqfun) ){
		gosolnp_idx = "c"
		gosolnp_R = 100
	}
	if(!is.null(ineqfun) && !is.null(eqfun) ){
		gosolnp_idx = "d"
		gosolnp_R = c(100,100)
	}
	gosolnp_rndpars = matrix(NA, ncol = N, nrow = n.sim * n.restarts)
	if(!is.null(fixed)) for(i in 1:length(fixed)) gosolnp_rndpars[,fixed[i]] = pars[fixed[i]]
	nf = 1:N
	if(!is.null(fixed)) nf = nf[-c(fixed)]
	gosolnp_m = length(nf)
	set.seed(rseed)
	for(i in 1:gosolnp_m){
		j = nf[i]
		gosolnp_rndpars[,j] = switch(distr[j],
				.distr1(LB[j], UB[j], n.restarts*n.sim),
				.distr2(LB[j], UB[j], n.restarts*n.sim, mean = distr.opt[[j]]$mean, sd = distr.opt[[j]]$sd),
				.distr3(n.restarts*n.sim, mean = distr.opt[[j]]$mean, sd = distr.opt[[j]]$sd)
		)
	}
	if( trace ) cat("ok!\n")
	# Barrier Function
	pclfn = function(x){
		z=x
		z[x<=0] = 0
		z[x>0] = (0.9+z[x>0])^2
		z
	}
	.lagrfun = function(pars, m, idx, fun, eqfun = NULL, eqB = 0, ineqfun = NULL, ineqLB = NULL, ineqUB = NULL, ...)
	{
		fn = switch(idx,
				"a" = fun(pars[1:m], ...),
				"b" = fun(pars[1:m], ...) + pars[m+1]* sum( pclfn( c(ineqLB - ineqfun(pars[1:m], ...), ineqfun(pars[1:m], ...) - ineqUB) ) ),
				"c" = fun(pars[1:m], ...) + sum( (eqfun(pars[1:m], ...) - eqB )^2 / pars[m+1]),
				"d" = fun(pars[1:m], ...) + sum( (eqfun(pars[1:m], ...) - eqB )^2 / pars[m+1]) + pars[m+2]* sum( pclfn( c(ineqLB - ineqfun(pars[1:m], ...), ineqfun(pars[1:m], ...) - ineqUB) ) ) )
		return(fn)
	}

	# evaluate function value
	if( trace ) cat("\nEvaluating Objective Function with Random Sampled Parameters...")
	if( !is.null(cluster) ){
		nx = dim(gosolnp_rndpars)[1]
		clusterExport(cluster, c("gosolnp_rndpars", "gosolnp_m", "gosolnp_idx",
						"gosolnp_R"), envir = environment())
		clusterExport(cluster, c("pclfn", ".lagrfun"), envir = environment())
		evfun = parallel::parLapply(cluster, as.list(1:nx), fun = function(i){
					.lagrfun(c(gosolnp_rndpars[i,], gosolnp_R), gosolnp_m,
							gosolnp_idx, fun, eqfun, eqB, ineqfun, ineqLB,
							ineqUB, ...)
				})
		evfun = as.numeric( unlist(evfun) )
	} else{
		evfun = apply(gosolnp_rndpars, 1, FUN = function(x){
					.lagrfun(c(x,gosolnp_R), gosolnp_m, gosolnp_idx, fun, eqfun,
							eqB, ineqfun, ineqLB, ineqUB, ...)})
	}
	if( trace ) cat("ok!\n")
	if( trace ) cat("\nSorting and Choosing Best Candidates for starting Solver...")
	z = sort.int(evfun, index.return = T)
	#distmat = dist(evfun, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
	ans = gosolnp_rndpars[z$ix[1:n.restarts],,drop = FALSE]
	prtable = cbind(ans, z$x[1:n.restarts])
	colnames(prtable) = c(paste("par", 1:N, sep = ""), "objf")
	if( trace ){
		cat("\nStarting Parameters and Starting Objective Function:\n")
		if(n.restarts == 1) print(t(prtable), digits = 4) else print(prtable, digits = 4)
	}
	return(prtable)
}


.distr1 = function(LB, UB, n)
{
	runif(n, min = LB, max = UB)
}

.distr2 = function(LB, UB, n, mean, sd)
{
	rtruncnorm(n, a = as.double(LB), b = as.double(UB), mean = as.double(mean), sd = as.double(sd))
}

.distr3 = function(n, mean, sd)
{
	rnorm(n, mean = mean, sd = sd)
}

.safefun = function(pars, fun, gosolnp_parnames, ...){
	# gosolnp_parnames = get("gosolnp_parnames", envir = .env)
	names(pars) = gosolnp_parnames
	v  = fun(pars, ...)
	if(is.na(v) | !is.finite(v) | is.nan(v)) {
		warning(paste("\ngosolnp-->warning: ", v , " detected in function call...check your function\n", sep = ""), immediate. = FALSE)
		v = 1e24
	}
	v
}
