#################################################################################
##
##   R package rmgarch by Alexios Ghalanos Copyright (C) 2008-2013.
##   This file is part of the R package rmgarch.
##
##   The R package rmgarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rmgarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################
.dcchessian = function(f, pars, arglist, fname)
{
	.eps = .Machine$double.eps
	cluster = arglist$cluster
	dccN = arglist$dccN
	arglist$returnType = "llh"
	fx = f(pars, arglist)
	n = length(pars)
	# Compute the stepsize (h)
	h = .eps ^ (1/3) * pmax( abs( pars ), 1 )
	xh = pars + h
	h  = xh - pars
	ee = as.matrix( diag( h ) )
	# Compute forward and backward steps
	g = vector(mode = "numeric", length = n)
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, require(rmgarch))
		clusterExport(cluster, c("pars", "ee", "arglist", "n", "fname"), envir = environment())
		tmp = parLapply(cluster, as.list(1:n), fun = function(i){
					tmpg =  eval(parse(text = paste(fname, "( pars = pars + ee[, i], arglist)", sep = "")))
					return( tmpg )
				})
		g = as.numeric( unlist(tmp) )
		H = h %*% t( h )
		clusterExport(cluster, c("H", "dccN", "g", "fx"), envir = environment())
		tmp = parLapply(cluster, as.list(1:n), fun = function(i){
					Htmp = H
					for(j in  (n - dccN + 1):n){
						if(i <= j){
							Htmp[i, j] = eval(parse(text = paste("(",fname, "( pars = pars + ee[, i] + ee[, j], arglist) - g[i] - g[j] + fx) / Htmp[i, j]", sep = "")))
							Htmp[j, i] = Htmp[i, j]
						}
					}
					return(Htmp)
				})
		
	} else{
		tmp = lapply(as.list(1:n), FUN = function(i){
					if(arglist$verbose) cat(paste("Evaluating StepValue ",i," out of ",n,"\n",sep=""))
					tmpg =  f( pars = pars + ee[, i], arglist)
					return( tmpg )
				})
		g = as.numeric( unlist(tmp) )
		H = h %*% t( h )
		tmp = lapply(as.list(1:n), FUN = function(i){
					Htmp = H
					for(j in  (n - dccN + 1):n){
						if(i <= j){
							Htmp[i, j] = (f( pars = pars + ee[, i] + ee[, j], arglist) - g[i] - g[j] + fx) / Htmp[i, j]
							Htmp[j, i] = Htmp[i, j]
						}
					}
					return(Htmp)
				})
	}
	for(i in 1:n){
		for(j in  (n - dccN + 1):n){
			if(i <= j){
				H[i, j] = tmp[[i]][i, j]
				H[j, i] = tmp[[i]][j, i]
			}
		}
	}	
	newH = H[(n - dccN + 1):n, ]
	H = newH
	return(H)
}


# 2-stage partitioned standard errors for DCC type models

.dccmakefitmodel = function(garchmodel, f, arglist, timer, message, fname)
{
	.eps = .Machine$double.eps
	mpars = arglist$mpars
	data = arglist$data
	cluster = arglist$cluster
	eval.se = arglist$eval.se
	fitlist = arglist$fitlist
	m = arglist$m
	midx = arglist$midx
	eidx = arglist$eidx
	dccN = arglist$dccN
	ipars = arglist$ipars
	estidx = arglist$estidx
	cnames = arglist$cnames
	model = arglist$model
	resids = residuals(fitlist)
	sigmas = sigma(fitlist)
	pars = mpars[which(eidx==1, arr.ind = TRUE)]
	arglist$returnType = "all"
	sol = f(pars, arglist)
	likelihoods 	= sol$lik
	loglikelihood 	= sol$llh
	Rtout = sol$Rt
	Qtout = sol$Qt
	N = dim(resids)[1]
	np = length(pars)
	Ht = array( 0, dim = c(m, m, N) )
	stdresid = matrix(0, nrow = N, ncol = m)
	
	if( !is.null(cluster) ){
		clusterExport(cluster, c("sigmas", "Rtout", "resids"), envir = environment())
		tmp = parLapply(cluster, as.list(1:N), fun = function(i){
					tmph = diag( sigmas[i, ] ) %*% Rtout[[i]] %*% diag( sigmas[i, ] )
					zz = eigen( tmph )
					sqrtzz = ( zz$vectors %*% diag( sqrt( zz$values ) ) %*% solve( zz$vectors ) )
					tmpz = as.numeric( resids[i, ] %*% solve( sqrtzz ) )
					return( list( H = tmph, Z = tmpz ) )
				})
		for(i in 1:N){
			Ht[,,i] = tmp[[i]]$H
			stdresid[i,] = tmp[[i]]$Z
		}
	} else{
		tmp = lapply(as.list(1:N), FUN = function(i){
					tmph = diag( sigmas[i, ] ) %*% Rtout[[i]] %*% diag( sigmas[i, ] )
					zz = eigen( tmph )
					sqrtzz = ( zz$vectors %*% diag( sqrt( zz$values ) ) %*% solve( zz$vectors ) )
					tmpz = as.numeric( resids[i, ] %*% solve( sqrtzz ) )
					return( list( H = tmph, Z = tmpz ) )
				})
		for(i in 1:N){
			Ht[,,i] = tmp[[i]]$H
			stdresid[i,] = tmp[[i]]$Z
		}
	}
	
	arglist$stdresid = stdresid
	arglist$Ht = Ht
	
	if(eval.se){
		A = zeros( np, np )
		tidx = 1
		for(i in 1:m){
			cvar = fitlist@fit[[i]]@fit$cvar
			workingsize = dim(cvar)[1]
			A[(tidx:(tidx + workingsize - 1)), (tidx:(tidx + workingsize - 1))] = solve(cvar)
			tidx = tidx + workingsize
		}
		if(arglist$verbose) cat("\n\nCalculating Standard Errors, this can take a while\n")
		otherA = .dcchessian(f = f, pars = pars, arglist, fname)
		A[(np - dccN + 1):np, ] = otherA
		jointscores = zeros(N, np)
		tidx = 1
		for(i in 1:m){
			cf = fitlist@fit[[i]]@model$pars[fitlist@fit[[i]]@model$pars[,4]==1,1]
			workingsize = length(cf)
			# head(fitlist@fit[[i]]@fit$scores, 22)
			scx = fitlist@fit[[i]]@fit$scores
			jointscores[,(tidx:(tidx + workingsize - 1))] = scx
			tidx = tidx + workingsize
		}
		h = pmax( abs( ipars[estidx,1]/2 ), 1e-2 ) * .eps^(1/3)
		hplus = ipars[estidx,1] + h
		hminus = ipars[estidx,1] - h
		likelihoodsplus =  zeros( N, dccN )
		likelihoodsminus = zeros( N, dccN )
		zparsplus = zparsminus = pars
		arglist$returnType = "lik"
		for(i in 1:dccN){
			hparameters1 = ipars[estidx,1]
			hparameters2 = ipars[estidx,1]
			hparameters1[i] = hplus[i]
			hparameters2[i] = hminus[i]
			# recombine
			zparsplus[(np-dccN+1):np] = hparameters1
			zparsminus[(np-dccN+1):np] = hparameters2
			LHT1 = f(pars = zparsplus, arglist)
			LHT2 = f(pars = zparsminus, arglist)
			likelihoodsplus[, i]  = LHT1
			likelihoodsminus[, i] = LHT2
		}
		
		sctemp = likelihoodsplus - likelihoodsminus
		DCCscores = matrix(NA, ncol = dim(sctemp)[2], nrow = dim(sctemp)[1])
		sdtemp = 2 * repmat( t( h ), N, 1 )
		for(i in 1:dim(sctemp)[2]){
			DCCscores[,i] = sctemp[,i] / sdtemp[,i]
		}
		jointscores[, (np-dccN+1):np] = DCCscores
		B = cov( jointscores )
		A = A/ (N) 
		dcccvar = ( solve( A ) %*% B %*% solve( A ) ) / N
		
		se.coef = sqrt(diag(abs(dcccvar)))
		tval = as.numeric( pars/se.coef )
		pval = 2* ( 1 - pnorm( abs( tval ) ) )
		matcoef = matrix(NA, nrow = length(pars), ncol = 4)
		matcoef[, 1] = pars
		matcoef[, 2] = se.coef
		matcoef[, 3] = tval
		matcoef[, 4] = pval
		allnames = NULL
		for(i in 1:m){
			allnames = c(allnames, paste("[",cnames[i],"].", rownames(eidx[eidx[,i]==1,i, drop = FALSE]), sep = ""))
		}
		garchnames = allnames
		dccnames = rownames(eidx[eidx[,m+1]==1,m+1, drop = FALSE])
		allnames = c(garchnames, paste("[Joint]", rownames(eidx[eidx[,m+1]==1,m+1, drop = FALSE]), sep = ""))
		dimnames(matcoef) = list(allnames, c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
		
	} else{

		se.coef = rep(NA, length(pars))
		tval = rep(NA, length(pars))
		pval = rep(NA, length(pars))
		matcoef = matrix(NA, nrow = length(pars), ncol = 4)
		matcoef[, 1] = pars
		allnames = NULL
		for(i in 1:m){
			allnames = c(allnames, paste("[",cnames[i],"].", rownames(eidx[eidx[,i]==1,i, drop = FALSE]), sep = ""))
		}
		garchnames = allnames
		dccnames = rownames(eidx[eidx[,m+1]==1,m+1, drop = FALSE])
		allnames = c(garchnames, paste("[Joint]", rownames(eidx[eidx[,m+1]==1,m+1, drop = FALSE]), sep = ""))
		dimnames(matcoef) = list(allnames, c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
		dcccvar = NULL
		jointscores = NULL
	}
	
	dccfit = list()
	dccfit$coef = pars
	names(dccfit$coef) = allnames
	dccfit$matcoef = matcoef
	dccfit$garchnames = garchnames
	dccfit$dccnames = dccnames
	dccfit$cvar = dcccvar
	dccfit$scores = jointscores
	dccfit$R = Rtout
	dccfit$H = Ht
	dccfit$Q = Qtout
	dccfit$stdresid = stdresid
	dccfit$llh = loglikelihood
	dccfit$log.likelihoods = likelihoods
	dccfit$timer = timer
	dccfit$convergence = 0
	dccfit$message = message
	return( dccfit )
}


.dccmakefiltermodel = function(garchmodel, f, arglist, timer, message, fname)
{
	.eps = .Machine$double.eps
	mpars = arglist$mpars
	data = arglist$data
	cluster = arglist$cluster
	m = arglist$m
	midx = arglist$midx
	eidx = arglist$eidx
	dccN = arglist$dccN
	ipars = arglist$ipars
	estidx = arglist$estidx
	cnames = arglist$cnames
	model = arglist$model
	
	filterlist = arglist$filterlist
	resids = residuals(filterlist)
	sigmas = sigma(filterlist)
	pars = mpars[which(midx==1, arr.ind = TRUE)]
	arglist$returnType = "all"
	sol = f(pars, arglist)
	likelihoods 	= sol$lik
	loglikelihood 	= sol$llh
	Rtout = sol$Rt
	Qtout = sol$Qt
	N = dim(resids)[1]
	np = length(pars)
	Ht = array( 0, dim = c(m, m, N) )
	stdresid = matrix(0, nrow = N, ncol = m)
	
	if( !is.null(cluster) ){
		clusterExport(cluster, c("sigmas", "Rtout", "resids"), envir = environment())
		tmp = parLapply(cluster, as.list(1:N), fun = function(i){
					tmph = diag( sigmas[i, ] ) %*% Rtout[[i]] %*% diag( sigmas[i, ] )
					zz = eigen( tmph )
					sqrtzz = ( zz$vectors %*% diag( sqrt( zz$values ) ) %*% solve( zz$vectors ) )
					tmpz = as.numeric( resids[i, ] %*% solve( sqrtzz ) )
					return( list( H = tmph, Z = tmpz ) )
				})
		for(i in 1:N){
			Ht[,,i] = tmp[[i]]$H
			stdresid[i,] = tmp[[i]]$Z
		}
	} else{
		tmp = lapply(as.list(1:N), FUN = function(i){
					tmph = diag( sigmas[i, ] ) %*% Rtout[[i]] %*% diag( sigmas[i, ] )
					zz = eigen( tmph )
					sqrtzz = ( zz$vectors %*% diag( sqrt( zz$values ) ) %*% solve( zz$vectors ) )
					tmpz = as.numeric( resids[i, ] %*% solve( sqrtzz ) )
					return( list( H = tmph, Z = tmpz ) )
				})
		for(i in 1:N){
			Ht[,,i] = tmp[[i]]$H
			stdresid[i,] = tmp[[i]]$Z
		}
	}
	
	arglist$stdresid = stdresid
	arglist$Ht = Ht
		
	se.coef = rep(NA, length(pars))
	tval = rep(NA, length(pars))
	pval = rep(NA, length(pars))
	matcoef = matrix(NA, nrow = length(pars), ncol = 4)
	matcoef[, 1] = pars
	allnames = NULL
	for(i in 1:m){
		allnames = c(allnames, paste("[",cnames[i],"].", rownames(eidx[eidx[,i]==1,i, drop = FALSE]), sep = ""))
	}
	garchnames = allnames
	dccnames = rownames(eidx[eidx[,m+1]==1,m+1, drop = FALSE])
	allnames = c(garchnames, paste("[Joint]", rownames(eidx[eidx[,m+1]==1,m+1, drop = FALSE]), sep = ""))
	dimnames(matcoef) = list(allnames, c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
	dcccvar = NULL
	jointscores = NULL
	
	dccfilter = list()
	dccfilter$coef = pars
	names(dccfilter$coef) = allnames
	dccfilter$garchnames = garchnames
	dccfilter$dccnames = dccnames
	dccfilter$dcccoef = pars
	dccfilter$matcoef = matcoef
	dccfilter$cvar = dcccvar
	dccfilter$scores = jointscores
	dccfilter$R = Rtout
	dccfilter$H = Ht
	dccfilter$Q = Qtout
	dccfilter$stdresid = stdresid
	dccfilter$llh = loglikelihood
	dccfilter$log.likelihoods = likelihoods
	dccfilter$timer = timer
	dccfilter$convergence = 0
	dccfilter$message = message
	return( dccfilter )
}


.sqrtsymmat = function( X )
{
	tmp = eigen( X )
	sqrttmp = ( tmp$vectors %*% diag( sqrt( tmp$values ) ) %*% solve( tmp$vectors ) )
	return( sqrttmp )
}