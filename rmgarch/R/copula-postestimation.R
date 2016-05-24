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
.cgarchhessian = function(f, pars, arglist, fname)
{
	cluster = arglist$cluster
	dccN = arglist$dccN
	arglist$returnType = "llh"
	fx = f(pars, arglist)
	.eps = .Machine$double.eps
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
		clusterExport(cluster, c("pars", "ee", "arglist", "n", "fname"), 
				envir = environment())
		tmp = parLapply(cluster, as.list(1:n), fun = function(i){
					tmpg =  eval(parse(text = paste(fname, "( pars = pars + ee[, i], arglist)", sep = "")))
					return( tmpg )
				})
			# sfStop() we'll stop it further down
		g = as.numeric( unlist(tmp) )
	} else{
		tmp = lapply(as.list(1:n), FUN = function(i){
					if(arglist$verbose) cat(paste("Evaluating StepValue ",i," out of ",n,"\n",sep=""))
					tmpg =  f( pars = pars + ee[, i], arglist)
					return( tmpg )
				})
		g = as.numeric( unlist(tmp) )
	}
	
	H = h %*% t( h )
	
	if( !is.null(cluster)  ){
		clusterEvalQ(cluster, require(rmgarch))
		clusterExport(cluster, c("pars", "H", "ee", "n", "dccN", 
						"g", "fx", "fname", "arglist"), envir = environment())
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
		for(i in 1:n){
			for(j in  (n - dccN + 1):n){
				if(i <= j){
					H[i, j] = tmp[[i]][i, j]
					H[j, i] = tmp[[i]][j, i]
				}
			}
		}
	} else{
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
		for(i in 1:n){
			for(j in  (n - dccN + 1):n){
				if(i <= j){
					H[i, j] = tmp[[i]][i, j]
					H[j, i] = H[i, j]
				}
			}
		}
	}
	newH = H[(n - dccN + 1):n, ]
	H = newH
	return(H)
}


.cgarchmakefitmodel1 = function(f, arglist, timer, message, fname)
{
	.eps = .Machine$double.eps
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
	mpars = arglist$mpars
	model = arglist$model
	maxgarchOrder = model$maxgarchOrder
	resids = residuals(fitlist)
	sigmas = sigma(fitlist)
	#if(maxgarchOrder>0){
	#	resids = resids[-c(1:maxgarchOrder), , drop = FALSE]
	#	sigmas = sigmas[-c(1:maxgarchOrder), , drop = FALSE]
	#}
	
	pars = mpars[which(eidx==1, arr.ind = TRUE)]
	arglist$returnType = "ALL"
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
		otherA = .cgarchhessian(f = f, pars = pars, arglist, fname)
		A[(np - dccN + 1):np, ] = otherA
		jointscores = zeros(N, np)
		tidx = 1
		for(i in 1:m){
			cf = fitlist@fit[[i]]@model$pars[fitlist@fit[[i]]@model$pars[,4]==1,1]
			workingsize = length(cf)
			# head(fitlist@fit[[i]]@fit$scores, 22)
			scx = fitlist@fit[[i]]@fit$scores
			#if(maxgarchOrder>0){
			#	scx = scx[-c(1:maxgarchOrder), , drop = FALSE]
			#}
			jointscores[,(tidx:(tidx + workingsize - 1))] = scx
			tidx = tidx + workingsize
		}
		## CONTINUE here
		#Now all we need to do is calculate the scores form the dcc estimator and we have everything
		h = pmax( abs( ipars[estidx,1]/2 ), 1e-2 ) * .eps^(1/3)
		hplus = ipars[estidx,1] + h
		hminus = ipars[estidx,1] - h
		likelihoodsplus =  zeros( N, dccN )
		likelihoodsminus = zeros( N, dccN )
		zparsplus = zparsminus = pars
		for(i in 1:dccN){
			hparameters1 = ipars[estidx,1]
			hparameters2 = ipars[estidx,1]
			hparameters1[i] = hplus[i]
			hparameters2[i] = hminus[i]
			# recombine
			zparsplus[(np-dccN+1):np] = hparameters1
			zparsminus[(np-dccN+1):np] = hparameters2
			arglist$returnType = "lik"
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
		if(!is.null(dccnames)){
			allnames = c(garchnames, paste("[Joint]", rownames(eidx[eidx[,m+1]==1,m+1, drop = FALSE]), sep = ""))
		} else{
			allnames = garchnames
		}
		allnames = c(garchnames, paste("[Joint]", rownames(eidx[eidx[,m+1]==1,m+1, drop = FALSE]), sep = ""))
		dimnames(matcoef) = list(allnames, c(" Estimate",
						" Std. Error", " t value", "Pr(>|t|)"))
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
		if(!is.null(dccnames)){
			allnames = c(garchnames, paste("[Joint]", rownames(eidx[eidx[,m+1]==1,m+1, drop = FALSE]), sep = ""))
		} else{
			allnames = garchnames
		}		
		allnames = c(garchnames, paste("[Joint]", rownames(eidx[eidx[,m+1]==1,m+1, drop = FALSE]), sep = ""))
		dimnames(matcoef) = list(allnames, c(" Estimate",
						" Std. Error", " t value", "Pr(>|t|)"))
		dcccvar = NULL
		jointscores = NULL
	}
	
	cfit = list()
	cfit$coef = pars
	names(cfit$coef) = allnames
	
	cfit$matcoef = matcoef
	cfit$garchnames = garchnames
	cfit$dccnames = dccnames
	cfit$cvar = dcccvar
	cfit$scores = jointscores
	cfit$H = Ht
	cfit$stdresid = stdresid
	cfit$timer = timer
	cfit$convergence = 0
	cfit$message = message
	return( cfit )
}

.cgarchmakefitmodel2 = function(f, arglist, timer, message, fname)
{
	.eps = .Machine$double.eps
	cluster = arglist$cluster
	eval.se = arglist$eval.se
	fitlist = arglist$fitlist
	m = arglist$m
	midx = arglist$midx
	eidx = arglist$eidx
	ipars = arglist$ipars
	estidx = arglist$estidx
	cnames = arglist$cnames
	mpars = arglist$mpars
	model = arglist$model
	method = model$modeldesc$cor.method
	maxgarchOrder = model$maxgarchOrder
	resids = residuals(fitlist)
	sigmas = sigma(fitlist)
	#if(maxgarchOrder>0){
	#	resids = resids[-c(1:maxgarchOrder), , drop = FALSE]
	#	sigmas = sigmas[-c(1:maxgarchOrder), , drop = FALSE]
	#}
	
	pars = mpars[which(eidx==1, arr.ind = TRUE)]
	arglist$returnType = "ALL"
	sol = f(pars, arglist)
	likelihoods 	= sol$lik
	loglikelihood 	= sol$llh
	Rt = sol$R
	N = dim(resids)[1]
	np = length(pars)
	Ht = array( 0, dim = c(m, m, N) )
	stdresid = matrix(0, nrow = N, ncol = m)
	
	if( !is.null(cluster) ){
		clusterExport(cluster, c("sigmas", "Rt", "resids"), envir = environment())
		tmp = parLapply(cluster, as.list(1:N), fun = function(i){
						tmph = diag( sigmas[i, ] ) %*% Rt %*% diag( sigmas[i, ] )
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
					tmph = diag( sigmas[i, ] ) %*% Rt %*% diag( sigmas[i, ] )
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
	#assign("stdresid", stdresid, envir = garchenv)
	#assign("Ht", Ht, envir = garchenv)
	if(eval.se){
		if(method == "ML" || model$modeldesc$distribution == "mvt"){
			dccN = arglist$npars
			arglist$dccN = dccN
			A = zeros( np, np )
			tidx = 1
			for(i in 1:m){
				cvar = fitlist@fit[[i]]@fit$cvar
				workingsize = dim(cvar)[1]
				A[(tidx:(tidx + workingsize - 1)), (tidx:(tidx + workingsize - 1))] = solve(cvar)
				tidx = tidx + workingsize
			}
			if(arglist$verbose) cat("\n\nCalculating Standard Errors, this can take a while\n")
			if(dccN>0){
				otherA = .cgarchhessian(f = f, pars = pars, arglist, fname)
				A[(np - dccN + 1):np, ] = otherA
				jointscores = zeros(N, np)
				tidx = 1
				for(i in 1:m){
					cf = fitlist@fit[[i]]@model$pars[fitlist@fit[[i]]@model$pars[,4]==1,1]
					workingsize = length(cf)
					# head(fitlist@fit[[i]]@fit$scores, 22)
					scx = fitlist@fit[[i]]@fit$scores
					#if(maxgarchOrder>0){
					#	scx = scx[-c(1:maxgarchOrder), , drop = FALSE]
					#}
					jointscores[,(tidx:(tidx + workingsize - 1))] = scx
					tidx = tidx + workingsize
				}
				h = pmax( abs( ipars[estidx,1]/2 ), 1e-2 ) * .eps^(1/3)
				hplus = ipars[estidx,1] + h
				hminus = ipars[estidx,1] - h
				likelihoodsplus =  zeros( N, dccN )
				likelihoodsminus = zeros( N, dccN )
				zparsplus = zparsminus = pars
				for(i in 1:dccN){
					hparameters1 = ipars[estidx,1]
					hparameters2 = ipars[estidx,1]
					hparameters1[i] = hplus[i]
					hparameters2[i] = hminus[i]
					# recombine
					zparsplus[(np-dccN+1):np] = hparameters1
					zparsminus[(np-dccN+1):np] = hparameters2
					arglist$returnType = "lik"
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
				if(!is.null(dccnames)){
					allnames = c(garchnames, paste("[Joint]", rownames(eidx[eidx[,m+1]==1,m+1, drop = FALSE]), sep = ""))
				} else{
					allnames = garchnames
				}
				dimnames(matcoef) = list(allnames, c(" Estimate"," Std. Error", " t value", "Pr(>|t|)"))
			} else{
				jointscores = zeros(N, np)
				tidx = 1
				for(i in 1:m){
					cf = fitlist@fit[[i]]@model$pars[fitlist@fit[[i]]@model$pars[,4]==1,1]
					workingsize = length(cf)
					# head(fitlist@fit[[i]]@fit$scores, 22)
					scx = fitlist@fit[[i]]@fit$scores
					#if(maxgarchOrder>0){
					#	scx = scx[-c(1:maxgarchOrder), , drop = FALSE]
					#}
					jointscores[,(tidx:(tidx + workingsize - 1))] = scx
					tidx = tidx + workingsize
				}
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
				if(!is.null(dccnames)){
					allnames = c(garchnames, paste("[Joint]", rownames(eidx[eidx[,m+1]==1,m+1, drop = FALSE]), sep = ""))
				} else{
					allnames = garchnames
				}
				dimnames(matcoef) = list(allnames, c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
			}
		} else{
			A = zeros( np, np )
			tidx = 1
			for(i in 1:m){
				cvar = fitlist@fit[[i]]@fit$cvar
				workingsize = dim(cvar)[1]
				A[(tidx:(tidx + workingsize - 1)), (tidx:(tidx + workingsize - 1))] = solve(cvar)
				tidx = tidx + workingsize
			}
			
			#cat("\n\nCalculating Standard Errors, this can take a while\n")
			#otherA = .cgarchhessian(f = f, pars = pars, garchenv = garchenv, fname)
			#A[(np - dccN + 1):np, ] = otherA
			jointscores = zeros(N, np)
			tidx = 1
			for(i in 1:m){
				cf = fitlist@fit[[i]]@model$pars[fitlist@fit[[i]]@model$pars[,4]==1,1]
				workingsize = length(cf)
				# head(fitlist@fit[[i]]@fit$scores, 22)
				scx = fitlist@fit[[i]]@fit$scores
				#if(maxgarchOrder>0){
				#	scx = scx[-c(1:maxgarchOrder), , drop = FALSE]
				#}
				jointscores[,(tidx:(tidx + workingsize - 1))] = scx
				tidx = tidx + workingsize
			}
			## CONTINUE here
			#Now all we need to do is calculate the scores form the dcc estimator and we have everything
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
			dccnames = NULL
			dimnames(matcoef) = list(allnames, c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
		}
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
		if(!is.null(dccnames)){
			allnames = c(garchnames, paste("[Joint]", rownames(eidx[eidx[,m+1]==1,m+1, drop = FALSE]), sep = ""))
		} else{
			allnames = garchnames
		}
		dimnames(matcoef) = list(allnames, c(" Estimate", " Std. Error", " t value", "Pr(>|t|)"))
		dcccvar = NULL
		jointscores = NULL
	}
	cgarchfit = list()
	cgarchfit$coef = pars
	names(cgarchfit$coef) = allnames
	cgarchfit$matcoef = matcoef
	cgarchfit$garchnames = garchnames
	cgarchfit$dccnames = dccnames
	cgarchfit$cvar = dcccvar
	cgarchfit$scores = jointscores
	cgarchfit$H = Ht
	cgarchfit$stdresid = stdresid
	cgarchfit$timer = timer
	cgarchfit$convergence = 0
	cgarchfit$message = message
	return( cgarchfit )
}