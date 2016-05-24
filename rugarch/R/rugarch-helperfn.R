#################################################################################
##
##   R package rugarch by Alexios Ghalanos Copyright (C) 2008-2015.
##   This file is part of the R package rugarch.
##
##   The R package rugarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rugarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################
# CONSTANTS used in package:
eps<-.Machine$double.eps
TinY = 1.0e-8

##################################################################################
# Helper Functions
.makearfimafitmodel = function(f, T, m, timer, convergence, message, hess, arglist, 
		numderiv.control = list(grad.eps=1e-4, grad.d=0.0001, grad.zero.tol=sqrt(.Machine$double.eps/7e-7),
				hess.eps=1e-4, hess.d=0.1, hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2))
{
	# Need to turn off stationarity check:
	fit.control = arglist$fit.control
	fit.control$stationarity = arglist$fit.control$stationarity = 0
	ipars = arglist$ipars
	model = arglist$model
	estidx = arglist$estidx
	idx = model$pidx
	data = arglist$data
	arglist$returnType = "llh"
	fit = list()
	if(is.null(hess)){
		fit$hessian = hessian(func = f, x = ipars[estidx, 1], method="Richardson", method.args=list(eps = numderiv.control$hess.eps, 
						d= numderiv.control$hess.d, zero.tol=numderiv.control$hess.zero.tol, r=numderiv.control$r, v=numderiv.control$v, 
						show.details=FALSE), arglist=arglist)
		#fit$hessian = .hessian2sided(f, ipars[estidx, 1], arglist = arglist)
	} else{
		fit$hessian = hess
	}
	fit$cvar = try(solve(fit$hessian), silent = TRUE)
	# (might also fail in which case user will see an error about the inversion failure)
	if(inherits(fit$cvar, "try-error")){
		zz = try(solve(.hessian2sided(f, ipars[estidx, 1], arglist = arglist)), silent=TRUE)
		if(inherits(zz, "try-error")) {
			fit$cvar = NULL
			warning("\nrugarch-->warning: failed to invert hessian\n")
		} else{
			fit$cvar = zz
		}
	}
	arglist$returnType = "all"
	temp = f(pars = ipars[estidx, 1], arglist = arglist)
	fit$z = temp$z
	fit$LLH = -temp$llh
	fit$log.likelihoods = temp$LHT
	fit$residuals = temp$res
	if(sum(ipars[,2])>0){
		pall = ipars[estidx | as.logical(ipars[,2]==1), 1]
		fixed = match(rownames(ipars[ipars[,2]==1, , drop = FALSE]), names(pall))
		fixedn = length(fixed)
		fNA = rep(NA, fixedn)
		nfixedn = length(pall) - fixedn
		fit$coef = pall
		if(is.null(fit$cvar)){
			fit$se.coef = rep(NA, nfixedn)
			fit$tval = rep(NA, nfixedn)
			# change here
			fit$matcoef = matrix(NA, ncol = 4, nrow = length(pall))
			fit$matcoef[-fixed,] = cbind(ipars[estidx, 1], fit$se.coef,
					fit$tval, rep(NA, nfixedn))
			fit$matcoef[fixed,] = cbind(fit$coef[fixed], fNA, fNA, fNA)
			fit$robust.se.coef = rep(NA, nfixedn)
			fit$robust.tval = rep(NA, nfixedn)
			fit$robust.matcoef = matrix(NA, ncol = 4, nrow = length(fit$coef))
			fit$robust.matcoef[-fixed,] = cbind(fit$coef[-fixed], fit$robust.se.coef, 
					fit$robust.tval, rep(NA, nfixedn))
			fit$robust.matcoef[fixed,] = cbind(fit$coef[fixed], fNA, fNA, fNA)
			fit$hessian.message = "failed to invert hessian"
		} else{
			arglist$returnType="LHT"
			tmp = robustvcv(fun = f, pars = ipars[estidx, 1], nlag = 0, hess = fit$hessian, n = T, arglist = arglist)
			fit$robust.cvar = tmp$vcv
			fit$scores = jacobian(func = f, x = ipars[estidx, 1], method="Richardson", method.args=list(eps = numderiv.control$grad.eps, 
							d = numderiv.control$grad.d, zero.tol=numderiv.control$grad.zero.tol, r=numderiv.control$r, v=numderiv.control$v, 
							show.details=FALSE), arglist=arglist)
			colnames(fit$scores) = names(ipars[estidx, 1])
			fit$se.coef = sqrt(diag(abs(fit$cvar)))
			fit$tval = fit$coef[-fixed]/fit$se.coef
			# change here
			fit$matcoef = matrix(NA, ncol = 4, nrow = length(fit$coef))
			fit$matcoef[-fixed,] = cbind(fit$coef[-fixed], fit$se.coef,
					fit$tval, 2*(1-pnorm(abs(fit$tval))))
			fit$matcoef[fixed,] = cbind(fit$coef[fixed], fNA,fNA, fNA)
			fit$robust.se.coef = sqrt(diag(fit$robust.cvar))
			fit$robust.tval = fit$coef[-fixed]/fit$robust.se.coef
			# change here
			fit$robust.matcoef = matrix(NA, ncol = 4, nrow = length(fit$coef))
			fit$robust.matcoef[-fixed,] = cbind(fit$coef[-fixed], fit$robust.se.coef,
					fit$robust.tval, 2*(1-pnorm(abs(fit$robust.tval))))
			fit$robust.matcoef[fixed,] = cbind(fit$coef[fixed], fNA, fNA, fNA)
			fit$hessian.message = NULL
		}
	} else{
		fit$coef = ipars[estidx, 1]
		if(is.null(fit$cvar)){
			fit$se.coef = rep(NA,length(fit$coef))
			fit$tval = rep(NA,length(fit$coef))
			# change here
			fit$matcoef = cbind(fit$coef, fit$se.coef,
					fit$tval, rep(NA,length(fit$coef)))
			fit$robust.se.coef = rep(NA,length(fit$coef))
			fit$robust.tval = rep(NA,length(fit$coef))
			# change here
			fit$robust.matcoef = cbind(fit$coef, fit$robust.se.coef,fit$robust.tval, 
					rep(NA,length(fit$coef)))
			fit$hessian.message = "failed to invert hessian"
		} else{
			nlag=min(floor(1.2*(T)^(1/3)),(T))
			# tmp = .robustSE(fun = f, pars = opt, hess = NULL, n = T-m, data = data, returnType = "LHT", garchenv)
			arglist$returnType = "LHT"
			tmp = robustvcv(fun = f, pars = ipars[estidx,1], nlag = nlag, hess = fit$hessian, n = T, arglist = arglist)
			fit$robust.cvar = tmp$vcv
			fit$scores = jacobian(func = f, x = ipars[estidx, 1], method="Richardson", method.args=list(eps = numderiv.control$grad.eps, 
							d = numderiv.control$grad.d, zero.tol=numderiv.control$grad.zero.tol, r=numderiv.control$r, v=numderiv.control$v, 
							show.details=FALSE), arglist=arglist)
			colnames(fit$scores) = names(ipars[estidx, 1])
			fit$se.coef = sqrt(diag(abs(fit$cvar)))
			fit$tval = fit$coef/fit$se.coef
			# change here
			fit$matcoef = cbind(fit$coef, fit$se.coef,
					fit$tval, 2*(1-pnorm(abs(fit$tval))))
			fit$robust.se.coef = sqrt(diag(fit$robust.cvar))
			fit$robust.tval = fit$coef/fit$robust.se.coef
			# change here
			fit$robust.matcoef = cbind(fit$coef, fit$robust.se.coef,fit$robust.tval, 
					2*(1-pnorm(abs(fit$robust.tval))))
			fit$hessian.message = NULL
		}
	}
	# return the correct indicators to the ipars (changed in the case of fixed pars and fixed.se = TRUE)
	ipars[,2:4] = model$pars[,2:4]
	dimnames(fit$matcoef) = list(names(fit$coef), c(" Estimate",
					" Std. Error", " t value", "Pr(>|t|)"))
	dimnames(fit$robust.matcoef) = list(names(fit$coef), c(" Estimate",
					" Std. Error", " t value", "Pr(>|t|)"))
	fit$fitted.values = data-fit$residuals
	fit$convergence = convergence
	fit$message = message
	fit$kappa = temp$kappa
	fit$persistence = temp$persistence
	fit$timer = timer
	fit$ipars = ipars
	return(fit)
}


.makefitmodel = function(garchmodel, f, T, m, timer, convergence, message, hess, arglist, 
		numderiv.control = list(grad.eps=1e-4, grad.d=0.0001, grad.zero.tol=sqrt(.Machine$double.eps/7e-7),
				hess.eps=1e-4, hess.d=0.1, hess.zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2))
{
	# Turn Stationarity Off for numerical derivative calculation
	fit.control = arglist$fit.control
	fit.control$stationarity = arglist$fit.control$stationarity = 0
	data = arglist$data
	ipars = arglist$ipars
	model = arglist$model
	estidx = arglist$estidx
	idx = model$pidx
	arglist$returnType = "llh"
	fit = vector(mode = "list")
	if(is.null(hess)){
		fit$hessian = hessian(func = f, x = ipars[estidx, 1], method="Richardson", method.args=list(eps = numderiv.control$hess.eps, 
						d= numderiv.control$hess.d, zero.tol=numderiv.control$hess.zero.tol, r=numderiv.control$r, v=numderiv.control$v, 
						show.details=FALSE), arglist=arglist)		
		# fit$hessian = .hessian2sided(f, ipars[estidx, 1], data = data, returnType = "llh", garchenv = garchenv)
		# fit$hessian = .hessian2sidedcpp(f, ipars[estidx, 1], arglist = arglist)
		E = eigen(fit$hessian)$values
		# approx. number of decimal places lost to roundoff/numerical estimation error
		if(any(E<0)){
			# This should silence the warnings
			condH = NaN
		} else{
			condH = log10(max(E)/min(E))			
		}
	} else{
		fit$hessian = hess
		E = eigen(fit$hessian)$values
		if(any(E<0)){
			# This should silence the warnings
			condH = NaN
		} else{
			condH = log10(max(E)/min(E))			
		}
	}
	fit$cvar = try(solve(fit$hessian), silent = TRUE)
	# (might also fail in which case user will see an error about the inversion failure)
	if(inherits(fit$cvar, "try-error")){
		zz = try(solve(.hessian2sided(f, ipars[estidx, 1], arglist = arglist)), silent=TRUE)
		if(inherits(zz, "try-error")) {
			fit$cvar = NULL
			warning("\nrugarch-->warning: failed to invert hessian\n")
		} else{
			fit$cvar = zz
		}
	}
	#fit$grad = grad(func = f, x = ipars[estidx, 1], method="Richardson", method.args=list(eps=1e-4, d=0.001, 
	#				zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE), data = data, returnType = "llh", garchenv = garchenv)
	# this is the sum of the scores: apply(fit$scores, 2, "sum")
	arglist$returnType = "all"
	temp = f(pars = ipars[estidx, 1],	arglist = arglist)
	if(garchmodel == "fGARCH" || garchmodel == "apARCH"){
		# apARCH and fGARCH are calculated in sigma rather than sigma^d (since d is a variable
		# to be calculated, unlike sGARCH where it is fixed at 2)
		fit$var = temp$h^2
		fit$sigma = temp$h
	} else{
		fit$var = abs(temp$h)
		fit$sigma = sqrt(abs(temp$h))
	}
	if(garchmodel == "csGARCH") fit$q = temp$q
	if(garchmodel == "realGARCH"){
		fit$u = temp$u
		fit$tau = temp$tau
		fit$partial.log.likelihoods = temp$LHT1P
	}
	fit$condH = condH
	fit$z = temp$z
	fit$LLH = -temp$llh
	fit$log.likelihoods = temp$LHT
	fit$residuals = temp$epsx
	if(sum(ipars[,2])>0){
		pall = ipars[estidx | as.logical(ipars[,2]==1), 1]
		fixed = match(rownames(ipars[ipars[,2]==1, , drop = FALSE]), names(pall))
		fixedn = length(fixed)
		fNA = rep(NA, fixedn)
		nfixedn = length(pall) - fixedn
		fit$coef = pall
		if(is.null(fit$cvar)){
			fit$se.coef = rep(NA, nfixedn)
			fit$tval = rep(NA, nfixedn)
			# change here
			fit$matcoef = matrix(NA, ncol = 4, nrow = length(pall))
			fit$matcoef[-fixed,] = cbind(ipars[estidx, 1], fit$se.coef,
					fit$tval, rep(NA, nfixedn))
			fit$matcoef[fixed,] = cbind(fit$coef[fixed], fNA, fNA, fNA)
			fit$robust.se.coef = rep(NA, nfixedn)
			fit$robust.tval = rep(NA, nfixedn)
			fit$robust.matcoef = matrix(NA, ncol = 4, nrow = length(fit$coef))
			fit$robust.matcoef[-fixed,] = cbind(fit$coef[-fixed], fit$robust.se.coef, 
					fit$robust.tval, rep(NA, nfixedn))
			fit$robust.matcoef[fixed,] = cbind(fit$coef[fixed], fNA, fNA, fNA)
			fit$hessian.message = "failed to invert hessian"
		} else{
			arglist$returnType = "LHT"
			tmp = robustvcv(fun = f, pars = ipars[estidx, 1], nlag = 0, hess = fit$hessian, n = T, arglist = arglist)
			fit$robust.cvar = tmp$vcv
			fit$scores = jacobian(func = f, x = ipars[estidx, 1], method="Richardson", method.args=list(eps = numderiv.control$grad.eps, 
							d = numderiv.control$grad.d, zero.tol=numderiv.control$grad.zero.tol, r=numderiv.control$r, v=numderiv.control$v, 
							show.details=FALSE), arglist=arglist)
			colnames(fit$scores) = names(ipars[estidx, 1])
			fit$se.coef = sqrt(diag(abs(fit$cvar)))
			fit$tval = fit$coef[-fixed]/fit$se.coef
			# change here
			fit$matcoef = matrix(NA, ncol = 4, nrow = length(fit$coef))
			fit$matcoef[-fixed,] = cbind(fit$coef[-fixed], fit$se.coef,
					fit$tval, 2*(1-pnorm(abs(fit$tval))))
			fit$matcoef[fixed,] = cbind(fit$coef[fixed], fNA,fNA, fNA)
			fit$robust.se.coef = sqrt(diag(fit$robust.cvar))
			fit$robust.tval = fit$coef[-fixed]/fit$robust.se.coef
			# change here
			fit$robust.matcoef = matrix(NA, ncol = 4, nrow = length(fit$coef))
			fit$robust.matcoef[-fixed,] = cbind(fit$coef[-fixed], fit$robust.se.coef,
					fit$robust.tval, 2*(1-pnorm(abs(fit$robust.tval))))
			fit$robust.matcoef[fixed,] = cbind(fit$coef[fixed], fNA, fNA, fNA)
			fit$hessian.message = NULL
		}
		if(model$modelinc[7] == 0){
			vtomega = ipars[idx["omega", 1], 1]
			names(vtomega) = "omega"
			fit$coef = c(fit$coef, vtomega)
			names(fit$coef)[length(fit$coef)] = "omega"
			fit$se.coef = c(fit$se.coef, NA)
			fit$tval = c(fit$tval, NA)
			# change here
			fit$matcoef = rbind(fit$matcoef, c(vtomega, NA, NA, NA))
			fit$robust.se.coef = c(fit$robust.se.coef, NA)
			fit$robust.tval = c(fit$robust.tval, NA)
			# change here
			fit$robust.matcoef = rbind(fit$robust.matcoef, c(vtomega, NA, NA, NA))
		}
	} else{
		fit$coef = ipars[estidx, 1]
		if(is.null(fit$cvar)){
			fit$se.coef = rep(NA,length(fit$coef))
			fit$tval = rep(NA,length(fit$coef))
			# change here
			fit$matcoef = cbind(fit$coef, fit$se.coef,
					fit$tval, rep(NA,length(fit$coef)))
			fit$robust.se.coef = rep(NA,length(fit$coef))
			fit$robust.tval = rep(NA,length(fit$coef))
			# change here
			fit$robust.matcoef = cbind(fit$coef, fit$robust.se.coef,fit$robust.tval, 
					rep(NA,length(fit$coef)))
			fit$hessian.message = "failed to invert hessian"
		} else{
			nlag=min(floor(1.2*(T)^(1/3)),(T))
			arglist$returnType = "LHT"
			tmp = robustvcv(fun = f, pars = ipars[estidx,1], nlag = nlag, hess = fit$hessian, n = T, arglist = arglist)
			fit$robust.cvar = tmp$vcv
			fit$scores = jacobian(func = f, x = ipars[estidx, 1], method="Richardson", method.args=list(eps = numderiv.control$grad.eps, 
							d = numderiv.control$grad.d, zero.tol=numderiv.control$grad.zero.tol, r=numderiv.control$r, v=numderiv.control$v, 
							show.details=FALSE), arglist=arglist)
			colnames(fit$scores) = names(ipars[estidx, 1])
			fit$se.coef = sqrt(diag(abs(fit$cvar)))
			fit$tval = fit$coef/fit$se.coef
			# change here
			fit$matcoef = cbind(fit$coef, fit$se.coef,
					fit$tval, 2*(1-pnorm(abs(fit$tval))))
			fit$robust.se.coef = sqrt(diag(fit$robust.cvar))
			fit$robust.tval = fit$coef/fit$robust.se.coef
			# change here
			fit$robust.matcoef = cbind(fit$coef, fit$robust.se.coef,fit$robust.tval, 
					2*(1-pnorm(abs(fit$robust.tval))))
			fit$hessian.message = NULL
		}
		# variance targeting case
		if(model$modelinc[7] == 0){
			vtomega = ipars[idx["omega", 1], 1]
			names(vtomega) = "omega"
			fit$coef = c(fit$coef, vtomega)
			names(fit$coef)[length(fit$coef)] = "omega"
			fit$se.coef = c(fit$se.coef, NA)
			fit$tval = c(fit$tval, NA)
			# change here
			fit$matcoef = rbind(fit$matcoef, c(vtomega, NA, NA, NA))
			
			fit$robust.se.coef = c(fit$robust.se.coef, NA)
			fit$robust.tval = c(fit$robust.tval, NA)
			# change here
			fit$robust.matcoef = rbind(fit$robust.matcoef, c(vtomega, NA, NA, NA))
		}
	}
	# return the correct indicators to the ipars (changed in the case of fixed pars and fixed.se = TRUE)
	ipars[,2:4] = model$pars[,2:4]
	dimnames(fit$matcoef) = list(names(fit$coef), c(" Estimate",
					" Std. Error", " t value", "Pr(>|t|)"))
	dimnames(fit$robust.matcoef) = list(names(fit$coef), c(" Estimate",
					" Std. Error", " t value", "Pr(>|t|)"))
	fit$fitted.values = data-fit$residuals
	fit$convergence = convergence
	fit$message = message
	fit$kappa = temp$kappa
	fit$persistence = temp$persistence
	fit$timer = timer
	fit$ipars = ipars
	return(fit)
}

.forcregressors = function(model, mregfor, vregfor, n.ahead, N, out.sample, n.roll)
{
	# N is the original length
	treq = n.ahead + n.roll
	mxn = model$modelinc[6]
	vxn = model$modelinc[15]
	if(mxn>0){
		if(!is.null(mregfor)){
			nmex = dim(as.matrix(mregfor))[1]
			mmex = dim(as.matrix(mregfor))[2]
		} else{
			nmex = 0
			mmex = 0
		}
		if(!is.null(mregfor) && mmex != mxn)
		{
			cat("\nugarchforecast-->error: Column dimension of external mean forecast matrix is wrong.")
			cat(paste("\nModel has ", mxn, " external regressors but forecast matrix has ", mmex, sep = ""))
			stop("\n...exiting\n")
		}
		
		if(!is.null(mregfor) && nmex < treq)
		{
			cat(paste("\nugarchforecast-->error: You requested ", treq ," actual forecasts 
(including the rolling periods) but external mean forecasts provided have only ", nmex, " rows",sep=""))
			cat(paste("\nA minimum of ", treq," rows are required", sep=""))
			stop("\n...exiting")
		}
		if(is.null(mregfor)){
			mxf = rbind(as.matrix(model$modeldata$mexdata)[1:(N-out.sample), ,drop = FALSE], matrix(0, ncol = mxn, nrow = treq))
		} else {
			mxf = rbind(as.matrix(model$modeldata$mexdata)[1:(N-out.sample), ,drop = FALSE], as.matrix(mregfor)[1:treq,,drop = FALSE])
		}
	} else{
		mxf = NULL
	}
	if(vxn>0){
		
		if(!is.null(vregfor)){
			nvex = dim(as.matrix(vregfor))[1]
			mvex = dim(as.matrix(vregfor))[2]
		} else{
			nvex = 0
			mvex = 0
		}
		if(!is.null(vregfor) && mvex != vxn)
		{
			cat("\nugarchforecast-->error: Column dimension of external variance forecast matrix is wrong.")
			cat(paste("\nModel has ",vxn," external regressors but forecast matrix has", mvex, sep = ""))
			stop("\n...exiting\n")
		}
		# N is the original length
		if(!is.null(vregfor) && nvex < treq)
		{
			cat(paste("\nugarchforecast-->error: You requested ", treq ," actual forecasts (including 
									the rolling periods) but external variance forecasts provided have only ",
							nvex, " rows",sep=""))
			cat(paste("\nA minimum of ", treq," rows are required", sep=""))
			stop("\n...exiting")
		}
		
		if(is.null(vregfor)){
			vxf = rbind(as.matrix(model$modeldata$vexdata)[1:(N-out.sample), ,drop = FALSE], matrix(0, ncol = vxn, nrow = treq))
		} else {
			vxf = rbind(as.matrix(model$modeldata$vexdata)[1:(N-out.sample), ,drop = FALSE], as.matrix(vregfor)[1:treq,,drop = FALSE])
		}
	} else{
		vxf = NULL
	}
	return(list(mxf = mxf, vxf = vxf))
}

.simregressors = function(model, mexsimdata, vexsimdata, N, n, m.sim, m)
{
	mxn = model$modelinc[6]
	vxn = model$modelinc[15]
	if(mxn>0){
		if(is.null(mexsimdata)){
			mexsimdata = vector(mode = "list", length = m.sim)
			for(i in 1:m.sim) mexsimdata[[i]] = matrix(0, ncol = mxn, nrow = n)
		}
		if(!is.null(mexsimdata))
		{
			if(!is.list(mexsimdata)) stop("\nugarchsim-->error: mexsimdata should be a list of length m.sim")
			if(length(mexsimdata) != m.sim){
				msd = vector(mode = "list", length = m.sim)
				for(i in 1:m.sim) msd[[i]] = as.matrix(mexsimdata[[1]])
				mexsimdata = msd
				warning("\nugarchsim-->warning: length of mexsimdata list not equal to m.sim...\nreplicating first list element m.sim times.\n")
			}
			for(i in 1:m.sim){
				if(dim(as.matrix(mexsimdata[[i]]))[2] != mxn ) 
					stop(paste("\nugarchsim-->error: mexsimdata ", i," has wrong no. of column", sep=""))
				if(dim(as.matrix(mexsimdata[[i]]))[1] != n )
					stop(paste("\nugarchsim-->error: mexsimdata ", i," has wrong no. of rows", sep=""))
			}		
		}
		mexsimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim){
			premexdata = model$modeldata$mexdata[(N-m+1):N,,drop=FALSE]
			mexsimlist[[i]] = matrix(rbind(premexdata, mexsimdata[[i]]), ncol = mxn)
		}
	} else{
		mexsimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim) mexsimlist[[i]]=0
	}
	if(vxn>0){
		if(is.null(vexsimdata)){
			vexsimdata = vector(mode = "list", length = m.sim)
			for(i in 1:m.sim) vexsimdata[[i]] = matrix(0, ncol = vxn, nrow = n)
		}
		if(!is.null(vexsimdata))
		{
			if(!is.list(vexsimdata)) 
				stop("\nugarchsim-->error: vexsimdata should be a list of length m.sim")
			if(length(vexsimdata) != m.sim){
				vsd = vector(mode = "list", length = m.sim)
				for(i in 1:m.sim) vsd[[i]] = as.matrix(vexsimdata[[1]])
				vexsimdata = vsd
				warning("\nugarchsim-->warning: length of vexsimdata list not equal to m.sim...\nreplicating first list element m.sim times.\n")
			}
			for(i in 1:m.sim){
				if(dim(as.matrix(vexsimdata[[i]]))[2] != vxn ) 
					stop(paste("\nugarchsim-->error: vexsimdata ", i," has wrong no. of column", sep=""))
				if(dim(as.matrix(vexsimdata[[i]]))[1] != n )
					stop(paste("\nugarchsim-->error: vexsimdata ", i," has wrong no. of rows", sep=""))
			}		
		}
		vexsimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim){
			prevexdata = model$modeldata$vexdata[(N-m+1):N,,drop=FALSE]
			vexsimlist[[i]] = matrix(rbind(prevexdata, vexsimdata[[i]]), ncol = vxn)
		}
	} else{
		vexsimlist = vector(mode = "list", length = m.sim)
		for(i in 1:m.sim) vexsimlist[[i]]=0
	}
	return(list(mexsimlist = mexsimlist, vexsimlist = vexsimlist))
}

.custzdist = function(custom.dist, zmatrix, m.sim, n)
{
	# ToDo: clean up (make use of rdist which is now vectorized)
	if(is.na(custom.dist$name) | is.na(custom.dist$distfit)[1]){
		z = matrix(apply(zmatrix, 1, FUN = function(x) .makeSample(as.character(x[1]), lambda = as.numeric(x[2]), 
		skew = as.numeric(x[3]), shape = as.numeric(x[4]), n = as.numeric(x[5]), seed = as.integer(x[6]))), n, m.sim)
	}
	if(!is.na(custom.dist$name) && !is.na(custom.dist$distfit)[1]){
		
		if(is.matrix(custom.dist$distfit))
		{
			if(dim(custom.dist$distfit)[2]!=m.sim) stop("column dimension of custom 
				innovations\n matrix must be equal to m.sim")
			if(dim(custom.dist$distfit)[1]!=n) stop("row dimension 
				of custom innovations\n matrix must be equal to n.sim+n.start")
			z = custom.dist$distfit
		} else{
			if(!is.character(custom.dist$name)) stop("custom distribution must be a 
				character string")
			temp = paste("r", custom.dist$name, sep="")
			if(is.null(custom.dist$distfit)) stop("custom distribution missing a 
				distfit object")
			# if this is often used we might consider using apply to
			# use a new seed for each m.sim
			.rdist = eval(parse(text=paste(temp)))
			tmp = .rdist(n*m.sim, custom.dist$distfit)
			z = matrix(as.numeric(tmp), ncol = m.sim, nrow = n, byrow=TRUE)
		}
	}
	return(z)
}
.stars = function(testvector,levels=c(0.01, 0.05, 0.1))
{
	N = length(testvector)
	ans = vector(mode="character",length=N)
	#recursive replacement
	z = which(testvector<levels[3])
	ans[z] = c("*")
	z = which(testvector<levels[2])
	ans[z] = c("**")
	z = which(testvector<levels[1])
	ans[z] = c("***")
	ans
}

# method from a spec and data object

.safefit = function(spec, data, out.sample, solver, fit.control, solver.control)
{
	ans = try(ugarchfit(spec = spec, data = data, out.sample = out.sample, fit.control = fit.control,
					solver = solver, solver.control = solver.control), silent = TRUE)
	if(inherits(ans, "try-error")) ans = NULL
	return(ans)
}

.safefitarfima = function(spec, data, out.sample, solver, fit.control, solver.control)
{
	ans = try(arfimafit(spec = spec, data = data, out.sample = out.sample, fit.control = fit.control,
					solver = solver, solver.control = solver.control), silent = TRUE)
	if(inherits(ans, "try-error")) ans = NULL
	return(ans)
}


.boxcoxtransform = function(x, lambda)
{
	if(lambda!=0) ret = (x^lambda - 1)/lambda else ret = log(x)
	return(ret)
}

repmat = function(a, n, m)
{
	kronecker(matrix(1, n, m), a)
}
size = function(x, n = NULL)
{
	x = as.matrix(x)
	if(missing(n)) sol = c(n = dim(x)[1], m = dim(x)[2]) else sol = dim(x)[n]
	return(sol)
}

zeros = function(n = 1, m = 1)
{
	if(missing(m)) m = n
	sol = matrix(0, nrow = n, ncol = m)
	return(sol)
}

ones = function(n = 1, m = 1)
{
	if(missing(m)) m = n
	sol = matrix(1, nrow = n, ncol = m)
	return(sol)
}

newlagmatrix = function(x,nlags,xc)
{
	nlags = nlags+1
	xt = size(x, 1);
	newX = rbind(x, zeros(nlags, 1))
	lagmatrix = repmat(newX, nlags, 1)
	lagmatrix = matrix(lagmatrix[1:(size(lagmatrix,1)-nlags)], nrow = (xt+nlags-1), ncol = nlags)
	lagmatrix = lagmatrix[nlags:xt,]
	y = lagmatrix[,1]
	x = lagmatrix[,2:nlags]
	if(xc == 1) x = cbind(ones(size(x,1), 1), x)
	return(data.frame(y = y, x = x))
}

.colorgradient = function(n = 50, low.col = 0.6, high.col=0.7, saturation = 0.8) {
	if (n < 2) stop("n must be greater than 2")
	n1 = n%/%2
	n2 = n - n1
	c(hsv(low.col, saturation, seq(1, 0.5, length = n1)),
			hsv(high.col, saturation, seq(0.5, 1, length = n2)))
}

.distinctcolors11 = function(){
	x = colors()
	c(x[24], x[31], x[32], x[41], x[46], x[51], x[62], x[90], x[139], x[146], x[132])
}

.simlayout = function(m)
{
	if(m == 1){
		nf = c(1, 1, 2, 2)
		nf = layout(matrix(nf, 2, 2, byrow = TRUE), respect = TRUE)
		middle.plot = 1
	}
	if(m == 2){
		nf = c(1, 1, 1, 1, 0, 2, 2, 0, 3, 3, 3, 3)
		nf = layout(matrix(nf, 3, 4, byrow = TRUE), respect = TRUE)
		middle.plot = 2
	}
	if(m == 3){
		nf = c(1, 0, 0, 2, 0, 3, 3, 0, 4, 0, 0, 5)
		nf = layout(matrix(nf, 3, 4, byrow = TRUE), respect = TRUE)
		middle.plot = 3
	}
	if(m == 12){
	nf = c(1, 2, 3, 4,
		   5, 6, 6, 7,
		   8, 6, 6, 9,
		  10, 11, 12, 13)		
	nf = layout(matrix(nf, 4, 4, byrow = TRUE), respect = TRUE)
	}
}

.sdigit = function(x){
	sid = as.numeric(strsplit(format(as.numeric(x), scientific=TRUE), "e")[[1]])[2]
	10^(-sid)
}

.embed<-function(data, k, by = 1, ascending = FALSE)
{
	# n = no of time points, k = number of columns
	# by = increment. normally =1 but if =b calc every b-th point
	# ascending If TRUE, points passed in ascending order else descending.
	# Note that embed(1:n,k) corresponds to embedX(n,k,by=1,rev=TRUE)
	# e.g. embedX(10,3)
	#if(is.null(dim(data)[1])) n<-length(data) else n<-dim(data)[1]
	data = matrix(data,ncol=1)
	n<-dim(data)[1]
	s <- seq(1,n-k+1,by)
	lens <- length(s)
	cols <- if (ascending) 1:k else k:1
	return(matrix(data[s + rep(cols,rep(lens,k))-1],lens))
}

.lagx = function(data, n.lag=1, removeNA = FALSE, pad = NA)
{
	# has NAs
	data = as.matrix(data)
	n = dim(data)[1]
	d = dim(data)[2]
	if(dim(data)[2]==1) data=matrix(data,ncol=1)
	z = apply(data,2,FUN=function(x) .embed(x,n.lag+1)[,n.lag+1])
	if(!removeNA) z = rbind(matrix(pad,ncol=d,nrow=n.lag),z)
	return(z)
}

.lagmatrix = function(data, n.lag = 1, pad = 0)
{
	n = length(as.numeric(data))
	z = matrix(NA, ncol = n.lag, nrow = n)
	for(i in 1:n.lag) z[,i] = .lagx(as.numeric(data), i, removeNA = FALSE, pad = pad)
	return(z)
}


.abind = function(x, y)
{
	m = dim(x)[1]
	if( is.matrix(x) ) {
		n1 = 1
		x = array(x, dim = c(m, m, 1))
	} else{
		n1 = dim(x)[3]
	}
	if( is.matrix(y) ) {
		n2 = 1
		y = array(y, dim = c(m, m, 1))
	} else{
		n2 = dim(y)[3]
	}
	nw = array(NA, dim = c(m, m, n1 + n2))
	nw[ , , 1:n1] = x[, , 1:n1]
	nw[ , , (n1+1):(n1+n2)] = y[, , 1:n2]
	nw
}

.checkrec = function(init, T){
	if(is.null(init)){
		type = 1
		n = T
	} else{
		if(is.character(init)){
			if(init == "all"){
				type = 1
				n = T
			} else{
				warning("\nugarchfit-->warning: unrecognized option in rec.init...using 'all' option instead.\n")
				type = 1
				n = T
			}
		} else{
			if(init>=1){
				type = 1
				if(init<=T){
					n = round(init, 0)
				} else{
					n = T
					warning("\nugarchfit-->warning: rec.init value > n.obs (less out.sample)...setting value to n.obs instead.\n")
				}
			} else{
				if(init>0){
					type = 2
					n = init
				} else{
					warning("\nugarchfit-->warning: unrecognized option in rec.init...using 'all' option instead \n")
					type = 1
					n = T
				}
			}
		}
	}
	return(list(type = type, n = n))
}

# Exponential smoothing backcast
backcastv = function(res, T, lambda, delta=2){
	s = mean(res^delta)
	v = (lambda^T)*s + (1-lambda)*sum(lambda^(0:(T-1))*(res^delta))
	return(v)
}