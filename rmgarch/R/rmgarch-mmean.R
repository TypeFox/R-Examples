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

# First Stage Conditional Mean Filter
# ggn  = goggarch numbering (+2 over copula and dcc models)
mvmean.varfit = function(model, data, VAR.fit = NULL, T, out.sample = 0, cluster = NULL, ggn = 0)
{
	if( model$modelinc[2+ggn]>0 ){
		ex = model$modeldata$mexdata
		if( dim(ex)[1] < dim(data)[1]  ) stop("\ngogarchfit-->error: external regressor matrix has less points than data!...", call. = FALSE)
	} else{
		ex = NULL
	}
	if( !is.null(VAR.fit) ){
		zdata = VAR.fit$xresiduals
		# should be N - p - n.start
		model$modelinc[1+ggn] = p = mp = VAR.fit$lag
		N = dim(zdata)[1]
		if(N != T) stop(paste("\ngogarchfit-->error: Residuals of VAR.fit object should be of length ", T, " (N - out.sample) ", sep=""))
		mu = VAR.fit$xfitted
		varcoef = VAR.fit$Bcoef
	} else{
		ic = model$varmodel$lag.criterion[1]
		if( !is.null(model$varmodel$lag.max) ){
			mp = .varxselect(y = data[1:T, ], lag.max = model$varmodel$lag.max, exogen = ex[1:T, , drop=FALSE])$selection
			mp = as.numeric( mp[which(substr(names(mp), 1, 2) == substr(ic, 1, 2))] )
		} else {
			mp = model$modelinc[1+ggn]
		}
		p = mp
		varmodel = varxfit(X = data[1:T, ,drop=FALSE], p = mp, constant = TRUE, exogen = ex[1:T, , drop=FALSE], 
				robust = model$varmodel$robust, gamma = model$varmodel$robust.control$gamma, 
				delta = model$varmodel$robust.control$delta, nc = model$varmodel$robust.control$nc, 
				ns = model$varmodel$robust.control$ns, postpad = "constant", cluster = cluster)
		zdata = varmodel$xresiduals
		# should be T
		N = dim(zdata)[1]
		mu = varmodel$xfitted
		model$modelinc[1+ggn] = mp
		varcoef = varmodel$Bcoef
	}
	# Filter Full Dataset now
	if(out.sample > 0){
		varmodel2 = varxfilter(X = data, p = mp, exogen = ex, Bcoef = varcoef, postpad = "constant")
		zdata = varmodel2$xresiduals
		# should be T - p  + out.sample
		N = dim(zdata)[1]
		mu = varmodel2$xfitted
	} else{
		zdata = zdata
	} 
	model$varcoef = varcoef
	return(list(model = model, zdata = zdata, mu = mu, varcoef = varcoef, p = p, N = N))
}
# arfit not used by dcc and copula models
mvmean.arfit = function(model, data, ARcoef = NULL, T, out.sample = 0, solver, 
		solver.control = list(), fit.control = list(), cluster = NULL)
{
	m = NCOL(data)
	if(is.null(ARcoef)){
		if( model$modelinc[4] > 0 ){
			ex = model$modeldata$mexdata
			if( dim(ex)[1] < dim(data)[1]  ) 
				stop("\ngogarchfit-->error: external regressor matrix has less points than data!...", call. = FALSE)
		} else{
			ex = NULL
		}
		p = model$modelinc[2]
		# Irrespsective of the Distribution, the mean filtration is always
		# Normal.
		arspec = arfimaspec(mean.model = list(armaOrder = c(p, 0), include.mean = TRUE, 
						arfima = FALSE, external.regressors = ex[1:T, , drop = FALSE]), 
				distribution.model = "norm")
		if( !is.null(cluster) ){
			clusterEvalQ(cluster, loadNamespace("rugarch"))
			clusterExport(cluster, c("T", "arspec", "data", "solver", 
							"solver.control", "fit.control"), envir = environment())
			arfit = parLapply(cluster, as.list(1:m), fun = function(i){
						rugarch::arfimafit(spec = arspec, data = data[1:T,i,drop=FALSE], 
								out.sample = 0, solver = solver, 
								solver.control = solver.control, 
								fit.control = fit.control)
					})
		} else{
			arfit = lapply(1:m, FUN = function(i) arfimafit(spec = arspec, 
								data = data[1:T,i,drop=FALSE], out.sample = 0, 
								solver = solver, solver.control = solver.control, 
								fit.control = fit.control))
		}
		
		
		zdata = sapply(arfit, FUN = function(x) residuals(x))
		#cat("done!\n")
		# rugarch always returns full length data
		#if( p > 0 ) zdata = zdata[-c(1:p), , drop = FALSE]
		# should be T - p
		N = NROW(zdata)
		mu = sapply(arfit, FUN = function(x) fitted(x))
		#if( p > 0 ) mu = mu[-c(1:p), , drop = FALSE]
		arcoef = sapply(arfit, FUN = function(x) coef(x))
	} else{
		if( model$modelinc[4] > 0 ){
			ex = model$modeldata$mexdata
			if( dim(ex)[1] < dim(data)[1]  ) 
				stop("\ngogarchfit-->error: external regressor matrix has less points than data!...", call. = FALSE)
		} else{
			ex = NULL
		}
		
		p = model$modelinc[2]
		# Irrespsective of the Distribution, the mean filtration is always Normal.
		if( !is.null(cluster) ){
			clusterEvalQ(cluster, loadNamespace("rugarch"))
			clusterExport(cluster, c("T", "ARcoef", "data", "solver", 
							"solver.control", "fit.control"), envir = environment())
			arfit = parLapply(cluster, as.list(1:m), fun = function(i){
							arspec = rugarch::arfimaspec(mean.model = list(armaOrder = c(p, 0), 
											include.mean = TRUE, arfima = FALSE, 
											external.regressors = ex[1:T, , drop = FALSE]), 
									distribution.model = "norm", fixed.pars = as.list(ARcoef[,i]));
							ans = rugarch::arfimafilter(spec = arspec, data = data[1:T,i,drop=FALSE], 
									out.sample = 0);
							return(ans)
						})
		} else{
			arfit = lapply(1:m, FUN = function(i){
						arspec = arfimaspec(mean.model = list(armaOrder = c(p, 0), include.mean = TRUE, 
										arfima = FALSE, external.regressors = ex[1:T, , drop = FALSE]), 
								distribution.model = "norm", fixed.pars = as.list(ARcoef[,i]));
						ans = arfimafilter(spec = arspec, data = data[1:T,i,drop=FALSE], out.sample = 0);
						return(ans)
					})
		}
		zdata = sapply(arfit, FUN = function(x) residuals(x))
		#cat("done!\n")
		# rugarch always returns full length data
		#if( p > 0 ) zdata = zdata[-c(1:p), , drop = FALSE]
		# should be T - p
		N = NROW(zdata)
		mu = sapply(arfit, FUN = function(x) fitted(x))
		#if( p > 0 ) mu = mu[-c(1:p), , drop = FALSE]
		arcoef = ARcoef
	}
	if(out.sample>0){
		if( !is.null(cluster) ){
			clusterExport(cluster, c("data", "p", "ex", "arcoef"), 
					envir = environment())			
			arfit2 = parLapply(cluster, as.list(1:m), fun = function(i) {
							specx = arfimaspec(mean.model = list(armaOrder = c(p, 0), 
											include.mean = TRUE, arfima = FALSE, 
											external.regressors = ex), 
									distribution.model = "norm",
									fixed.pars = as.list(arcoef[,i]));
							ans = arfimafilter(spec = specx, data = data[,i, drop=FALSE], out.sample = 0);
							return(ans)
						})
		} else{
			arfit2 = lapply(1:m, FUN = function(i){
						specx = arfimaspec(mean.model = list(armaOrder = c(p, 0), include.mean = TRUE, 
										arfima = FALSE, external.regressors = ex), distribution.model = "norm",
								fixed.pars = as.list(arcoef[,i]));
						ans = arfimafilter(spec = specx, data = data[,i, drop = FALSE], out.sample = 0);
						return(ans)
					})
		}
		zdata = sapply(arfit2, FUN = function(x) residuals(x))
		#if( p > 0 ) zdata = zdata[-c(1:p), , drop = FALSE]
		N = NROW(zdata)
	} else{
		zdata = zdata
	}
	model$arcoef = arcoef
	return(list(model = model, zdata = zdata, mu = mu, arcoef = arcoef, p = p, N = N))
	
}


mvmean.varfilter = function(model, data, varcoef, T, out.sample = 0, ggn = 0)
{
	if( model$modelinc[2+ggn]>0 ){
		ex = model$modeldata$mexdata
		if( dim(ex)[1] < dim(data)[1]  ) 
			stop("\ngogarchfilter-->error: external regressor matrix has less points than data!...", call. = FALSE)
	} else{
		ex = NULL
	}
	if(is.null(varcoef)) 
		stop("\ncgarchfilter-->error: VAR model requires varcoef in specification to be supplied.")
	p = model$modelinc[1+ggn]
	vrmodel = varxfilter(X = data[1:T, , drop=FALSE], p = p, exogen = ex[1:T, , drop=FALSE], Bcoef = varcoef, postpad = "constant")
	#cat("done!\n")
	zdata = vrmodel$xresiduals
	# should be N - p
	N = NROW(data)
	mu = vrmodel$xfitted
	# if we have out of sample data we need to filter
	if(out.sample > 0){
		varmodel2 = varxfilter(X = data, p = p, exogen = ex, Bcoef = varcoef, postpad = "constant")
		zdata = varmodel2$xresiduals
		N = NROW(zdata)
		mu = varmodel2$xfitted
	} else{
		zdata = zdata
	}
	model$varcoef = varcoef
	model$residuals = zdata
	return(list(model = model, zdata = zdata, mu = mu, p = p, N = N))
}


mvmean.arfilter = function(model, data, arcoef, T, out.sample = 0, cluster = NULL)
{
	if( model$modelinc[4] > 0 ){
		ex = model$modeldata$mexdata
		if( dim(ex)[1] < dim(data)[1]  ) 
			stop("\ngogarchfilter-->error: external regressor matrix has less points than data!...", call. = FALSE)
	} else{
		ex = NULL
	}
	m = NCOL(data)
	p = model$modelinc[2]
	if( !is.null(cluster) ){
			clusterEvalQ(cluster, loadNamespace("rugarch"))
			clusterExport(cluster, c("T", "data", "p", "ex", "arcoef"), 
					envir = environment())
			armodel = parLapply(cluster, as.list(1:m), fun = function(i) {
						specx = rugarch::arfimaspec(mean.model = list(armaOrder = c(p, 0), 
										include.mean = TRUE, arfima = FALSE, 
										external.regressors = ex[1:T, ]), 
								distribution.model = "norm",
								fixed.pars = as.list(arcoef[, i]))
						ans = rugarch::arfimafilter(spec = specx, data = data[1:T, i], out.sample = 0)
						return(ans)
					})

	} else{
		armodel = lapply(1:m, FUN = function(i){
					specx = arfimaspec(mean.model = list(armaOrder = c(p, 0), include.mean = TRUE, 
									arfima = FALSE, external.regressors = ex[1:T, ]), distribution.model = "norm",
							fixed.pars = as.list(arcoef[, i]))
					ans = arfimafilter(spec = specx, data = data[1:T, i], out.sample = 0)
					return(ans)
				})
	}
	#cat("done!\n")
	zdata = sapply(armodel, FUN = function(x) residuals(x))
	#if( p > 0 ) zdata = zdata[-c(1:p), , drop = FALSE]
	# should be N - p
	N = NROW(data)
	mu = sapply(armodel, FUN = function(x) fitted(x))
	#if( p > 0 ) mu = mu[-c(1:p), , drop = FALSE]		
	# if we have out of sample data we need to filter
	if(out.sample > 0){
		if( !is.null(cluster) ){
			clusterExport(cluster, c("data", "p", "ex", "arcoef"), 
					envir = environment())			
			armodel2 = parLapply(cluster, as.list(1:m), fun = function(i){
						specx = arfimaspec(mean.model = list(armaOrder = c(p, 0), 
										include.mean = TRUE, arfima = FALSE, 
										external.regressors = ex), 
								distribution.model = "norm",
								fixed.pars = as.list(arcoef[,i]))
						ans = arfimafilter(spec = specx, data = data[,i], out.sample = 0)
						return(ans)
					})
		} else{
			armodel2 = lapply(1:m, FUN = function(i){
						specx = arfimaspec(mean.model = list(armaOrder = c(p, 0), 
										include.mean = TRUE, arfima = FALSE, 
										external.regressors = ex), 
								distribution.model = "norm",
								fixed.pars = as.list(arcoef[,i]))
						ans = arfimafilter(spec = specx, data = data[,i], out.sample = 0)
						return(ans)
					})
		}			
		zdata = sapply(armodel2, FUN = function(x) residuals(x))
		#if( p > 0 ) zdata = zdata[-c(1:p), , drop = FALSE]
	} else{
		zdata = zdata
	}
	return(list(zdata = zdata, mu = mu, p = p, N = N))
}


mvmean.varforecast = function(model, mregfor, varcoef, n.ahead, n.roll, ns, ggn = 0){
	m = NCOL(model$modeldata$data)
	tf = n.ahead + n.roll
	if( model$modelinc[2+ggn] > 0 ){
		if( is.null( mregfor ) ){
			warning("\nExternal Regressor Forecasts Matrix NULL...setting to zero...\n")
			mregfor = matrix(0, ncol = model$modelinc[2+ggn], nrow = (n.roll + n.ahead) )
		} else{
			if( dim(mregfor)[2] != model$modelinc[2+ggn] ) 
				stop("\ngogarchforecast-->error: wrong number of external regressors!...", call. = FALSE)
			if( dim(mregfor)[1] < (n.roll + n.ahead) ) 
				stop("\ngogarchforecast-->error: external regressor matrix has less points than requested forecast length (1+n.roll) x n.ahead!...", call. = FALSE)
		}
	} else{
		mregfor = NULL
	}
	if( n.ahead == 1 && (n.roll > ns) ) 
		stop("\ngogarchforecast-->error: n.roll greater than out.sample!", call. = FALSE)
	
	# varxforecast returns array --> change to matrix since we do not allow
	# mixed n.roll and n.ahead
	Mu = varxforecast(X = model$modeldata$data, Bcoef = varcoef, p = model$modelinc[1+ggn], out.sample = ns, n.ahead, n.roll, mregfor)
	return(Mu)
}

mvmean.arforecast = function(model, mregfor, arcoef, n.ahead, n.roll, ns, cluster = NULL)
{
	m = NCOL(model$modeldata$data)
	tf = n.ahead + n.roll
	if( !is.null(cluster) ){
		clusterEvalQ(cluster, loadNamespace("rugarch"))
		clusterExport(cluster, c("model", "arcoef", "mregfor", "ns", 
						"n.ahead", "n.roll"), envir = environment())
		arfx = parLapply(cluster, as.list(1:m), fun = function(i){
					specx = rugarch::arfimaspec(mean.model = list(armaOrder = c(model$modelinc[2], 0), 
									include.mean = TRUE, arfima = FALSE, 
									external.regressors = model$modeldata$mexdata), 
							distribution.model = "norm", fixed.pars = as.list(arcoef[,i]))
					ans = rugarch::arfimaforecast(fitORspec = specx, data = model$modeldata$data[,i], 
							out.sample = ns, n.ahead = n.ahead, n.roll = n.roll, 
							external.forecasts = list(mregfor = mregfor))
					return(ans)
				})
	} else{
		arfx = lapply(1:m, FUN = function(i){
					specx = arfimaspec(mean.model = list(armaOrder = c(model$modelinc[2], 0), 
									include.mean = TRUE, arfima = FALSE, 
									external.regressors = model$modeldata$mexdata), 
							distribution.model = "norm",
							fixed.pars = as.list(arcoef[,i]))
					ans = arfimaforecast(fitORspec = specx, data = model$modeldata$data[,i], 
							out.sample = ns, n.ahead = n.ahead, n.roll = n.roll, 
							external.forecasts = list(mregfor = mregfor))
					return(ans)
				})
	}
	Mu = array(NA, dim=c(n.ahead, m, n.roll+1))
	for(i in 1:(n.roll+1)){
		Mu[,,i] = matrix(sapply(arfx, FUN = function(x) fitted(x)[,i], simplify=TRUE), ncol = m)
	}
	return(Mu)
}


mvmean.varsim = function(model, Data, res, mexsimdata, prereturns, m.sim, 
		n.sim, n.start, startMethod = "sample", cluster = NULL, ggn = 0)
{
	seriesSim = vector(mode = "list", length = m.sim)
	m = NCOL(Data)
	varcoef = model$varcoef
	if( model$modelinc[2+ggn] > 0 ){
		if( !is.null(mexsimdata) ){
			if( !is.list(mexsimdata) | length(mexsimdata) != m.sim ) 
				stop("\ngogarchsim-->error: mexsimdata must be a list of length equal to m.sim...", call. = FALSE)
			for(j in 1:m.sim){
				if( !is.matrix(mexsimdata[[j]]) ) 
					stop("\ngogarchsim-->error: mexsimdata list submatrix must be a matrix...", call. = FALSE)
				if( dim(mexsimdata[[j]])[2] != model$modelinc[2+ggn] ) 
					stop("\ngogarchsim-->error: wrong number of external regressors in list submatrix!...", call. = FALSE)
				if( dim(mexsimdata[[j]])[1] < (n.sim + n.start) ) 
					stop("\ngogarchsim-->error: external regressor list submatrix has less points than requested simulation length (n.sim + n.start)...", call. = FALSE)
			}
		} else{
			mexsimdata = vector(mode = "list", length = m.sim)
			for(j in 1:m.sim) mexsimdata[[j]] = matrix(0, ncol = model$modelinc[2+ggn], nrow = (n.sim + n.start))
		}
	} else{
		mexsimdata = NULL
	}
	p = model$modelinc[1+ggn]
	if(startMethod == "sample"){
		if( !is.null(prereturns) ){
			if(!is.matrix(prereturns)) 
				stop("\nprereturns must be a matrix\n")
			if(dim(prereturns)[1] < p)
				stop("\nwrong number of rows for prereturns matrix")
			if(dim(prereturns)[2] != m) 
				stop("\nwrong number of cols for prereturns matrix")
			prereturns = matrix(tail(prereturns, p), nrow = p, ncol = m)
		} else{
			prereturns = as.matrix( tail(Data, p) )
		}
	} else{
		prereturns = matrix(0, nrow = p, ncol = m, byrow = TRUE)
	}
	if( n.sim == 1 && n.start == 0 ){
		# SPECIAL ROUTINE TO SPEED UP 1-AHEAD ESTIMATION
		tmp = varxsimXX(X = Data, varcoef, p = p, m.sim = m.sim, prereturns, resids = t(sapply(res, FUN = function(x) x)), if(!is.null(mexsimdata)) t(sapply(mexsimdata, FUN = function(x) x)) else NULL)
		for(i in 1:m.sim) seriesSim[[i]] = tmp[i, , drop = FALSE]
	} else{
		if( !is.null(cluster) ){
			# when the GARCH model is without mean equation, simX == simulated residuals
			clusterEvalQ(cluster, require(rmgarch))
			clusterExport(cluster, c("Data", "varcoef", "p", "n.sim", 
							"n.start", "prereturns", "res", "mexsimdata"), envir = environment())
			seriesSim = parLapply(cluster, as.list(1:m.sim), fun = function(i){
						varxsim(X = Data, varcoef, p = p, n.sim = n.sim, 
								n.start = n.start, prereturns = prereturns, 
								resids = res[[i]], mexsimdata = mexsimdata[[i]])
					})
		} else{
			seriesSim = lapply(as.list(1:m.sim), FUN = function(i){
						varxsim(X = Data, varcoef, p = p, n.sim = n.sim, 
								n.start = n.start, prereturns = prereturns, 
								resids = res[[i]], mexsimdata = mexsimdata[[i]])
					})
		}
	}
	return( seriesSim )
}

mvmean.arsim = function(model, Data, res, arcoef, mxn, mexdata, mexsimdata, 
		prereturns, preresiduals, m.sim, n.sim, n.start, startMethod = "sample", 
		cluster = NULL)
{
	seriesSim = vector(mode = "list", length = m.sim)
	m = NCOL(Data)
	if(model$modelinc[4]>0){
		if( !is.null(mexsimdata) ){
			if( !is.list(mexsimdata) | length(mexsimdata) != m.sim ) 
				stop("\ngogarchsim-->error: mexsimdata must be a list of length equal to m.sim...", call. = FALSE)
			for(j in 1:m.sim){
				if( !is.matrix(mexsimdata[[j]]) ) 
					stop("\ngogarchsim-->error: mexsimdata list submatrix must be a matrix...", call. = FALSE)
				if( dim(mexsimdata[[j]])[2] != model$modelinc[4] ) 
					stop("\ngogarchsim-->error: wrong number of external regressors in list submatrix!...", call. = FALSE)
				if( dim(mexsimdata[[j]])[1] < (n.sim + n.start) ) 
					stop("\ngogarchsim-->error: external regressor list submatrix has less points than requested simulation length (n.sim + n.start)...", call. = FALSE)
			}
		} else{
			mexsimdata = vector(mode = "list", length = m.sim)
			for(j in 1:m.sim) mexsimdata[[j]] = matrix(0, ncol = model$modelinc[4], nrow = (n.sim + n.start))
		}
	} else{
		mexsimdata = NULL
	}
	p = model$modelinc[2]
	if(!is.na(prereturns)[1]){
		if(!is.matrix(prereturns)) 
			stop("\nprereturns must be a matrix\n")
		if(p>0) prereturns = tail(prereturns, p)
		#if(dim(prereturns)[1] != p) stop("\nwrong number of rows for prereturns matrix")
		if(p>0) if(NCOL(prereturns) != m) 
			stop("\nwrong number of cols for prereturns matrix")
	} else{
		if(startMethod == "sample"){
			if(p>0) prereturns = as.matrix( tail(Data, p) )
		} else{
			if(p>0) prereturns = matrix(0, nrow = p, ncol = m, byrow = TRUE)
		}
	}
	
	if(!is.na(preresiduals)[1]){
		if(!is.matrix(preresiduals)) 
			stop("\npreresiduals must be a matrix\n")
		if(p>0) preresiduals = tail(preresiduals, p)
		if(p>0) if(NCOL(preresiduals) != m) 
				stop("\nwrong number of cols for prereturns matrix")
	} else{
		if(startMethod != "sample") 
			if(p>0) preresiduals = matrix(0, nrow = p, ncol = m, byrow = TRUE)
	}
	
	if(model$modelinc[4]>0){
		ex = mexdata[1:T, , drop = FALSE]
	} else{
		ex = NULL
	}
	if( !is.null(cluster) ){
			tres = matrix(tail(preresiduals, p), ncol = m)
			clusterEvalQ(cluster, loadNamespace("rugarch"))
			clusterExport(cluster, c("p", "ex", "arcoef", "mexsimdata", 
							"res", "n.start", "tres", "n.sim", "m.sim", 
							"prereturns"), envir = environment())
			arxsim = parLapply(cluster, as.list(1:m), fun = function(i){
						specx = rugarch::arfimaspec(mean.model = list(armaOrder = c(p, 0), 
										include.mean = TRUE, arfima = FALSE, 
										external.regressors = ex), 
								distribution.model = "norm",
								fixed.pars = as.list(arcoef[,i]))
						ans = rugarch::arfimapath(spec = specx, n.sim = n.sim, 
								n.start = n.start, m.sim = m.sim, 
								prereturns = if(p>0) tail(prereturns[,i],p) else NA,
								preresiduals = if(p>0) tres[,i] else NA, 
								custom.dist = list(name = "sample", 
										distfit = matrix(sapply(res, FUN = function(x) x[,i]), 
												ncol = m.sim, nrow = n.sim, byrow = TRUE), 
										type = "res"), 
								mexsimdata = mexsimdata)
						return(ans)
					})
	} else{
		arxsim = lapply(1:m, FUN = function(i){
					specx = arfimaspec(mean.model = list(armaOrder = c(p, 0), 
									include.mean = TRUE, arfima = FALSE, 
									external.regressors = ex), 
							distribution.model = "norm", 
							fixed.pars = as.list(arcoef[,i]))
					ans = arfimapath(spec = specx, n.sim = n.sim, n.start = n.start, 
							m.sim = m.sim, prereturns = if(p>0) tail(prereturns[,i], p) else NA,
							preresiduals = if(p>0) tail(preresiduals[,i], p) else NA, 
							custom.dist = list(name = "sample", 
									distfit = matrix(sapply(res, FUN = function(x) x[,i]), 
											ncol = m.sim, nrow = n.sim, byrow = TRUE), 
									type = "res"), 
							mexsimdata = mexsimdata)
					return(ans)
				})
	}
	for(i in 1:m.sim) seriesSim[[i]] = matrix(sapply(arxsim, FUN = function(x) x@path$seriesSim[,i]), ncol = m, nrow = n.sim, byrow = TRUE)
	
	return( seriesSim )
}