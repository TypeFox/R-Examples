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

# DCC CGARCH and GOGARCH
.fmoments = function(spec, Data, n.ahead = 1, roll = 0, solver = "solnp", 
		solver.control = list(), fit.control = list(eval.se=FALSE), 
		cluster = NULL, save.output = FALSE, save.dir = getwd(), 
		save.name = paste("M", sample(1:1000, 1), sep= ""), ...)
{
	if(class(spec)[1]!= "goGARCHspec" & class(spec)[1]!="DCCspec") stop("\nonly gogarch and dcc models supported.")
	if(n.ahead>1 && roll>0) stop("\nn.ahead must be equal to 1 when using roll")
	asset.names = colnames(Data)
	m = dim(Data)[2]
	xmodel = list(asset.names = asset.names, assets = m, 
			model = ifelse(class(spec)[1]=="goGARCHspec", "gogarch", "dcc"), 
			n.ahead = n.ahead, roll = roll, save.output = save.output, save.dir = save.dir,
			save.name = save.name)
	# ... in mom.dcc are to pass the realizedVol matrix whereas in the GOGARCH model
	# (which does not use realizedVol), they are related to extra arguments passed to
	# the fastica routine
	mom = switch(class(spec)[1],
			goGARCHspec = mom.gogarch(Data = Data, spec = spec, n.ahead = n.ahead, 
					roll  = roll, solver = solver, solver.control = solver.control, 
					fit.control = fit.control, cluster = cluster,
					save.output = save.output, save.dir = save.dir, 
					save.name = save.name, ...),
			DCCspec = mom.dcc(Data = Data, spec = spec, n.ahead = n.ahead, 
					roll  = roll, solver = solver, solver.control = solver.control, 
					fit.control = fit.control, cluster = cluster,
					save.output = save.output, save.dir = save.dir, 
					save.name = save.name, ...))
	
	if(!save.output){
		sol = new("fMoments",
				moments = mom,
				model = xmodel)
	} else{
		sol = new("fMoments",
				moments = list(),
				model = xmodel)
	}
	return(sol)
}

mom.gogarch = function(Data, spec, n.ahead = 1, roll  = 0, 
		solver = "solnp", solver.control = list(), fit.control = list(), 
		cluster = NULL, save.output = FALSE, save.dir = getwd(), 
		save.name = NULL, ...)
{
	fit = gogarchfit(spec, data = Data, out.sample = roll, solver = solver, 
			cluster = cluster, fit.control = fit.control, solver.control = solver.control, ...)
	A = fit@mfit$A
	K = fit@mfit$K	
	T = fit@model$modeldata$T
	T0 = fit@model$modeldata$index[(T):(T+roll)]
	indexf = c(T0[-1], generatefwd(tail(T0, 1), length.out=1, by = fit@model$modeldata$period))
	nam = fit@model$modeldata$asset.names
	m  = NCOL(Data)
	# sim will have an analytical forecast and simulated scenario
	forc = gogarchforecast(fit, n.ahead = n.ahead, n.roll = roll, cluster = cluster)
	forecastMu = matrix(forc@mforecast$mu, ncol = m, nrow = roll+1, byrow=TRUE)
	forecastCov = rcov(forc)
	colnames(forecastMu) = nam
	rownames(forecastMu) = as.character(indexf)
	forecastCov = array(unlist(rcov(forc)), dim=c(m,m,roll+1), dimnames=list(nam, nam, as.character(indexf)))
	
	if(spec@model$modeldesc$distribution != "mvnorm"){
		forecastM3 = rcoskew(forc, standardize  = FALSE, roll="all")
		forecastM4 = rcokurt(forc, standardize  = FALSE, roll="all")
		dimnames(forecastM3) = list(NULL, NULL, as.character(indexf))
		dimnames(forecastM4) = list(NULL, NULL, as.character(indexf))
	} else{
		forecastM3 = NULL
		forecastM4 = NULL
	}
	if(save.output){
		olddir = getwd()
		setwd(save.dir)
		mom = list(forecastMu = forecastMu, forecastCov = forecastCov, 
				forecastM3 = forecastM3, forecastM4 = forecastM4,
				gogarch.options = list(A = A, K = K), n.ahead = n.ahead, 
				roll = roll, model = "gogarch", index = indexf)
		eval(parse(text = paste("save(mom, file = '", save.name,".rda')",sep="")))
		setwd(olddir)
		return(1)
	} else{
		return(list(forecastMu = forecastMu, forecastCov = forecastCov, 
						forecastM3 = forecastM3, forecastM4 = forecastM4,
						gogarch.options = list(A = A, K = K), index = indexf))
	}
}

mom.dcc = function(Data, spec, n.ahead = 1, roll  = 0, 
		solver = "solnp", solver.control = list(), fit.control = list(), 
		cluster = NULL, save.output = FALSE, save.dir = getwd(), 
		save.name = NULL, ...)
{
	fit = dccfit(spec, data = Data, out.sample = roll, solver = solver, 
			cluster = cluster, fit.control = fit.control, 
			solver.control = solver.control, ...)
	T0 = fit@model$modeldata$index[(T):(T+roll)]
	indexf = c(T0[-1], generatefwd(tail(T0, 1), length.out=1, by = fit@model$modeldata$period))
	nam = fit@model$modeldata$asset.names
	m  = NCOL(Data)
	# sim will have an analytical forecast and simulated scenario
	forc = dccforecast(fit, n.ahead = n.ahead, n.roll = roll, cluster = cluster)	
	forecastMu = matrix(forc@mforecast$mu, ncol = m, nrow = roll+1, byrow=TRUE)
	forecastCov = rcov(forc)
	colnames(forecastMu) = nam
	rownames(forecastMu) = as.character(indexf)
	forecastCov = array(unlist(rcov(forc)), dim=c(m,m,roll+1), dimnames=list(nam, nam, as.character(indexf)))	
	if(spec@model$modeldesc$distribution == "mvt"){
		skewness = 0
		kurtosis = unname(6/(coef(fit)["[Joint]mshape"] - 4) + 4)
	} else{
		kurtosis = NULL
		skewness = NULL
	}
	if(save.output){
		olddir = getwd()
		setwd(save.dir)
		mom = list(forecastMu = forecastMu, forecastCov = forecastCov, 
				skewness = skewness, kurtosis = kurtosis,
				n.ahead = n.ahead, roll = roll, model = "dcc", index = indexf)
		eval(parse(text = paste("save(mom, file = '",save.name,".rda')",sep="")))
		setwd(olddir)
		return(1)
	} else{
		return(list(forecastMu = forecastMu, forecastCov = forecastCov, 
						skewness = skewness, kurtosis = kurtosis, index = indexf))		
	}
}

####----------------------------------------------------------------------------
.fscenario = function(Data, sim = 1000, roll = 0, 
		model = c("gogarch", "dcc", "cgarch", "var", "mdist"), spec = NULL,
		var.model = list(
				lag = 1, lag.max = NULL, 
				lag.criterion = c("AIC", "HQ", "SC", "FPE"), 
				robust = FALSE, 
				robust.control = list("gamma" = 0.25, "delta" = 0.01, "nc" = 10, "ns" = 500)),
		mdist.model = list(distribution = c("mvn", "mvt", "manig"), AR = TRUE, lag = 1),
		spd.control = list(lower = 0.1, upper = 0.9, type = "pwm", kernel = "epanech"),
		cov.method = c("ML", "LW", "EWMA", "MVE", "MCD", "MVT", "BS"),
		cov.options = list(shrinkage=-1, lambda = 0.96),
		solver = "solnp", solver.control = list(), fit.control = list(eval.se = 1), 
		cluster = NULL, save.output = FALSE, save.dir = getwd(), 
		save.name = paste("S", sample(1:1000, 1), sep= ""), rseed  = NULL, ...)
{
	# ... additional parameters passed for example to gogarchfit routine (gfun...)
	tmp = match(tolower(model[1]), c("gogarch", "dcc", "cgarch", "var", "mdist"))
	model = model[1]
	sim = max(1, as.integer(sim))
	roll = max(0, as.integer(roll))
	m = NCOL(Data)
	if(is.na(tmp)) stop("\nunrecongnized model choice.")
	if(tolower(model) == "gogarch"){
		if(!is(spec, "goGARCHspec")) stop("\ngogarch chosen but spec not of goGARCHspec class")
		asset.names = colnames(Data)
		Scen = scenario.gogarch(Data = Data, spec = spec, sim = sim, 
				roll = roll, solver = solver, solver.control = solver.control,
				fit.control = fit.control, cluster = cluster,
				rseed = rseed, save.output = save.output, save.dir = save.dir, 
				save.name = save.name, ...)
		xmodel = list(asset.names = asset.names, assets = m, model = "gogarch", 
				sim = sim, roll = roll, save.output = save.output, save.dir = save.dir,
				save.name = save.name)
		if(!save.output){
			sol = new("fScenario",
					scenario = Scen,
					model = xmodel)
		} else{
			sol = new("fScenario",
					scenario = list(),
					model = xmodel)
		}
	} else if(tolower(model) == "dcc"){
		if(!is(spec, "DCCspec")) stop("\ndcc chosen but spec not of DCCspec class")
		asset.names = colnames(Data)
		Scen = scenario.dcc(Data = Data, spec = spec, sim = sim, 
				roll = roll, solver = solver, solver.control = solver.control,
				fit.control = fit.control, cluster = cluster,
				rseed = rseed, save.output = save.output, save.dir = save.dir, 
				save.name = save.name, ...)
		xmodel = list(asset.names = asset.names, assets = m, model = "dcc", 
				sim = sim, roll = roll, save.output = save.output, save.dir = save.dir,
				save.name = save.name)
		if(!save.output){
			sol = new("fScenario",
					scenario = Scen,
					model = xmodel)
		} else{
			sol = new("fScenario",
					scenario = list(),
					model = xmodel)
		}
	} else if(tolower(model) == "cgarch"){
		if(!is(spec, "cGARCHspec")) stop("\ncgarch chosen but spec not of cGARCHspec class")
		m = NCOL(Data)[2]
		asset.names = colnames(Data)
		Scen = scenario.cgarch(Data = Data, spec = spec, 
				spd.control = spd.control, sim = sim, roll = roll, 
				solver = solver, solver.control = solver.control,
				fit.control = fit.control, cluster = cluster,
				rseed = rseed, save.output = save.output, save.dir = save.dir, 
				save.name = save.name, ...)
		xmodel = list(asset.names = asset.names, assets = m, model = "cgarch", 
				sim = sim, roll = roll, save.output = save.output, save.dir = save.dir,
				save.name = save.name)
		if(!save.output){
			sol = new("fScenario",
					scenario = Scen,
					model = xmodel)
		} else{
			sol = new("fScenario",
					scenario = list(),
					model = xmodel)
		}
	} else if(tolower(model) == "var"){
		asset.names = colnames(Data)
		Scen = scenario.var(Data = Data, sim = sim, roll = roll, model = var.model, 
				cov.method = cov.method, cov.options = cov.options, rseed = rseed, 
				save.output = save.output, save.dir = save.dir, save.name = save.name)
		xmodel = list(asset.names = asset.names, assets = m, model = "var", 
				sim = sim, roll = roll, save.output = save.output, save.dir = save.dir,
				save.name = save.name)
		if(!save.output){
			sol = new("fScenario",
					scenario = Scen,
					model = xmodel)
		} else{
			sol = new("fScenario",
					scenario = list(),
					model = xmodel)
		}
	} else if(tolower(model) == "mdist"){
		asset.names = colnames(Data)
		Scen = scenario.mdist(Data = Data, sim = sim, roll = roll, model = mdist.model, 
				rseed = rseed, save.output = save.output, save.dir = save.dir,
				save.name = save.name)
		if(!save.output){
			sol = new("fScenario",
					scenario = Scen,
					model = xmodel)
		} else{
			sol = new("fScenario",
					scenario = list(),
					model = xmodel)
		}
	} else{
		sol = NULL
	}
	return(sol)
}

scenario.gogarch = function(Data, spec, sim = 1000, roll = 0, solver = "solnp", 
		solver.control = list(), fit.control = list(stationarity = 1), 
		cluster = NULL, rseed = NULL, save.output = FALSE, save.dir = getwd(), 
		save.name = NULL, ...)
{
	fit = gogarchfit(spec, data = Data, out.sample = roll, solver = solver, 
			cluster = cluster, fit.control = fit.control, 
			solver.control = solver.control, ...)
	T = fit@model$modeldata$T
	T0 = fit@model$modeldata$index[(T):(T+roll)]
	indexf = c(T0[-1], generatefwd(tail(T0, 1), length.out=1, by = fit@model$modeldata$period))
	nam = fit@model$modeldata$asset.names
	m = NCOL(Data)
	#xfit<<-fit
	A = fit@mfit$A
	K = fit@mfit$K
	W = fit@mfit$W
	fsim = vector(mode="list", length = roll)
	# sim will have an analytical forecast and simulated scenario
	forc = gogarchforecast(fit, n.ahead = 1, n.roll = roll, cluster = cluster)
	forecastMu = forecastMu = matrix(forc@mforecast$mu, ncol = m, nrow = roll+1, byrow=TRUE)
	forecastCov = rcov(forc)
	colnames(forecastMu) = nam
	rownames(forecastMu) = as.character(indexf)
	forecastCov = array(unlist(rcov(forc)), dim=c(m,m,roll+1), dimnames=list(nam, nam, as.character(indexf)))
	# max lag for mean.model = 2
	p = sum(fit@model$modelinc[1:3])	
	xseed = simRes = simMu = vector(mode = "list", length = roll+1)
	filt = gogarchfilter(fit, data = Data, out.sample = 0, n.old = T, 
			cluster = cluster)
	# use parallel here:
	if(!is.null(rseed)){
		rseed = as.numeric(rseed)
		if(length(rseed)==(roll+1)) zseed = rseed
		if(length(rseed)>(roll+1)) zseed = rseed[1:(roll+1)]
		if(length(rseed)<(roll+1)) zseed = rseed[1]+seq_len(roll)
	} else{
		zseed = as.integer( runif( roll+1, 0, as.integer(Sys.time()) ) )
	}
	for(i in 1:(roll+1)){
		preres = matrix(fit@mfit$Y[(T-p+i):(T-1+i),,drop=FALSE], ncol = NCOL(fit@mfit$Y), byrow = TRUE)
		presigma = matrix(filt@mfilter$factor.sigmas[(T-p+i):(T-1+i),,drop=FALSE], ncol = NCOL(fit@mfit$Y), byrow = TRUE)
		prereturns = as.matrix(fit@model$modeldata$data[(T-p+i):(T-1+i),,drop=FALSE])
		fsim = gogarchsim(fit, n.sim = 1, m.sim = sim, startMethod = "sample",
				preresiduals = preres, presigma = presigma, 
				prereturns = prereturns, rseed = zseed[i], cluster = cluster)
		simMu[[i]] = t(sapply(fsim@msim$seriesSim, FUN = function(x) x))
		simRes[[i]] = t(sapply(fsim@msim$residSim, FUN = function(x) x))
		xseed[[i]] = fsim@msim$rseed
	}
	if(save.output){
		olddir = getwd()
		setwd(save.dir)
		scenario = list(simMu = simMu, simRes = simRes, forecastMu = forecastMu, 
				forecastCov = forecastCov, gogarch.options = list(A = A, K = K), 
				rseed = xseed, sim = sim, roll = roll, model = "gogarch", index = indexf)
		eval(parse(text = paste("save(scenario, file = '",save.name,".rda')",sep="")))
		setwd(olddir)
		return(1)
	} else{
		return(list(simMu = simMu, simRes = simRes, forecastMu = forecastMu, 
						forecastCov = forecastCov, 
						gogarch.options = list(A = A, K = K, W = W), rseed = xseed,
						index = indexf))	
	}
}

scenario.dcc = function(Data, spec, sim, roll, solver = "solnp", 
		fit.control = list(eval.se = FALSE), solver.control = list(), 
		cluster = NULL, rseed = NULL, save.output = FALSE, save.dir = getwd(), 
		save.name = NULL, debug = FALSE, ...)
{
	pars = list(...)
	modv = spec@spec[[1]]@model$modeldesc$vmodel
	if(is.null(pars$fit)) xfit = NULL else xfit = pars$fit
	if(is.null(pars$VAR.fit)) VAR.fit = NULL else VAR.fit = pars$VAR.fit
	if(modv=="realGARCH"){
		if(is.null(pars$realizedVol)){
			stop("\nfscenario-->error: realizedVol required for realGARCH model.\n")
		} else{
			realizedVol = pars$realizedVol
		}
	} else{
		realizedVol = NULL
	}
	fit = dccfit(spec, data =  Data, out.sample = roll, solver = solver, 
			solver.control = list(), fit.control = fit.control, cluster = cluster,
			fit = xfit, VAR.fit = VAR.fit, realizedVol = realizedVol)
	
	T = fit@model$modeldata$T
	T0 = fit@model$modeldata$index[(T):(T+roll)]
	indexf = c(T0[-1], generatefwd(tail(T0, 1), length.out=1, by = fit@model$modeldata$period))
	nam = fit@model$modeldata$asset.names
	m  = NCOL(Data)
	
	fsim = vector(mode="list", length = roll)
	p = spec@model$maxgarchOrder
	# sim will have an analytical forecast and simulated scenario
	forc = dccforecast(fit, n.ahead = 1, n.roll = roll, cluster = cluster)
	forecastMu = t(apply(forc@mforecast$mu, 3, FUN = function(x) x))
	colnames(forecastMu) = nam
	rownames(forecastMu) = as.character(indexf)
	forecastCov = array(unlist(rcov(forc)), dim=c(m,m,roll+1), dimnames=list(nam, nam, as.character(indexf)))
	xseed = simRes = simMu = vector(mode = "list", length = roll+1)
	rcovfilt = rcorfilt = rcovsim = rcorsim = array(NA, dim = c(m,m,roll+1))
	if(spec@model$modelinc[1]>0) varcoef = fit@model$varcoef else varcoef = NULL
	model = fit@model
	umodel  = model$umodel
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, 
			umodel$modeldesc$vsubmodel, umodel$modeldata$mexdata, 
			umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, NULL)
	
	midx = fit@model$midx
	mpars = fit@model$mpars
	for(i in 1:m){
		setfixed(mspec@spec[[i]])<-as.list(mpars[midx[,i]==1, i])
	}
	dccfix = as.list(coef(fit, "dcc"))
	specf = dccspec(uspec = mspec,  VAR = ifelse(spec@model$modelinc[1]>0, TRUE, FALSE), 
			lag = spec@model$modelinc[1], dccOrder = model$modelinc[3:4], 
			model = spec@model$DCC, groups = spec@model$fdccindex, 
			distribution = model$modeldesc$distribution, fixed.pars = dccfix)
	filt = dccfilter(specf, data = Data, out.sample = 0, filter.control = list(n.old = T),
			cluster = cluster, varcoef = varcoef, realizedVol = realizedVol)
	sig = sigma(filt)
	res = residuals(filt)
	C = rcor(filt, type = "Q")
	Z = filt@mfilter$stdresid
	if(!is.null(rseed)){
		rseed = as.numeric(rseed)
		if(length(rseed) == (roll+1)) zseed = rseed
		if(length(rseed) >  (roll+1)) zseed = rseed[1:(roll+1)]
		if(length(rseed) <  (roll+1)) zseed = rseed[1]+seq_len(roll)
	} else{
		zseed = as.integer( runif( roll+1, 0, as.integer(Sys.time()) ) )
	}
	for(i in 1:(roll+1)){
		preQ = C[,,(T+i-1)]
		presigma = sig[(T-p+i):(T-1+i), , drop=FALSE]
		if(modv=="realGARCH") prerealized = coredata(realizedVol)[(T-p+i):(T-1+i), , drop=FALSE] else prerealized = NULL
		preresiduals = res[(T-p+i):(T-1+i), , drop=FALSE]
		prereturns = as.matrix(Data[(T+i-p):(T+i-1),], ncol = m)
		preZ = matrix(Z[(T-p+i):(T-1+i),], ncol = m)
		fsim = dccsim(fitORspec = fit, n.sim = 1, n.start = 0, m.sim = sim, startMethod="sample",
				presigma = presigma, preresiduals = preresiduals, 
				prereturns = prereturns, preQ = preQ, preZ = preZ, Qbar = fit@mfit$Qbar, 
				rseed = zseed[i], cluster = cluster, prerealized = prerealized)
		simMu[[i]] = t(sapply(fsim@msim$simX, FUN = function(x) x))
		simZ = sapply(fsim@msim$simZ, FUN = function(x) tail(x,1))		
		HH = array(unlist(fsim@msim$simH), dim=c(m,m,sim))
		simRes[[i]] = .Call("Cov2Res", HH, simZ, as.integer(c(sim, m)), PACKAGE="rmgarch")
		xseed[[i]] = fsim@msim$rseed
		# simMu - simRes == forecastMu
	}
	if(save.output){
		olddir = getwd()
		setwd(save.dir)
		scenario = list(simMu = simMu, simRes = simRes, forecastMu = forecastMu, 
				forecastCov = forecastCov, rseed = xseed, sim = sim, roll = roll, 
				model = "dcc", index = indexf)
		eval(parse(text = paste("save(scenario, file = '",save.name,".rda')",sep="")))		
		setwd(olddir)
		return(1)
	} else{
		return(list(simMu = simMu, simRes = simRes, forecastMu = forecastMu, 
						forecastCov = forecastCov, rseed = xseed, index = indexf))
	}
}

scenario.cgarch = function(Data, spec, spd.control, sim, roll, solver = "solnp", 
		fit.control = list(eval.se = FALSE), solver.control = list(), 
		cluster = NULL, rseed = NULL, save.output = FALSE, save.dir = getwd(), 
		save.name = NULL, debug = FALSE, ...)
{
	pars = list(...)
	modv = spec@spec[[1]]@model$modeldesc$vmodel
	if(is.null(pars$fit)) xfit = NULL else xfit = pars$fit
	if(is.null(pars$VAR.fit)) VAR.fit = NULL else VAR.fit = pars$VAR.fit
	if(modv=="realGARCH"){
		if(is.null(pars$realizedVol)){
			stop("\nfscenario-->error: realizedVol required for realGARCH model.\n")
		} else{
			realizedVol = pars$realizedVol
		}
	} else{
		realizedVol = NULL
	}
	fit = cgarchfit(spec, data = Data, spd.control = spd.control, 
			out.sample = roll, solver = solver, 
			solver.control = solver.control, fit.control = fit.control, 
			cluster = cluster, fit = xfit, VAR.fit = VAR.fit, realizedVol = realizedVol, ...)
	fsim = vector(mode="list", length = roll)
	p = spec@model$maxgarchOrder
	T = fit@model$modeldata$T
	T0 = fit@model$modeldata$index[(T):(T+roll)]
	indexf = c(T0[-1], generatefwd(tail(T0, 1), length.out=1, by = fit@model$modeldata$period))
	xseed = simRes = simMu = vector(mode = "list", length = roll+1)
	m  = NCOL(Data)
	rcovfilt = rcorfilt = rcovsim = rcorsim = array(NA, dim = c(m,m,roll+1))
	if(spec@model$modelinc[1]>0) varcoef = fit@model$varcoef else varcoef = NULL
	
	model = fit@model
	umodel  = model$umodel
	mspec = .makemultispec(umodel$modelinc, umodel$modeldesc$vmodel, umodel$modeldesc$vsubmodel, 
			umodel$modeldata$mexdata, umodel$modeldata$vexdata, umodel$start.pars, 
			umodel$fixed.pars, NULL)
	
	midx = fit@model$midx
	mpars = fit@model$mpars
	for(i in 1:m){
		setfixed(mspec@spec[[i]])<-as.list(mpars[midx[,i]==1, i])
	}
	dccfix = as.list(coef(fit, "dcc"))
	specf = cgarchspec(uspec = mspec,  VAR = ifelse(spec@model$modelinc[1]>0, TRUE, FALSE), 
			lag = spec@model$modelinc[1], dccOrder = model$modelinc[4:5], 
			asymmetric = ifelse(spec@model$modelinc[6]>0, TRUE, FALSE), 
			distribution.model = list(copula = model$modeldesc$distribution, 
					method = model$modeldesc$cor.method, 
					time.varying = model$modeldesc$timecopula,
					transformation = model$modeldesc$transformation), 
			fixed.pars = dccfix)
	# we really only need to use this if roll>0
	filt = cgarchfilter(specf, data = Data, out.sample = 0, filter.control = list(n.old = T),
			spd.control = spd.control, cluster = cluster, varcoef = varcoef,
			realizedVol = realizedVol)
	if(roll>0){
		forecastMu = rbind(tail(filt@model$mu, roll), matrix(NA, ncol = m))
	} else{
		forecastMu = matrix(NA, ncol = m)
	}
	sig = coredata(sigma(filt))
	res = coredata(residuals(filt))
	R = rcor(filt)
	Z = filt@mfilter$stdresid
	Q = filt@mfilter$Qt
	if(!is.null(rseed)){
		rseed = as.integer(rseed)
		if(length(rseed)==(roll+1)) zseed = rseed
		if(length(rseed)>(roll+1)) zseed = rseed[1:(roll+1)]
		if(length(rseed)<(roll+1)) zseed = rseed[1]+seq_len(roll)
	} else{
		zseed = as.integer( runif( roll+1, 0, as.integer(Sys.time()) ) )
	}
	# distinguish between static and time varying copula
	if(spec@model$modeldesc$timecopula){
		for(i in 1:(roll+1)){
			preR = R[,,(T+i-1)]
			# extra check
			diag(preR) = 1
			preQ = Q[[(T+i-1)]]
			presigma = sig[(T-p+i):(T-1+i),, drop=FALSE]
			if(modv=="realGARCH") prerealized = coredata(realizedVol)[(T-p+i):(T-1+i), , drop=FALSE] else prerealized = NULL
			preresiduals = res[(T-p+i):(T-1+i),, drop=FALSE]
			prereturns = as.matrix(Data[(T-p+i):(T-1+i),], ncol = m)
			preZ = matrix(Z[(T-p+i):(T-1+i),], ncol = m)
			fsim = cgarchsim(fit, n.sim = 1, n.start = 0, m.sim = sim, 
					startMethod = "sample", presigma = presigma, 
					preresiduals = preresiduals, prereturns = prereturns, 
					preR = preR, preQ = preQ, preZ = preZ, mexsimdata = NULL, 
					vexsimdata = NULL, cluster = cluster, rseed = zseed[i],
					prerealized = prerealized)
			simMu[[i]] = t(sapply(fsim@msim$simX, FUN = function(x) x))
			simZ = t(apply(fsim@msim$simZ, 3, FUN = function(x) tail(x, 1)))
			HH = array(unlist(fsim@msim$simH), dim=c(m,m,sim))
			simRes[[i]] = .Call("Cov2Res", HH, simZ, as.integer(c(sim, m)), PACKAGE="rmgarch")
			xseed[[i]] = fsim@msim$rseed
		}
	} else{
		for(i in 1:(roll+1)){
			preR = R
			# extra check
			diag(preR) = 1
			presigma = sig[(T-p+i):(T-1+i),, drop=FALSE]
			if(modv=="realGARCH") prerealized = coredata(realizedVol)[(T-p+i):(T-1+i), , drop=FALSE] else prerealized = NULL
			
			preresiduals = res[(T-p+i):(T-1+i),, drop=FALSE]
			prereturns = as.matrix(Data[(T-p+i):(T-1+i),], ncol = m)
			preZ = matrix(Z[(T-p+i):(T-1+i),], ncol = m)
			fsim = cgarchsim(fit, n.sim = 1, n.start = 0, m.sim = sim, 
					startMethod = "sample", presigma = presigma, 
					preresiduals = preresiduals, prereturns = prereturns, 
					preR = preR, preQ = preQ, preZ = preZ, mexsimdata = NULL, 
					vexsimdata = NULL, cluster = cluster, rseed = zseed[i],
					prerealized = prerealized)
			simMu[[i]] = t(sapply(fsim@msim$simX, FUN = function(x) x))
			simZ = t(apply(fsim@msim$simZ, 3, FUN = function(x) tail(x, 1)))
			HH = array(unlist(fsim@msim$simH), dim=c(m,m,sim))
			simRes[[i]] = .Call("Cov2Res", HH, simZ, as.integer(c(sim, m)), PACKAGE="rmgarch")
			xseed[[i]] = fsim@msim$rseed
		}
	}
	if(save.output){
		olddir = getwd()
		setwd(save.dir)
		scenario = list(simMu = simMu, simRes = simRes, forecastMu = forecastMu, 
				rseed = xseed, sim = sim, roll = roll, model = "cgarch", index = indexf)
		eval(parse(text = paste("save(scenario, file = '",save.name,".rda')",sep="")))		
		setwd(olddir)
		return(1)
	} else{
		return(list(simMu = simMu, simRes = simRes, forecastMu = forecastMu,
						rseed = xseed, index = indexf))
	}
}


scenario.var = function(Data, sim, roll, model = list(lag = 1, lag.max = NULL, 
				lag.criterion = c("AIC", "HQ", "SC", "FPE"), robust = FALSE, 
				robust.control = list("gamma" = 0.25, "delta" = 0.01, "nc" = 10, "ns" = 500)),
		cov.method = c("ML", "LW", "EWMA", "MVE", "MCD", "MVT", "BS"),
		cov.options = list(shrinkage=-1, lambda = 0.96),
		rseed = NULL, save.output = FALSE, save.dir = getwd(), save.name = NULL)
{
	if(is.null(rseed)) rseed = as.integer(Sys.time())
	rseed = as.numeric(rseed)
	Data = as.matrix(Data)
	n = dim(Data)[1]
	m = dim(Data)[2]
	T = n - roll
	X = Data[1:T, , drop = FALSE]
	if( !is.null(model$lag.max) ){
		ic = model$lag.criterion[1]
		mp = .varxselect(y = X, lag.max = model$lag.max, exogen = NULL)$selection
		mp = as.numeric( mp[which(substr(names(mp), 1, 2) == substr(ic, 1, 2))] )
	} else {
		mp = max(1, model$lag)
	}
	rb = model$robust.control
	fit = varxfit(X, p = mp, exogen = NULL, robust = model$robust, 
			gamma = rb$gamma, delta = rb$delta, nc = rb$nc, ns = rb$ns, 
			postpad = c("none", "constant", "zero", "NA")[2])
	varcoef = fit$Bcoef
	forc = varxforecast(Data, p = mp, Bcoef = fit$Bcoef, out.sample = roll, 
			n.ahead = 1, n.roll = roll, mregfor = NULL)
	filt = varxfilter(Data, p = mp, Bcoef = fit$Bcoef, postpad = "constant", 
			exogen = NULL)
	
	forecastMu = t(apply(forc, 3, FUN = function(x) x))
	forecastCov = array(NA, dim = c(m,m,roll+1))
	res = fit$xresiduals
	preret = tail(X, mp)
	preres = tail(res, mp)
	xseed = simRes = simMu = vector(mode = "list", length = roll+1)
	if(!is.null(rseed)){
		rseed = as.numeric(rseed)
		if(length(rseed)==(roll+1)) zseed = rseed
		if(length(rseed)>(roll+1)) zseed = rseed[1:(roll+1)]
		if(length(rseed)<(roll+1)) zseed = rseed[1]+seq_len(roll)
	} else{
		zseed = as.integer( runif( roll+1, 0, as.integer(Sys.time()) ) )
	}
	for(i in 1:(roll+1)){
		if(cov.method[1]!="ML"){
			fS = fun.cov(Data = res, method = cov.method, 
					shrinkage=cov.options$shrinkage, 
					lambda = cov.options$lambda, 
					demean = FALSE)$S
		} else{
			fS = (t(res) %*% (res))/T
		}
		set.seed(zseed[i])
		xseed[[i]] = zseed
		Z = .rmvnorm(sim, mean = rep(0, m), sigma = fS)
		tmp = varxsimXX(X = Data[1:(T+1-i), ], Bcoef = varcoef, p = mp, 
				m.sim = sim, prereturns = preret, resids = Z, NULL)
		simMu[[i]] = tmp
		simRes[[i]] = Z
		forecastCov[,,i] = fS
		if(i<(roll+1)){ 
			preres = filt$xresiduals[(T+i-mp+1):(T+i), , drop = FALSE]
			preret = Data[(T+i-mp+1):(T+i), , drop = FALSE]
			res = filt$xresiduals[1:(T+i), ]
		}
	}
	if(save.output){
		olddir = getwd()
		setwd(save.dir)
		scenario = list(simMu = simMu, simRes = simRes, forecastMu = forecastMu, 
				forecastCov = forecastCov, rseed = xseed, sim = sim, roll = roll, 
				model = "var")
		eval(parse(text = paste("save(scenario, file = '",save.name,".rda')",sep="")))
		setwd(olddir)
		return(1)
	} else{
		return(list(simMu = simMu, simRes = simRes, forecastMu = forecastMu, 
						rseed = xseed))		
	}
}

scenario.mdist = function(Data, sim, roll, model, rseed, save.output = FALSE, 
		save.dir = getwd(), save.name = NULL)
{
	stop("\nNot yet implemented")
}


goget = function(object, arg){
	UseMethod("goget")
}

.gogetscen <- function(object, arg){
	valid.choices = "scenario"
	tmp = match.arg(tolower(arg[1]), valid.choices)
	scenario = object@scenario$simMu
	ans = NA
	eval(parse(text = paste("ans=",tmp,sep="")))
	return(ans)
}

setMethod("goget", signature(object = "fScenario"), definition = .gogetscen)



goload  = function(object, ...)
{
	UseMethod("goload")
}

.goload1 = function(object){
	if(object@model$save.output){
		oldir = getwd()
		scenario = list()
		setwd(object@model$save.dir)
		eval(parse(text = paste("load(file='",object@model$save.name, ".rda')", sep = "")))
		
		object@model$sim = scenario$sim
		object@model$roll = scenario$roll
		object@model$model = scenario$model
		scenario$model = NULL
		scenario$roll = NULL
		scenario$sim = NULL
		object@scenario = scenario
		object@model$save.output = FALSE
		object@model$assets = dim(scenario$simMu[[1]])[2]
		if(object@model$asset.names[1]=="") object@model$asset.names = paste("A_",1:object@model$assets,sep="")
		setwd(oldir)
	}
	return(object)
}
setMethod("goload", signature(object = "fScenario"), .goload1)



.goload2 = function(object){
	if(object@model$save.output){
		oldir = getwd()
		mom = list()
		setwd(object@model$save.dir)
		eval(parse(text = paste("load(file='",object@model$save.name, ".rda')", sep = "")))
		object@model$n.ahead = mom$n.ahead
		object@model$roll = mom$roll
		object@model$model = mom$model
		mom$model = NULL
		mom$roll = NULL
		mom$n.ahead = NULL
		object@moments = mom
		object@model$assets = dim(mom[[1]])[2]
		object@model$save.output = FALSE
		if(object@model$asset.names[1]=="") object@model$asset.names = paste("A_",1:object@model$assets,sep="")
		setwd(oldir)
	}
	return(object)
}
setMethod("goload", signature(object = "fMoments"), .goload2)

# covariance functions
fun.cov = function(Data, method = c("LW", "EWMA", "MVE", "MCD", "MVT", "BS"), 
		shrinkage=-1, lambda = 0.96, demean = FALSE){
	ans = switch(tolower(method),
			lw = lw.cov(Data, shrinkage, demean),
			ewma = ewma.cov(Data, lambda, demean),
			mve = rob.cov(Data, "mve", demean),
			mcd = rob.cov(Data, "mcd", demean),
			mvt = rob.cov(Data, "mcd", demean),
			bs = bs.cov(Data))
	return( ans )
}

# Ledoit and Wolfe Shrinkage Estimator
lw.cov = function(X, shrink = -1, demean = TRUE){
	n = dim(X)[1]
	m = dim(X)[2]
	meanx = colMeans(X)
	if(demean) X = scale(X, scale = FALSE)
	# compute sample covariance matrix
	samplecov = (t(X)%*%X)/n
	# compute prior
	meanvar = mean(diag(samplecov))
	prior = meanvar * diag(m)
	if(shrink == -1){
		# compute shrinkage parameters
		# p in paper 
		y = X^2
		phiMat = (t(y)%*%y)/n - 2*(t(X)%*%X)*samplecov/n + samplecov^2
		phi = sum(apply(phiMat, 1, "sum"))
		# c in paper
		cgamma = norm(samplecov-prior,'F')^2
		# shrinkage constant
		kappa = phi/cgamma
		shrinkage = max(0, min(1, kappa/n))
	} else{
		shrinkage = shrink
	}
	sigma = shrinkage * prior + (1 - shrinkage)*samplecov
	return(list(S = sigma, par = shrinkage))
}

ewma.cov = function(X, lambda = 0.96, demean = TRUE)
{
	n = dim(X)[1]
	i = (0:(n-1))
	ewma.wt = lambda^i
	ewma.wt = ewma.wt/sum(ewma.wt)
	covm = cov.wt(X, wt = rev(ewma.wt), center = demean)$cov
	return(list(S = covm, par = lambda))
}

rob.cov = function(X, method = c("mve", "mcd", "mvt"), demean = TRUE){
	S = switch(tolower(method), 
			mve = cov.rob(x = X, method = "mve"),
			mcd = cov.rob(x = X, method = "mcd"),
			mvt = cov.trob(x = X, center = demean))
	return(list(S = S$cov, par = NA))
}

bs.cov = function(X)
{
	# Bayes Stein estimator
	# Alexios Ghalanos 2008
	# This function encapsulates an example of shrinking the returns 
	#   and covariance using Bayes-Stein shrinkage as described in 
	#   Jorion, 1986
	mu = as.matrix(apply(X, 2, FUN = function(x) mean(x)))
	S = cov(X)
	m = dim(X)[2]
	n = dim(X)[1]
	one = as.matrix(rep(1, m))
	a = solve(S, one)
	# Constant non informative prior
	mu.prior = one * as.numeric(t(mu) %*% a/t(one) %*% a)
	S.inv = solve(S)
	d = t(mu - mu.prior) %*% S.inv %*% (mu - mu.prior)
	d = as.numeric(d)
	lambda = (m + 2) / d
	w = lambda / (n+lambda)
	mu.pred = (1 - w) * mu + w * mu.prior
	
	wc1 = 1 / (n + lambda)
	wc2 = lambda*(n-1) / (n*(n+lambda)*(n - m - 2))
	wc2 = wc2 / as.numeric(t(one) %*% a)
	V.post = wc1 * S + wc2 * one %*% t(one)
	V.pred = S + V.post
	sigma.post = sqrt(diag(V.post))
	sigma.pred = sqrt(diag(V.pred))
	
	return(list(S = V.pred, par = mu.pred[,1]))
}

###############################################################################
fmoments = function(spec, Data, n.ahead = 1, roll  = 0, 
		solver = "solnp", solver.control = list(), 
		fit.control = list(eval.se=FALSE), 
		cluster = NULL, save.output = FALSE, save.dir = getwd(), 
		save.name = paste("M", sample(1:1000, 1), sep= ""), ...)
{
	UseMethod("fmoments")
}

setMethod("fmoments", definition = .fmoments)

fscenario = function(Data, sim = 1000, roll = 0, 
		model = c("gogarch", "dcc", "cgarch", "var", "mdist"), spec = NULL,
		var.model = list(
				lag = 1, lag.max = NULL, 
				lag.criterion = c("AIC", "HQ", "SC", "FPE"), 
				robust = FALSE, 
				robust.control = list("gamma" = 0.25, "delta" = 0.01, 
						"nc" = 10, "ns" = 500)),
		mdist.model = list(distribution = c("mvn", "mvt", "manig"), AR = TRUE, 
				lag = 1),
		spd.control = list(lower = 0.1, upper = 0.9, type = "pwm", 
				kernel = "epanech"),
		cov.method = c("ML", "LW", "EWMA", "MVE", "MCD", "MVT", "BS"),
		cov.options = list(shrinkage=-1, lambda = 0.96),
		solver = "solnp", solver.control = list(), 
		fit.control = list(eval.se=FALSE), cluster =  NULL, 
		save.output = FALSE, save.dir = getwd(), 
		save.name = paste("S", sample(1:1000, 1), sep = ""), 
		rseed  = NULL, ...)
{
	UseMethod("fscenario")
}

setMethod("fscenario", definition = .fscenario)


.fitted.fscenario = function(object){
	if(length(object@scenario)>0){
		roll = object@model$roll+1
		m = NCOL(object@scenario$simMu[[1]])
		n = object@model$sim
		index = as.character(object@scenario$index)
		nam = object@model$asset.names
		ans = array(data = unlist(object@scenario$simMu), dim=c(n, m, roll), dimnames=list(NULL, nam, index))
		#if(object@model$model!="cgarch") attr(ans, "forecast")<-object@scenario$forecastMu
		return(ans)
	} else{
		stop("\nrmgarch-->error: scenario slot not populated with data.")
	}
}
setMethod("fitted", signature(object = "fScenario"), .fitted.fscenario)

.fitted.fmoments = function(object){
	if(length(object@moments)>0){
		return(object@moments$forecastMu)
	} else{
		stop("\nrmgarch-->error: moments slot not populated with data.")
	}
}
setMethod("fitted", signature(object = "fMoments"), .fitted.fmoments)

.rcov.fmoments = function(object){
	if(length(object@moments)>0){
		return(object@moments$forecastCov)
	} else{
		stop("\nrmgarch-->error: moments slot not populated with data.")
	}
}
setMethod("rcov", signature(object = "fMoments"), .rcov.fmoments)

.rcoskew.fmoments = function(object){
	if(length(object@moments)>0){
		if(!is.null(object@moments$forecastM3)){
			return(object@moments$forecastM3)
		} else{
			stop("\nrmgarch-->error: model does not return a coskewness matrix.")
		}
	} else{
		stop("\nrmgarch-->error: moments slot not populated with data.")
	}
}
setMethod("rcoskew", signature(object = "fMoments"), .rcoskew.fmoments)

.rcokurt.fmoments = function(object){
	if(length(object@moments)>0){
		if(!is.null(object@moments$forecastM4)){
			return(object@moments$forecastM4)
		} else{
			stop("\nrmgarch-->error: model does not return a cokurtosis matrix.")
		}
	} else{
		stop("\nrmgarch-->error: moments slot not populated with data.")
	}
}
setMethod("rcokurt", signature(object = "fMoments"), .rcokurt.fmoments)