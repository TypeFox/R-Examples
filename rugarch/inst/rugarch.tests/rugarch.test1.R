#################################################################################
##
##   R package rugarch by Alexios Ghalanos Copyright (C) 2008-2013.
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


rugarch.test1a = function(cluster = NULL){
	tic = Sys.time()
	# simulated parameter distribution
	spec = arfimaspec( mean.model = list(armaOrder = c(2,2), include.mean = TRUE, arfima = FALSE), 
			distribution.model = "norm", fixed.pars = list(ar1=0.6, ar2=0.21, ma1=-0.7, ma2=0.3, mu = 0.02,
					sigma = 0.02))
	dist = arfimadistribution(spec, n.sim = 2000, n.start = 100, m.sim = 100, recursive = TRUE, 
			recursive.length = 10000, recursive.window = 1000, cluster = cluster)
	
	save(dist, file = "test1a.rda")
	options(width=150)
	zz <- file("test1a.txt", open="wt")
	sink(zz)
	# slots:
	slotNames(dist)
	# methods:
	# summary
	show(dist)
	# as.data.frame(...., window, which=c("rmse", "stats", "coef", "coefse"))
	# default
	as.data.frame(dist)
	as.data.frame(dist, window = 1, which = "rmse")
	as.data.frame(dist, window = 1, which = "stats")
	as.data.frame(dist, window = 1, which = "coef")
	as.data.frame(dist, window = 1, which = "coefse")
	as.data.frame(dist, window = 8, which = "rmse")
	as.data.frame(dist, window = 8, which = "stats")
	as.data.frame(dist, window = 8, which = "coef")
	as.data.frame(dist, window = 8, which = "coefse")
	sink(type="message")
	sink()
	close(zz)
	
	# create some plots
	nwindows = dist@dist$details$nwindows
	# 2000/3000/4000/5000/6000/7000/8000/9000/10000
	
	# expected reduction factor in RMSE for sqrt(N) consistency
	expexcted.rmsegr = sqrt(2000/seq(3000,10000,by=1000))
	
	# actual RMSE reduction
	actual.rmsegr = matrix(NA, ncol = 8, nrow = 6)
	rownames(actual.rmsegr) = c("mu", "ar1", "ar2", "ma2", "ma2", "sigma")
	# start at 2000 (window 1)
	rmse.start = as.data.frame(dist, window = 1, which = "rmse")
	for(i in 2:nwindows) actual.rmsegr[,i-1] = as.numeric(as.data.frame(dist, window = i, which = "rmse")/rmse.start)
	postscript("test1a.eps", bg = "white", width = 800, height = 800)
	par(mfrow = c(2,3))
	for(i in 1:6){
		plot(seq(3000,10000,by=1000), actual.rmsegr[i,], type = "l", lty = 2, ylab = "RMSE Reduction", xlab = "N (sim)", 
				main = rownames(actual.rmsegr)[i])
		lines(seq(3000,10000,by=1000), expexcted.rmsegr, col = 2)
		legend("topright", legend = c("Actual", "Expected"), col = 1:2, bty = "m", lty = c(2,1))
	}
	dev.off()
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}


rugarch.test1b = function(cluster=NULL){
	# fit/filter
	tic = Sys.time()
	data(sp500ret)	
	fit = vector(mode = "list", length = 9)
	dist = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig", "ghyp", "jsu")
	for(i in 1:9){
		spec = arfimaspec( mean.model = list(armaOrder = c(1,1), include.mean = TRUE, arfima = FALSE), 
				distribution.model = dist[i])
		fit[[i]] = arfimafit(spec = spec, data = sp500ret, solver = "solnp", fit.control = list(scale = 1))
	}
	cfmatrix = matrix(NA, nrow = 9, ncol = 7)
	colnames(cfmatrix) = c("mu", "ar1", "ma1", "sigma", "skew", "shape", "ghlambda")
	rownames(cfmatrix) = dist
	
	for(i in 1:9){
		cf = coef(fit[[i]])
		cfmatrix[i, match(names(cf), colnames(cfmatrix))] =  cf
	}
	sk = ku = rep(0, 9)
	for(i in 1:9){
		cf = coef(fit[[i]])
		if(fit[[i]]@model$modelinc[16]>0) sk[i] = dskewness(distribution = dist[i], 
					skew = cf["skew"], shape = cf["shape"], lambda = cf["ghlambda"])		
		if(fit[[i]]@model$modelinc[17]>0) ku[i] = dkurtosis(distribution = dist[i], 
					skew = cf["skew"], shape = cf["shape"], lambda = cf["ghlambda"])
	}
	hq = sapply(fit, FUN = function(x) infocriteria(x)[4])
	cfmatrix = cbind(cfmatrix, sk, ku, hq)
	colnames(cfmatrix) = c(colnames(cfmatrix[,1:7]), "skewness", "ex.kurtosis","HQIC")
	
	
	# filter the data to check results
	filt = vector(mode = "list", length = 9)
	for(i in 1:9){
		spec = arfimaspec( mean.model = list(armaOrder = c(1,1), include.mean = TRUE, arfima = FALSE), 
				distribution.model = dist[i])
		setfixed(spec) = as.list(coef(fit[[i]]))
		filt[[i]] = arfimafilter(spec = spec, data = sp500ret)
	}
	
	options(width = 120)
	zz <- file("test1b.txt", open="wt")
	sink(zz)
	print(cfmatrix, digits = 4)
	cat("\nARFIMAfit and ARFIMAfilter residuals check:\n")
	print(head(sapply(filt, FUN = function(x) residuals(x))) == head(sapply(fit, FUN = function(x) residuals(x))))
	cat("\ncoef method:\n")
	print(cbind(coef(filt[[1]]), coef(fit[[1]])))
	cat("\nfitted method:\n")
	print(cbind(head(fitted(filt[[1]])), head(fitted(fit[[1]]))))
	cat("\ninfocriteria method:\n")
	# For filter, it assumes estimation of parameters else does not make sense!
	print(cbind(infocriteria(filt[[1]]), infocriteria(fit[[1]])))
	cat("\nlikelihood method:\n")
	print(cbind(likelihood(filt[[1]]), likelihood(fit[[1]])))
	cat("\nresiduals method:\n")
	print(cbind(head(residuals(filt[[1]])), head(residuals(fit[[1]]))))
	cat("\nuncmean method:\n")
	print(cbind(uncmean(filt[[1]]), uncmean(fit[[1]])))
	cat("\nuncmean method (by simulation):\n")
	# For spec and fit
	spec = arfimaspec( mean.model = list(armaOrder = c(1,1), include.mean = TRUE, 
					arfima = FALSE), distribution.model = dist[1])
	setfixed(spec) = as.list(coef(fit[[1]]))
	print(cbind(uncmean(spec, method = "simulation", n.sim = 100000, rseed = 100), 
					uncmean(fit[[1]], method = "simulation", n.sim = 100000, rseed = 100)))
	cat("\nsummary method:\n")
	print(show(filt[[1]]))
	print(show(fit[[1]]))
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

rugarch.test1c = function(cluster=NULL){
	# unconditional forecasting
	tic = Sys.time()
	require(xts)
	data(sp500ret)	
	fit = vector(mode = "list", length = 9)
	dist = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig", "ghyp", "jsu")
	for(i in 1:9){
		spec = arfimaspec( mean.model = list(armaOrder = c(1,1), include.mean = TRUE, arfima = FALSE), 
				distribution.model = dist[i])
		fit[[i]] = arfimafit(spec = spec, data = sp500ret, solver = "solnp", fit.control = list(scale = 1))
	}
	cfmatrix = matrix(NA, nrow = 9, ncol = 7)
	colnames(cfmatrix) = c("mu", "ar1", "ma1", "sigma", "skew", "shape", "ghlambda")
	rownames(cfmatrix) = dist
	
	for(i in 1:9){
		cf = coef(fit[[i]])
		cfmatrix[i, match(names(cf), colnames(cfmatrix))] =  cf
	}
	
	umean = rep(0, 9)
	for(i in 1:9){
		umean[i] = uncmean(fit[[i]])
	}
	
	forc = vector(mode = "list", length = 9)
	for(i in 1:9){
		forc[[i]] = arfimaforecast(fit[[i]], n.ahead = 100)
	}
	
	lmean40 = sapply(forc, FUN = function(x) as.numeric(fitted(x)[40,1]))
	cfmatrix1 = cbind(cfmatrix, umean, lmean40)
	colnames(cfmatrix1) = c(colnames(cfmatrix1[,1:7]), "uncmean", "forecast40")
	
	# forecast with spec to check results
	forc2 = vector(mode = "list", length = 9)
	for(i in 1:9){
		spec = arfimaspec( mean.model = list(armaOrder = c(1,1), include.mean = TRUE, arfima = FALSE), 
				distribution.model = dist[i])
		setfixed(spec) = as.list(coef(fit[[i]]))
		forc2[[i]] = arfimaforecast(spec, data = sp500ret, n.ahead = 100)
	}
	lmean240 = sapply(forc2, FUN = function(x) as.numeric(fitted(x)[40,1]))
	cfmatrix2 = cbind(cfmatrix, umean, lmean240)
	colnames(cfmatrix2) = c(colnames(cfmatrix2[,1:7]), "uncmean", "forecast40")
	
	# Test Methods  on object
	
	options(width = 120)
	zz <- file("test1c.txt", open="wt")
	sink(zz)
	cat("\nARFIMAforecast from ARFIMAfit and ARFIMAspec check:")
	cat("\nFit\n")	
	print(cfmatrix1, digits = 4)
	cat("\nSpec\n")	
	print(cfmatrix2, digits = 4)
	slotNames(forc[[1]])
	# summary
	print(show(forc[[1]]))
	sink(type="message")
	sink()
	close(zz)
	
	nforc = sapply(forc, FUN = function(x) t(as.numeric(fitted(x))))
	postscript("test1c.eps", width = 12, height = 5)
	# generate FWD dates:
	dx = as.POSIXct(tail(rownames(sp500ret),50)) 
	df = generatefwd(tail(dx, 1), length.out = 100+1, by = forc[[1]]@model$modeldata$period)[-1]		
	dd = c(dx, df)
	clrs = rainbow(9, alpha = 1, start = 0.4, end = 0.95)
	plot(xts(c(tail(sp500ret[,1], 50), nforc[,1]), dd), col = "lightgrey", ylim=c(-0.02, 0.02),
			ylab = "", xlab = "", main = "100-ahead Unconditional Forecasts")
	addLegend("topright", dist, col = clrs, fill = clrs, on=1)
	for(i in 1:9){
		tmp = c(tail(sp500ret[,1], 50), rep(NA, 100))
		tmp[51:150] = nforc[1:100,i]
		lines(xts(c(rep(NA, 50), tmp[-(1:50)]),dd), col = clrs[i], on=1)
	}
	dev.off()
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

rugarch.test1d = function(cluster=NULL){
	# rolling forecast
	tic = Sys.time()
	
	data(sp500ret)
	fit = vector(mode = "list", length = 9)
	dist = c("norm", "snorm", "std", "sstd", "ged", "sged", "nig", "ghyp", "jsu")
	for(i in 1:9){
		spec = arfimaspec( mean.model = list(armaOrder = c(1,1), include.mean = TRUE, arfima = FALSE), 
				distribution.model = dist[i])
		fit[[i]] = arfimafit(spec = spec, data = sp500ret, solver = "solnp", 
				out.sample = 1000, fit.control = list(scale = 1))
	}
	cfmatrix = matrix(NA, nrow = 9, ncol = 7)
	colnames(cfmatrix) = c("mu", "ar1", "ma1", "sigma", "skew", "shape", "ghlambda")
	rownames(cfmatrix) = dist
	
	for(i in 1:9){
		cf = coef(fit[[i]])
		cfmatrix[i, match(names(cf), colnames(cfmatrix))] =  cf
	}
	
	
	forc = vector(mode = "list", length = 9)
	for(i in 1:9){
		forc[[i]] = arfimaforecast(fit[[i]], n.ahead = 1, n.roll = 999)
	}
	rollforc = sapply(forc, FUN = function(x) t(fitted(x)))
	
	# forecast performance measures:
	fpmlist = vector(mode = "list", length = 9)
	for(i in 1:9){
		fpmlist[[i]] = fpm(forc[[i]], summary = FALSE)
	}
	
	postscript("test1d.eps", width = 16, height = 5)
	dd = as.POSIXct(tail(rownames(sp500ret), 1250)) 
	clrs = rainbow(9, alpha = 1, start = 0.4, end = 0.95)
	plot(xts(tail(sp500ret[,1], 1250), dd), type = "l", ylim = c(-0.02, 0.02), 
			col = "lightgrey", ylab = "", xlab = "", 
			main = "Rolling 1-ahead Forecasts\nvs Actual", minor.ticks=FALSE, 
			auto.grid=FALSE)
	addLegend("topleft", dist, col = clrs, fill = clrs, on=1)
	for(i in 1:9){
		tmp = tail(sp500ret[,1], 1250)
		tmp[251:1250] = rollforc[1:1000,i]
		lines(xts(c(rep(NA, 250), tmp[-(1:250)]), dd), col = clrs[i], on=1)
	}
	
	# plot deviation measures and range
	tmp = vector(mode = "list", length = 9)
	for(i in 1:9){
		tmp[[i]] = fpmlist[[i]][,"AE"]
		names(tmp[[i]]) = dist[i]
	}
	boxplot(tmp, col = clrs, names = dist, range  = 6, notch = TRUE, 
			main = "Rolling 1-ahead Forecasts\nAbsolute Deviation Loss")
	dev.off()
	
	# fpm comparison
	compm = matrix(NA, nrow = 3, ncol = 9)
	compm = sapply(fpmlist, FUN = function(x) c(mean(x[,"SE"]), mean(x[,"AE"]), mean(x[,"DAC"])))
	colnames(compm) = dist
	rownames(compm) = c("MSE", "MAD", "DAC")
	
	zz <- file("test1d.txt", open="wt")
	sink(zz)
	cat("\nRolling Forecast FPM\n")
	print(compm, digits = 4)
	cat("\nMethods Check\n")
	print(fitted(forc[[1]])[,1:10,drop=FALSE])
	print(fpm(forc[[1]], summary = TRUE))
	print(show(forc[[1]]))
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

rugarch.test1e = function(cluster=NULL){
	# Multi-Methods
	tic = Sys.time()
	
	data(dji30ret)
	Dat = dji30ret[, 1:3, drop = FALSE]
	
	#------------------------------------------------
	# Unequal Spec
	# Fit
	spec1 = arfimaspec(mean.model = list(armaOrder = c(2,1)))
	spec2 = arfimaspec(mean.model = list(armaOrder = c(2,2)))
	spec3 = arfimaspec(mean.model = list(armaOrder = c(1,1)), 
			distribution.model = "sstd")
	speclist = as.list(c(spec1, spec2, spec3))
	mspec = multispec( speclist )
	mfit1 = multifit(multispec = mspec, data = Dat, fit.control = list(stationarity=1), 
			cluster = cluster)
	# Filter
	fspec = vector(mode = "list", length = 3)
	fspec[[1]] = spec1
	fspec[[2]] = spec2
	fspec[[3]] = spec3
	for(i in 1:3){
		setfixed(fspec[[i]])<-as.list(coef(mfit1)[[i]])
	}
	mspec1 = multispec( fspec )	
	mfilt1 = multifilter(multifitORspec = mspec1, data = Dat, cluster = cluster)
	# Forecast from Fit
	mforc1 = multiforecast(mfit1, n.ahead = 10, cluster = cluster)
	# Forecast from Spec
	mforc11 = multiforecast(mspec1, data = Dat, n.ahead = 10, cluster = cluster)	
	#------------------------------------------------
	
	#------------------------------------------------
	# Equal Spec
	# Fit
	spec1 = arfimaspec(mean.model = list(armaOrder = c(1,1)))
	mspec = multispec( replicate(3, spec1) )
	mfit2 = multifit(multispec = mspec, data = Dat, cluster = cluster)
	# Filter
	fspec = vector(mode = "list", length = 3)
	fspec = replicate(3, spec1)
	for(i in 1:3){
		setfixed(fspec[[i]])<-as.list(coef(mfit2)[,i])
	}
	mspec2 = multispec( fspec )
	mfilt2 = multifilter(multifitORspec = mspec2, data = Dat, cluster = cluster)
	# Forecast From Fit
	mforc2 = multiforecast(mfit2, n.ahead = 10)
	# Forecast From Spec
	mforc21 = multiforecast(mspec2, data = Dat, n.ahead = 10, cluster = cluster)
	#------------------------------------------------

	#------------------------------------------------
	# Equal Spec/Same Data
	# Fit
	spec1 = arfimaspec(mean.model = list(armaOrder = c(1,1)))
	spec2 = arfimaspec(mean.model = list(armaOrder = c(2,1)))
	spec3 = arfimaspec(mean.model = list(armaOrder = c(3,1)))
	speclist = as.list(c(spec1, spec2, spec3))
	mspec = multispec( speclist )
	mfit3 = multifit(multispec = mspec, data = cbind(Dat[,1], Dat[,1], Dat[,1]),
			cluster = cluster)
	# Forecast
	mforc3 = multiforecast(mfit3, n.ahead = 10, cluster = cluster)
	#------------------------------------------------
	
	zz <- file("test1e.txt", open="wt")
	sink(zz)
	cat("\nMultifit Evaluation\n")
	cat("\nUnequal Spec\n")
	print(mfit1)
	print(likelihood(mfit1))
	print(coef(mfit1))
	print(head(fitted(mfit1)))
	print(head(residuals(mfit1)))
	print(mfilt1)
	print(likelihood(mfilt1))
	print(coef(mfilt1))
	print(head(fitted(mfilt1)))
	print(head(residuals(mfilt1)))
	print(mforc1)
	print(fitted(mforc1))
	print(mforc11)
	print(fitted(mforc11))
	cat("\nEqual Spec\n")
	print(mfit2)
	print(likelihood(mfit2))
	print(coef(mfit2))
	print(head(fitted(mfit2)))
	print(head(residuals(mfit2)))
	print(mfilt2)
	print(likelihood(mfilt2))
	print(coef(mfilt2))
	print(head(fitted(mfilt2)))
	print(head(residuals(mfilt2)))
	print(mforc2)
	print(fitted(mforc2))
	print(mforc21)
	print(fitted(mforc21))
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

rugarch.test1f = function(cluster=NULL){
	# rolling fit/forecast
	tic = Sys.time()
	
	data(sp500ret)
	spec = arfimaspec()
	roll1 = arfimaroll(spec,  data = sp500ret, n.ahead = 1, forecast.length = 500, 
			refit.every = 25, refit.window = "moving", cluster = cluster,
			solver = "hybrid", fit.control = list(), solver.control = list() ,
			calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05))
	
	# as.ARFIMAforecast
	# as.data.frame
	
	zz <- file("test1f.txt", open="wt")
	sink(zz)
	cat("\nForecast Evaluation\n")
	report(roll1, "VaR")
	report(roll1, "fpm")
	# Extractor Functions:
	# default:
	print(head(as.data.frame(roll1, which = "density"), 25))
	print(tail(as.data.frame(roll1, which = "density"), 25))
	print(head(as.data.frame(roll1, which = "VaR"), 25))
	print(tail(as.data.frame(roll1, which = "VaR"), 25))
	print(coef(roll1)[[1]])
	print(coef(roll1)[[20]])
	print(head(fpm(roll1, summary=FALSE)))
	sink(type="message")
	sink()
	close(zz)
	
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

rugarch.test1g = function(cluster=NULL){
	# simulation
	tic = Sys.time()
	require(fracdiff)
	
	spec1 = arfimaspec( mean.model = list(armaOrder = c(2,1), include.mean = TRUE, arfima = TRUE), 
			distribution.model = "std", fixed.pars = list(mu = 0.02, ar1 = 0.6, ar2 = 0.01, ma1 = -0.7, arfima = 0,
					shape = 5, sigma = 0.0123))
	spec2 = arfimaspec( mean.model = list(armaOrder = c(2,1), include.mean = TRUE, arfima = FALSE), 
			distribution.model = "std", fixed.pars = list(mu = 0.02, ar1 = 0.6, ar2 = 0.01, ma1 = -0.7,
					shape = 5, sigma = 0.0123))
	
	sim1 = arfimapath(spec1, n.sim = 1000, m.sim = 1, rseed = 100, preresiduals = c(0,0), prereturns = c(0,0),
			n.start=1)
	sim2 = arfimapath(spec2, n.sim = 1000, m.sim = 1, rseed = 100)
	
	sim1 = arfimapath(spec1, n.sim = 1000, m.sim = 1, rseed = 100, n.start=1)
	sim2 = arfimapath(spec2, n.sim = 1000, m.sim = 1, rseed = 100, n.start=1)
	
	zz <- file("test1g-1.txt", open="wt")
	sink(zz)
	cat("\nARFIMA and ARMA simulation tests:\n")
	print(tail(fitted(sim1)), digits = 5)
	print(tail(fitted(sim2)), digits = 5)
	sink(type="message")
	sink()
	close(zz)
	# Now the rugarch simulation of ARFIMA/ARMA with arima.sim of R
	# Note that arima.sim simulates the residuals (i.e no mean):
	#  ARMA(2,2)
	set.seed(33)
	inn = rdist("std", 1000, mu = 0, sigma = 1, lambda = 0, skew = 0, shape = 5)
	
	spec1 = arfimaspec( mean.model = list(armaOrder = c(2,2), include.mean = FALSE, arfima = TRUE), 
			distribution.model = "std", fixed.pars = list(ar1 = 0.6, ar2 = 0.21, arfima = 0, ma1 = -0.7,
					ma2 = 0.3, shape = 5, sigma = 0.0123))
	spec2 = arfimaspec( mean.model = list(armaOrder = c(2,2), include.mean = FALSE, arfima = FALSE), 
			distribution.model = "std", fixed.pars = list(ar1 = 0.6, ar2 = 0.21, ma1 = -0.7,
					ma2 = 0.3, shape = 5,sigma = 0.0123))
	# Notice the warning...it would be an error had we not added 2 extra zeros to the custom distribution
	# equal to the MA order since n.start >= MA order in arfima model
	sim1 = arfimapath(spec1, n.sim = 1000, m.sim = 1, rseed = 100, preresiduals = c(0,0), prereturns = c(0,0), 
			custom.dist = list(name = "sample",
					distfit = matrix(c(0,0,inn), ncol = 1), type = "z"))
	sim2 = arfimapath(spec2, n.sim = 1000, m.sim = 1, rseed = 100, preresiduals = c(0,0), prereturns = c(0,0), 
			custom.dist = list(name = "sample",
					distfit = matrix(inn, ncol = 1), type = "z"))
	
	
	# Test with a GARCH specification as well (with alpha=beta=0)
	specx = ugarchspec( mean.model = list(armaOrder = c(2,2), include.mean = FALSE, arfima = TRUE), 
			distribution.model = "std", fixed.pars = list(ar1 = 0.6, ar2 = 0.21, ma1 = -0.7,
					ma2 = 0.3, arfima=0, shape = 5, omega = 0.0123^2, alpha1 = 0, beta1=0))
	
	simx  = ugarchpath(specx, n.sim = 1000, m.sim = 1, rseed = 100, preresiduals = c(0,0), prereturns = c(0,0), 
			presigma = c(0,0), custom.dist = list(name = "sample",
					distfit = matrix(c(0,0,inn), ncol = 1), type = "z"))
	
	
	# Note that we pass the non-standardized innovations to arima.sim (i.e. multiply by sigma)
	sim3 = arima.sim(model = list(ar = c(0.6, 0.21), ma = c(-0.7, 0.3)), n = 1000,
			n.start = 4, start.innov = c(0,0,0,0), innov = inn*0.0123)
	
	# set fracdiff setting to n.start=0 and allow.0.nstart=TRUE
	sim4 = fracdiff.sim(n=1000, ar = c(0.6, 0.21), ma = c(0.7, -0.3), d = 0,
			innov = c(0,0,inn*0.0123), n.start = 0, backComp = TRUE, allow.0.nstart = TRUE,
			 mu = 0)
	
	tst1 = cbind(head(fitted(sim1)), head(fitted(sim2)), head(sim3), head(sim4$series), head(fitted(simx)))
	tst2 = cbind(tail(fitted(sim1)), tail(fitted(sim2)), tail(sim3), tail(sim4$series), tail(fitted(simx)))
	
	colnames(tst1) = colnames(tst2) = c("ARFIMA(d = 0)", "ARMA", "arima.sim", "fracdiff", "GARCH(0,0)")
	
	zz <- file("test1g-2.txt", open="wt")
	sink(zz)
	cat("\nARFIMA, ARMA arima.sim simulation tests:\n")
	print(tst1, digits = 6)
	print(tst2, digits = 6)
	
	sink(type="message")
	sink()
	close(zz)
	#  ARMA(2,1)
	set.seed(33)
	inn = rdist("std", 1000, mu = 0, sigma = 1, lambda = 0, skew = 0, shape = 5)
	
	spec1 = arfimaspec( mean.model = list(armaOrder = c(2,1), include.mean = FALSE, arfima = TRUE), 
			distribution.model = "std", fixed.pars = list(ar1 = 0.6, ar2 = 0.21, arfima = 0, ma1 = -0.7,
					shape = 5, sigma = 0.0123))
	spec2 = arfimaspec( mean.model = list(armaOrder = c(2,1), include.mean = FALSE, arfima = FALSE), 
			distribution.model = "std", fixed.pars = list(ar1 = 0.6, ar2 = 0.21, ma1 = -0.7,
					shape = 5,sigma = 0.0123))

	
	sim1 = arfimapath(spec1, n.sim = 1000, m.sim = 1, rseed = 100, preresiduals = c(0,0), prereturns = c(0,0), 
			custom.dist = list(name = "sample",
					distfit = matrix(c(0,inn), ncol = 1), type = "z"))
	sim2 = arfimapath(spec2, n.sim = 1000, m.sim = 1, rseed = 100, preresiduals = c(0,0), prereturns = c(0,0), 
			custom.dist = list(name = "sample",
					distfit = matrix(inn, ncol = 1), type = "z"))
	
	
	# Test with a GARCH specification as well (with alpha=beta=0)
	specx = ugarchspec( mean.model = list(armaOrder = c(2,1), include.mean = FALSE, arfima = TRUE), 
			distribution.model = "std", fixed.pars = list(ar1 = 0.6, ar2 = 0.21, ma1 = -0.7,
					arfima=0, shape = 5, omega = 0.0123^2, alpha1 = 0, beta1=0))
	
	simx  = ugarchpath(specx, n.sim = 1000, m.sim = 1, rseed = 100, preresiduals = c(0,0), prereturns = c(0,0), 
			presigma = c(0,0), custom.dist = list(name = "sample",
					distfit = matrix(c(0,inn), ncol = 1), type = "z"))
	# Note that we pass the non-standardized innovations to arima.sim (i.e. multiply by sigma)
	sim3 = arima.sim(model = list(ar = c(0.6, 0.21), ma = c(-0.7)), n = 1000,
			n.start = 3, start.innov = c(0,0,0), innov = inn*0.0123)
	
	tst1 = cbind(head(fitted(sim1)), head(fitted(sim2)), head(sim3), head(fitted(simx)))
	tst2 = cbind(tail(fitted(sim1)), tail(fitted(sim2)), tail(sim3), tail(fitted(simx)))
	
	colnames(tst1) = colnames(tst2) = c("ARFIMA(d = 0)", "ARMA", "arima.sim", "GARCH(0,0)")
	
	zz <- file("test1g-3.txt", open="wt")
	sink(zz)
	cat("\nARFIMA, ARMA arima.sim simulation tests:\n")
	print(tst1, digits = 6)
	print(tst2, digits = 6)
	
	sink(type="message")
	sink()
	close(zz)
	
	#  Pure AR
	set.seed(33)
	inn = rdist("std", 1000, mu = 0, sigma = 1, lambda = 0, skew = 0, shape = 5)
	
	spec1 = arfimaspec( mean.model = list(armaOrder = c(2,0), include.mean = FALSE, arfima = TRUE), 
			distribution.model = "std", fixed.pars = list(ar1 = 0.6, ar2 = 0.21, arfima = 0, ma1 = -0.7,
					shape = 5, sigma = 0.0123))
	spec2 = arfimaspec( mean.model = list(armaOrder = c(2,0), include.mean = FALSE, arfima = FALSE), 
			distribution.model = "std", fixed.pars = list(ar1 = 0.6, ar2 = 0.21, ma1 = -0.7,
					shape = 5,sigma = 0.0123))
	
	sim1 = arfimapath(spec1, n.sim = 1000, m.sim = 1, rseed = 100, preresiduals = c(0,0), prereturns = c(0,0), 
			custom.dist = list(name = "sample",
					distfit = matrix(inn, ncol = 1), type = "z"))
	sim2 = arfimapath(spec2, n.sim = 1000, m.sim = 1, rseed = 100, preresiduals = c(0,0), prereturns = c(0,0), 
			custom.dist = list(name = "sample",
					distfit = matrix(inn, ncol = 1), type = "z"))
	
	specx = ugarchspec( mean.model = list(armaOrder = c(2,0), include.mean = FALSE, arfima = TRUE), 
			distribution.model = "std", fixed.pars = list(ar1 = 0.6, ar2 = 0.21, 
					arfima=0, shape = 5, omega = 0.0123^2, alpha1 = 0, beta1=0))
	
	simx  = ugarchpath(specx, n.sim = 1000, m.sim = 1, rseed = 100, preresiduals = c(0,0), prereturns = c(0,0), 
			presigma = c(0,0), custom.dist = list(name = "sample",
					distfit = matrix(c(inn), ncol = 1), type = "z"))
	
	# Note that we pass the non-standardized innovations to arima.sim (i.e. multiply by sigma)
	sim3 = arima.sim(model = list(ar = c(0.6, 0.21), ma = NULL), n = 1000,
			n.start = 2, start.innov = c(0,0), innov = inn*0.0123)
	
	tst1 = cbind(head(fitted(sim1)), head(fitted(sim2)), head(sim3), head(fitted(simx)))
	tst2 = cbind(tail(fitted(sim1)), tail(fitted(sim2)), tail(sim3), tail(fitted(simx)))
	
	colnames(tst1) = colnames(tst2) = c("ARFIMA(d = 0)", "ARMA", "arima.sim", "GARCH(0,0)")
	
	zz <- file("test1g-4.txt", open="wt")
	sink(zz)
	cat("\nARFIMA, ARMA arima.sim simulation tests:\n")
	print(tst1, digits = 6)
	print(tst2, digits = 6)
	
	sink(type="message")
	sink()
	close(zz)
	
	#  Pure MA
	set.seed(33)
	inn = rdist("std", 1000, mu = 0, sigma = 1, lambda = 0, skew = 0, shape = 5)
	
	spec1 = arfimaspec( mean.model = list(armaOrder = c(0,2), include.mean = FALSE, arfima = TRUE), 
			distribution.model = "std", fixed.pars = list(ma1 = 0.6, ma2 = -0.21, arfima = 0,
					shape = 5, sigma = 0.0123))
	spec2 = arfimaspec( mean.model = list(armaOrder = c(0,2), include.mean = FALSE, arfima = FALSE), 
			distribution.model = "std", fixed.pars = list(ma1 = 0.6, ma2 = -0.21,
					shape = 5,sigma = 0.0123))
	
	sim1 = arfimapath(spec = spec1, n.sim = 1000, m.sim = 1, rseed = 100, preresiduals = c(0,0), prereturns = c(0,0), 
			custom.dist = list(name = "sample",
					distfit = matrix(c(0,0,inn), ncol = 1), type = "z"))
	sim2 = arfimapath(spec2, n.sim = 1000, m.sim = 1, rseed = 100, preresiduals = c(0,0), prereturns = c(0,0), 
			custom.dist = list(name = "sample",
					distfit = matrix(inn, ncol = 1), type = "z"))
	
	specx = ugarchspec( mean.model = list(armaOrder = c(0,2), include.mean = FALSE, arfima = TRUE), 
			distribution.model = "std", fixed.pars = list(ma1 = 0.6, ma2 = -0.21, 
					arfima=0, shape = 5, omega = 0.0123^2, alpha1 = 0, beta1=0))
	
	simx  = ugarchpath(specx, n.sim = 1000, m.sim = 1, rseed = 100, preresiduals = c(0,0), prereturns = c(0,0), 
			presigma = c(0,0), custom.dist = list(name = "sample",
					distfit = matrix(c(0,0,inn), ncol = 1), type = "z"))
	
	# Note that we pass the non-standardized innovations to arima.sim (i.e. multiply by sigma)
	set.seed(33)
	inn = rdist("std", 1000, mu = 0, sigma = 1, lambda = 0, skew = 0, shape = 5)

	sim3 = arima.sim(model = list(ar = NULL, ma = c(0.6, -0.21)), n = 1000,
			n.start = 2, start.innov = c(0,0), innov = inn*0.0123)
	
	tst1 = cbind(head(fitted(sim1)), head(fitted(sim2)), head(sim3), head(fitted(simx)))
	tst2 = cbind(tail(fitted(sim1)), tail(fitted(sim2)), tail(sim3), tail(fitted(simx)))
	colnames(tst1) = colnames(tst2) = c("ARFIMA(d = 0)", "ARMA", "arima.sim", "GARCH(0,0)")
	
	zz <- file("test1g-5.txt", open="wt")
	sink(zz)
	cat("\nARFIMA, ARMA arima.sim simulation tests:\n")
	print(tst1, digits = 6)
	print(tst2, digits = 6)
	sink(type="message")
	sink()
	close(zz)
	# arfimasim + exogenous regressors + custom innovations
	data(dji30ret)
	Dat = dji30ret[,1, drop = FALSE]
	T = dim(Dat)[1]
	Bench = as.matrix(cbind(apply(dji30ret[,2:10], 1, "mean"), apply(dji30ret[,11:20], 1, "mean")))	

	spec = arfimaspec( mean.model = list(armaOrder = c(1,1), include.mean = TRUE, arfima = FALSE, 
					external.regressors = Bench), distribution.model = "std")
	fit = arfimafit(spec = spec, data = Dat, solver = "solnp", out.sample = 500)
	
	# lag1 Benchmark
	BenchF = Bench[(T-500):(T-500+9), , drop = FALSE]
	
	exsim = vector(mode = "list", length = 10000)
	for(i in 1:10000) exsim[[i]] = as.matrix(BenchF)
	# simulated residuals
	res = residuals(fit)
	ressim = matrix(NA, ncol = 10000, nrow = 10)
	set.seed(10000)
	for(i in 1:10000) ressim[,i] = sample(res, 10, replace = TRUE)
	
	sim = arfimasim(fit, n.sim = 10, m.sim = 10000, startMethod="sample", 
			custom.dist = list(name = "sample", distfit = ressim, type = "res"), mexsimdata = exsim)
	forc = fitted(arfimaforecast(fit, n.ahead = 10, external.forecasts = list(mregfor = BenchF)))
	simx = fitted(sim)
	actual10 = Dat[(T-500+1):(T-500+10), 1, drop = FALSE]
	
	simm = apply(simx, 1 ,"mean")
	simsd = apply(simx, 1 ,"sd")
	
	zz <- file("test1g-6.txt", open="wt")
	sink(zz)
	print(round(cbind(actual10, forc, simm, simsd),5), digits = 4)
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

# ARFIMA benchmark tests
rugarch.test1h = function(cluster=NULL){
	tic = Sys.time()	
	# ARFIMA(2,d,1)
	require(fracdiff)
	truecoef1 = list(mu = 0.005, ar1 = 0.6, ar2 = 0.01, ma1 = -0.7, arfima = 0.3, sigma = 0.0123)
	spec1 = arfimaspec( 
	mean.model = list(armaOrder = c(2,1), include.mean = TRUE, arfima = TRUE), 
	distribution.model = "norm", 
	fixed.pars = truecoef1)
	sim1 = arfimapath(spec1, n.sim = 5000, n.start = 100, m.sim = 1, rseed = 101)
	data1 = fitted(sim1)
	#write.csv(data1[,1], file = "D:/temp1.csv")
	spec1 = arfimaspec( 
	mean.model = list(armaOrder = c(2,1), include.mean = TRUE, arfima = TRUE), 
	distribution.model = "norm")
	fit1 = arfimafit(spec1, data = data1)
	fit1.fd = fracdiff(as.numeric(data1[,1])-coef(fit1)["mu"], nar = 2, nma = 1)
	
	# Commercial Implementation Program Fit (NLS-with imposed stationarity):
	commcheck1 = c(0.00488381, 0.537045,  0.0319251,  -0.721266, 0.348604, 0.0122415)
	fdcheck1 = c(NA, coef(fit1.fd)[2:3], -coef(fit1.fd)[4], coef(fit1.fd)[1], fit1.fd$sigma)
	
	chk1 = cbind(coef(fit1), commcheck1, fdcheck1, unlist(truecoef1))
	colnames(chk1) = c("rugarch", "commercial", "fracdiff", "true")
	chk1lik = c(likelihood(fit1),   14920.4279, fit1.fd$log.likelihood)
	
	# ARFIMA(2,d,0)
	truecoef2 = list(mu = 0.005, ar1 = 0.6, ar2 = 0.01, arfima = 0.1, sigma = 0.0123)
	spec2 = arfimaspec( 
	mean.model = list(armaOrder = c(2,0), include.mean = TRUE, arfima = TRUE), 
	distribution.model = "norm", 
	fixed.pars = truecoef2)
	sim2 = arfimapath(spec2, n.sim = 5000, n.start = 100, m.sim = 1,  rseed = 102)	
	data2 = fitted(sim2)
	#write.csv(data2[,1], file = "D:/temp2.csv")
	spec2 = arfimaspec( 
	mean.model = list(armaOrder = c(2,0), include.mean = TRUE, arfima = TRUE), 
	distribution.model = "norm")
	fit2 = arfimafit(spec2, data = data2)
	fit2.fd = fracdiff(as.numeric(data2[,1])-coef(fit2)["mu"], nar = 2, nma = 0)
	fdcheck2 = c(NA, coef(fit2.fd)[2:3], coef(fit2.fd)[1], fit2.fd$sigma)
	
	commcheck2 = c( 0.00585040, 0.692693, 0.000108778,0.00466664,0.0122636)
	
	chk2 = cbind(coef(fit2), commcheck2, fdcheck2, unlist(truecoef2))
	colnames(chk2) = c("rugarch", "commercial", "fracdiff", "true")
	chk2lik = c(likelihood(fit2), 14954.5702, fit2.fd$log.likelihood)
	

	# ARFIMA(0,d,2)
	truecoef3 = list(mu = 0.005, ma1 = 0.3, ma2 = 0.2, arfima = 0.1, sigma = 0.0123)
	spec3 = arfimaspec( 
	mean.model = list(armaOrder = c(0,2), include.mean = TRUE, arfima = TRUE), 
	distribution.model = "norm", 
	fixed.pars = truecoef3)
	sim3 = arfimapath(spec3, n.sim = 5000, n.start = 100, m.sim = 1, rseed = 103)
	data3 = fitted(sim3)
	#write.csv(data3[,1], file = "D:/temp3.csv")
	spec3 = arfimaspec( 
	mean.model = list(armaOrder = c(0,2), include.mean = TRUE, arfima = TRUE), 
	distribution.model = "norm")
	fit3 = arfimafit(spec3, data = data3, solver="hybrid")
	
	fit3.fd = fracdiff(as.numeric(data3[,1])-coef(fit3)["mu"], nar = 0, nma = 2)
	fdcheck3 = c(NA, -coef(fit3.fd)[2:3], coef(fit3.fd)[1], fit3.fd$sigma)
	
	commcheck3 = c( 0.00580941, 0.320205, 0.206786, 0.0546052, 0.0120114)
	
	chk3 = cbind(coef(fit3), commcheck3, fdcheck3, unlist(truecoef3))
	colnames(chk3) = c("rugarch", "commercial", "fracdiff", "true")
	chk3lik = c(likelihood(fit3), 15015.2957, fit3.fd$log.likelihood)

	
	# ARFIMA(2,d,1) simulation (using rugarch path)
	truecoef = list(mu = 0.005, ar1 = 0.6, ar2 = 0.01, ma1 = -0.7, arfima = 0.45, sigma = 0.0123)
	spec = arfimaspec( 
	mean.model = list(armaOrder = c(2,1), include.mean = TRUE, arfima = TRUE), 
	distribution.model = "norm", fixed.pars = truecoef)
	sim = arfimapath(spec, n.sim = 5000, n.start = 100, m.sim = 50, rseed = 1:50)
	Data = fitted(sim)
	spec = arfimaspec( 
	mean.model = list(armaOrder = c(2,1), include.mean = TRUE, arfima = TRUE), 
	distribution.model = "norm")
	coefx = matrix(NA, ncol = 6, nrow = 50)
	coefy = matrix(NA, ncol = 6, nrow = 50)
	if(!is.null(cluster)){
		parallel::clusterEvalQ(cluster, require(rugarch))
		parallel::clusterEvalQ(cluster, require(fracdiff))
		parallel::clusterExport(cluster, c("Data", "spec"), envir = environment())
		sol = parallel::parLapply(cluster, as.list(1:50), fun = function(i){
					fit = arfimafit(spec, data = Data[,i], solver="hybrid")
					if(fit@fit$convergence == 0) coefx = coef(fit) else coefx = rep(NA, 6)
					if(fit@fit$convergence == 0){ 
						fit = fracdiff(as.numeric(Data[,i]) - coef(fit)["mu"], nar = 2, nma = 1)
					} else{
						fit = fracdiff(scale(as.numeric(Data[,i]), scale=F), nar = 2, nma = 1)
					}
					coefy = c(NA, coef(fit)[2:3], -coef(fit)[4], coef(fit)[1], fit$sigma)
					return(list(coefx = coefx, coefy = coefy))
				})
		coefx = t(sapply(sol, FUN = function(x) x$coefx))
		coefy = t(sapply(sol, FUN = function(x) x$coefy))
	} else{
		for(i in 1:50){
			fit = arfimafit(spec, data = Data[,i], solver="hybrid")
			if(fit@fit$convergence == 0) coefx[i,] = coef(fit)
			fit = fracdiff(scale(as.numeric(Data[,i]), scale=F), nar = 2, nma = 1)
			coefy[i,] = c(NA, coef(fit)[2:3], -coef(fit)[4], coef(fit)[1], fit$sigma)
		}
	}
	
	zz <- file("test1h-1.txt", open="wt")
	sink(zz)
	cat("\nARFIMA(2,d,1)\n")
	print(chk1)
	print(chk1lik)
	cat("\nARFIMA(2,d,0)\n")
	print(chk2)
	print(chk2lik)
	cat("\nARFIMA(0,d,2)\n")
	print(chk3)
	print(chk3lik)
	cat("\nARFIMA(2,d,1) mini-simulation/fit\n")
	# small sample/simulation also use median:
	cat("\nMedian (rugarch, fracdiff)\n")
	print(	data.frame(rugarch=round(apply(coefx, 2, "median"),5), fracdiff = round(apply(coefy, 2, "median"),5), true=unlist(truecoef) ) )
	cat("\nMean (rugarch, fracdiff)\n")
	print(	data.frame(rugarch=round(apply(coefx, 2, "mean"),5), fracdiff = round(apply(coefy, 2, "mean"),5), true=unlist(truecoef) ) )
	print(	data.frame(rugarch.sd =round(apply(coefx, 2, "sd"),5), fracdiff.sd = round(apply(coefy, 2, "sd"),5) ) )
	sink(type="message")
	sink()
	close(zz)
	
	
	# ARFIMA(2,d,1) simulation (using fracdiff path)
	truecoef = list(mu = 0.005, ar1 = 0.6, ar2 = 0.01, ma1 = -0.7, arfima = 0.45, sigma = 0.0123)
	Data = matrix(NA, ncol = 50, nrow = 5000)
	for(i in 1:50){
		set.seed(i)
		sim = fracdiff.sim(n=5000, ar = c(0.6, 0.01), ma = c(0.7), d = 0.45,
				rand.gen = rnorm, n.start = 100, backComp = TRUE, sd = 0.0123, mu = 0.005)
		Data[,i] = sim$series
	}
	spec = arfimaspec( 
			mean.model = list(armaOrder = c(2,1), include.mean = TRUE, arfima = TRUE), 
			distribution.model = "norm")
	coefx = matrix(NA, ncol = 6, nrow = 50)
	coefy = matrix(NA, ncol = 6, nrow = 50)
	if(!is.null(cluster)){
		parallel::clusterEvalQ(cluster, require(rugarch))
		parallel::clusterEvalQ(cluster, require(fracdiff))
		parallel::clusterExport(cluster, c("Data", "spec"), envir = environment())
		sol = parallel::parLapply(cluster, as.list(1:50), fun = function(i){
					fit = arfimafit(spec, data = Data[,i], solver="hybrid")
					if(fit@fit$convergence == 0) coefx = coef(fit) else coefx = rep(NA, 6)
					if(fit@fit$convergence == 0){ 
						fit = fracdiff(as.numeric(Data[,i]) - coef(fit)["mu"], nar = 2, nma = 1)
					} else{
						fit = fracdiff(scale(as.numeric(Data[,i]), scale=F), nar = 2, nma = 1)
					}
					coefy = c(NA, coef(fit)[2:3], -coef(fit)[4], coef(fit)[1], fit$sigma)
					return(list(coefx = coefx, coefy = coefy))
				})
		coefx = t(sapply(sol, FUN = function(x) x$coefx))
		coefy = t(sapply(sol, FUN = function(x) x$coefy))
	} else{
		for(i in 1:50){
			fit = arfimafit(spec, data = Data[,i], solver="hybrid")
			if(fit@fit$convergence == 0) coefx[i,] = coef(fit)
			fit = fracdiff(scale(as.numeric(Data[,i]), scale=F), nar = 2, nma = 1)
			coefy[i,] = c(NA, coef(fit)[2:3], -coef(fit)[4], coef(fit)[1], fit$sigma)
		}
	}
	
	zz <- file("test1h-2.txt", open="wt")
	sink(zz)
	cat("\nARFIMA(2,d,1) mini-simulation/fit2 (simulation from fracdiff.sim)\n")
	# small sample/simulation also use median:
	cat("\nMedian (rugarch, fracdiff)\n")
	print(	data.frame(rugarch=round(apply(coefx, 2, "median"),5), fracdiff = round(apply(coefy, 2, "median"),5), true=unlist(truecoef) ) )
	cat("\nMean (rugarch, fracdiff)\n")
	print(	data.frame(rugarch=round(apply(coefx, 2, "mean"),5), fracdiff = round(apply(coefy, 2, "mean"),5), true=unlist(truecoef) ) )
	print(	data.frame(rugarch.sd =round(apply(coefx, 2, "sd"),5), fracdiff.sd = round(apply(coefy, 2, "sd"),5) ) )
	sink(type="message")
	sink()
	close(zz)
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}