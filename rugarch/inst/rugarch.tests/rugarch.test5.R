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


#################################################################################
# Model Simulation Tests
#################################################################################
rugarch.seeds = function(n, test)
{
	rseed = round(runif(n,100,10000),0)
	options(width=120)
	zz <- file(paste(test,"-seeds.txt", sep=""), open="wt")
	sink(zz)
	print(rseed)
	sink(type="message")
	sink()
	close(zz)
	return(rseed)
}

rugarch.test5a = function(cluster=NULL)
{
	#cat("\nrugarch-->test5-1: Simulation Test (sGARCH)\n")
	tic = Sys.time()
	
	data(dji30ret)
	
	# sGARCH (with variations in the simulation inputs/methods)
	# ---------------------------------------------------------------------------------
	# basic simulation (no external regressors in simulation)
	# use external regressors
	dates = rownames(dji30ret[,"AA", drop = FALSE])
	monday = rugarch:::.WeekDayDummy(dates, date.format = "%Y-%m-%d", weekday = "Monday")
	# convert to matrix which is what the specification expects
	monday = matrix(monday, ncol = 1)
	# create a dummy day-of-week variable for the variance regression (Friday)
	friday = rugarch:::.WeekDayDummy(dates, date.format="%Y-%m-%d", weekday = "Friday")
	# convert to matrix which is what the specification expects
	friday = matrix(friday, ncol = 1)
	# define the specification
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1), 
					external.regressors = friday), 
			mean.model = list(armaOrder = c(0,0), include.mean = TRUE, 
					external.regressors = monday), 
			distribution.model = "std")
	
	sgarch.fit = ugarchfit(data = dji30ret[,"AA",drop=FALSE], spec = spec, 
			solver = "solnp")
	
	sgarch.sim1 = ugarchsim(fit = sgarch.fit, n.sim = 1200, n.start = 0, 
			m.sim = 1, startMethod = "unconditional", rseed = 100)
		
	sgarch.sim2 = ugarchsim(fit = sgarch.fit, n.sim = 1200, n.start = 0, 
			m.sim = 1, startMethod = "sample", rseed = 100)
	
	# use pre-sample data (we use last points of data set so should agree with sim2)
	sgarch.sim3 = ugarchsim(fit = sgarch.fit, n.sim = 1200, n.start = 0, 
			m.sim = 1, startMethod = "sample", 
			presigma = sgarch.fit@fit$sigma[5521], 
			prereturns = sgarch.fit@model$modeldata$data[5521], 
			preresiduals = NA, rseed = 100)
	
	# use external regressors
	fwd1200 = rugarch:::.ForwardDates(dates, n.ahead = 1200, date.format="%Y-%m-%d", 
			periodicity="days")
	
	# create a dummy vector for those forward days which are Mondays and Fridays
	fwdMonday = rugarch:::.WeekDayDummy(as.character(fwd1200), date.format="%Y-%m-%d", 
			weekday="Monday")
	fwdFriday = rugarch:::.WeekDayDummy(as.character(fwd1200), date.format="%Y-%m-%d", 
			weekday="Friday")
	
	sgarch.sim4 = ugarchsim(sgarch.fit, n.start = 0, n.sim = 1200, m.sim = 1, 
			startMethod = "sample", mexsimdata = list(matrix(fwdMonday,ncol=1)), 
			vexsimdata = list(matrix(fwdFriday,ncol=1)), rseed = 100)
	
	# simulate with 10 different paths
	# (note the warning about providing only 1 simulated path for the external regressors,
	# in which case it will be reused in all simulations)
	sgarch.sim5 = ugarchsim(sgarch.fit, n.start = 0, n.sim = 1200, m.sim = 10, 
			startMethod = "sample", mexsimdata = list(matrix(fwdMonday,ncol=1)), 
			vexsimdata = list(matrix(fwdFriday,ncol=1)), rseed = 100:109)
	
	options(width=120)
	zz <- file("test5a.txt", open="wt")
	sink(zz)
	cat("\nsGARCH simulation 1:\n")
	show(sgarch.sim1)
	cat("\nsGARCH simulation 2:\n")
	show(sgarch.sim2)
	cat("\nsGARCH simulation 3:\n")
	show(sgarch.sim3)
	cat("\nsGARCH simulation 4:\n")
	show(sgarch.sim4)
	cat("\nsGARCH simulation 5:\n")
	show(sgarch.sim5)
	sink(type="message")
	sink()
	close(zz)
	
	postscript("test5a.eps", width = 12, height = 8)
	par(mfrow = c(3,3))
	for(i in 1:9) plot(sgarch.sim5, which = 3, m.sim = i)
	dev.off()
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

rugarch.test5b = function(cluster=NULL)
{
	#cat("\nrugarch-->test5-2: Simulation Test (eGARCH)\n")
	tic = Sys.time()
	
	data(dji30ret)
	# eGARCH (with variations in the simulation inputs/methods)
	# ---------------------------------------------------------------------------------
	
	# basic simulation (no external regressors in simulation)
	# use external regressors
	dates = rownames(dji30ret[,"AA", drop = FALSE])
	monday = rugarch:::.WeekDayDummy(dates, date.format = "%Y-%m-%d", weekday = "Monday")
	# convert to matrix which is what the specification expects
	monday = matrix(monday, ncol = 1)
	# create a dummy day-of-week variable for the variance regression (Friday)
	friday = rugarch:::.WeekDayDummy(dates, date.format = "%Y-%m-%d", weekday = "Friday")
	# convert to matrix which is what the specification expects
	friday = matrix(friday, ncol = 1)
	# define the specification
	spec = ugarchspec(
			variance.model = list(model = "eGARCH", garchOrder=c(1,1), 
					external.regressors=friday), 
			mean.model = list(armaOrder = c(0,0), 
					include.mean = TRUE, external.regressors = monday), 
			distribution.model = "std")
	
	egarch.fit = ugarchfit(data = dji30ret[,"AA",drop=FALSE], spec = spec, 
			solver = "solnp", fit.control = list(stationarity = 0), 
			solver.control = list(tol = 1e-9, delta = 1e-8))
	
	egarch.sim1 = ugarchsim(fit = egarch.fit, n.sim = 1200, n.start = 0, 
			m.sim = 1, startMethod = "unconditional", rseed = 100)
	
	
	# startMethod = sample
	# if it is not unconditional, it defaults to sample (whatever else is entered is
	# ignored)
	
	egarch.sim2 = ugarchsim(fit = egarch.fit, n.sim = 1200, n.start = 0, m.sim = 1, 
			startMethod = "sample", rseed=100)
	
	# use pre-sample data (we use last points of data set so should agree with sim2)
	egarch.sim3 = ugarchsim(fit = egarch.fit, n.sim = 1200, n.start = 0, 
			m.sim = 1, startMethod = "sample", presigma = egarch.fit@fit$sigma[5521], 
			prereturns = egarch.fit@model$modeldata$data[5521], preresiduals = NA, 
			rseed = 100)
	
	# use external regressors
	fwd1200 = rugarch:::.ForwardDates(dates, n.ahead = 1200, date.format = "%Y-%m-%d", 
			periodicity = "days")
	
	# create a dummy vector for those forward days which are Mondays and Fridays
	fwdMonday = rugarch:::.WeekDayDummy(as.character(fwd1200), date.format="%Y-%m-%d", 
			weekday="Monday")
	fwdFriday = rugarch:::.WeekDayDummy(as.character(fwd1200), date.format="%Y-%m-%d", 
			weekday="Friday")
	
	egarch.sim4 = ugarchsim(egarch.fit, n.start = 0, n.sim = 1200, m.sim = 1, 
			startMethod = "sample", 
			mexsimdata = list(matrix(fwdMonday,ncol=1)), 
			vexsimdata = list(matrix(fwdFriday,ncol=1)), rseed = 100)
	
	# simulate with 10 different paths
	# (note the warning about providing only 1 simulated path for the external regressors,
	# in which case it will be reused in all simulations)
	egarch.sim5 = ugarchsim(egarch.fit, n.start = 0, n.sim = 1200, m.sim = 10, 
			startMethod = "sample", 
			mexsimdata = list(matrix(fwdMonday,ncol=1)), 
			vexsimdata = list(matrix(fwdFriday,ncol=1)), rseed = 100:109)
	
	options(width=120)
	zz <- file("test5b.txt", open="wt")
	sink(zz)
	cat("\neGARCH simulation 1:\n")
	show(egarch.sim1)
	cat("\neGARCH simulation 2:\n")
	show(egarch.sim2)
	cat("\neGARCH simulation 3:\n")
	show(egarch.sim3)
	cat("\neGARCH simulation 4:\n")
	show(egarch.sim4)
	cat("\neGARCH simulation 5:\n")
	show(egarch.sim5)
	sink(type="message")
	sink()
	close(zz)
	
	postscript("test5b.eps", width = 12, height = 8)
	par(mfrow = c(3,3))
	for(i in 1:9) plot(egarch.sim5, which = 3, m.sim = i)
	dev.off()
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

rugarch.test5c = function(cluster=NULL)
{
	#cat("\nrugarch-->test5-3: Simulation Test (sGARCH w/th ARFIMA)\n")
	tic = Sys.time()
	
	# simulation of arfima process, and check via refitting:
	rseed = rugarch.seeds(100, "test5c")
	
	data(sp500ret)
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder=c(1,1), arfima = TRUE, 
					include.mean = TRUE, 
					archm = FALSE, archpow = 2), 
			distribution.model="norm")
	
	fit = ugarchfit(data = sp500ret[, 1, drop = FALSE], spec = spec, 
			solver = "solnp")
	tc = coef(fit)
	sim = ugarchsim(fit, n.sim = 5000, n.start = 100, m.sim = 100)
	simdf = fitted(sim)
	setstart(spec) <- as.list(coef(fit))
	if( !is.null(cluster) ){
		parallel::clusterEvalQ(cluster, require(rugarch))
		parallel::clusterExport(cluster, c("spec", "simdf"), envir = environment())
		fitlist = parallel::parLapply(cluster, as.list(1:100), fun = function(i){
					try(ugarchfit(spec = spec, data = simdf[,i]))
				})
	} else{
		fitlist = apply(simdf, 2, FUN = function(x) try(ugarchfit(data = x, spec = spec)))
	}
	
	coefl = t(sapply(fitlist, FUN = function(x) if(is.null(coef(x))) rep(NA, 7) else coef(x)))
	
	options(width=150)
	zz <- file("test5c1.txt", open="wt")
	sink(zz)
	print(as.data.frame(coefl), digits = 6)
	sink(type="message")
	sink()
	close(zz)
	
	exc = which(is.na(coefl[,1]))
	if(length(exc) > 0) exc = exc else exc = -1:-100
	postscript("test5c1.eps", width = 12, height = 8)
	par(mfrow = c(3,3))
	plot(density(coefl[-exc,1]), 
			main = paste("ARFIMA: mu parameter\n true parameter=", 
							round(uncmean(fit),5), sep=""))
	plot(density(coefl[-exc,2]), 
			main = paste("ARFIMA: ar parameter\n true parameter=", 
					round(tc["ar1"],3), sep=""))
	plot(density(coefl[-exc,3]), 
			main = paste("ARFIMA: ma parameter\n true parameter=", 
					round(tc["ma1"],3), sep=""))
	plot(density(coefl[-exc,4]), 
			main = paste("ARFIMA: arfima parameter\n true parameter=", 
					round(tc["arfima"],4), sep=""))
	plot(density(coefl[-exc,5]), 
			main = paste("sGARCH: omega parameter\n true parameter=", 
					round(tc["omega"],5), sep=""))
	plot(density(coefl[-exc,6]), 
			main = paste("sGARCH: alpha parameter\n true parameter=",
					round(tc["alpha1"],3), sep=""))
	plot(density(coefl[-exc,7]), 
			main = paste("sGARCH: beta parameter\n true parameter=", 
					round(tc["beta1"],3), sep=""))
	dev.off()
	
	
	spec2 = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder=c(2,1), arfima = TRUE, 
					include.mean = TRUE, archm = FALSE, archpow = 2), 
			distribution.model="norm")
	
	fit2 = ugarchfit(data = sp500ret[,1,drop=FALSE], spec = spec2, 
			solver = "solnp", fit.control=list(scale=1))
	
	tc = coef(fit2)
	sim2 = ugarchsim(fit2, n.sim = 5000, n.start = 100, m.sim = 100, 
			rseed = rseed)
	simdf2 = fitted(sim2)
	setstart(spec2)<-as.list(coef(fit2))
	
	
	if( !is.null(cluster) ){
		# already evaluated, but in case it is not run as one file.
		parallel::clusterEvalQ(cluster,require(rugarch))
		parallel::clusterExport(cluster, c("spec2", "simdf2"), envir = environment())
		fitlist2 = parallel::parLapply(cluster, as.list(1:100), fun = function(i){
					try(ugarchfit(spec = spec2, data = simdf2[,i], fit.control=list(scale=1)))
				})
	} else{
		fitlist2 = apply(simdf2, 2, FUN = function(x) try(ugarchfit(data = x, spec = spec2)))
	}
	
	
	coefl2 = matrix(NA, ncol = 8, nrow = 100)
	for(i in 1:100) if(fitlist2[[i]]@fit$convergence==0) coefl2[i,] = coef(fitlist2[[i]])
	
	options(width = 150)
	zz <- file("test5c2.txt", open="wt")
	sink(zz)
	print(as.data.frame(coefl2), digits = 6)
	sink(type="message")
	sink()
	close(zz)
	
	exc = which(is.na(coefl2[,1]))
	if(length(exc) > 0) exc = exc else exc = -1:-100
	postscript("test5c2.eps", width = 12, height = 8)
	par(mfrow = c(3,3))
	plot(density(coefl2[-exc,1]), 
			main = paste("ARFIMA: mu parameter\n true parameter=", 
					round(uncmean(fit),5), sep=""))
	plot(density(coefl2[-exc,2]), 
			main = paste("ARFIMA: ar1 parameter\n true parameter=", 
					round(tc["ar1"],3), sep=""))
	plot(density(coefl2[-exc,3]), 
			main = paste("ARFIMA: ar2 parameter\n true parameter=", 
					round(tc["ar2"],3), sep=""))
	plot(density(coefl2[-exc,4]), 
			main = paste("ARFIMA: ma parameter\n true parameter=", 
					round(tc["ma1"],3), sep=""))
	plot(density(coefl2[-exc,5]), 
			main = paste("ARFIMA: arfima parameter\n true parameter=", 
					round(tc["arfima"],4), sep=""))
	plot(density(coefl2[-exc,6]), 
			main = paste("sGARCH: omega parameter\n true parameter=", 
					round(tc["omega"],5), sep=""))
	plot(density(coefl2[-exc,7]), 
			main = paste("sGARCH: alpha parameter\n true parameter=", 
					round(tc["alpha1"],3), sep=""))
	plot(density(coefl2[-exc,8]), 
			main = paste("sGARCH: beta parameter\n true parameter=", 
					round(tc["beta1"],3), sep=""))
	dev.off()
	
	spec3 = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,2), arfima = TRUE, 
					include.mean = TRUE, 
					archm = FALSE, archpow = 2), 
			distribution.model="norm")
	
	fit3 = ugarchfit(data = sp500ret[, 1, drop = FALSE], spec = spec3, 
			solver="solnp", solver.control=list(tol=1e-10, trace=1))
	sim3 = ugarchsim(fit3, n.sim = 5000, n.start = 100, m.sim = 100, 
			rseed = rseed)
	simdf3 = fitted(sim3)
	setstart(spec3)<-as.list(coef(fit3))
	tc = coef(fit3)
	
	# Add parallel
	if( !is.null(cluster) ){
		# already evaluated, but in case it is not run as one file.
		parallel::clusterEvalQ(cluster,require(rugarch))
		parallel::clusterExport(cluster, c("spec3", "simdf3"), envir = environment())
		fitlist3 = parallel::parLapply(cluster, as.list(1:100), fun = function(i){
					ugarchfit(spec = spec3, data = simdf3[,i])
				})
	} else{
		fitlist3 = apply(simdf3, 2, FUN = function(x) ugarchfit(data = x, spec = spec3))
	}

	coefl3 = matrix(NA, ncol = 8, nrow = 100)
	for(i in 1:100) if(fitlist3[[i]]@fit$convergence==0) 
			coefl3[i,] = coef(fitlist3[[i]])
	
	options(width=150)
	zz <- file("test5c3.txt", open="wt")
	sink(zz)
	print(as.data.frame(coefl3), digits = 6)
	sink(type="message")
	sink()
	close(zz)
	
	exc = which(is.na(coefl3[,1]))
	if(length(exc) > 0) exc = exc else exc = -1:-100
	postscript("test5c3.eps", width = 12, height = 8)
	par(mfrow = c(3,3))
	plot(density(coefl3[-exc,1]), 
			main = paste("ARFIMA: mu parameter\n true parameter=", 
					round(uncmean(fit3),5), sep=""))
	plot(density(coefl3[-exc,2]), 
			main = paste("ARFIMA: ar1 parameter\n true parameter=", 
					round(tc["ar1"],3), sep=""))
	plot(density(coefl3[-exc,3]), 
			main = paste("ARFIMA: ma1 parameter\n true parameter=", 
					round(tc["ma1"],3), sep=""))
	plot(density(coefl3[-exc,4]), 
			main = paste("ARFIMA: ma2 parameter\n true parameter=", 
					round(tc["ma2"],3), sep=""))
	plot(density(coefl3[-exc,5]), 
			main = paste("ARFIMA: arfima parameter\n true parameter=", 
					round(tc["arfima"],4), sep=""))
	plot(density(coefl3[-exc,6]), 
			main = paste("sGARCH: omega parameter\n true parameter=", 
					round(tc["omega"],5), sep=""))
	plot(density(coefl3[-exc,7]), 
			main = paste("sGARCH: alpha parameter\n true parameter=", 
					round(tc["alpha1"],3), sep=""))
	plot(density(coefl3[-exc,8]), 
			main = paste("sGARCH: beta parameter\n true parameter=", 
					round(tc["beta1"],3), sep=""))
	dev.off()
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}