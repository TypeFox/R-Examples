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


ugarchbench = function( benchmark = c("commercial", "published") )
{
	switch(benchmark[1],
			commercial = .commbench(),
			published = .publishbench())
	return(invisible(benchmark))
}

.publishbench = function()
{
	#Example: GARCH w/th Bollerslev-Ghysels Benchmark
	# load the Deutschemark-Sterling Benchmark Returns Data
	# Log Relative Error Test measures the number of digits of accuracy of rugarch versus benchmark
	newv<-new.env(hash = TRUE, parent = parent.frame(), size = 29L)
	data("dmbp", envir=newv)
	dmbp<-get("dmbp", envir = newv)
	dmbp = as.matrix(dmbp)
	spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(0,0), include.mean = TRUE), 
			distribution.model = "norm")
	# fit the model
	fit1 = ugarchfit(data = dmbp[,1,drop=FALSE], spec = spec, solver = "solnp", solver.control = list(trace=0))
	# test against the benchmark using log relative error measure:
	benchmark.pars = c(-0.00619041, 0.0107613, 0.153134, 0.805974)
	benchmark.se = c(0.00846212, 0.00285271, 0.0265228, 0.0335527)
	rugarch.pars = coef(fit1)
	rugarch.se = fit1@fit$matcoef[,2]
	LRE.vars = -log(abs(rugarch.pars-benchmark.pars)/abs(benchmark.pars), base = 10)
	LRE.se = -log(abs(rugarch.se-benchmark.se)/abs(benchmark.se), base = 10)
	test = cbind(LRE.vars, LRE.se)
	colnames(test) = c("coef", "st.error")
	rownames(test)<-c("mu", "omega", "alpha", "beta")
	cat("\nBollerslev-Ghysels Benchmark 1/2: mu-ARMA(0,0)-sGARCH(1,1)-norm\n")
	cat("\nparameters:\n")
	tmp = t(cbind(rugarch.pars, benchmark.pars))
	rownames(tmp) = c("rugarch", "benchmark")
	print(tmp)
	cat("\nstandard errors:\n")
	tmp = t(cbind(rugarch.se, benchmark.se))
	rownames(tmp) = c("rugarch", "benchmark")
	print(tmp)
	cat("\nLog Relative Error Test:\n")
	print(test, digits=4)
	
	spec = ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(0,0), include.mean = TRUE), 
			distribution.model = "norm")
	# fit the model
	fit2 = ugarchfit(data = dmbp[,1,drop=FALSE], spec = spec, solver = "solnp")
	benchmark.pars = c(-.01167873487, -.12633933747, -.03845788444, .91265373928, .33305592776)
	benchmark.se = c(.00886, .0285, .0192, .0168, .0406)
	rugarch.pars = coef(fit2)
	rugarch.se = fit2@fit$matcoef[,2]
	LRE.vars = -log(abs(rugarch.pars-benchmark.pars)/abs(benchmark.pars), base = 10)
	LRE.se = -log(abs(rugarch.se-benchmark.se)/abs(benchmark.se), base = 10)
	test = cbind(LRE.vars, LRE.se)
	colnames(test) = c("coef", "st.error")
	rownames(test)<-c("mu", "omega", "alpha", "gamma", "beta")
	cat("\nBollerslev-Ghysels Benchmark 2/2: mu-ARMA(0,0)-eGARCH(1,1)-norm\n")
	cat("\nparameters:\n")
	tmp = t(cbind(rugarch.pars, benchmark.pars))
	rownames(tmp) = c("rugarch", "benchmark")
	print(tmp)
	cat("\nstandard errors:\n")
	tmp = t(cbind(rugarch.se, benchmark.se))
	rownames(tmp) = c("rugarch", "benchmark")
	print(tmp)
	cat("\nLog Relative Error Test:\n")
	print(test, digits=4)
	return(invisible(0))
}
.commbench = function()
{
	newv<-new.env(hash = TRUE, parent = parent.frame(), size = 29L)
	data("dji30ret", envir=newv)
	dji30ret<-get("dji30ret", envir = newv)	
	dji30ret = as.matrix(dji30ret)
	ctrl = list(rho = 1, delta = 1e-9, outer.iter = 100, inner.iter = 650, tol = 1e-9, trace = 0)
	benchmark = vector(mode="list")
	benchmark$sgarch = vector(mode="list", length = 3)
	benchmark$aparch = vector(mode="list", length = 1)
	benchmark$gjrgarch = vector(mode="list", length = 1)
	benchmark$egarch = vector(mode="list", length = 1)
	cat("\n...calculating rugarch models...please wait...\n")
	
	spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(2,2), include.mean = TRUE), 
			distribution.model = "std")
	sgarch.1102210002 = ugarchfit(data =  dji30ret[,"AA",drop=FALSE], spec = spec, 
			solver = "solnp", solver.control = ctrl)
	
	benchmark$sgarch[[1]]$model = 1102210002
	benchmark$sgarch[[1]]$pars = cbind(coef(sgarch.1102210002), 
			c( 0.000373, 0.418913, 0.519832, -0.369890, -0.587810, 2.8405e-06, 0.044669, 0.949917, 6.727497))
	colnames(benchmark$sgarch[[1]]$pars)=c("rugarch", "benchmark")
	benchmark$sgarch[[1]]$se = cbind(sgarch.1102210002@fit$matcoef[,2],
			c(0.00015827, 0.078308, 0.074981, 0.073396, 0.070594, 1.0829e-06, 0.0072783, 0.0087306, 0.65459))
	colnames(benchmark$sgarch[[1]]$se)=c("rugarch", "benchmark")
	benchmark$sgarch[[1]]$llh = c(likelihood(sgarch.1102210002), 14004.118)
	names(benchmark$sgarch[[1]]$llh)=c("rugarch", "benchmark")
	benchmark$sgarch[[1]]$time = c(round(as.numeric(sgarch.1102210002@fit$timer),3), 2.877)
	names(benchmark$sgarch[[1]]$time)=c("rugarch", "benchmark")
	
	cat("\nBenchmark Test 1/5: mu-ARMA(2,2)-sGARCH(1,1)-std\n")
	cat("\nparameters:\n")
	print(t(benchmark$sgarch[[1]]$pars), digits=4)
	cat("\nstandard errors:\n")
	print(t(benchmark$sgarch[[1]]$se), digits=4)
	cat("\nLogLikelihood:\n")
	print((benchmark$sgarch[[1]]$llh), digits=4)
	cat("\nTimings:\n")
	print((benchmark$sgarch[[1]]$time), digits=4)
	
	spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(2,2), include.mean = TRUE, archm = TRUE, 
					archpow = 2), distribution.model = "std")
	sgarch.1102211202 = ugarchfit(data =  dji30ret[,"AA",drop=FALSE], spec = spec, 
			solver = "solnp", solver.control = ctrl)
	benchmark$sgarch[[2]]$model = 1102211202
	benchmark$sgarch[[2]]$pars = cbind(coef(sgarch.1102211202), 
			c( 0.000608, 0.420190, 0.520847,-0.371267,-0.588981, -0.698302, 2.8905e-06, 0.044969, 0.949498, 6.724380))
	colnames(benchmark$sgarch[[2]]$pars)=c("rugarch", "benchmark")
	benchmark$sgarch[[2]]$se = cbind(sgarch.1102211202@fit$matcoef[,2],
			c(0.00026683, 0.077325, 0.074642, 0.072494, 0.070277, 0.65505, 1.0995E-06, 0.0073556, 0.0088467, 0.6534))
	colnames(benchmark$sgarch[[2]]$se)=c("rugarch", "benchmark")
	benchmark$sgarch[[2]]$llh = c(likelihood(sgarch.1102211202), 14004.667)
	names(benchmark$sgarch[[2]]$llh)=c("rugarch", "benchmark")
	benchmark$sgarch[[2]]$time = c(round(as.numeric(sgarch.1102211202@fit$timer),3), 4.016)
	names(benchmark$sgarch[[2]]$time)=c("rugarch", "benchmark")
	
	cat("\nBenchmark Test 2/5: mu-ARMA(2,2)-inmean(2)-sGARCH(1,1)-std\n")
	cat("\nparameters:\n")
	print(t(benchmark$sgarch[[2]]$pars), digits=4)
	cat("\nstandard errors:\n")
	print(t(benchmark$sgarch[[2]]$se), digits=4)
	cat("\nLogLikelihood:\n")
	print((benchmark$sgarch[[2]]$llh), digits=4)
	cat("\nTimings:\n")
	print((benchmark$sgarch[[2]]$time), digits=4)
	
	spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(2,2), include.mean = TRUE, archm = TRUE, 
					archpow = 1), distribution.model = "std")
	
	sgarch.1102211102 = ugarchfit(data =  dji30ret[,"AA",drop=FALSE], spec = spec, 
			solver = "solnp", solver.control = ctrl)
	
	benchmark$sgarch[[3]]$model = 1102211102
	benchmark$sgarch[[3]]$pars = cbind(coef(sgarch.1102211102), 
			c(0.000966, 0.41983, 0.520097, -0.370886, -0.588035, -0.033642, 2.9003E-06, 0.045041, 0.949407, 6.721048))
	colnames(benchmark$sgarch[[3]]$pars)=c("rugarch", "benchmark")
	benchmark$sgarch[[3]]$se = cbind(sgarch.1102211102@fit$matcoef[,2],
			c(0.00065137, 0.077789, 0.074867, 0.072935, 0.07049, 0.035944, 1.1066E-06, 0.0073928, 0.0089013, 0.65312))
	colnames(benchmark$sgarch[[3]]$se)=c("rugarch", "benchmark")
	benchmark$sgarch[[3]]$llh = c(likelihood(sgarch.1102211102),  14004.571)
	names(benchmark$sgarch[[3]]$llh)=c("rugarch", "benchmark")
	benchmark$sgarch[[3]]$time = c(round(as.numeric(sgarch.1102211102@fit$timer),3) , 3.728)
	names(benchmark$sgarch[[3]]$time)=c("rugarch", "benchmark")
	
	cat("\nBenchmark Test 3/5: mu-ARMA(2,2)-inmean(1)-sGARCH(1,1)-std\n")
	cat("\nparameters:\n")
	print(t(benchmark$sgarch[[3]]$pars), digits=4)
	cat("\nstandard errors:\n")
	print(t(benchmark$sgarch[[3]]$se), digits=4)
	cat("\nLogLikelihood:\n")
	print((benchmark$sgarch[[3]]$llh), digits=4)
	cat("\nTimings:\n")
	print((benchmark$sgarch[[3]]$time), digits=4)
	
	spec = ugarchspec(variance.model = list(model = "apARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std")
	
	aparch.1101110002 = ugarchfit(data = dji30ret[,"AA",drop=FALSE], spec = spec, 
			solver = "solnp", solver.control = list(tol=1e-8, trace=0))
	
	benchmark$aparch[[1]]$model = 1101110002
	benchmark$aparch[[1]]$pars = cbind(coef(aparch.1101110002), 
			c( 0.000175, -0.585745, 0.6407, 6.42972E-05, 0.051403,  0.952363, 0.292031, 1.224256, 7.118846))
	colnames(benchmark$aparch[[1]]$pars)=c("rugarch", "benchmark")
	benchmark$aparch[[1]]$se = cbind(aparch.1101110002@fit$matcoef[,2],
			c(0.00024266, 0.073315, 0.068083, 0.000037743, 0.0075942, 0.0081283, 0.077745, 0.1373, 0.72147))
	colnames(benchmark$aparch[[1]]$se)=c("rugarch", "benchmark")
	benchmark$aparch[[1]]$llh = c(likelihood(aparch.1101110002), 14011.7)
	names(benchmark$aparch[[1]]$llh)=c("rugarch", "benchmark")
	benchmark$aparch[[1]]$time = c(round(as.numeric(aparch.1101110002@fit$timer),3), 4.185)
	names(benchmark$aparch[[1]]$time)=c("rugarch", "benchmark")
	
	cat("\nBenchmark Test 4/5: mu-ARMA(1,1)-apARCH(1,1)-std\n")
	cat("\nparameters:\n")
	print(t(benchmark$aparch[[1]]$pars), digits=4)
	cat("\nstandard errors:\n")
	print(t(benchmark$aparch[[1]]$se), digits=4)
	cat("\nLogLikelihood:\n")
	print((benchmark$aparch[[1]]$llh), digits=4)
	cat("\nTimings:\n")
	print((benchmark$aparch[[1]]$time), digits=4)
	
	spec = ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std")
	
	gjrgarch.1101110002 = ugarchfit(data = dji30ret[,"AA",drop=FALSE], spec = spec, 
			solver = "solnp", solver.control = list(tol=1e-8, trace=0))
	
	
	benchmark$gjrgarch[[1]]$model = 1101110002
	benchmark$gjrgarch[[1]]$pars = cbind(coef(gjrgarch.1101110002), 
			c( 0.000228, -0.585394, 0.639836, 3.6631E-06, 0.03235,  0.945559, 0.030108,6.925304))
	colnames(benchmark$gjrgarch[[1]]$pars)=c("rugarch", "benchmark")
	benchmark$gjrgarch[[1]]$se = cbind(gjrgarch.1101110002@fit$matcoef[,2],
			c(0.00023834, 0.07242, 0.067309, 1.3806E-06, 0.0061025, 0.010113, 0.012071, 0.67969 ))
	colnames(benchmark$gjrgarch[[1]]$se)=c("rugarch", "benchmark")
	benchmark$gjrgarch[[1]]$llh = c(likelihood(gjrgarch.1101110002), 14002)
	names(benchmark$gjrgarch[[1]]$llh)=c("rugarch", "benchmark")
	benchmark$gjrgarch[[1]]$time = c(round(as.numeric(gjrgarch.1101110002@fit$timer),3), 2.324)
	names(benchmark$gjrgarch[[1]]$time)=c("rugarch", "benchmark")
	
	cat("\nBenchmark Test 5/5: mu-ARMA(1,1)-gjrGARCH(1,1)-std\n")
	cat("\nparameters:\n")
	print(t(benchmark$gjrgarch[[1]]$pars), digits=4)
	cat("\nstandard errors:\n")
	print(t(benchmark$gjrgarch[[1]]$se), digits=4)
	cat("\nLogLikelihood:\n")
	print((benchmark$gjrgarch[[1]]$llh), digits=4)
	cat("\nTimings:\n")
	print((benchmark$gjrgarch[[1]]$time), digits=4)
	return(invisible(0))
	
}
