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
# Model path Tests
#################################################################################

# ToDo: need to update it to take advantage of parallel...

# ---------------------------------------------------------------------------------

# Generating GARCH model paths and testing the parameter distributions.
# The example will demonstrate how to use the ugarchpath function to simulate garch 
# models with fixed parameters and then test the distribution of those parameters by 
# fitting the simulations with the ugarchfit function. 

# Generating GARCH model paths and testing the parameter distributions.
# generate common seeds for all models:

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
# ---------------------------------------------------------------------------------
# The sGARCH model
# ---------------------------------------------------------------------------------

rugarch.test6a = function(cluster=NULL)
{
	#cat("\nrugarch-->test6-1: Path Test (sGARCH)\n")
	tic = Sys.time()
	
	data(dji30ret)
	rseed = rugarch.seeds(100, "test6a")
	fixpars = list(mu = 0.001,omega = 0.00001, ar1 = 0.7, ma1 = -0.3, alpha1 = 0.05, 
			beta1 = 0.80, skew = 2, shape = 4)
	
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE, 
					archm = FALSE, archpow = 2), 
			distribution.model = "sstd", fixed.pars = fixpars)
	
	# generate the various paths
	sgarch.path = ugarchpath(spec, n.sim = 4000, n.start = 1, m.sim = 100, 
			rseed = rseed)
	
	# create a specification for fitting (same as before without fixed parameters)
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE, 
					archm =FALSE, archpow = 2), 
			distribution.model = "sstd")
	
	#setup coefficient matrix
	# for this example we will use the option of scaling the data as it leads
	# to better results with the solnp solver. The nlminb does not require this.
	fit = vector(mode = "list", length = 100)
	simcoef.s = matrix(NA, ncol = 8, nrow = 100)
	path.df = fitted(sgarch.path)
	
	fit = apply(path.df, 2, FUN = function(x){
				ugarchfit(data = x, spec = spec, solver = "solnp", 
						fit.control = list(scale = 1))
			})
		
	for(i in 1:100){
		if(fit[[i]]@fit$convergence == 0) simcoef.s[i,] = coef(fit[[i]])
	}
	
	# remove those which did not converge
	exc = which(is.na(simcoef.s[,1]))
	if(length(exc) > 0) exc = exc else exc = -1:-100
	postscript("test6a.eps", width = 12, height = 8)
	par(mfrow=c(3,3))

	plot(density(simcoef.s[-exc,1]), 
			main = "ARMA: mu parameter\n true parameter=001")
	plot(density(simcoef.s[-exc,2]), 
			main = "ARMA: ar parameter\n true parameter=0.7")
	plot(density(simcoef.s[-exc,3]), 
			main = "ARMA: ma parameter\n true parameter=-0.3")
	plot(density(simcoef.s[-exc,4]), 
			main = "sGARCH: omega parameter\n true parameter=0.00001")
	plot(density(simcoef.s[-exc,5]), 
			main = "sGARCH: alpha parameter\n true parameter=0.05")
	plot(density(simcoef.s[-exc,6]), 
			main = "sGARCH: beta parameter\n true parameter=0.8")
	plot(density(simcoef.s[-exc,7]), 
			main = "Skew Student Distribution: skew parameter\n true parameter=2")
	plot(density(simcoef.s[-exc,8]), 
			main = "Skew Student Distribution: shape parameter\n true parameter=4")
	dev.off()
	
	colnames(simcoef.s) = c("mu", "ar1", "ma1", "omega", "alpha1", "beta1", 
			"skew", "shape")
	rownames(simcoef.s) = 1:100
	zz <- file("test6a.txt", open="wt")
	sink(zz)
	print(round(simcoef.s,5))
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

# ---------------------------------------------------------------------------------
# The gjrGARCH model
rugarch.test6b = function(cluster=NULL)
{
	#cat("\nrugarch-->test6-2: Path Test (gjrGARCH)\n")
	tic = Sys.time()
	
	rseed = rugarch.seeds(100, "test6b")
	fixpars = list(mu = 0.001,omega = 0.00001, ar1 = 0.7, ma1 = -0.3, 
			alpha1 = 0.05, gamma1 = 0.2, beta1 = 0.80, skew = 2, shape = 4)
	
	spec = ugarchspec(
			variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "sstd", fixed.pars = fixpars)
	
	# generate the various paths
	gjrgarch.path = ugarchpath(spec, n.sim = 4000, n.start = 1, m.sim = 100, 
			rseed = rseed)
	
	# create a specification for fitting (same as before without fixed parameters)
	spec = ugarchspec(
			variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "sstd")
	
	#setup coefficient matrix
	simcoef.gjr = matrix(NA, ncol = 9, nrow = 100)
	path.df = fitted(gjrgarch.path)
	if(!is.null(cluster)){
		parallel::clusterEvalQ(cluster, require(rugarch))
		parallel::clusterExport(cluster, c("path.df", "spec"), envir = environment())
		fit = parallel::parLapply(cluster, as.list(1:100), fun = function(i){
					ugarchfit(data = path.df[,i], spec = spec, solver = "solnp",
							fit.control = list(scale = 1))
				})
	} else{
		fit = apply(path.df, 2, FUN = function(x){
				ugarchfit(data = x, spec = spec, solver = "solnp",
						fit.control = list(scale = 1))
			})
	}
	#as.numeric(sapply(fit, FUN=function(x) x@fit$convergence))
	
	for(i in 1:100){
		if(fit[[i]]@fit$convergence == 0) simcoef.gjr[i,] = coef(fit[[i]])
	}
	
	# remove those which did not converge
	exc = which(is.na(simcoef.gjr[,1]))
	if(length(exc) > 0) exc = exc else exc = -1:-100
	postscript("test6b.eps", width = 12, height = 8)
	par(mfrow=c(3,3))
	plot(density(simcoef.gjr[-exc,1]), 
			main = "ARMA: mu parameter\n true parameter=0.001")
	plot(density(simcoef.gjr[-exc,2]), 
			main = "ARMA: ar parameter\n true parameter=0.7")
	plot(density(simcoef.gjr[-exc,3]), 
			main = "ARMA: ma parameter\n true parameter=-0.3")
	plot(density(simcoef.gjr[-exc,4]), 
			main = "gjrGARCH: omega parameter\n true parameter=0.00001")
	plot(density(simcoef.gjr[-exc,5]), 
			main = "gjrGARCH: alpha parameter\n true parameter=0.05")
	plot(density(simcoef.gjr[-exc,6]), 
			main = "gjrGARCH: beta parameter\n true parameter=0.8")
	plot(density(simcoef.gjr[-exc,7]), 
			main = "gjrGARCH: gamma parameter\n true parameter=0.2")
	plot(density(simcoef.gjr[-exc,8]), 
			main = "Student Distribution: skew parameter\n true parameter=2")
	plot(density(simcoef.gjr[-exc,9]), 
			main = "Student Distribution: shape parameter\n true parameter=4")
	dev.off()
	
	colnames(simcoef.gjr) = c("mu", "ar1", "ma1", "omega", "alpha1", 
			"beta1", "gamma1", "skew", "shape")
	rownames(simcoef.gjr) = 1:100
	zz <- file("test6b.txt", open="wt")
	sink(zz)
	print(round(simcoef.gjr,5))
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

# ---------------------------------------------------------------------------------
# The eGARCH model
# we generate a realistic parameter set by fitting from a dataset
rugarch.test6c = function(cluster=NULL)
{
	#cat("\nrugarch-->test6-3: Path Test (eGARCH)\n")
	tic = Sys.time()
	
	rseed = rugarch.seeds(100, "test6c")
	data(dmbp)
	spec = ugarchspec(
			variance.model = list(model = "eGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "sstd")
	fit = ugarchfit(data = dmbp[,1], spec = spec, solver = "solnp")
	
	fixpars = as.list(coef(fit))
	um = uncmean(fit)
	
	spec = ugarchspec(
			variance.model = list(model = "eGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "sstd", fixed.pars = fixpars)
	
	# generate the various paths
	egarch.path = ugarchpath(spec, n.sim = 4000, n.start = 1, m.sim = 100, 
			rseed = rseed)
	
	# create a specification for fitting (same as before without fixed parameters)
	spec = ugarchspec(
			variance.model = list(model = "eGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "sstd")
	
	#setup coefficient matrix
	simcoef.e = matrix(NA, ncol = 9, nrow = 100)
	path.df = fitted(egarch.path)
	
	if(!is.null(cluster)){
		parallel::clusterEvalQ(cluster, require(rugarch))
		parallel::clusterExport(cluster, c("path.df", "spec"), envir = environment())
		fit = parallel::parLapply(cluster, as.list(1:100), fun = function(i){
					ugarchfit(data = path.df[,i], spec = spec, solver = "solnp",
							fit.control = list(scale = 1))
				})
	} else{
		fit = apply(path.df, 2, FUN = function(x){
					ugarchfit(data = x, spec = spec, solver = "solnp",
							fit.control = list(scale = 1))
				})
	}	
	for(i in 1:100){
		if(fit[[i]]@fit$convergence == 0) simcoef.e[i,] = coef(fit[[i]])
	}
	
	# remove those which did not converge
	exc = which(is.na(simcoef.e[,1]))
	if(length(exc) > 0) exc = exc else exc = -1:-100
	postscript("test6c.eps", width = 12, height = 8)
	par(mfrow=c(3,3))
	plot(density(simcoef.e[-exc,1]), 
			main = paste("ARMA: mu parameter\n true parameter=", round(um, 4),sep=""))
	plot(density(simcoef.e[-exc,2]), 
			main = "ARMA: ar parameter\n true parameter=-0.548")
	plot(density(simcoef.e[-exc,3]), 
			main = "ARMA: ma parameter\n true parameter=0.581")
	plot(density(simcoef.e[-exc,4]), 
			main = "eGARCH: omega parameter\n true parameter=-0.039")
	plot(density(simcoef.e[-exc,5]), 
			main = "eGARCH: alpha parameter\n true parameter= -0.041")
	plot(density(simcoef.e[-exc,6]), 
			main = "eGARCH: beta parameter\n true parameter=0.976")
	plot(density(simcoef.e[-exc,7]), 
			main = "eGARCH: gamma parameter\n true parameter=0.258")
	plot(density(simcoef.e[-exc,8]), 
			main = "Student Distribution: skew parameter\n true parameter=0.909")
	plot(density(simcoef.e[-exc,9]), 
			main = "Student Distribution: shape parameter\n true parameter= 4.21")
	dev.off()
	
	colnames(simcoef.e) = c("mu", "ar1", "ma1", "omega", "alpha1", "beta1", 
			"gamma1", "skew", "shape")
	rownames(simcoef.e) = 1:100
	zz <- file("test6c.txt", open="wt")
	sink(zz)
	print(round(simcoef.e,5))
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

# ---------------------------------------------------------------------------------
# The apARCH model
rugarch.test6d = function(cluster=NULL)
{
	#cat("\nrugarch-->test6-4: Path Test (apARCH)\n")
	tic = Sys.time()
	
	rseed = rugarch.seeds(100, "test6d")
	fixpars = list(mu = 0.001,omega = 0.0008, ar1 = 0.7, ma1 = -0.3, 
			alpha1 = 0.1, gamma1 = 0.2, delta = 1.5, beta1 = 0.80, skew = 2, 
			shape = 4)
	
	spec = ugarchspec(
			variance.model = list(model = "apARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "sstd", fixed.pars = fixpars)
	
	# generate the various paths
	aparch.path = ugarchpath(spec, n.sim = 4000, n.start = 1, m.sim = 100, 
			rseed = rseed)
	
	spec = ugarchspec(
			variance.model = list(model = "apARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "sstd")
	
	#setup coefficient matrix
	simcoef.a = matrix(NA, ncol = 10, nrow = 100)
	path.df = fitted(aparch.path)
	
	if(!is.null(cluster)){
		parallel::clusterEvalQ(cluster, require(rugarch))
		parallel::clusterExport(cluster, c("path.df", "spec"), envir = environment())
		fit = parallel::parLapply(cluster, as.list(1:100), fun = function(i){
					ugarchfit(data = path.df[,i], spec = spec, solver = "solnp",
							fit.control = list(scale = 1))
				})
	} else{
		fit = apply(path.df, 2, FUN = function(x){
					ugarchfit(data = x, spec = spec, solver = "solnp",
							fit.control = list(scale = 1))
				})
	}
	
	for(i in 1:100){
		if(fit[[i]]@fit$convergence == 0) simcoef.a[i,] = coef(fit[[i]])
	}
	# remove those which did not converge
	exc = which(is.na(simcoef.a[,1]))
	if(length(exc) > 0) exc = exc else exc = -1:-100
	postscript("test6d.eps", width = 12, height = 8)
	par(mfrow=c(3,4))
	plot(density(simcoef.a[-exc,1]), 
			main = "ARMA: mu parameter\n true parameter=0.001")
	plot(density(simcoef.a[-exc,2]), 
			main = "ARMA: ar parameter\n true parameter=0.7")
	plot(density(simcoef.a[-exc,3]), 
			main = "ARMA: ma parameter\n true parameter=-0.3")
	plot(density(simcoef.a[-exc,4]), 
			main = "apARCH: omega parameter\n true parameter=0.0008")
	plot(density(simcoef.a[-exc,5]), 
			main = "apARCH: alpha parameter\n true parameter=0.1")
	plot(density(simcoef.a[-exc,6]), 
			main = "apARCH: beta parameter\n true parameter=0.8")
	plot(density(simcoef.a[-exc,7]), 
			main = "apARCH: gamma parameter\n true parameter=0.2")
	plot(density(simcoef.a[-exc,8]), 
			main = "apARCH: delta parameter\n true parameter=1.5")
	plot(density(simcoef.a[-exc,9]), 
			main = "Student Distribution: skew parameter\n true parameter=2")
	plot(density(simcoef.a[-exc,10]), 
			main = "Student Distribution: shape parameter\n true parameter=4")
	dev.off()
	options(width=100)
	colnames(simcoef.a) = c("mu", "ar1", "ma1", "omega", "alpha1", "beta1", 
			"gamma1", "delta", "skew", "shape")
	rownames(simcoef.a) = 1:100
	zz <- file("test6d.txt", open="wt")
	sink(zz)
	print(round(simcoef.a,5))
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

# ---------------------------------------------------------------------------------
# The fGARCH model - NAGARCH
rugarch.test6e = function(cluster=NULL)
{
	#cat("\nrugarch-->test6-5: Path Test (fGARCH/NAGARCH)\n")
	tic = Sys.time()
	
	rseed = rugarch.seeds(100, "test6e")
	fixpars = list(mu = 0.001,omega = 0.00001, ar1 = 0.7, ma1 = -0.3, 
			alpha1 = 0.04, eta21 = 0.2, beta1 = 0.80, shape = 3)
	
	spec = ugarchspec(
			variance.model = list(model = "fGARCH", garchOrder = c(1,1), 
					submodel="NAGARCH"), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std", fixed.pars = fixpars)
	# generate the various paths
	fgarch.path = ugarchpath(spec, n.sim = 4000, n.start = 1, m.sim = 100, 
			rseed = rseed)
	# create a specification for fitting (same as before without fixed parameters)
	spec = ugarchspec(
			variance.model = list(model = "fGARCH", garchOrder = c(1,1), 
					submodel="NAGARCH"), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std")
	
	#setup coef-matrix
	simcoef.fna = matrix(NA, ncol = 8, nrow = 100)
	path.df = fitted(fgarch.path)
	
	if(!is.null(cluster)){
		parallel::clusterEvalQ(cluster, require(rugarch))
		parallel::clusterExport(cluster, c("path.df", "spec"), envir = environment())
		fit = parallel::parLapply(cluster, as.list(1:100), fun = function(i){
					ugarchfit(data = path.df[,i], spec = spec, solver = "solnp",
							fit.control = list(scale = 1))
				})
	} else{
		fit = apply(path.df, 2, FUN = function(x){
					ugarchfit(data = x, spec = spec, solver = "solnp",
							fit.control = list(scale = 1))
				})
	}
	
	for(i in 1:100){
		if(fit[[i]]@fit$convergence == 0) simcoef.fna[i,] = coef(fit[[i]])
	}
	# remove those which did not converge
	exc = which(is.na(simcoef.fna[,1]))
	if(length(exc) > 0) exc = exc else exc = -1:-100
	postscript("test6e.eps", width = 12, height = 8)
	par(mfrow=c(3,3))
	plot(density(simcoef.fna[-exc,1]), 
			main = "ARMA: mu parameter\n true parameter=0.001")
	plot(density(simcoef.fna[-exc,2]), 
			main = "ARMA: ar parameter\n true parameter=0.7")
	plot(density(simcoef.fna[-exc,3]), 
			main = "ARMA: ma parameter\n true parameter=-0.3")
	plot(density(simcoef.fna[-exc,4]), 
			main = "NAGARCH: omega parameter\n true parameter=0.00001")
	plot(density(simcoef.fna[-exc,5]), 
			main = "NAGARCH: alpha parameter\n true parameter=0.04")
	plot(density(simcoef.fna[-exc,6]), 
			main = "NAGARCH: beta parameter\n true parameter=0.8")
	plot(density(simcoef.fna[-exc,7]), 
			main = "NAGARCH: (asymmetry-2) parameter\n true parameter=0.2")
	plot(density(simcoef.fna[-exc,8]), 
			main = "Student Distribution: shape parameter\n true parameter=3")
	dev.off()
	options(width=100)
	colnames(simcoef.fna) = c("mu", "ar1", "ma1", "omega", "alpha1", "beta1", 
			"eta21", "shape")
	rownames(simcoef.fna) = 1:100
	zz = file("test6e.txt", open="wt")
	sink(zz)
	print(round(simcoef.fna,5))
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

# ---------------------------------------------------------------------------------
# The fGARCH model - NGARCH
rugarch.test6f = function(cluster=NULL)
{
	#cat("\nrugarch-->test6-6: Path Test (fGARCH/NGARCH)\n")
	tic = Sys.time()
	
	rseed = rugarch.seeds(100, "test6f")
	fixpars = list(mu = 0.001,omega = 0.0003, ar1 = 0.7, ma1 = -0.3, 
			alpha1 = 0.1, lambda = 1.5, beta1 = 0.80, shape = 3)
	
	spec = ugarchspec(
			variance.model = list(model = "fGARCH", garchOrder = c(1,1), 
					submodel = "NGARCH"), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std", fixed.pars = fixpars)
	
	# generate the various paths
	fgarch.path = ugarchpath(spec, n.sim = 4000, n.start = 1, m.sim = 100, 
			rseed = rseed)
	
	# create a specification for fitting (same as before without fixed parameters)
	spec = ugarchspec(
			variance.model = list(model = "fGARCH", garchOrder = c(1,1), 
					submodel = "NGARCH"), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std")
	
	#setup coef-matrix
	simcoef.fn = matrix(NA, ncol = 8, nrow = 100)
	path.df = fitted(fgarch.path)
	
	if(!is.null(cluster)){
		parallel::clusterEvalQ(cluster, require(rugarch))
		parallel::clusterExport(cluster, c("path.df", "spec"), envir = environment())
		fit = parallel::parLapply(cluster, as.list(1:100), fun = function(i){
					ugarchfit(data = path.df[,i], spec = spec, solver = "solnp",
							fit.control = list(scale = 1))
				})
	} else{
		fit = apply(path.df, 2, FUN = function(x){
					ugarchfit(data = x, spec = spec, solver = "solnp",
							fit.control = list(scale = 1))
				})
	}
	
	for(i in 1:100){
		if(fit[[i]]@fit$convergence == 0) simcoef.fn[i,] = coef(fit[[i]])
	}
	# remove those which did not converge
	exc = which(is.na(simcoef.fn[,1]))
	if(length(exc) > 0) exc = exc else exc = -1:-100
	postscript("test6f.eps", width = 12, height = 8)
	par(mfrow=c(3,3))
	plot(density(simcoef.fn[-exc,1]), 
			main = "ARMA: mu parameter\n true parameter=0.001")
	plot(density(simcoef.fn[-exc,2]), 
			main = "ARMA: ar parameter\n true parameter=0.7")
	plot(density(simcoef.fn[-exc,3]), 
			main = "ARMA: ma parameter\n true parameter=-0.3")
	plot(density(simcoef.fn[-exc,4]), 
			main = "NGARCH: omega parameter\n true parameter=0.0003")
	plot(density(simcoef.fn[-exc,5]), 
			main = "NGARCH: alpha parameter\n true parameter=0.1")
	plot(density(simcoef.fn[-exc,6]), 
			main = "NGARCH: beta parameter\n true parameter=0.8")
	plot(density(simcoef.fn[-exc,7]), 
			main = "NGARCH: lambda (nonlinear) parameter\n true parameter=1.5")
	plot(density(simcoef.fn[-exc,8]), 
			main = "Student Distribution: shape parameter\n true parameter=3")
	dev.off()
	options(width=100)
	colnames(simcoef.fn) = c("mu", "ar1", "ma1", "omega", "alpha1", "beta1", 
			"lambda", "shape")
	rownames(simcoef.fn) = 1:100
	zz = file("test6f.txt", open="wt")
	sink(zz)
	print(round(simcoef.fn,5))
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

# ---------------------------------------------------------------------------------
# The fGARCH model - AVGARCH
rugarch.test6g = function(cluster=NULL)
{
	#cat("\nrugarch-->test6-7: Path Test (fGARCH/AVGARCH)\n")
	tic = Sys.time()
	
	rseed = rugarch.seeds(100, "test6g")
	fixpars = list(mu = 0.001,omega = 0.003, ar1 = 0.7, ma1 = -0.3, 
			alpha1 = 0.1, eta21 = -0.2, eta11 = 0.4, beta1 = 0.80, shape = 3)
	
	spec = ugarchspec(
			variance.model = list(model = "fGARCH", garchOrder = c(1,1), 
					submodel = "AVGARCH"), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std", fixed.pars = fixpars)
	# generate the various paths
	fgarch.path = ugarchpath(spec, n.sim = 4000, n.start = 1, m.sim = 100, 
			rseed = rseed)
	
	# create a specification for fitting (same as before without fixed parameters)
	spec = ugarchspec(
			variance.model = list(model = "fGARCH", garchOrder = c(1,1), 
					submodel = "AVGARCH"), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std")
	
	#setup coef-matrix
	simcoef.av = matrix(NA, ncol = 9, nrow = 100)
	path.df = fitted(fgarch.path)
	
	if(!is.null(cluster)){
		parallel::clusterEvalQ(cluster, require(rugarch))
		parallel::clusterExport(cluster, c("path.df", "spec"), envir = environment())
		fit = parallel::parLapply(cluster, as.list(1:100), fun = function(i){
					ugarchfit(data = path.df[,i], spec = spec, solver = "solnp",
							fit.control = list(scale = 1))
				})
	} else{
		fit = apply(path.df, 2, FUN = function(x){
					ugarchfit(data = x, spec = spec, solver = "solnp",
							fit.control = list(scale = 1))
				})
	}
	
	for(i in 1:100){
		if(fit[[i]]@fit$convergence == 0) simcoef.av[i,] = coef(fit[[i]])
	}
	
	# remove those which did not converge
	exc = which(is.na(simcoef.av[,1]))
	if(length(exc) > 0) exc = exc else exc = -1:-100
	postscript("test6g.eps", width = 12, height = 8)
	par(mfrow=c(3,3))
	plot(density(simcoef.av[-exc,1]), 
			main = "ARMA: mu parameter\n true parameter=0.001")
	plot(density(simcoef.av[-exc,2]), 
			main = "ARMA: ar parameter\n true parameter=0.7")
	plot(density(simcoef.av[-exc,3]), 
			main = "ARMA: ma parameter\n true parameter=-0.3")
	plot(density(simcoef.av[-exc,4]), 
			main = "AVGARCH: omega parameter\n true parameter=0.003")
	plot(density(simcoef.av[-exc,5]), 
			main = "AVGARCH: alpha parameter\n true parameter=0.1")
	plot(density(simcoef.av[-exc,6]), 
			main = "AVGARCH: beta parameter\n true parameter=0.8")
	plot(density(simcoef.av[-exc,7]), 
			main = "AVGARCH: (assymetry-1) parameter\n true parameter=0.4")
	plot(density(simcoef.av[-exc,8]), 
			main = "AVGARCH: (assymetry-1) parameter\n true parameter=-0.2")
	plot(density(simcoef.av[-exc,9]), 
			main = "Student Distribution: shape parameter\n true parameter=3")
	dev.off()
	options(width=100)
	colnames(simcoef.av) = c("mu", "ar1", "ma1", "omega", "alpha1", "beta1", 
			"eta11", "eta21", "shape")
	rownames(simcoef.av) = 1:100
	zz = file("test6g.txt", open="wt")
	sink(zz)
	print(round(simcoef.av,5))
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}