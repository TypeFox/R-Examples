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
# GARCH archex tests (new functionality in mean.model)
#################################################################################

# Fit and Forecast
rugarch.test9a = function(cluster=NULL)
{
	tic = Sys.time()	
	data(dji30ret)
	# Individual Stock
	X = dji30ret[1:1000,2,drop=FALSE]
	# The 'Market'
	M = apply(dji30ret[1:1000,], 1, "mean")
	# Add Random Regressor, and scale M
	Y = cbind( matrix(rnorm(1000, 0, 1), ncol = 1), M/sd(M))
		
	spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(0,0), include.mean = TRUE, external.regressors = as.matrix(Y), 
					archex = 1), distribution.model = "std")
	fit = ugarchfit(spec, X, out.sample = 10)
	
	Mforc = apply(dji30ret[1001:1010,], 1, "mean")
	# Add Random Regressor, and scale M
	Yforc = cbind( matrix(rnorm(10, 0, 1), ncol = 1), Mforc/sd(Mforc))
	
	forc1 = ugarchforecast(fit, n.ahead = 1, n.roll = 9, external.forecasts = list(mregfor = Yforc))
	
	# check:
	Ycheck1 = rep(0,10)
	sigmafor = sigma(forc1)
	for(i in 1:10){
		Ycheck1[i] = coef(fit)["mu"] + coef(fit)["mxreg1"]*Yforc[i,1] + coef(fit)["mxreg2"]*Yforc[i,2]*sigmafor[i]
	}
	
	# Individual Stock
	X = dji30ret[1:1000,1,drop=FALSE]
	# The 'Market'
	M = apply(dji30ret[1:1000,], 1, "mean")
	# Add 2 Random Regressor, and scale M
	Y = cbind( matrix(rnorm(2000, 0, 1), ncol = 2), M/sd(M))
	
	spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(0,0), include.mean = TRUE, external.regressors = as.matrix(Y), 
					archex = 2), distribution.model = "std")
	fit = ugarchfit(spec, X, out.sample = 10)
	
	Mforc = apply(dji30ret[1001:1010,], 1, "mean")
	# Add Random Regressor, and scale M
	Yforc = cbind( matrix(rnorm(20, 0, 1), ncol = 2), Mforc/sd(Mforc))
	
	forc2 = ugarchforecast(fit, n.ahead = 1, n.roll = 9, external.forecasts = list(mregfor = Yforc))
	
	# check:
	Ycheck2 = rep(0,10)
	sigmafor = sigma(forc2)
	for(i in 1:10){
		Ycheck2[i] = coef(fit)["mu"] + coef(fit)["mxreg1"]*Yforc[i,1] + coef(fit)["mxreg2"]*Yforc[i,2]*sigmafor[i] + coef(fit)["mxreg3"]*Yforc[i,3]*sigmafor[i]
	}	
	
	# Individual Stock
	X = dji30ret[1:1000,3,drop=FALSE]
	# The 'Market'
	M = apply(dji30ret[1:1000,], 1, "mean")
	# Add 2 Random Regressor, and scale M
	Y = cbind( M/sd(M) )
	spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(0,0), include.mean = TRUE, external.regressors = as.matrix(Y), 
					archex = 1), distribution.model = "std")
	fit = ugarchfit(spec, X, out.sample = 10)
	
	Mforc = apply(dji30ret[1001:1010,], 1, "mean")
	# Add Random Regressor, and scale M
	Yforc = cbind(Mforc/sd(Mforc))
	
	forc3 = ugarchforecast(fit, n.ahead = 1, n.roll = 9, external.forecasts = list(mregfor = Yforc))
	
	# check:
	Ycheck3 = rep(0,10)
	sigmafor = sigma(forc3)
	for(i in 1:10){
		Ycheck3[i] = coef(fit)["mu"] + coef(fit)["mxreg1"]*Yforc[i,1]*sigmafor[i]
	}

	z4 <- file("test9a.txt", open="wt")
	sink(z4)
	cat("\n[1] Forecast Check with archex.\n Pass:\n")
	print( all.equal(Ycheck1, as.numeric(fitted(forc1))))
	cat("\n[2] Forecast Check with archex.\n Pass:\n")
	print( all.equal(Ycheck2,  as.numeric(fitted(forc2))))
	cat("\n[3] Forecast Check with archex.\n Pass:\n")
	print( all.equal(Ycheck3,  as.numeric(fitted(forc3))))
	sink(type="message")
	sink()
	close(z4)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}


# Filter and Forecast from Spec
rugarch.test9b = function(cluster=NULL)
{
	tic = Sys.time()	
	data(dji30ret)
	# Individual Stock
	X = dji30ret[1:1000,2,drop=FALSE]
	# The 'Market'
	M = apply(dji30ret[1:1000,], 1, "mean")
	# Add Random Regressor, and scale M
	Y = cbind( matrix(rnorm(1000, 0, 1), ncol = 1), M/sd(M))
	
	spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(0,0), include.mean = TRUE, external.regressors = as.matrix(Y), 
					archex = 1), distribution.model = "std")
	fit = ugarchfit(spec, X, out.sample = 10)
	spec2 = spec
	setfixed(spec2)<-as.list(coef(fit))
	filt = ugarchfilter(spec2, X, out.sample = 10)
	
	Mforc = apply(dji30ret[1001:1010,], 1, "mean")
	# Add Random Regressor, and scale M
	Yforc = cbind( matrix(rnorm(10, 0, 1), ncol = 1), Mforc/sd(Mforc))
	
	forc1 = ugarchforecast(fit, n.ahead = 1, n.roll = 9, external.forecasts = list(mregfor = Yforc))
	forc2 = ugarchforecast(spec2, data = X, out.sample = 10, n.ahead = 1, n.roll = 9, external.forecasts = list(mregfor = Yforc))
	
	z4 <- file("test9b.txt", open="wt")
	sink(z4)
	cat("\nFilter Test.\nPass:\n")
	print(all.equal(fitted(fit), fitted(filt)))
	cat("\nForecast Test.\nPass:\n")
	print(all.equal(as.numeric(fitted(forc1)), as.numeric(fitted(forc2))))
	sink(type="message")
	sink()
	close(z4)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

# Simulation1
rugarch.test9c = function(cluster=NULL)
{
	tic = Sys.time()	
	N <- 2000
	Ncoefs <- matrix(rep(NA, N * 4), ncol = 4)
	Nstder <- matrix(rep(NA, N * 4), ncol = 4)
	Y = matrix(NA, ncol = 2000, nrow = 1000)
	spec = vector(mode="list", length = N)
	# first do a "manual simulation":
	for(i in 1:N) {
		H <- 1000
		residual <- rnorm(H,0,1)
		vcoef <- c(0.001,0.4,0.55)
		mcoef <- c(1)
		exreg <- rnorm(H,0,1)
		variance <- rep(NA,H)
		variance[1] <- vcoef[1] / (1 - vcoef[2] - vcoef[3])
		Y[1,i] <- mcoef[1] * sqrt(variance[1]) * exreg[1] + sqrt(variance[1]) * residual[1]
		for(j in 2:H) {
			variance[j] <- vcoef[1] + vcoef[2] * variance[j-1] + vcoef[3] * variance[j-1] * residual[j-1]^2
			Y[j,i] <- mcoef[1] * sqrt(variance[j]) * exreg[j] + sqrt(variance[j]) * residual[j]
		}
		spec[[i]] <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1))
				, mean.model = list(armaOrder = c(0,0), include.mean = FALSE
						, external.regressors = as.matrix(exreg), archex = 1)
				, distribution.model = "norm")
	}
	mspec = multispec(spec)
	mfit = multifit(mspec, data = Y, cluster = cluster)
	for(i in 1:N){
		if(convergence(mfit@fit[[i]]) == 0){
			Ncoefs[i, ] <- coef(mfit@fit[[i]])
			Nstder[i, ] <- sqrt(diag(abs(vcov(mfit@fit[[i]]))))
		}
	}
	
	z4 <- file("test9c.txt", open="wt")
	sink(z4)
	cat("\nCoefs:\n")
	tmp1 = cbind(c(1, 0.001,0.55,0.4), colMeans(Ncoefs, na.rm = TRUE))
	tmp2 = matrix(colMeans(Nstder, na.rm = TRUE), ncol = 1)
	colnames(tmp2) = c("Simulated s.e.")
	colnames(tmp1) = c("True", "Simulated")
	rownames(tmp1) = rownames(tmp2) = c("mxreg1", "omega", "alpha1", "beta1")
	cat("\nPath Simulation\n")
	print(tmp1)
	cat("\n")
	print(tmp2)	
	sink(type="message")
	sink()
	close(z4)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)	
}


# Simulation2
rugarch.test9d = function(cluster=NULL)
{
	tic = Sys.time()	
	N <- 2000
	Ncoefs <- matrix(rep(NA, N * 4), ncol = 4)
	Nstder <- matrix(rep(NA, N * 4), ncol = 4)
	Y = matrix(NA, ncol = 2000, nrow = 1000)
	spec = vector(mode="list", length = N)
	for(i in 1:N) {
		H <- 1000
		exreg <- rnorm(H,0,1)
		spec[[i]] <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1))
				, mean.model = list(armaOrder = c(0,0), include.mean = FALSE
						, external.regressors = as.matrix(exreg), archex = 1),
						fixed.pars = list(mxreg1 = 1, omega = 0.001, alpha1 = 0.55, beta1 = 0.4)
						, distribution.model = "norm")
		sim = ugarchpath(spec[[i]], n.sim = H, mexsimdata = list(as.matrix(exreg)))
		Y[,i] = as.numeric( fitted(sim)[,1] )
		rm(sim)
		spec[[i]] <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1))
				, mean.model = list(armaOrder = c(0,0), include.mean = FALSE
						, external.regressors = as.matrix(exreg), archex = 1),
				, distribution.model = "norm")
	}
	mspec = multispec(spec)
	mfit = multifit(mspec, data = Y, cluster = cluster)
	for(i in 1:N){
		if(convergence(mfit@fit[[i]]) == 0){
			Ncoefs[i, ] <- coef(mfit@fit[[i]])
			Nstder[i, ] <- sqrt(diag(abs(vcov(mfit@fit[[i]]))))
		}
	}
	
	z4 <- file("test9d.txt", open="wt")
	sink(z4)
	cat("\nCoefs:\n")
	tmp1 = cbind(c(1, 0.001,0.55,0.4), colMeans(Ncoefs, na.rm = TRUE))
	tmp2 = matrix(colMeans(Nstder, na.rm = TRUE), ncol = 1)
	colnames(tmp2) = c("Simulated s.e.")
	colnames(tmp1) = c("True", "Simulated")
	rownames(tmp1) = rownames(tmp2) = c("mxreg1", "omega", "alpha1", "beta1")
	cat("\nPath Simulation\n")
	print(tmp1)
	cat("\n")
	print(tmp2)	
	sink(type="message")
	sink()
	close(z4)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}
