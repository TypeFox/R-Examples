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
# Fixed and Starting Parameter Testing
#################################################################################

rugarch.test2a = function(cluster=NULL)
{
	tic = Sys.time()
	
	#cat("\nrugarch-->test2-1: Fixed and Starting Parameter Test (fGARCH/ALLGARCH)\n")
	data(sp500ret)
	# No Parameters Fixed, Starting Parameters
	spars = list(mu = 1.9917e-04, ar1 = -1.7519e-02, omega = 6.5805e-05, 
			alpha1 = 6.0165e-02, beta1 = 9.3376e-01, lambda = 1.1702e+00,  
			eta21 = 4.2051e-02, eta11 = 7.9775e-01, skew = -1.4682e-01, 
			shape = 2.0338e+00)
	
	spec = ugarchspec(
			variance.model = list(model = "fGARCH", garchOrder = c(1,1), 
					submodel = "ALLGARCH"), 
			mean.model = list(armaOrder = c(1,0), include.mean = TRUE, 
					archm = FALSE, archpow = 2), 
			distribution.model = "nig")
	setstart(spec)<-spars
	
	fgarch.fit1 = ugarchfit(data = sp500ret, spec = spec, solver = "solnp")
	
	# Some Parameters Fixed (we obtain from fgarch.fit1)
	fpars = as.list(coef(fgarch.fit1)[1:6])
	
	spec = ugarchspec(
			variance.model = list(model = "fGARCH", garchOrder = c(1,1), 
					submodel = "ALLGARCH"), 
			mean.model = list(armaOrder = c(1,0), include.mean = TRUE, 
					archm = FALSE, archpow = 2), distribution.model = "nig", 
			fixed.pars = fpars)
	# alternative use setfixed(spec)<-fpars
	fgarch.fit2 = ugarchfit(data = sp500ret, spec = spec, solver = "solnp", 
			fit.control=list(scale=1))
	
	# notice that the standard errors of the fixed parameters is NA (they are fixed!).
	# However, we can calculate those post estimation with the fixed.se argument:
	
	fgarch.fit3 = ugarchfit(data = sp500ret, spec = spec, solver = "solnp", 
			fit.control = list(fixed.se = TRUE))
	
	# All Parameters Fixed (we obtain from fgarch.fit1)
	fpars = as.list(coef(fgarch.fit1))
	spec = ugarchspec(
			variance.model = list(model = "fGARCH", garchOrder = c(1,1), 
					submodel = "ALLGARCH"), 
			mean.model = list(armaOrder = c(1,0), include.mean = TRUE, 
					archm = FALSE, archpow = 2), distribution.model = "nig", 
			fixed.pars = fpars)
	
	fgarch.fit4 = ugarchfit(data = sp500ret,spec = spec, solver = "solnp", 
			fit.control = list(fixed.se = TRUE))
	
	# compare LLH, coefficients and std.errors:
	fgarch.lik = cbind(
			likelihood(fgarch.fit1), 
			likelihood(fgarch.fit2), 
			likelihood(fgarch.fit3), 
			likelihood(fgarch.fit4))
	fgarch.coefs = round(cbind(
					coef(fgarch.fit1), 
					coef(fgarch.fit2), 
					coef(fgarch.fit3), 
					coef(fgarch.fit4)),5)
	
	fgarch.se = round(cbind(
					fgarch.fit1@fit$matcoef[,2], 
					fgarch.fit2@fit$matcoef[,2], 
					fgarch.fit3@fit$matcoef[,2], 
					fgarch.fit4@fit$matcoef[,2]),5)
	
	rownames(fgarch.lik) = "likelihood"
	colnames(fgarch.lik) = paste("test", 1:4, sep=".")
	colnames(fgarch.coefs) = paste("test", 1:4, sep=".")
	colnames(fgarch.se) = paste("test", 1:4, sep=".")
	
	options(width=100)
	zz <- file("test2a.txt", open="wt")
	sink(zz)
	print(fgarch.lik)
	cat("\nparameters:\n")
	print(fgarch.coefs)
	cat("\nstd.errors:\n")
	print(fgarch.se)
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}


rugarch.test2b = function(cluster=NULL)
{
	tic = Sys.time()
	
	#cat("\nrugarch-->test2-2: Fixed and Starting Parameter Test (eGARCH)\n")	
	data(sp500ret)
	spars = list(mu = 1.9917e-04, ar1 = -1.7519e-02, omega = 6.5805e-05, 
			alpha1 = 6.0165e-02, beta1 = 9.3376e-01)
	
	spec = ugarchspec(
			variance.model = list(model = "eGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,0), include.mean = TRUE), 
			distribution.model = "nig", start.pars = spars)
	
	egarch.fit1 = ugarchfit(data = sp500ret, spec = spec, solver = "solnp")
	
	# Some Parameters Fixed (we obtain from egarch.fit1)
	fpars = as.list(coef(egarch.fit1)[1:6])
	
	spec = ugarchspec(
			variance.model = list(model = "eGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,0), include.mean = TRUE), 
			distribution.model = "nig", fixed.pars = fpars)
	
	egarch.fit2 = ugarchfit(data = sp500ret, spec = spec, solver = "hybrid")
	
	# notice that the standard errors of the fixed parameters is NA (they are fixed!).
	# However, we can calculate those post estimation with the fixed.se argument:
	
	egarch.fit3 = ugarchfit(data = sp500ret,spec = spec, solver = "hybrid", 
			fit.control = list(fixed.se = TRUE))
	
	# All Parameters Fixed (we obtain from egarch.fit1)
	fpars = as.list(coef(egarch.fit1))
	spec = ugarchspec(
			variance.model = list(model = "eGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,0), include.mean = TRUE), 
			distribution.model = "nig", fixed.pars = fpars)
	
	egarch.fit4 = ugarchfit(data = sp500ret,spec = spec, solver = "hybrid", 
			fit.control = list(fixed.se = TRUE))
	
	# compare LLH, coefficients and std.errors:
	egarch.lik = cbind(
			likelihood(egarch.fit1), 
			likelihood(egarch.fit2), 
			likelihood(egarch.fit3), 
			likelihood(egarch.fit4))
	
	egarch.coefs = round(cbind(
					coef(egarch.fit1), 
					coef(egarch.fit2), 
					coef(egarch.fit3), 
					coef(egarch.fit4)),5)
	
	egarch.se = round(cbind(
					egarch.fit1@fit$matcoef[,2], 
					egarch.fit2@fit$matcoef[,2], 
					egarch.fit3@fit$matcoef[,2], 
					egarch.fit4@fit$matcoef[,2]),5)
	
	rownames(egarch.lik) = "likelihood"
	colnames(egarch.lik) = paste("test", 1:4, sep=".")
	colnames(egarch.coefs) = paste("test", 1:4, sep=".")
	colnames(egarch.se) = paste("test", 1:4, sep=".")
	
	options(width=100)
	zz <- file("test2b.txt", open="wt")
	sink(zz)
	print(egarch.lik)
	cat("\nparameters:\n")
	print(egarch.coefs)
	cat("\nstd.errors:\n")
	print(egarch.se)
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

rugarch.test2c = function(cluster=NULL)
{
	tic = Sys.time()
	
	#cat("\nrugarch-->test2-3: Fixed and Starting Parameter Test (apARCH)\n")
	data(sp500ret)
	spars = list(mu = 1.9917e-04, ar1 = -1.7519e-02, omega = 6.5805e-05, 
			alpha1 = 6.0165e-02, beta1 = 9.3376e-01)
	
	spec = ugarchspec(
			variance.model = list(model = "apARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,0), include.mean = TRUE), 
			distribution.model = "nig", start.pars = spars)
	
	aparch.fit1 = ugarchfit(data = sp500ret, spec = spec, solver = "solnp")
	likelihood(aparch.fit1)
	coef(aparch.fit1)
	# Some Parameters Fixed (we obtain from aparch.fit1)
	fpars = as.list(coef(aparch.fit1)[1:6])
	
	spec = ugarchspec(
			variance.model = list(model = "apARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,0), include.mean = TRUE), 
			distribution.model = "nig", fixed.pars = fpars)
	
	aparch.fit2 = ugarchfit(data = sp500ret,spec = spec, solver = "solnp")
	
	# notice that the standard errors of the fixed parameters is NA (they are fixed!).
	# However, we can calculate those post estimation with the fixed.se argument:
	
	aparch.fit3 = ugarchfit(data = sp500ret,spec = spec, solver = "solnp", 
			fit.control = list(fixed.se = TRUE))
	
	# All Parameters Fixed (we obtain from aparch.fit1)
	fpars = as.list(coef(aparch.fit1))
	spec = ugarchspec(
			variance.model = list(model = "apARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,0), include.mean = TRUE), 
			distribution.model = "nig", fixed.pars = fpars)
	
	aparch.fit4 = ugarchfit(data = sp500ret,spec = spec, solver = "solnp", 
			fit.control = list(fixed.se = TRUE))
	
	# compare LLH, coefficients and std.errors:
	aparch.lik = cbind(
			likelihood(aparch.fit1), 
			likelihood(aparch.fit2), 
			likelihood(aparch.fit3), 
			likelihood(aparch.fit4))
	
	aparch.coefs = round(cbind(
					coef(aparch.fit1), 
					coef(aparch.fit2), 
					coef(aparch.fit3), 
					coef(aparch.fit4)),8)
	
	aparch.se = round(cbind(
					aparch.fit1@fit$matcoef[,2], 
					aparch.fit2@fit$matcoef[,2], 
					aparch.fit3@fit$matcoef[,2], 
					aparch.fit4@fit$matcoef[,2]),5)
	
	rownames(aparch.lik) = "likelihood"
	colnames(aparch.lik) = paste("test", 1:4, sep=".")
	colnames(aparch.coefs) = paste("test", 1:4, sep=".")
	colnames(aparch.se) = paste("test", 1:4, sep=".")
	
	options(width=100)
	zz <- file("test2c.txt", open="wt")
	sink(zz)
	print(aparch.lik)
	cat("\nparameters:\n")
	print(aparch.coefs)
	cat("\nstd.errors:\n")
	print(aparch.se)
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}


rugarch.test2d = function(cluster=NULL)
{
	tic = Sys.time()
	
	#cat("\nrugarch-->test1-4: Fixed and Starting Parameter Test (gjrGARCH)\n")
	
	data(sp500ret)
	spars = list(mu = 1.9917e-04, ar1 = -1.7519e-02, omega = 6.5805e-05, 
			alpha1 = 6.0165e-02, beta1 = 9.3376e-01)
	
	spec = ugarchspec(
			variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,0), include.mean = TRUE), 
			distribution.model = "nig", start.pars = spars)
	
	gjrgarch.fit1 = ugarchfit(data = sp500ret, spec = spec, solver = "solnp")
	
	# Some Parameters Fixed (we obtain from gjrgarch.fit1)
	fpars = as.list(coef(gjrgarch.fit1)[1:6])
	
	spec = ugarchspec(
			variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,0), include.mean = TRUE), 
			distribution.model = "nig",fixed.pars = fpars)
	
	gjrgarch.fit2 = ugarchfit(data = sp500ret,spec = spec, solver = "solnp")
	
	# Well, fixing all but 2 parameters leads to a slightly higher likelihood
	# as the skew and shape move slightly...
	
	# notice that the standard errors of the fixed parameters is NA (they are fixed!).
	# However, we can calculate those post estimation with the fixed.se argument:
	
	gjrgarch.fit3 = ugarchfit(data = sp500ret,spec = spec, solver = "solnp", 
			fit.control = list(fixed.se = TRUE))
	
	# All Parameters Fixed (we obtain from gjrgarch.fit1)
	fpars = as.list(coef(gjrgarch.fit1))
	spec = ugarchspec(
			variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,0), include.mean = TRUE), 
			distribution.model = "nig", fixed.pars = fpars)
	
	gjrgarch.fit4 = ugarchfit(data = sp500ret,spec = spec, solver = "solnp", 
			fit.control = list(fixed.se = TRUE))
	
	# compare LLH, coefficients and std.errors:
	gjrgarch.lik = cbind(
			likelihood(gjrgarch.fit1), 
			likelihood(gjrgarch.fit2), 
			likelihood(gjrgarch.fit3), 
			likelihood(gjrgarch.fit4))
	
	gjrgarch.coefs = round(cbind(
					coef(gjrgarch.fit1), 
					coef(gjrgarch.fit2), 
					coef(gjrgarch.fit3), 
					coef(gjrgarch.fit4)),5)
	
	gjrgarch.se = round(cbind(
					gjrgarch.fit1@fit$matcoef[,2], 
					gjrgarch.fit2@fit$matcoef[,2], 
					gjrgarch.fit3@fit$matcoef[,2], 
					gjrgarch.fit4@fit$matcoef[,2]),5)
	
	rownames(gjrgarch.lik) = "likelihood"
	colnames(gjrgarch.lik) = paste("test", 1:4, sep=".")
	colnames(gjrgarch.coefs) = paste("test", 1:4, sep=".")
	colnames(gjrgarch.se) = paste("test", 1:4, sep=".")
	
	options(width=100)
	zz <- file("test2d.txt", open="wt")
	sink(zz)
	print(gjrgarch.lik)
	cat("\nparameters:\n")
	print(gjrgarch.coefs)
	cat("\nstd.errors:\n")
	print(gjrgarch.se)
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

rugarch.test2e = function(cluster=NULL)
{
	tic = Sys.time()
	
	#cat("\nrugarch-->test1-5: Fixed and Starting Parameter Test (sGARCH)\n")
	
	data(sp500ret)
	spars = list(mu = 1.9917e-04, ar1 = -1.7519e-02, omega = 6.5805e-05, 
			alpha1 = 6.0165e-02, beta1 = 9.3376e-01)
	
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,0), include.mean = TRUE), 
			distribution.model = "nig", start.pars = spars)
	
	sgarch.fit1 = ugarchfit(data = sp500ret, spec = spec, solver = "solnp")
	
	# Some Parameters Fixed (we obtain from sgarch.fit1)
	fpars = as.list(coef(sgarch.fit1)[1:3])
	
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,0), include.mean = TRUE), 
			distribution.model = "nig", fixed.pars = fpars)
	
	sgarch.fit2 = ugarchfit(data = sp500ret,spec = spec, solver = "solnp")
	
	# notice that the standard errors of the fixed parameters is NA (they are fixed!).
	# However, we can calculate those post estimation with the fixed.se argument:
	
	sgarch.fit3 = ugarchfit(data = sp500ret,spec = spec, solver = "solnp", 
			fit.control = list(fixed.se = TRUE))
	
	# All Parameters Fixed (we obtain from sgarch.fit1)
	fpars = as.list(coef(sgarch.fit1))
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,0), include.mean = TRUE), 
			distribution.model = "nig", fixed.pars = fpars)
	
	sgarch.fit4 = ugarchfit(data = sp500ret,spec = spec, solver = "solnp", 
			fit.control = list(fixed.se = TRUE))
	
	# compare LLH, coefficients and std.errors:
	sgarch.lik = cbind(
			likelihood(sgarch.fit1), 
			likelihood(sgarch.fit2), 
			likelihood(sgarch.fit3), 
			likelihood(sgarch.fit4))
	
	sgarch.coefs = round(cbind(
					coef(sgarch.fit1), 
					coef(sgarch.fit2), 
					coef(sgarch.fit3), 
					coef(sgarch.fit4)),5)
	
	sgarch.se = round(cbind(
					sgarch.fit1@fit$matcoef[,2], 
					sgarch.fit2@fit$matcoef[,2], 
					sgarch.fit3@fit$matcoef[,2], 
					sgarch.fit4@fit$matcoef[,2]),5)
	
	rownames(sgarch.lik) = "likelihood"
	colnames(sgarch.lik) = paste("test", 1:4, sep=".")
	colnames(sgarch.coefs) = paste("test", 1:4, sep=".")
	colnames(sgarch.se) = paste("test", 1:4, sep=".")
	
	options(width=100)
	zz <- file("test2e.txt", open="wt")
	sink(zz)
	print(sgarch.lik)
	cat("\nparameters:\n")
	print(sgarch.coefs)
	cat("\nstd.errors:\n")
	print(sgarch.se)
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}