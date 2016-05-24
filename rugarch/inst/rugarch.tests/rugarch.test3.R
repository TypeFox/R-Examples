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
# Fitting (estimation) tests
#################################################################################

rugarch.test3a = function(cluster=NULL)
{
	tic = Sys.time()
	
	#cat("\nrugarch-->test3-1: Fit Test (fGARCH)\n")
	
	data(dji30ret)
	# 1. GARCH
	spec = ugarchspec(
			variance.model = list(model = "fGARCH", garchOrder = c(1,1), 
					submodel = "GARCH"), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std")
	
	fgarch.s = ugarchfit(data = dji30ret[,"AA",drop=FALSE], spec = spec, solver = "solnp")
	
	# 2. TGARCH
	spec = ugarchspec(
			variance.model = list(model = "fGARCH", garchOrder = c(1,1), 
					submodel = "TGARCH"), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std")
	
	fgarch.t = ugarchfit(data = dji30ret[,"AA",drop=FALSE], spec = spec, solver = "solnp")
	
	# 3. NGARCH
	spec = ugarchspec(
			variance.model = list(model = "fGARCH", garchOrder = c(1,1),  
					submodel = "NGARCH"), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std")
	
	fgarch.n = ugarchfit(data = dji30ret[,"AA",drop=FALSE], spec = spec, solver = "solnp")
	
	
	# 4. NAGARCH
	spec = ugarchspec(
			variance.model = list(model = "fGARCH", garchOrder = c(1,1),  
					submodel = "NAGARCH"), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std")
	
	fgarch.na = ugarchfit(data = dji30ret[,"AA",drop=FALSE], spec = spec, solver = "solnp")
	
	# 5. AVGARCH
	spec = ugarchspec(
			variance.model = list(model = "fGARCH", garchOrder = c(1,1), 
					submodel = "AVGARCH"), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std")
	
	fgarch.av = ugarchfit(data = dji30ret[,"AA",drop=FALSE], spec = spec, solver = "solnp")
	
	# 6. APARCH
	spec = ugarchspec(
			variance.model = list(model = "fGARCH", garchOrder = c(1,1),  
					submodel = "APARCH"), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std")
	
	fgarch.ap = ugarchfit(data = dji30ret[,"AA",drop=FALSE], spec = spec, solver = "solnp")
	
	
	# 7. GJR
	spec = ugarchspec(
			variance.model = list(model = "fGARCH", garchOrder = c(1,1),  
					submodel = "GJRGARCH"), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std")
	
	fgarch.gjr = ugarchfit(data = dji30ret[,"AA",drop=FALSE], spec = spec, solver = "solnp")
	
	
	# 8. ALLGARCH (needs some starting values, we use aparch fit)
	spec = ugarchspec(
			variance.model = list(model = "fGARCH", garchOrder = c(1,1), 
					submodel = "ALLGARCH"), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std", start.pars = as.list(coef(fgarch.ap)))
	
	fgarch.al = ugarchfit(data = dji30ret[,"AA",drop=FALSE], spec = spec, solver = "solnp")
	
	# compare coefficients and likelihoods
	fgarch.lik = c(
			likelihood(fgarch.s), 
			likelihood(fgarch.t), 
			likelihood(fgarch.n), 
			likelihood(fgarch.na), 
			likelihood(fgarch.av), 
			likelihood(fgarch.ap),
			likelihood(fgarch.gjr),
			likelihood(fgarch.al))
	
	fgarch.bic = c(
			infocriteria(fgarch.s)[2], 
			infocriteria(fgarch.t)[2], 
			infocriteria(fgarch.n)[2], 
			infocriteria(fgarch.na)[2], 
			infocriteria(fgarch.av)[2], 
			infocriteria(fgarch.ap)[2],
			infocriteria(fgarch.gjr)[2],
			infocriteria(fgarch.al)[2])
	
	fgarch.cnames = c("mu", "ar1", "ma1", "omega", "alpha1", "eta11", "eta21", 
			"beta1", "lambda", "shape")
	
	fgarch.models = matrix(NA, ncol = 8, nrow = 10)
	fgarch.models[,1] = coef(fgarch.s)[fgarch.cnames]
	fgarch.models[,2] = coef(fgarch.t)[fgarch.cnames]
	fgarch.models[,3] = coef(fgarch.n)[fgarch.cnames]
	fgarch.models[,4] = coef(fgarch.na)[fgarch.cnames]
	fgarch.models[,5] = coef(fgarch.av)[fgarch.cnames]
	fgarch.models[,6] = coef(fgarch.ap)[fgarch.cnames]
	fgarch.models[,7] = coef(fgarch.gjr)[fgarch.cnames]
	fgarch.models[,8] = coef(fgarch.al)[fgarch.cnames]
	fgarch.models = rbind(fgarch.models, fgarch.lik)
	fgarch.models = rbind(fgarch.models, fgarch.bic)
	
	colnames(fgarch.models) = c("GARCH", "TGARCH", "NGARCH", "NAGARCH", 
			"AVGARCH", "APARCH", "GJRGARCH", "ALLGARCH")
	
	rownames(fgarch.models) = c(fgarch.cnames, "LLH", "BIC")
	fgarch.models = fgarch.models[,order(fgarch.models[12,], decreasing = TRUE)]
	options(width=175)
	zz <- file("test3a.txt", open="wt")
	sink(zz)
	print(round(fgarch.models,5), digits = 6)
	sink(type="message")
	sink()
	close(zz)
	
	# do some newsimpact plots
	fgarch.al.ni = newsimpact(fgarch.al)
	zall = fgarch.al.ni$zx
	nilist = vector(mode = "list", length = 8)
	nilist[[1]] = newsimpact(fgarch.s, z = zall)
	nilist[[2]] = newsimpact(fgarch.t, z = zall)
	nilist[[3]] = newsimpact(fgarch.n, z = zall)
	nilist[[4]] = newsimpact(fgarch.na, z = zall)
	nilist[[5]] = newsimpact(fgarch.av, z = zall)
	nilist[[6]] = newsimpact(fgarch.ap, z = zall)
	nilist[[7]] = newsimpact(fgarch.gjr, z = zall)
	
	postscript("test3a.eps", bg = "white", width = 800, height = 800)
	colx = rugarch:::.distinctcolors11()
	plot(zall, fgarch.al.ni$zy, ylab = fgarch.al.ni$yexpr, xlab = fgarch.al.ni$xexpr, type = "l", lwd=2,
			col = 1, main="News Impact Curves\nfGARCH Family")
	for(i in 1:7) lines(zall, nilist[[i]]$zy, col = colx[i+1], lty=i)
	legend("bottomleft", legend = c("ALLGARCH", "GARCH", "TGARCH", "NGARCH", "NAGARCH", 
			"AVGARCH", "APARCH", "GJRGARCH"), col = c(1, colx[2:8]), lty=c(1,1:7))
	dev.off()
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}


# ---------------------------------------------------------------------------------
# Models with variation in mean model choices and external regressors
# in mean and variance equation (fully parametrized model - we follow an
# incremental building strategy using the coefficients from previous fits
# as starting parameters into more complicated models)

rugarch.test3b = function(cluster=NULL)
{
	tic = Sys.time()
	# As of versions 1.01-3. WeekDayDummy no longer exported
	#cat("\nrugarch-->test3-2: Fit Test (sGARCH)\n")
	data(dji30ret)
	dates = rownames(dji30ret[,"AA", drop = FALSE])
	monday = rugarch:::.WeekDayDummy(dates, date.format = "%Y-%m-%d", weekday = "Monday")
	# convert to matrix which is what the specification expects
	monday = matrix(monday, ncol = 1)
	# create a dummy day-of-week variable for the variance regression (Friday)
	friday = rugarch:::.WeekDayDummy(dates, date.format = "%Y-%m-%d", weekday = "Friday")
	# convert to matrix which is what the specification expects
	friday = matrix(friday, ncol = 1)
	# sGARCH(1,1)
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(0,0), include.mean = FALSE), 
			distribution.model = "std")
	
	sgarch.fit1 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], spec = spec, 
			solver = "solnp", fit.control = list(scale=1))
	
	# sGARCH(1,1) + MU
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(0,0), include.mean = TRUE), 
			distribution.model = "std")
	
	sgarch.fit2 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], spec = spec, 
			solver = "solnp")
	
	# sGARCH(1,1) + MU + ARMA(1,1)
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std")
	
	sgarch.fit3 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], spec = spec, 
			solver = "solnp")
	
	# sGARCH(1,1) + XREG(V) + ARMA(1,1) + XREG(M)
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1), 
					external.regressors = friday), 
			mean.model = list(armaOrder = c(1,1), include.mean = FALSE, 
					external.regressors = monday), distribution.model = "std", 
			start.pars = as.list(coef(sgarch.fit3)))
	
	sgarch.fit4 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], spec = spec, 
			solver = "solnp")
	
	# sGARCH(1,1) + XREG(V) + MU + ARMA(1,1) + XREG(M)
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1), 
					external.regressors = friday), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE, 
					external.regressors = monday), 
			distribution.model = "std", 
			start.pars = as.list(coef(sgarch.fit4 )))
	
	sgarch.fit5 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], spec = spec, 
			solver = "solnp")
	
	# sGARCH(1,1) + XREG(V) + ARMA(1,1) + XREG(M) + INMEAN(1)
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1), 
					external.regressors = friday), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE, 
					archm = TRUE, archpow = 1, external.regressors = monday), 
			distribution.model = "std", 
			start.pars = as.list(coef(sgarch.fit5)))
	
	sgarch.fit6 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], spec = spec, 
			solver = "solnp")
	
	# compare coefficients and likelihoods
	sgarch.lik = c(
			likelihood(sgarch.fit1), 
			likelihood(sgarch.fit2) ,
			likelihood(sgarch.fit3), 
			likelihood(sgarch.fit4), 
			likelihood(sgarch.fit5), 
			likelihood(sgarch.fit6))
	
	sgarch.bic = c(
			infocriteria(sgarch.fit1)[2], 
			infocriteria(sgarch.fit2)[2],
			infocriteria(sgarch.fit3)[2], 
			infocriteria(sgarch.fit4)[2], 
			infocriteria(sgarch.fit5)[2], 
			infocriteria(sgarch.fit6)[2])
	
	sgarch.cnames = c("mu", "ar1", "ma1", "archm", "mxreg1", "omega", "alpha1", 
			"beta1", "vxreg1", "shape")
	
	sgarch.models = matrix(NA, ncol = 6, nrow = 10)
	sgarch.models[,1] = coef(sgarch.fit1)[sgarch.cnames]
	sgarch.models[,2] = coef(sgarch.fit2)[sgarch.cnames]
	sgarch.models[,3] = coef(sgarch.fit3)[sgarch.cnames]
	sgarch.models[,4] = coef(sgarch.fit4)[sgarch.cnames]
	sgarch.models[,5] = coef(sgarch.fit5)[sgarch.cnames]
	sgarch.models[,6] = coef(sgarch.fit6)[sgarch.cnames]
	sgarch.models = rbind(sgarch.models, sgarch.lik)
	sgarch.models = rbind(sgarch.models, sgarch.bic)
	
	colnames(sgarch.models) = c("model1", "model2", "model3", "model4", "model5", 
			"model6")
	rownames(sgarch.models) = c(sgarch.cnames, "LLH", "BIC")
	options(width = 150)
	zz <- file("test3b.txt", open="wt")
	sink(zz)
	print(round(sgarch.models[,order(sgarch.models[12,], decreasing = TRUE)],5), digits = 6)
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}


#test2-2.txt


# ---------------------------------------------------------------------------------
# apARCH
# ---------------------------------------------------------------------------------
rugarch.test3c = function(cluster=NULL)
{
	tic = Sys.time()
	
	#cat("\nrugarch-->test3-3: Fit Test (apARCH)\n")
	
	data(dji30ret)
	dates = rownames(dji30ret[,"AA", drop = FALSE])
	monday = rugarch:::.WeekDayDummy(dates, date.format = "%Y-%m-%d", weekday = "Monday")
	# convert to matrix which is what the specification expects
	monday = matrix(monday, ncol = 1)
	# create a dummy day-of-week variable for the variance regression (Friday)
	friday = rugarch:::.WeekDayDummy(dates, date.format = "%Y-%m-%d", weekday = "Friday")
	# convert to matrix which is what the specification expects
	friday = matrix(friday, ncol = 1)
	# apARCH(1,1)
	spec = ugarchspec(
			variance.model = list(model = "apARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(0,0), include.mean = FALSE), 
			distribution.model = "std")
	
	aparch.fit1 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], spec = spec, 
			solver = "solnp")
	
	# apARCH(1,1) + MU
	spec = ugarchspec(
			variance.model = list(model = "apARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(0,0), include.mean = TRUE), 
			distribution.model = "std")
	
	aparch.fit2 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], spec = spec, 
			solver = "solnp")
	
	# apARCH(1,1) + MU + ARMA(1,1)
	spec = ugarchspec(
			variance.model = list(model = "apARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std")
	
	aparch.fit3 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], spec = spec, 
			solver = "solnp")
	
	# apARCH(1,1) + XREG(V) + ARMA(1,1) + XREG(M)
	spec = ugarchspec(
			variance.model = list(model = "apARCH", garchOrder = c(1,1), 
					external.regressors = friday), 
			mean.model = list(armaOrder = c(1,1), include.mean = FALSE, 
					external.regressors = monday), 
			distribution.model = "std", 
			start.pars = as.list(coef(aparch.fit3)))
	
	aparch.fit4 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], spec = spec, 
			solver = "solnp")
	
	# apARCH(1,1) + XREG(V) + MU + ARMA(1,1) + XREG(M)
	spec = ugarchspec(
			variance.model = list(model = "apARCH", garchOrder = c(1,1), 
					external.regressors = friday), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE, 
					external.regressors = monday), 
			distribution.model = "std", 
			start.pars = as.list(coef(aparch.fit4)))
	
	aparch.fit5 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], spec = spec, 
			solver = "solnp")
	
	# apARCH(1,1) + XREG(V) + ARMA(1,1) + XREG(M) + INMEAN(1)
	spec = ugarchspec(
			variance.model = list(model = "apARCH", garchOrder = c(1,1), 
					external.regressors = friday), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE, 
					archm = TRUE, archpow = 1, external.regressors = monday), 
			distribution.model = "std", 
			start.pars = as.list(coef(aparch.fit5)))
	
	aparch.fit6 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], spec = spec, 
			solver = "solnp")
	
	# compare coefficients and likelihoods
	aparch.lik = c(
			likelihood(aparch.fit1), 
			likelihood(aparch.fit2),
			likelihood(aparch.fit3), 
			likelihood(aparch.fit4), 
			likelihood(aparch.fit5), 
			likelihood(aparch.fit6))
	
	aparch.bic = c(
			infocriteria(aparch.fit1)[2], 
			infocriteria(aparch.fit2)[2],
			infocriteria(aparch.fit3)[2], 
			infocriteria(aparch.fit4)[2], 
			infocriteria(aparch.fit5)[2], 
			infocriteria(aparch.fit6)[2])
	
	aparch.cnames = c("mu", "ar1", "ma1", "inmean", "mxreg1", "omega", 
			"alpha1", "beta1", "gamma1", "delta", "vxreg1", "shape")
	
	aparch.models = matrix(NA, ncol = 6, nrow = 12)
	aparch.models[,1] = coef(aparch.fit1)[aparch.cnames]
	aparch.models[,2] = coef(aparch.fit2)[aparch.cnames]
	aparch.models[,3] = coef(aparch.fit3)[aparch.cnames]
	aparch.models[,4] = coef(aparch.fit4)[aparch.cnames]
	aparch.models[,5] = coef(aparch.fit5)[aparch.cnames]
	aparch.models[,6] = coef(aparch.fit6)[aparch.cnames]
	aparch.models = rbind(aparch.models, aparch.lik)
	aparch.models = rbind(aparch.models, aparch.bic)
	
	colnames(aparch.models) = c("model1", "model2", "model3", "model4", 
			"model5", "model6")
	rownames(aparch.models) = c(aparch.cnames, "LLH", "BIC")
	
	options(width=150)
	zz <- file("test3c.txt", open="wt")
	sink(zz)
	print(round(aparch.models[,order(aparch.models[14,], decreasing = TRUE)],5), digits = 6)
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

rugarch.test3d = function(cluster=NULL)
{
	tic = Sys.time()
	
	#cat("\nrugarch-->test3-4: Filter Test (sGARCH, iGARCH & apARCH)\n")
	
	data(dji30ret)
	# sGARCH Model
	# ---------------------------------------------------------------------------------
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std")
	
	sgarch.fit = ugarchfit(data = dji30ret[,"AA",drop=FALSE], spec = spec, 
			solver = "solnp")
	
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std", fixed.pars = as.list(coef(sgarch.fit)))
	sgarch.filter = ugarchfilter(data = dji30ret[,"AA",drop=FALSE], spec = spec)	
	
	# iGARCH Model
	# ---------------------------------------------------------------------------------
	spec = ugarchspec(
			variance.model = list(model = "iGARCH", garchOrder = c(2,2)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std")
	
	igarch.fit = ugarchfit(data = dji30ret[,"AA",drop=FALSE], spec = spec, 
			solver = "solnp")
	
	spec = ugarchspec(
			variance.model = list(model = "iGARCH", garchOrder = c(2,2)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std", fixed.pars = as.list(coef(igarch.fit)))
	igarch.filter = ugarchfilter(data = dji30ret[,"AA",drop=FALSE], spec = spec)
		
	# apARCH Model
	# ---------------------------------------------------------------------------------
	spec = ugarchspec(
			variance.model = list(model = "apARCH", garchOrder = c(1,2)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std")
	
	aparch.fit = ugarchfit(data = dji30ret[,"AA",drop=FALSE], spec = spec, 
			solver = "solnp", fit.control = list(scale = 1))
	
	spec = ugarchspec(
			variance.model = list(model = "apARCH", garchOrder = c(1,2)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std", fixed.pars = as.list(coef(aparch.fit)))
	aparch.filter = ugarchfilter(data = dji30ret[,"AA",drop=FALSE], spec = spec)
	
	x1 = cbind(head(sigma(sgarch.fit),10) , head(sigma(sgarch.filter),10))
	x2 = cbind(head(sigma(igarch.fit),10) , head(sigma(igarch.filter),10))
	x3 = cbind(head(sigma(aparch.fit),10) , head(sigma(aparch.filter),10))
	colnames(x1) = c("fit", "filter")
	colnames(x2) = c("fit", "filter")
	colnames(x3) = c("fit", "filter")
	zz <- file("test3d.txt", open="wt")
	sink(zz)
	cat("\nsGARCH sigma\n")
	print(x1, digits=8)
	cat("\niGARCH sigma\n")
	print(x2, digits=8)
	cat("\napARCH sigma\n")
	print(x3, digits=8)
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}