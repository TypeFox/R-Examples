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

# The realGARCH model
rugarch.test13a = function(cluster=NULL)
{
	# replicate results of paper for the SPY dataset
	data(spyreal)
	spec = ugarchspec(mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
			variance.model = list(model = "realGARCH", garchOrder = c(2,1), 
					submodel = NULL, external.regressors = NULL, variance.targeting = FALSE),
			distribution="norm")
	
	fit = ugarchfit(spec, spyreal[,1]*100, out.sample = 0, solver = "hybrid",
			realizedVol = spyreal[,2]*100)
	cf = coef(fit)
	se = fit@fit$matcoef[,2]
	names(se) = names(cf)
	# paper:
	# They divide the period into in-(T=1492) and out- of sample, even though they are estimating
	# the model using the entire dataset (T=1659)!
	benchmark.LL = c("logL" = -2388.75, "pLogL" = -1710.34)
	benchmark.pars = c("logh0" = -0.289213754, "omega" = 0.041246175, "alpha1" =0.450676773, "alpha2" = -0.176042911,
			"beta1" = 0.701210369, "xi" = -0.179994885, "delta" = 1.037491982, "eta1" = -0.067809509,
			"eta2" = 0.07015778, "lambda" = sqrt(0.145369801))
	benchmark.se  = c("logh0" = 0.278787604, "omega" = 0.015515016, "alpha1" = 0.038027571, "alpha2" = 0.06335417,
			"beta1" = 0.055974942, "xi" = 0.043611245, "delta" = 0.057726499, "eta1" = 0.010309694,
			"eta2" = 0.006620828, "lambda" = 0.00594321)
	# rugarch does not estimate h0, instead uses either mean(residuals^2), else a choice of variance targeting
	# with options
	rugarch.LL = c("logL" =	sum(-fit@fit$log.likelihoods[1:1492]), "pLogL" = sum(-fit@fit$partial.log.likelihoods[1:1492]))
	rugarch.pars = c("logh0" = NA, "omega" = cf["omega"], "alpha1" = cf["alpha1"], "alpha2" = cf["alpha2"],
			"beta1" = cf["beta1"], "xi" = cf["xi"], "delta" = cf["delta"], "eta1" = cf["eta11"],
			"eta2" = cf["eta21"], "lambda" = cf["lambda"])
	rugarch.se = c("logh0" = NA, "omega" = se["omega"], "alpha1" = se["alpha1"], "alpha2" = se["alpha2"],
			"beta1" = se["beta1"], "xi" = se["xi"], "delta" = se["delta"], "eta1" = se["eta11"],
			"eta2" = se["eta21"], "lambda" = se["lambda"])
	names(rugarch.pars) = names(rugarch.se) = c("logh0","omega","alpha1","alpha2","beta1","xi","delta","eta1","eta2","lambda")
	parsdf = cbind(benchmark.pars, rugarch.pars, benchmark.pars-rugarch.pars)
	sedf = cbind(benchmark.se, rugarch.se, benchmark.se-rugarch.se)
	
	LRE.vars = -log(abs(rugarch.pars-benchmark.pars)/abs(benchmark.pars), base = 10)
	LRE.se = -log(abs(rugarch.se-benchmark.se)/abs(benchmark.se), base = 10)
	test = cbind(LRE.vars, LRE.se)
	
	options(width=120)
	zz <- file("test13a-1.txt", open="wt")
	sink(zz)
	cat("\nRealized GARCH model benchmark:\n")
	cat("\nparameters:\n")
	tmp = t(cbind(rugarch=rugarch.pars, benchmark=benchmark.pars))
	print(round(tmp, 4))
	cat("\nstandard errors:\n")
	tmp = t(cbind(rugarch=rugarch.se, benchmark=benchmark.se))
	print(round(tmp, 4))
	cat("\nLog Relative Error Test:\n")
	print(round(t(test), 4))
	sink(type="message")
	sink()
	close(zz)
	
}
	

rugarch.test13b = function(cluster=NULL)
{
	# fit/filter and forecast
	require(xts)
	data(spyreal)
	spec = ugarchspec(mean.model=list(armaOrder=c(1,1), include.mean=TRUE),
			variance.model = list(model = "realGARCH", garchOrder = c(1,1), 
					submodel = NULL, external.regressors = NULL, variance.targeting = FALSE),
			distribution="norm")
	
	fit = ugarchfit(spec, spyreal[,1]*100, out.sample = 25, solver = "hybrid",
			realizedVol = spyreal[,2]*100)
	specf = spec
	setfixed(specf)<-as.list(coef(fit))
	filt = ugarchfilter(specf, data = spyreal[,1]*100, n.old = nrow(spyreal)-25,
			realizedVol = spyreal[,2]*100)
	
	options(width=120)
	zz <- file("test13b-1.txt", open="wt")
	sink(zz)
	print(all.equal(head(sigma(fit)), head(sigma(filt))))
	print(all.equal(head(fitted(fit)), head(fitted(filt))))
	sink(type="message")
	sink()
	close(zz)
	
	# forecast from fit
	forc1 = ugarchforecast(fit, n.ahead=1, n.roll = 25)
	# forecast from spec
	forc2 = ugarchforecast(specf, n.ahead=1, n.roll = 25, data = spyreal[,1]*100, out.sample = 25, realizedVol = spyreal[,2]*100)
	
	filts = tail(sigma(filt), 25)
	colnames(filts)="filter"
	forcs1 = xts(sigma(forc1)[1,], move(as.Date(names(sigma(forc1)[1,])), by=1))
	forcs2 = xts(sigma(forc2)[1,], move(as.Date(names(sigma(forc2)[1,])), by=1))
	colnames(forcs1)  = "fit2forecast"
	colnames(forcs2)  = "spec2forecast"
	ftest = cbind(filts, forcs1, forcs2)
	
	# last forecast is completely out of sample, so not available from the filter method (which filters given T-1)
	options(width=120)
	zz <- file("test13b-2.txt", open="wt")
	sink(zz)
	print(round(ftest,5))
	sink(type="message")
	sink()
	close(zz)
	
	set.seed(55)
	forc3 = ugarchforecast(fit, n.ahead=400, n.sim=5000)
	sqrt(uncvariance(fit))
	
	postscript("test13b-1.eps")
	plot(sigma(forc3), type="l", main="realGARCH long-run forecast")
	abline(h=sqrt(uncvariance(fit)), col=2)
	legend("topright", c("long-run forecast", "unconditional value"), col=1:2, lty=c(1,1), bty="n")
	dev.off()
	
}
	
rugarch.test13c = function(cluster=NULL)
{
	# simulation
	data(spyreal)
	spec = ugarchspec(mean.model=list(armaOrder=c(1,1), include.mean=TRUE),
			variance.model = list(model = "realGARCH", garchOrder = c(1,1), 
					submodel = NULL, external.regressors = NULL, variance.targeting = FALSE),
			distribution="norm")
	
	fit = ugarchfit(spec, spyreal[,1]*100, out.sample = 25, solver = "hybrid",
			realizedVol = spyreal[,2]*100)
	specf = spec
	setfixed(specf)<-as.list(coef(fit))
	
	sim = ugarchsim(fit, n.sim = 10000, m.sim=400, n.start = 500, startMethod = "sample")	
	msim = apply(log(sim@simulation$sigmaSim^2), 2, "mean")
	exp(mean(msim))
	uncvariance(fit)
	
	postscript("test13c-1.eps")
	hist(msim, main="simulated average log-sigma \n(400 paths x 10000 horizon)")
	abline(v=log(uncvariance(fit)), col=2)
	legend("topright", c(paste("mean (density): ", round(mean(msim),3), sep=""), 
					paste("mean (analytical): ", round(log(uncvariance(fit)), 3))), 
			col=1:2, lty=c(1,1), bty="n")
	dev.off()
	
	# pre-values
	T = nrow(spyreal)-25
	sim1 = ugarchsim(fit, n.sim = 1000, m.sim=1, n.start = 0, startMethod = "sample", rseed=12)
	sim2 = ugarchsim(fit, n.sim = 1000, m.sim=1, n.start = 0, startMethod = "sample", rseed=12,
			prereturns = as.numeric(tail(spyreal[1:T,1]*100,1)), 
			presigma = as.numeric(tail(sigma(fit),1)), 
			preresiduals = as.numeric(tail(residuals(fit),1)),
			prerealized = as.numeric(tail(spyreal[1:T,2]*100, 1)))
	sim3 = ugarchpath(specf, n.sim = 1000, m.sim=1, n.start = 0, rseed=12,
			prereturns = as.numeric(tail(spyreal[1:T,1]*100,1)), 
			presigma = as.numeric(tail(sigma(fit),1)), 
			preresiduals = as.numeric(tail(residuals(fit),1)),
			prerealized = as.numeric(tail(spyreal[1:T,2]*100, 1)))
	
	s1 = sigma(sim1)
	s2 = sigma(sim2)
	s3 = sigma(sim3)
	
	f1 = fitted(sim1)
	f2 = fitted(sim2)
	f3 = fitted(sim3)
	
	options(width=120)
	zz <- file("test13c-1.txt", open="wt")
	sink(zz)
	print(all.equal(s1,s2,s3))
	print(all.equal(f1,f2,f3))
	sink(type="message")
	sink()
	close(zz)
}

rugarch.test13d = function(cluster=NULL)
{
	# rolling backtest
	data(spyreal)
	spec1 = ugarchspec(mean.model=list(armaOrder=c(1,1), include.mean=TRUE),
			variance.model = list(model = "realGARCH", garchOrder = c(1,1), 
					submodel = NULL, external.regressors = NULL, variance.targeting = FALSE),
			distribution="norm")
	spec2 = ugarchspec(mean.model=list(armaOrder=c(1,1), include.mean=TRUE),
			variance.model = list(model = "eGARCH", garchOrder = c(1,1), 
					submodel = NULL, external.regressors = NULL, variance.targeting = FALSE),
			distribution="norm")
	roll1 = ugarchroll(spec1, spyreal[,1]*100, forecast.length = 500, solver="hybrid",
			refit.every=25, refit.window="recursive", realizedVol = spyreal[,2]*100,
			cluster = cluster)
	
	roll2 = ugarchroll(spec2, spyreal[,1]*100, forecast.length = 500, 
			refit.every=25, refit.window="recursive", cluster = cluster)
	
	# compare likelihoods:
	LL = data.frame(realGARCH = roll1@model$partial.loglik, eGARCH = roll2@model$loglik)
	options(width=120)
	zz <- file("test13d-1.txt", open="wt")
	sink(zz)
	cat("\nLog-Likelihood Comparison:\n\n")
	print(LL)
	cat("\nVaR Reports:\n\n")
	report(roll1)
	cat("\n\n")
	report(roll2)
	sink(type="message")
	sink()
	close(zz)
	
	postscript("test13d-1.eps")
	plot(as.xts(as.data.frame(roll1)[,"Sigma",drop=FALSE]), main = "realGARCH vs eGARCH\n(out-of-sample volatility forecast)",
			auto.grid = FALSE, minor.ticks = FALSE)
	lines(as.xts(as.data.frame(roll2)[,"Sigma",drop=FALSE]), col=2)
	legend("topleft", c("realGARCH","eGARCH"), col=1:2, lty=c(1,1), bty="n")
	grid()
	dev.off()
	
}

rugarch.test13e = function(cluster=NULL)
{
	# forecast density of sigma and realvol
}