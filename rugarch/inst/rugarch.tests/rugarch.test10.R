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
# Distribution Checks
rugarch.test10a = function(dist = "ghst", skew = -0.5, shape = 8.1, cluster=NULL)
{
	tic = Sys.time()
	#cat("\nrugarch-->test10-1: Fixed and Starting Parameter Test (fGARCH/ALLGARCH)\n")
	data(sp500ret)
	# No Parameters Fixed, Starting Parameters
	spars = list(mu = 1.9917e-04, ar1 = -1.7519e-02, omega = 6.5805e-05, 
			alpha1 = 6.0165e-02, beta1 = 9.3376e-01, lambda = 1.1702e+00,  
			eta21 = 4.2051e-02, eta11 = 7.9775e-01, skew = skew, 
			shape = shape)
	
	spec = ugarchspec(
			variance.model = list(model = "fGARCH", garchOrder = c(1,1), 
					submodel = "ALLGARCH"), 
			mean.model = list(armaOrder = c(1,0), include.mean = TRUE, 
					archm = FALSE, archpow = 2), 
			distribution.model = dist)
	setstart(spec)<-spars
	
	fgarch.fit1 = ugarchfit(data = sp500ret, spec = spec, solver = "solnp", 
			solver.control = list(trace=1, tol=1e-10, delta=1e-9,rho=1))
	
	# Some Parameters Fixed (we obtain from fgarch.fit1)
	fpars = as.list(coef(fgarch.fit1)[1:6])
	
	spec = ugarchspec(
			variance.model = list(model = "fGARCH", garchOrder = c(1,1), 
					submodel = "ALLGARCH"), 
			mean.model = list(armaOrder = c(1,0), include.mean = TRUE, 
					archm = FALSE, archpow = 2), distribution.model = dist, 
			fixed.pars = fpars, start.pars = list(skew = -0.8))
	# alternative use setfixed(spec)<-fpars
	fgarch.fit2 = ugarchfit(data = sp500ret, spec = spec, solver = "solnp", 
			fit.control=list(scale=0))	
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
					archm = FALSE, archpow = 2), distribution.model = dist, 
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
	zz <- file(paste("test10a-",dist,".txt",sep=""), open="wt")
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


rugarch.test10b = function(dist = "ghst", skew = -0.5, shape = 8.1, cluster=NULL)
{
	tic = Sys.time()
	
	data(dji30ret)
	rseed = rugarch.seeds(300, "test10b")
	spec = ugarchspec(
			variance.model = list(model = "sGARCH"), 
			distribution.model = dist)
	# set lower bounds on the ghst distribution for existence of kurtosis
	if(dist=="ghst"){
		setbounds(spec)<-list(shape=c(8.1, 25))
	# best to set a starting parameter above the lower bound else the solver
	# might fail:
		setstart(spec)<-list(shape = 8.1)
	}
	fit = ugarchfit(spec, data = dji30ret[,"XOM"], fit.control = list(scale = 1))	
	distx = ugarchdistribution(fit, n.sim = 2000, n.start = 50, m.sim = 300, 
			solver = "solnp", fit.control = list(scale = 1), rseed  = rseed, 
			cluster = cluster)
	
	postscript(paste("test10b1-",dist,".eps",sep=""), width = 12, height = 8)
	plot(distx, which = 1)
	dev.off()
	
	postscript(paste("test10bs-",dist,".eps",sep=""), width = 12, height = 8)
	plot(distx, which = 2)
	dev.off()
	
	postscript(paste("test10b3-",dist,".eps",sep=""), width = 12, height = 8)
	plot(distx, which = 3)
	dev.off()
	
	z1 <- file(paste("test10b1-",dist,".txt",sep=""), open="wt")
	sink(z1)
	print(as.data.frame(distx, which = "coef"))
	sink(type="message")
	sink()
	close(z1)
	
	z2 <- file(paste("test10b2-",dist,".txt",sep=""), open="wt")
	sink(z2)
	print(as.data.frame(distx, which = "rmse"))
	sink(type="message")
	sink()
	close(z2)
	
	
	z3 <- file(paste("test10b3-",dist,".txt",sep=""), open="wt")
	sink(z3)
	print(as.data.frame(distx, which = "stats"))
	sink(type="message")
	sink()
	close(z3)
	
	z4 <- file(paste("test10b4-",dist,".txt",sep=""), open="wt")
	sink(z4)
	print(as.data.frame(distx, which = "coefse"))
	sink(type="message")
	sink()
	close(z4)
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}


rugarch.test10c = function(dist = "ghst", skew = -0.5, shape = 8.1, cluster=NULL)
{
	tic = Sys.time()
	
	#cat("\nrugarch-->test3-4: Filter Test (sGARCH, iGARCH & apARCH)\n")
	
	data(dji30ret)
	# sGARCH Model
	# ---------------------------------------------------------------------------------
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = dist)
	
	sgarch.fit = ugarchfit(data = dji30ret[,"AA",drop=FALSE], spec = spec, 
			solver = "solnp")
	
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = dist, fixed.pars = as.list(coef(sgarch.fit)))
	sgarch.filter = ugarchfilter(data = dji30ret[,"AA",drop=FALSE], spec = spec)	
	
	# iGARCH Model
	# ---------------------------------------------------------------------------------
	spec = ugarchspec(
			variance.model = list(model = "iGARCH", garchOrder = c(2,2)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = dist)
	
	igarch.fit = ugarchfit(data = dji30ret[,"AA",drop=FALSE], spec = spec, 
			solver = "solnp")
	
	spec = ugarchspec(
			variance.model = list(model = "iGARCH", garchOrder = c(2,2)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = dist, fixed.pars = as.list(coef(igarch.fit)))
	igarch.filter = ugarchfilter(data = dji30ret[,"AA",drop=FALSE], spec = spec)
	
	# apARCH Model
	# ---------------------------------------------------------------------------------
	spec = ugarchspec(
			variance.model = list(model = "apARCH", garchOrder = c(1,2)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = dist)
	
	aparch.fit = ugarchfit(data = dji30ret[,"AA",drop=FALSE], spec = spec, 
			solver = "solnp", fit.control = list(scale = 1))
	
	spec = ugarchspec(
			variance.model = list(model = "apARCH", garchOrder = c(1,2)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = dist, fixed.pars = as.list(coef(aparch.fit)))
	aparch.filter = ugarchfilter(data = dji30ret[,"AA",drop=FALSE], spec = spec)
	
	x1 = cbind(head(sigma(sgarch.fit),10) , head(sigma(sgarch.filter),10))
	x2 = cbind(head(sigma(igarch.fit),10) , head(sigma(igarch.filter),10))
	x3 = cbind(head(sigma(aparch.fit),10) , head(sigma(aparch.filter),10))
	colnames(x1) = c("fit", "filter")
	colnames(x2) = c("fit", "filter")
	colnames(x3) = c("fit", "filter")
	zz <- file(paste("test10c-",dist,".txt",sep=""), open="wt")
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


rugarch.test10d = function(dist = "ghst", skew = -0.5, shape = 8.1, cluster=NULL)
{
	#cat("\nrugarch-->test4-2: Forecast Test (sGARCH)\n")
	tic = Sys.time()
	
	# Forecasting from the sGARCH model with variations (in sample example)
	# ---------------------------------------------------------------------------------
	data(dji30ret)
	# create weekday dummies for external regressors
	dates = rownames(dji30ret[,"MSFT", drop = FALSE])
	monday = rugarch:::.WeekDayDummy(dates, date.format = "%Y-%m-%d", weekday = "Monday")
	# convert to matrix which is what the specification expects
	monday = matrix(monday, ncol = 1)
	# create a dummy day-of-week variable for the variance regression (Friday)
	friday = rugarch:::.WeekDayDummy(dates, date.format = "%Y-%m-%d", weekday = "Friday")
	# convert to matrix which is what the specification expects
	friday = matrix(friday, ncol = 1)
	
	# ---------------------------------------------------------------------------------
	# Test 2sGARCH
	# ---------------------------------------------------------------------------------
	
	# sGARCH(1,1)
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(0,0), include.mean = FALSE), 
			distribution.model = dist)
	
	sgarch.fit1 = ugarchfit(data=dji30ret[,"MSFT", drop = FALSE], 
			out.sample = 200, spec = spec, solver = "nlminb")
	
	
	sgarch.pred1 = ugarchforecast(sgarch.fit1, n.ahead = 50, n.roll = 10)
	
	# sGARCH(1,1) + MU
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(0,0), include.mean = TRUE), 
			distribution.model = dist)
	
	sgarch.fit2 = ugarchfit(data=dji30ret[,"MSFT", drop = FALSE], 
			out.sample = 200, spec = spec, solver = "nlminb")
	
	sgarch.pred2 = ugarchforecast(sgarch.fit2, n.ahead = 50, n.roll = 10)
	
	
	# sGARCH(1,1) + MU + ARMA(1,1)
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = dist)
	
	sgarch.fit3 = ugarchfit(data=dji30ret[,"MSFT", drop = FALSE], 
			out.sample = 200, spec = spec, solver = "nlminb")
	
	sgarch.pred3 = ugarchforecast(sgarch.fit3, n.ahead = 50, n.roll = 10)
	
	# sGARCH(1,1) + XREG(V) + ARMA(1,1) + XREG(M)
	
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1), 
					external.regressors = friday), 
			mean.model = list(armaOrder = c(1,1), include.mean = FALSE, 
					external.regressors = monday), 
			distribution.model = dist, start.pars = as.list(coef(sgarch.fit3)))
	
	sgarch.fit4 = ugarchfit(data=dji30ret[,"MSFT", drop = FALSE], 
			out.sample = 200, spec = spec, solver = "gosolnp")
	
	xregm = matrix(monday[5322:(5322+50+10 - 1)], ncol=1)
	xregv = matrix(friday[5322:(5322+50+10 - 1)], ncol=1)
	sgarch.pred4 = ugarchforecast(sgarch.fit4, n.ahead = 50, n.roll = 10, 
			external.forecasts = list(mregfor = xregm, vregfor = xregv)) 
	
	# sGARCH(1,1) + XREG(V) + MU + ARMA(1,1) + XREG(M)
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1), 
					external.regressors = friday), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE, 
					external.regressors = monday), 
			distribution.model = dist, start.pars = as.list(coef(sgarch.fit4)))
	
	sgarch.fit5 = ugarchfit(data=dji30ret[,"MSFT", drop = FALSE], 
			out.sample = 200, spec = spec, solver = "gosolnp")
	
	xregm = matrix(monday[5322:(5322+50+10 - 1)], ncol=1)
	xregv = matrix(friday[5322:(5322+50+10 - 1)], ncol=1)
	sgarch.pred5 = ugarchforecast(sgarch.fit5, n.ahead = 50, n.roll = 10, 
			external.forecasts = list(mregfor = xregm, vregfor = xregv)) 
	
	# sGARCH(1,1) + XREG(V) + ARMA(1,1) + XREG(M) + INMEAN(1)
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1), 
					external.regressors = friday), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE, 
					archm = TRUE, archpow = 1, external.regressors = monday), 
			distribution.model = dist, start.pars = as.list(coef(sgarch.fit5)))
	
	sgarch.fit6 = ugarchfit(data=dji30ret[,"MSFT", drop = FALSE], 
			out.sample = 200, spec = spec, solver = "gosolnp")
	
	xregm = matrix(monday[5322:(5322+50+10 - 1)], ncol=1)
	xregv = matrix(friday[5322:(5322+50+10 - 1)], ncol=1)
	sgarch.pred6 = ugarchforecast(sgarch.fit6, n.ahead = 50, n.roll = 10, 
			external.forecasts = list(mregfor = xregm, vregfor = xregv)) 
	
	# compare the 1-step ahead rolling forecasts (10 rolls)
	sgarch.fmu  = cbind(
			fitted(sgarch.pred1)[1,], 
			fitted(sgarch.pred2)[1,],
			fitted(sgarch.pred3)[1,],
			fitted(sgarch.pred4)[1,],
			fitted(sgarch.pred5)[1,], 
			fitted(sgarch.pred6)[1,])
	
	sgarch.fsigma = cbind(
			sigma(sgarch.pred1)[1,], 
			sigma(sgarch.pred2)[1,],
			sigma(sgarch.pred3)[1,],
			sigma(sgarch.pred4)[1,],
			sigma(sgarch.pred5)[1,], 
			sigma(sgarch.pred6)[1,])
	
	sgarch.shape = cbind(
			rep(coef(sgarch.fit1)["shape"],10), 
			rep(coef(sgarch.fit2)["shape"],10),
			rep(coef(sgarch.fit3)["shape"],10), 
			rep(coef(sgarch.fit4)["shape"],10),
			rep(coef(sgarch.fit5)["shape"],10), 
			rep(coef(sgarch.fit6)["shape"],10))
	
	sgarch.skew = cbind(
			rep(coef(sgarch.fit1)["skew"],10), 
			rep(coef(sgarch.fit2)["skew"],10),
			rep(coef(sgarch.fit3)["skew"],10), 
			rep(coef(sgarch.fit4)["skew"],10),
			rep(coef(sgarch.fit5)["skew"],10), 
			rep(coef(sgarch.fit6)["skew"],10))
	# plot the forecast 1-step rolling density
	postscript(paste("test10d-",dist,".eps",sep=""), width = 12, height = 8)
	zseq = seq(-0.2, 0.2, length.out = 1000)
	colr = heat.colors(10, alpha = 1)
	par(mfrow = c(2,3))
	for(i in 1:6){
		plot(zseq, ddist(distribution = dist, y = zseq, mu = sgarch.fmu[1,i], 
						sigma = sgarch.fsigma[1,i], skew = sgarch.skew[1,i],
						shape = sgarch.shape[1,i]), 
				main = "", xlab="", ylab="", ylim=c(0,24))
		title(paste("model-", i, sep=""), line = 0.4, cex = 0.9)
		for(j in 2:10){
			lines(zseq, ddist(distribution = dist, y = zseq, mu = sgarch.fmu[j,i], 
							sigma = sgarch.fsigma[j,i], skew = sgarch.skew[1,i],
							shape = sgarch.shape[j,i]), col = colr[j])
		}
	}
	title(main = list(paste("Rolling 1-ahead Forecast Densities\nstudent distribution",sep=""), 
					cex = 1.2, col = "steelblue", font=2), outer = TRUE, line = -2)
	dev.off()
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}



rugarch.test10e = function(dist = "ghst", skew = -0.5, shape = 8.1, cluster=NULL)
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
			distribution.model = dist)
	
	sgarch.fit = ugarchfit(data = dji30ret[,"AA",drop=FALSE], spec = spec, 
			solver = "nlminb", solver.control = list(trace=0))
	
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
	zz <- file(paste("test10e-",dist,".txt",sep=""), open="wt")
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
	
	postscript(paste("test10e-",dist,".eps",sep=""), width = 12, height = 8)
	par(mfrow = c(3,3))
	for(i in 1:9) plot(sgarch.sim5, which = 3, m.sim = i)
	dev.off()
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}


rugarch.test10f= function(dist = "ghst", skew = -1, shape = 8.2, cluster=NULL)
{
	tic = Sys.time()
	set.seed(100)
	x = rdist(dist, n = 1000, skew = skew, shape = shape)
	p = pdist(dist, x, skew = skew, shape = shape)
	q = qdist(dist, p, skew = skew, shape = shape)
	
	postscript(paste("test10fa-",dist,".eps",sep=""), width = 12, height = 8)
	par(mfrow=c(2,2))
	plot(x, pch = 3, pty = 2, col = "steelblue", main = "quantile function check")
	points(q, col = "tomato1")
	legend("bottomleft", c("x", "q*(p)"), col=c("steelblue", "tomato1"), pch=c(3,1), bty="n")
	hist( abs(x-q), main = "Absolute Error (quantile function)", col = "steelblue")
	hist(x, col = "steelblue", freq  = FALSE, ylim = c(0, 0.5), main="Density of x", ylab="f(x)")
	curve(ddist(dist, x, skew = skew, shape = shape), lty = 2, lwd = 2, add = TRUE, col = "tomato1")
	plot(sort(x), ppoints(x), col  = "steelblue", xlab="x", ylab="F(x)", main = "Distribution of x")
	lines(sort(x), pdist(dist, sort(x), skew = skew, shape = shape), col = "tomato1", lwd = 2)
	dev.off()
	
	# fit test
	if(!is.null(cluster)){
		set.seed(10121)
		mu = 0.015
		sigma = 0.035
		x = matrix(rdist(dist, n = 100000, mu = mu, sigma = sigma, skew = skew, shape = shape), ncol = 100)
		parallel::clusterEvalQ(cluster, require(rugarch))
		parallel::clusterExport(cluster, c("x","dist"), envir = environment())
		pars = parallel::parLapply(cluster, as.list(1:100), fun = function(i){
					fit = try(fitdist(dist, x[,i], control=list(trace=0, tol=1e-10)))
					if(!inherits(fit, 'try-error') && fit$convergence==0){
						ans = fit$par
					} else{
						ans = rep(NA, 4)
					}
					return(ans)
				})
	} else{
		set.seed(10121)
		mu = 0.015
		sigma = 0.035
		x = matrix(rdist(dist, n = 100000, mu = mu, sigma = sigma, skew = skew, shape = shape), ncol = 100)
		pars = lapply(as.list(1:100), FUN = function(i){
					fit = try(fitdist(dist, x[,i], control=list(trace=0, tol=1e-10)))
					if(!inherits(fit, 'try-error') && fit$convergence==0){
						ans = fit$par
					} else{
						ans = rep(NA, 4)
					}
					return(ans)
				})
	}
	pars = matrix(unlist(pars), ncol = 4, byrow = TRUE)
	retmat = rbind(c(mu, sigma, skew, shape), apply(pars, 2, "mean"), apply(pars, 2, "sd"))
	colnames(retmat) = c("mu", "sigma", "skew", "shape")
	rownames(retmat) = c("real", "esimated", "s.d.")
	
	options(width=120)
	zz <- file(paste("test10fb-",dist,".txt",sep=""), open="wt")
	sink(zz)
	print(retmat)
	sink(type="message")
	sink()
	close(zz)
	
	postscript(paste("test10fb-",dist,".eps",sep=""), width = 12, height = 8)
	par(mfrow=c(2,2))
	for(i in 1:4){
		boxplot(pars[,i], main = colnames(retmat)[i])
		abline(h = retmat[1,i], col = "grey")
	}
	dev.off()
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}