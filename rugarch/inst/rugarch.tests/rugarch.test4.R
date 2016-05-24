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
# Model Forecast Tests
#################################################################################


rugarch.test4a = function(cluster=NULL)
{
	#cat("\nrugarch-->test4-1: Forecast Test (sGARCH)\n")
	tic = Sys.time()
	
	data(dji30ret)
	# ---------------------------------------------------------------------------------
	# Test 1 : Various Tests on the forecast function
	# sGARCH(1,1) + MU + XREG(V) + ARMA(1,1) + XREG(M)
	# create weekday dummies for external regressors
	dates = rownames(dji30ret[,"AA", drop = FALSE])
	monday = rugarch:::.WeekDayDummy(dates, date.format = "%Y-%m-%d", weekday = "Monday")
	# convert to matrix which is what the specification expects
	monday = matrix(monday, ncol = 1)
	# create a dummy day-of-week variable for the variance regression (Friday)
	friday = rugarch:::.WeekDayDummy(dates, date.format = "%Y-%m-%d", weekday = "Friday")
	# convert to matrix which is what the specification expects
	friday = matrix(friday, ncol = 1)
	
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1), 
					external.regressors = friday), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE, 
					external.regressors = monday), 
			distribution.model = "std")
	
	sgarch.fit1 = ugarchfit(data = dji30ret[,"AA", drop = FALSE], 
			out.sample = 200, spec = spec, solver = "solnp")
	# first predit without providing external forecasts
	# the model will proceed as if they were zero
	sgarch.pred1 = ugarchforecast(sgarch.fit1, n.ahead = 50, n.roll = 10)
	# no add some forecasts
	xregm = matrix(monday[5322:(5322+50+10 - 1)], ncol=1)
	xregv = matrix(friday[5322:(5322+50+10 - 1)], ncol=1)
	sgarch.pred2 = ugarchforecast(sgarch.fit1, n.ahead = 50, n.roll = 10, 
			external.forecasts = list(mregfor = xregm, vregfor = xregv)) 
	# we confirm that the first forecast is the same for the 2 preditions
	# while the second starts to differ for the variance forecast (Friday)
	# and the third start to differ for the mean forecast(Monday).
	#head(cbind(xregm, xregv), 5)
	#head(cbind(as.data.frame(sgarch.pred1), as.data.frame(sgarch.pred2)), 5)
	
	# now we overide with some other external forecasts
	# note that they must be equal to n.ahead + n.roll
	xregm = matrix(rep(0, 60), ncol=1)
	xregv = matrix(rep(0, 60), ncol=1)
	sgarch.pred3 = ugarchforecast(sgarch.fit1, n.ahead = 50, n.roll = 10, 
			external.forecasts = list(mregfor = xregm, vregfor = xregv)) 
	# they shoud be the same
	#head(cbind(as.data.frame(sgarch.pred1), as.data.frame(sgarch.pred3)), 5)
	
	# we might also like to see some of the rolling forecasts:
	# which we do by calling the as.data.frame extractor with
	# rollframe="all". Read the documentation for a description of how to read the
	# matrix
	
	options(width=150)
	zz <- file("test4a.txt", open="wt")
	sink(zz)
	print(head(sigma(sgarch.pred3)))
	print(head(fitted(sgarch.pred3)))
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

rugarch.test4b = function(cluster=NULL)
{
	#cat("\nrugarch-->test4-2: Forecast Test (sGARCH)\n")
	tic = Sys.time()
	
	# Forecasting from the sGARCH model with variations (in sample example)
	# ---------------------------------------------------------------------------------
	data(dji30ret)
	# create weekday dummies for external regressors
	dates = rownames(dji30ret[,"AA", drop = FALSE])
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
			distribution.model = "std")
	
	sgarch.fit1 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], 
			out.sample = 200, spec = spec, solver = "solnp")
	
	
	sgarch.pred1 = ugarchforecast(sgarch.fit1, n.ahead = 50, n.roll = 10)
	
	# sGARCH(1,1) + MU
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(0,0), include.mean = TRUE), 
			distribution.model = "std")
	
	sgarch.fit2 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], 
			out.sample = 200, spec = spec, solver = "solnp")
	
	sgarch.pred2 = ugarchforecast(sgarch.fit2, n.ahead = 50, n.roll = 10)
	
	
	# sGARCH(1,1) + MU + ARMA(1,1)
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std")
	
	sgarch.fit3 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], 
			out.sample = 200, spec = spec, solver = "solnp")
	
	sgarch.pred3 = ugarchforecast(sgarch.fit3, n.ahead = 50, n.roll = 10)
	
	# sGARCH(1,1) + XREG(V) + ARMA(1,1) + XREG(M)
	
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1), 
					external.regressors = friday), 
			mean.model = list(armaOrder = c(1,1), include.mean = FALSE, 
					external.regressors = monday), 
			distribution.model = "std", start.pars = as.list(coef(sgarch.fit3)))
	
	sgarch.fit4 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], 
			out.sample = 200, spec = spec, solver = "solnp")
	
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
			distribution.model = "std", start.pars = as.list(coef(sgarch.fit4)))
	
	sgarch.fit5 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], 
			out.sample = 200, spec = spec, solver = "solnp")
	
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
			distribution.model = "std", start.pars = as.list(coef(sgarch.fit5)))
	
	sgarch.fit6 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], 
			out.sample = 200, spec = spec, solver = "solnp")
	
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
	
	
	# plot the forecast 1-step rolling density
	postscript("test4b.eps", width = 12, height = 8)
	zseq = seq(-0.2, 0.2, length.out = 1000)
	colr = heat.colors(10, alpha = 1)
	par(mfrow = c(2,3))
	for(i in 1:6){
		plot(zseq, ddist(distribution = "std", y = zseq, mu = sgarch.fmu[1,i], 
						sigma = sgarch.fsigma[1,i], 
						shape = sgarch.shape[1,i]), 
				main = "", xlab="", ylab="", ylim=c(0,24))
		title(paste("model-", i, sep=""), line = 0.4, cex = 0.9)
		for(j in 2:10){
			lines(zseq, ddist(distribution = "std", y = zseq, mu = sgarch.fmu[j,i], 
							sigma = sgarch.fsigma[j,i], 
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


rugarch.test4c = function(cluster=NULL)
{
	#cat("\nrugarch-->test4-3: Forecast Test (gjrGARCH)\n")
	tic = Sys.time()
	
	data(dji30ret)
	# create weekday dummies for external regressors
	dates = rownames(dji30ret[,"AA", drop = FALSE])
	monday = rugarch:::.WeekDayDummy(dates, date.format = "%Y-%m-%d", weekday = "Monday")
	# convert to matrix which is what the specification expects
	monday = matrix(monday, ncol = 1)
	# create a dummy day-of-week variable for the variance regression (Friday)
	friday = rugarch:::.WeekDayDummy(dates, date.format = "%Y-%m-%d", weekday = "Friday")
	# convert to matrix which is what the specification expects
	friday = matrix(friday, ncol = 1)
	# gjrGARCH(1,1)
	spec = ugarchspec(
			variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(0,0), include.mean = FALSE), 
			distribution.model = "sstd")
	
	gjrgarch.fit1 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], 
			out.sample = 200, spec = spec, solver = "solnp")
	
	gjrgarch.pred1 = ugarchforecast(gjrgarch.fit1, n.ahead = 50, n.roll = 10)
	
	# gjrGARCH(1,1) + MU
	spec = ugarchspec(
			variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(0,0), include.mean = TRUE), 
			distribution.model = "sstd")
	
	gjrgarch.fit2 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], 
			out.sample = 200, spec = spec, solver = "solnp")
	
	gjrgarch.pred2 = ugarchforecast(gjrgarch.fit2, n.ahead = 50, n.roll = 10)
	
	
	# gjrGARCH(1,1) + MU + ARMA(1,1)
	spec = ugarchspec(
			variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "sstd")
	
	gjrgarch.fit3 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], 
			out.sample = 200, spec = spec, solver = "solnp")
	
	gjrgarch.pred3 = ugarchforecast(gjrgarch.fit3, n.ahead = 50, n.roll = 10)
	
	# gjrGARCH(1,1) + XREG(V) + ARMA(1,1) + XREG(M)
	spec = ugarchspec(
			variance.model = list(model = "gjrGARCH", garchOrder = c(1,1), 
					external.regressors = friday), 
			mean.model = list(armaOrder = c(1,1), include.mean = FALSE, 
					external.regressors = monday), 
			distribution.model = "sstd", 
			start.pars = as.list(coef(gjrgarch.fit3)))
	
	gjrgarch.fit4 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], 
			out.sample = 200, spec = spec, solver = "solnp")
	
	xregm = matrix(monday[5322:(5322+50+10 - 1)], ncol=1)
	xregv = matrix(friday[5322:(5322+50+10 - 1)], ncol=1)
	gjrgarch.pred4 = ugarchforecast(gjrgarch.fit4, n.ahead = 50, n.roll = 10, 
			external.forecasts = list(mregfor = xregm, vregfor = xregv)) 
	
	# gjrGARCH(1,1) + XREG(V) + MU + ARMA(1,1) + XREG(M)
	spec = ugarchspec(
			variance.model = list(model = "gjrGARCH", garchOrder = c(1,1), 
					external.regressors = friday), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE, 
					external.regressors = monday), 
			distribution.model = "sstd", 
			start.pars = as.list(coef(gjrgarch.fit4)))
	
	gjrgarch.fit5 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], 
			out.sample = 200, spec = spec, solver = "solnp")
	
	xregm = matrix(monday[5322:(5322+50+10 - 1)], ncol=1)
	xregv = matrix(friday[5322:(5322+50+10 - 1)], ncol=1)
	gjrgarch.pred5 = ugarchforecast(gjrgarch.fit5, n.ahead = 50, n.roll = 10, 
			external.forecasts = list(mregfor = xregm, vregfor = xregv)) 
	
	# gjrGARCH(1,1) + XREG(V) + ARMA(1,1) + XREG(M) + INMEAN(1)
	spec = ugarchspec(
			variance.model = list(model = "gjrGARCH", garchOrder = c(1,1), 
					external.regressors = friday), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE, 
					archm = TRUE, archpow = 1, external.regressors = monday), 
			distribution.model = "sstd", 
			start.pars = as.list(coef(gjrgarch.fit5)))
	
	gjrgarch.fit6 = ugarchfit(data=dji30ret[,"AA", drop = FALSE], 
			out.sample = 200, spec = spec, solver = "solnp")
	
	xregm = matrix(monday[5322:(5322+50+10 - 1)], ncol=1)
	xregv = matrix(friday[5322:(5322+50+10 - 1)], ncol=1)
	gjrgarch.pred6 = ugarchforecast(gjrgarch.fit6, n.ahead = 50, n.roll = 10, 
			external.forecasts = list(mregfor = xregm, vregfor = xregv)) 
		
	gjrgarch.fmu  = cbind(
			fitted(gjrgarch.pred1)[1,], 
			fitted(gjrgarch.pred2)[1,],
			fitted(gjrgarch.pred3)[1,], 
			fitted(gjrgarch.pred4)[1,], 
			fitted(gjrgarch.pred5)[1,], 
			fitted(gjrgarch.pred6)[1,])
	
	gjrgarch.fsigma = cbind(
			sigma(gjrgarch.pred1)[1,], 
			sigma(gjrgarch.pred2)[1,],
			sigma(gjrgarch.pred3)[1,], 
			sigma(gjrgarch.pred4)[1,], 
			sigma(gjrgarch.pred5)[1,], 
			sigma(gjrgarch.pred6)[1,])
	
	gjrgarch.skew = cbind(
			rep(coef(gjrgarch.fit1)["skew"],10), 
			rep(coef(gjrgarch.fit2)["skew"],10), 
			rep(coef(gjrgarch.fit3)["skew"],10), 
			rep(coef(gjrgarch.fit4)["skew"],10), 
			rep(coef(gjrgarch.fit5)["skew"],10), 
			rep(coef(gjrgarch.fit6)["skew"],10))
	
	gjrgarch.shape = cbind(
			rep(coef(gjrgarch.fit1)["shape"],10), 
			rep(coef(gjrgarch.fit2)["shape"],10), 
			rep(coef(gjrgarch.fit3)["shape"],10), 
			rep(coef(gjrgarch.fit4)["shape"],10), 
			rep(coef(gjrgarch.fit5)["shape"],10), 
			rep(coef(gjrgarch.fit6)["shape"],10))
	
	# plot the forecast 1-step rolling density
	postscript("test4c.eps", width = 12, height = 8)
	zseq = seq(-0.2, 0.2, length.out = 1000)
	colr = heat.colors(10, alpha = 1)
	par(mfrow = c(2,3))
	for(i in 1:6){
		plot(zseq, ddist(distribution = "sstd", y = zseq, mu = gjrgarch.fmu[1,i], 
						sigma= gjrgarch.fsigma[1,i], skew = gjrgarch.skew[1,i], 
						shape = gjrgarch.shape[1,i]), 
				main = "", xlab="", ylab="", ylim=c(0,24))
		title(paste("model-", i, sep=""), line = 0.4, cex = 0.9)
		for(j in 2:10){
			lines(zseq, ddist( distribution = "sstd", y = zseq, mu = gjrgarch.fmu[j,i], 
							sigma = gjrgarch.fsigma[j,i], shape = gjrgarch.shape[j,i], 
							skew = gjrgarch.skew[j,i]), col = colr[j])
		}
	}
	title(main = list(paste("Rolling 1-ahead Forecast Densities\nskew student distribution", 
							sep = ""), cex = 1.2, col = "steelblue", font=2), 
			outer = TRUE, line = -2)
	dev.off()
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

rugarch.test4d = function(cluster=NULL)
{
	#cat("\nrugarch-->test4-4: Forecast Performance Measures Test (sGARCH)\n")
	tic = Sys.time()
	
	# fpm tests
	data(dmbp)
	spec = ugarchspec(
			variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
			mean.model = list(armaOrder = c(1,1), include.mean = TRUE), 
			distribution.model = "std")
	
	fit1 = ugarchfit(data = dmbp[,1,drop=FALSE], out.sample = 50, spec = spec)
	
	pred1 = ugarchforecast(fit1, n.ahead = 50, n.roll = 2)
	pred2 = ugarchforecast(fit1, n.ahead = 1, n.roll = 50)
	
	options(width=150)
	zz <- file("test4d.txt", open="wt")
	sink(zz)
	print(fpm(pred1, summary = TRUE))
	print(fpm(pred1, summary = FALSE))
	print(fpm(pred2, summary = TRUE))
	print(fpm(pred2, summary = FALSE))
	sink(type="message")
	sink()
	close(zz)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}