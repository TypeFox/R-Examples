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
# GARCH distribution tests
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
# non recursive, simple test

rugarch.test8a = function(cluster=NULL)
{
	#cat("\nRuGARCH-->test8-1: Parameter Distribution Test - Non recursive (fGARCH/NAGARCH)\n")
	tic = Sys.time()
	
	data(dji30ret)
	rseed = rugarch.seeds(100, "test8a")
	spec = ugarchspec(
			variance.model = list(model = "fGARCH", submodel = "NAGARCH"), 
			distribution.model = "std")
	fit = ugarchfit(spec, data = dji30ret[,"XOM",drop=FALSE], fit.control = list(scale = 1))
	dist = ugarchdistribution(fit, n.sim = 2000, n.start = 50, m.sim = 100, 
			solver = "solnp", fit.control = list(scale = 1), rseed  = rseed, 
			cluster = cluster)
	#.dist1 <<- dist
	
	png("test8a1.png", width = 1000, height = 800)
	plot(dist, which = 1)
	dev.off()
	
	png("test8a2.png", width = 1000, height = 800)
	plot(dist, which = 2)
	dev.off()
	
	png("test8a3.png", width = 1000, height = 800)
	plot(dist, which = 3)
	dev.off()
	
	z1 <- file("test8a1.txt", open="wt")
	sink(z1)
	print(as.data.frame(dist, which = "coef"))
	sink(type="message")
	sink()
	close(z1)
	
	z2 <- file("test8a2.txt", open="wt")
	sink(z2)
	print(as.data.frame(dist, which = "rmse"))
	sink(type="message")
	sink()
	close(z2)
	
	
	z3 <- file("test8a3.txt", open="wt")
	sink(z3)
	print(as.data.frame(dist, which = "stats"))
	sink(type="message")
	sink()
	close(z3)
	
	z4 <- file("test8a4.txt", open="wt")
	sink(z4)
	print(as.data.frame(dist, which = "coefse"))
	sink(type="message")
	sink()
	close(z4)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

# recursive test and asymptotic consistency
rugarch.test8b = function(cluster=NULL)
{
	#cat("\nrgarch-->test8-2: Parameter Distribution Test - recursive (fGARCH/ALLGARCH)\n")
	tic = Sys.time()
	
	data(dji30ret)
	spec = ugarchspec(
			variance.model = list(model = "fGARCH", submodel = "ALLGARCH"), 
			distribution.model = "std")
	fit = ugarchfit(spec, data = dji30ret[,"XOM"], solver.control = list(scale = 0))
	dist = ugarchdistribution(fit, n.sim = 1000, n.start = 50, m.sim = 100, 
			recursive = TRUE, recursive.length = 8000, recursive.window = 1000,
			solver.control = list(tol = 1e-8, delta = 1e-7), fit.control = list(scale = 0), 
			cluster = cluster)
	# save(dist, file = "dist.rda")
	
	png("test8b1.png", width = 1000, height = 800)
	plot(dist, which = 1, window = 8)
	dev.off()
	
	png("test8b2.png", width = 1000, height = 800)
	plot(dist, which = 2, window = 6)
	dev.off()
	
	png("test8b3.png", width = 1000, height = 800)
	plot(dist, which = 3)
	dev.off()
	
	png("test8b4.png", width = 1000, height = 800)
	plot(dist, which = 4)
	dev.off()
	
	zz1 <- file("test8b2.txt", open="wt")
	sink(zz1)
	print(as.data.frame(dist, which = "coef", window = 6))
	sink(type="message")
	sink()
	close(zz1)
	
	z2 <- file("test8b2.txt", open="wt")
	sink(z2)
	print(rbind(as.data.frame(dist, which = "rmse", window = 1), 
			as.data.frame(dist, which = "rmse", window = 2),
			as.data.frame(dist, which = "rmse", window = 3),
			as.data.frame(dist, which = "rmse", window = 4),
			as.data.frame(dist, which = "rmse", window = 5),
			as.data.frame(dist, which = "rmse", window = 6)))
	sink(type="message")
	sink()
	close(z2)
	
	
	z3 <- file("test8b3.txt", open="wt")
	sink(z3)
	print(as.data.frame(dist, which = "stats", window = 6))
	sink(type="message")
	sink()
	close(z3)
	
	z4 <- file("test8b4.txt", open="wt")
	sink(z4)
	print(as.data.frame(dist, which = "coefse", window = 6))
	sink(type="message")
	sink()
	close(z4)
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}



# ARFIMA - GARCH benchmark
rugarch.test8c = function(cluster=NULL)
{
	tic = Sys.time()
	# change this for testing other models and change the truecoef
	vmodel = "sGARCH"
	vsubmodel = NULL
	# Check Simulation of C code versus C++
	# sim1 is C code (n.sim<20 && m.sim>100)
	truecoef = list(mu = 0.005, ar1 = 0.6, ar2 = 0.01, ma1 = -0.7, arfima = 0.3, 
			omega = 2.5e-5, alpha1  = 0.03, beta1 = 0.93)
	spec = ugarchspec( mean.model = list(armaOrder = c(2,1), include.mean = TRUE, 
					arfima = TRUE), 
			variance.model = list(model = vmodel, submodel = vsubmodel),
			distribution.model = "norm", 
			fixed.pars = truecoef)
	
	sim1 = ugarchpath(spec, n.sim = 15, n.start = 0, m.sim = 10, rseed = c(100:109))
	# sim2 is C++ code
	sim2 = ugarchpath(spec, n.sim = 15, n.start = 0, m.sim = 101, rseed = c(100:200))
	
	options(width=100)
	zz <- file("test8c1.txt", open="wt")
	sink(zz)
	print(all.equal(sim1@path$seriesSim[1:10,1:10],  sim2@path$seriesSim[1:10,1:10]))
	print(all.equal(sim1@path$seriesSim[11:15,1:10], sim2@path$seriesSim[11:15,1:10]))
	sink(type="message")
	sink()
	close(zz)
	
	# ARFIMA(2,d,1)-sGARCH(1,1)
	truecoef1 = list(mu = 0.005, ar1 = 0.5, ar2 = 0.1, ma1 = -0.7, arfima = 0.3, 
			omega = 2.5e-6, alpha1 = 0.03, beta1 = 0.94)
	spec1 = ugarchspec( 
			mean.model = list(armaOrder = c(2,1), include.mean = TRUE, arfima = TRUE), 
			variance.model = list(model = vmodel, submodel = vsubmodel),
			distribution.model = "norm", fixed.pars = truecoef1)
	sim1 = ugarchpath(spec1, n.sim = 5000, n.start = 100, m.sim = 1, rseed = 125)
	data1 = fitted(sim1)
	spec1 = ugarchspec( 
			mean.model = list(armaOrder = c(2,1), include.mean = TRUE, arfima = TRUE), 
			variance.model = list(model = vmodel, submodel = vsubmodel),
			distribution.model = "norm", start.pars = truecoef1)
	fit1 = ugarchfit(spec1, data = data1, solver = "solnp")
	chk1 = cbind(round(coef(fit1),4), round(unlist(truecoef1),4))
	colnames(chk1) = c("rugarch", "true")
	
	# ARFIMA(2,d,0)-sGARCH(1,1)
	truecoef2 = list(mu = 0.005, ar1 = 0.5, ar2 = 0.1, arfima = 0.3, 
			omega = 2.5e-6, alpha1 = 0.03, beta1 = 0.94)
	spec2 = ugarchspec(
			mean.model = list(armaOrder = c(2,0), include.mean = TRUE, arfima = TRUE), 
			variance.model = list(model = vmodel, submodel = vsubmodel),
			distribution.model = "norm", fixed.pars = truecoef2)
	sim2 = ugarchpath(spec2, n.sim = 5000, n.start = 100, m.sim = 1, rseed = 125)
	data2 = fitted(sim2)
	spec2 = ugarchspec( 
			mean.model = list(armaOrder = c(2,0), include.mean = TRUE, arfima = TRUE), 
			variance.model = list(model = vmodel, submodel = vsubmodel),
			distribution.model = "norm", start.pars = truecoef2)
	fit2 = ugarchfit(spec2, data = data2, solver = "solnp")
	chk2 = cbind(round(coef(fit2),4), round(unlist(truecoef2),4))
	colnames(chk2) = c("rugarch", "true")
	
	# ARFIMA(0,d,2)-sGARCH(1,1)
	truecoef3 = list(mu = 0.005, ma1 = 0.5, ma2 = 0.1, arfima = 0.3, 
			omega = 2.5e-6, alpha1 = 0.03, beta1 = 0.94)
	spec3 = ugarchspec(
			mean.model = list(armaOrder = c(0,2), include.mean = TRUE, arfima = TRUE), 
			variance.model = list(model = vmodel, submodel = vsubmodel),
			distribution.model = "norm", fixed.pars = truecoef3)
	sim3 = ugarchpath(spec3, n.sim = 5000, n.start = 100, m.sim = 1, rseed = 125)
	data3 = fitted(sim3)
	spec3 = ugarchspec( 
			mean.model = list(armaOrder = c(0,2), include.mean = TRUE, arfima = TRUE), 
			variance.model = list(model = vmodel, submodel = vsubmodel),
			distribution.model = "norm", start.pars = truecoef3)
	fit3 = ugarchfit(spec3, data = data3, solver = "solnp")
	chk3 = cbind(round(coef(fit3),4), round(unlist(truecoef3),4))
	colnames(chk3) = c("rugarch", "true")
	
	
	
	options(width=100)
	zz <- file("test8c2.txt", open="wt")
	sink(zz)
	cat(paste("\nARFIMA(2,d,1)-",vmodel,"(1,1)\n", sep = ""))
	print(chk1)
	cat(paste("\nARFIMA(2,d,0)-",vmodel,"(1,1)\n", sep = ""))
	print(chk2)
	cat(paste("\nARFIMA(0,d,2)-",vmodel,"(1,1)\n", sep = ""))
	print(chk3)
	sink(type="message")
	sink()
	close(zz)
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}