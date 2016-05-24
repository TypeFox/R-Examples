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
# GARCH rolling test
#################################################################################


rugarch.test7a = function(cluster = NULL)
{
	#cat("\nrugarch-->test7-1: Roll Test (apARCH)\n")
	tic = Sys.time()
	
	data(sp500ret)
	spec = ugarchspec(
			variance.model = list(model = "eGARCH"), distribution.model = "jsu")
	roll = ugarchroll(spec,  data = sp500ret, n.ahead = 1, n.start = 1000,
			refit.every = 100, refit.window = "recursive", 
			cluster = cluster,
			solver = "hybrid", fit.control = list(scale = 1), 
			solver.control = list(tol = 1e-5, delta = 1e-6, trace=0),
			calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05))
	
	postscript("test7a1.eps", width = 12, height = 8)
	plot(roll, which = "all")
	dev.off()
	
	postscript("test7a2.eps", width = 12, height = 8)
	plot(roll, which = 5)
	dev.off()
	
	z1 <- file("test7a1.txt", open="wt")
	sink(z1)
	report(roll, type = "VaR", VaR.alpha = 0.01, conf.level = 0.95)
	sink(type="message")
	sink()
	close(z1)
	
	z2 <- file("test7a2.txt", open="wt")
	sink(z2)
	report(roll, type = "fpm")
	sink(type="message")
	sink()
	close(z2)
	
	z3 <- file("test7a3.txt", open="wt")
	sink(z3)
	print(coef(roll)[[1]])
	sink(type="message")
	sink()
	close(z3)
	
	z4 <- file("test7a4.txt", open="wt")
	sink(z4)
	print(head(as.data.frame(roll, which = "density")))
	print(tail(as.data.frame(roll, which = "density")))
	sink(type="message")
	sink()
	close(z4)
	
	z5 <- file("test7a5.txt", open="wt")
	sink(z5)
	print(head(as.data.frame(roll, which = "VaR")))
	print(tail(as.data.frame(roll, which = "VaR")))
	sink(type="message")
	sink()
	close(z5)
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}

rugarch.test7b = function(cluster = NULL)
{
	#cat("\nrugarch-->test7-2: Roll and Resume Test\n")
	tic = Sys.time()
	
	data(sp500ret)
	spec = ugarchspec(
			variance.model = list(model = "sGARCH"), distribution.model = "jsu")
	# create a solver control which is bound to fail
	roll = ugarchroll(spec,  data = sp500ret, n.ahead = 1, n.start = 1000,
			refit.every = 250, refit.window = "moving", 
			cluster = cluster,
			solver = "lbfgs", fit.control = list(scale = 1), 
			solver.control = list(trace=1),
			calculate.VaR = TRUE, VaR.alpha = c(0.01, 0.05))
	
	# see the warnings printed
	roll2 = resume(roll, cluster = cluster,
			solver = "hybrid", fit.control = list(scale = 1), 
			solver.control = list(trace=0))
	# once more to fix the one remaining window with problems
	roll2 = resume(roll2, cluster = cluster,
			solver = "solnp", solver.control = list(trace=0, 
					rho=0.5, tol=1e-6, delta=1e-5))
	
	z5 <- file("test7b.txt", open="wt")
	sink(z5)
	print(convergence(roll))
	show(roll2)
	print(convergence(roll2))
	sink(type="message")
	sink()
	close(z5)
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}