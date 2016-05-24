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

# Check the csGARCH model
rugarch.test11a = function(cluster=NULL)
{
	tic = Sys.time()
	#cat("\nrugarch-->test11-1: Component GARCH parameter distribution\n")
	data(sp500ret)
	
	spec = ugarchspec(variance.model = list(model="csGARCH"))
	fit = ugarchfit(spec, sp500ret)
	
	dist = ugarchdistribution(fit, n.sim = 1000, n.start = 1, 
			m.sim = 200,  recursive = TRUE, recursive.length = 6000, recursive.window = 1000,
			fit.control = list(), solver = "solnp", solver.control = list(), cluster = cluster)
	

	options(width=120)
	zz <- file("test11a.txt", open="wt")
	sink(zz)
	show(dist)
	sink(type="message")
	sink()
	close(zz)
	
	postscript("test11a-1.eps")
	plot(dist, which = 1, window=6)
	dev.off()
	
	png("test11a-1.png", width=800, height=800)
	plot(dist, which = 2, window = 6)
	dev.off()
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
}
