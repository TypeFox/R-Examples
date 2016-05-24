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

# bootstrap

rugarch.test12a = function(cluster=NULL)
{
	tic = Sys.time()
	#cat("\nrugarch-->test12-1: Partial Bootstrap\n")
	data(sp500ret)
	spec = ugarchspec(variance.model=list(model="csGARCH"), distribution="std")
	fit = ugarchfit(spec, sp500ret, out.sample=500)
	bootp  = ugarchboot(fit, method = c("Partial", "Full")[1], 
			n.ahead = 500, n.bootpred = 500,
			cluster = cluster, verbose = TRUE)
	
	options(width=120)
	zz <- file("test12a-1.txt", open="wt")
	sink(zz)
	show(bootp)
	sink(type="message")
	sink()
	close(zz)
	
	options(width=120)
	zz <- file("test12a-2.txt", open="wt")
	sink(zz)
	# n.bootpred by n.ahead matrix:
	cat("\nSeries Forecasts:\n")
	print(head(ser<-as.data.frame(bootp, which = "series", type = "raw")))
	cat("\nSeries Summary:\n")
	print(as.data.frame(bootp, which = "series", type = "summary"))
	cat("\nQuantile Summary:\n")
	print(as.data.frame(bootp, which = "series", type = "q", qtile= c(0.25, 0.75)))
	sink(type="message")
	sink()
	close(zz)
	
	options(width=120)
	zz <- file("test12a-3.txt", open="wt")
	sink(zz)
	# n.bootpred by n.ahead matrix:
	cat("\nSigma Forecasts:\n")
	print(head(sig<-as.data.frame(bootp, which = "sigma", type = "raw")))
	cat("\nSigma Summary:\n")
	print(as.data.frame(bootp, which = "sigma", type = "summary"))
	sink(type="message")
	sink()
	close(zz)
	
	postscript("test12a-1.eps")
	plot(bootp, which = "all")
	dev.off()
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
	
}


rugarch.test12b = function(cluster=NULL)
{
	tic = Sys.time()
	#cat("\nrugarch-->test12-2: Full Bootstrap\n")
	data(sp500ret)
	spec = ugarchspec(variance.model=list(model="csGARCH"), distribution="std")
	fit = ugarchfit(spec, sp500ret, out.sample=10)
	bootp  = ugarchboot(fit, method = c("Partial", "Full")[2], 
			n.ahead = 10, n.bootpred = 500, n.bootfit = 100,
			cluster = cluster)
	
	options(width=120)
	zz <- file("test12b-1.txt", open="wt")
	sink(zz)
	show(bootp)
	sink(type="message")
	sink()
	close(zz)
	
	options(width=120)
	zz <- file("test12b-2.txt", open="wt")
	sink(zz)
	# n.bootpred by n.ahead matrix:
	cat("\nSeries Forecasts:\n")
	print(head(ser<-as.data.frame(bootp, which = "series", type = "raw")))
	cat("\nSeries Summary:\n")
	print(as.data.frame(bootp, which = "series", type = "summary"))
	cat("\nQuantile Summary:\n")
	print(as.data.frame(bootp, which = "series", type = "q", qtile= c(0.25, 0.75)))
	sink(type="message")
	sink()
	close(zz)
	
	options(width=120)
	zz <- file("test12b-3.txt", open="wt")
	sink(zz)
	# n.bootpred by n.ahead matrix:
	cat("\nSigma Forecasts:\n")
	print(head(sig<-as.data.frame(bootp, which = "sigma", type = "raw")))
	cat("\nSigma Summary:\n")
	print(as.data.frame(bootp, which = "sigma", type = "summary"))
	sink(type="message")
	sink()
	close(zz)
	
	postscript("test12b-1.eps")
	plot(bootp, which = "all")
	dev.off()
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
	
}


rugarch.test12c = function(cluster=NULL)
{
	tic = Sys.time()
	#cat("\nrugarch-->test12-3: Partial Bootstrap/SPD\n")
	data(sp500ret)
	library(spd)
	spec = ugarchspec(variance.model=list(model="csGARCH"), distribution="std")
	fit = ugarchfit(spec, sp500ret, out.sample=125)
	bootp  = ugarchboot(fit, method = c("Partial", "Full")[1], 
			sampling = "spd", n.ahead = 125, n.bootpred = 500,
			cluster = cluster)
	
	options(width=120)
	zz <- file("test12c-1.txt", open="wt")
	sink(zz)
	show(bootp)
	sink(type="message")
	sink()
	close(zz)
	
	options(width=120)
	zz <- file("test12c-2.txt", open="wt")
	sink(zz)
	# n.bootpred by n.ahead matrix:
	cat("\nSeries Forecasts:\n")
	print(head(ser<-as.data.frame(bootp, which = "series", type = "raw")))
	cat("\nSeries Summary:\n")
	print(as.data.frame(bootp, which = "series", type = "summary"))
	cat("\nQuantile Summary:\n")
	print(as.data.frame(bootp, which = "series", type = "q", qtile= c(0.25, 0.75)))
	sink(type="message")
	sink()
	close(zz)
	
	options(width=120)
	zz <- file("test12c-3.txt", open="wt")
	sink(zz)
	# n.bootpred by n.ahead matrix:
	cat("\nSigma Forecasts:\n")
	print(head(sig<-as.data.frame(bootp, which = "sigma", type = "raw")))
	cat("\nSigma Summary:\n")
	print(as.data.frame(bootp, which = "sigma", type = "summary"))
	sink(type="message")
	sink()
	close(zz)
	
	postscript("test12c-1.eps")
	plot(bootp, which = "all")
	dev.off()
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
	
}


rugarch.test12d = function(cluster=NULL)
{
	tic = Sys.time()
	#cat("\nrugarch-->test12-3: Partial Bootstrap/Kernel\n")
	data(sp500ret)
	library(ks)
	spec = ugarchspec(variance.model=list(model="csGARCH"), distribution="std")
	fit = ugarchfit(spec, sp500ret, out.sample=125)
	bootp  = ugarchboot(fit, method = c("Partial", "Full")[1], 
			sampling = "kernel", n.ahead = 125, n.bootpred = 500,
			cluster = cluster)
	
	options(width=120)
	zz <- file("test12d-1.txt", open="wt")
	sink(zz)
	show(bootp)
	sink(type="message")
	sink()
	close(zz)
	
	options(width=120)
	zz <- file("test12d-2.txt", open="wt")
	sink(zz)
	# n.bootpred by n.ahead matrix:
	cat("\nSeries Forecasts:\n")
	print(head(ser<-as.data.frame(bootp, which = "series", type = "raw")))
	cat("\nSeries Summary:\n")
	print(as.data.frame(bootp, which = "series", type = "summary"))
	cat("\nQuantile Summary:\n")
	print(as.data.frame(bootp, which = "series", type = "q", qtile= c(0.25, 0.75)))
	sink(type="message")
	sink()
	close(zz)
	
	options(width=120)
	zz <- file("test12d-3.txt", open="wt")
	sink(zz)
	# n.bootpred by n.ahead matrix:
	cat("\nSigma Forecasts:\n")
	print(head(sig<-as.data.frame(bootp, which = "sigma", type = "raw")))
	cat("\nSigma Summary:\n")
	print(as.data.frame(bootp, which = "sigma", type = "summary"))
	sink(type="message")
	sink()
	close(zz)
	
	postscript("test12d-1.eps")
	plot(bootp, which = "all")
	dev.off()
	
	toc = Sys.time()-tic
	cat("Elapsed:", toc, "\n")
	return(toc)
	
}