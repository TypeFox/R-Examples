#################################################################################
##
##   R package rmgarch by Alexios Ghalanos Copyright (C) 2008-2013.
##   This file is part of the R package rmgarch.
##
##   The R package rmgarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rmgarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################


.plotdccfit = function(x, which = "ask", series = c(1, 2), ...)
{
	n = dim(x@model$umodel$modelinc)[2]
	xseries = unique( as.integer( series ) )
	xseries = xseries[xseries > 0]
	if( any(xseries > n) ) stop( "\nrmgarch-->error: series index out of bounds in DCCfit plot.")
	if( length(xseries) < 2 ) stop( "\nrmgarch-->error: series length must be 2 or more.")
	
	series = xseries
	choices = c(
			"Conditional Mean (vs Realized Returns)",
			"Conditional Sigma (vs Realized Absolute Returns)",
			"Conditional Covariance",
			"Conditional Correlation",
			"EW Portfolio Plot with conditional density VaR limits")
	.interdccfitPlot(x, choices = choices, plotFUN = paste(".plot.dccfit", 1:5, sep = "."), which = which, series = series, ...)
	# Return Value:
	invisible(x)
}

.interdccfitPlot = function(x, choices, plotFUN, which, series = c(1, 2), ...)
{
	if (is.numeric(which)) {
		if(which>length(choices)) stop("Not a valid choice. Plots choices are 1-5.\n",call. = FALSE)
		FUN = match.fun(plotFUN[which])
		FUN(x, series = series, ...)
	}
	if(is.character(which))
	{
		if( which!="ask" ) stop("Not a valid choice.\n",call. = FALSE)
		.multdccfitPlot(x, choices, series = series, ...)
	}
	
	invisible(x)
}

.multdccfitPlot = function(x, choices, series = c(1, 2), ...)
{
	
	pick = 1
	while (pick > 0) {
		pick = menu (
				choices = paste(" ", choices),
				title = "\nMake a plot selection (or 0 to exit):")
		switch (pick,
				.plot.dccfit.1(x, series, ...),  .plot.dccfit.2(x, series, ...),  
				.plot.dccfit.3(x, series, ...),  .plot.dccfit.4(x, series, ...),
				.plot.dccfit.5(x, ...))
	}
}

.plot.dccfit.1 = function(x, series = c(1, 2), ...)
{
	ops = list(...)
	n = length( series )
	cnames = x@model$modeldata$asset.names
	p = x@model$maxgarchOrder
	T = x@model$modeldata$T
	y = fitted(x)
	if( n > 16 ){
		scr = floor( n/16 )
		z	= n - 16 * floor( n/16 )
		start	= dev.next( which = dev.cur() )
		xn	= 0
		for( j in 1:scr ){
			dev.new( start + j )
			par( mfrow = c(4, 4) )
			for(i in 1:16){
				plot(  y[, series[i + xn]], col = colors()[16], 
						main = cnames[series[i + xn]], ylab = "", xlab = "",
						minor.ticks = FALSE, auto.grid=FALSE)
				lines( y[, series[i + xn]], col = colors()[132])
				mtext(paste("rmgarch  : DCC model fit"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
			}
			title( paste( "DCC Conditional Mean\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			grid()
			xn = xn + 16
		}
		if( z != 0 ){
			dev.new( dev.next( which = dev.cur() ) + 1 )
			par( mfrow = c(4, 4) )
			for(i in 1:z){
				plot(  y[, series[i + xn]], type = "l", col = colors()[16], 
						main = cnames[series[i + xn]], ylab = "", xlab = "",
						minor.ticks = FALSE, auto.grid=FALSE)
				lines( y[, series[i + xn]], type = "l", col = colors()[132])
				mtext(paste("rmgarch  : DCC model fit"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
			}
			title( paste( "DCC Conditional Mean\n(page...", j, ")", sep = ""), outer =TRUE, line = -1.5, cex = 0.75)
			grid()
		}
	} else{
		d = .divisortable( n )
		par( mfrow= c(d[1], d[2]) )
		for(i in 1:n){
			plot(  y[, series[i]], type = "l", col = colors()[16], 
					main = cnames[series[i]], ylab = "", xlab = "",
					minor.ticks = FALSE, auto.grid=FALSE)
			lines( y[, series[i]], type = "l", col = colors()[132])
			mtext(paste("rmgarch  : DCC model fit"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
			title( paste( "DCC Conditional Mean", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			grid()
		}
	}
	invisible()
}


.plot.dccfit.2 = function(x, series = c(1, 2), ...)
{
	ops = list(...)
	
	ops = list(...)
	n = length( series )
	cnames = x@model$modeldata$asset.names
	p = x@model$maxgarchOrder
	T = x@model$modeldata$T
	actual = xts(x@model$modeldata$data[1:T, ], x@model$modeldata$index[1:T])
	y = sigma(x)	
	if( n > 16 ){
		scr = floor( n/16 )
		z	= n - 16 * floor( n/16 )
		start	= dev.next( which = dev.cur() )
		xn	= 0
		for( j in 1:scr ){
			dev.new( start + j )
			par( mfrow = c(4, 4) )
			for(i in 1:16){
				plot(  abs(actual[, series[i + xn]]), type = "l", col = colors()[16], 
						main = cnames[series[i + xn]], ylab = "", xlab = "",
						minor.ticks = FALSE, auto.grid=FALSE)
				lines( y[, series[i + xn]], type = "l", col = colors()[132])
				mtext(paste("rmgarch  : DCC model fit"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
			}
			title( paste( "DCC Conditional Sigma vs |returns|\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			grid()
			xn = xn + 16
		}
		if( z != 0 ){
			dev.new( dev.next( which = dev.cur() ) + 1 )
			par( mfrow = c(4, 4) )
			for(i in 1:z){
				plot(  abs(actual[, series[i + xn]]), type = "l", col = colors()[16], 
						main = cnames[series[i + xn]], ylab = "", xlab = "",
						minor.ticks = FALSE, auto.grid=FALSE)
				lines( y[, series[i + xn]], type = "l", col = colors()[132])
				mtext(paste("rmgarch  : DCC model fit"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
			}
			title( paste( "DCC Conditional Sigma vs |returns|\n(page...", j, ")", sep = ""), outer =TRUE, line = -1.5, cex = 0.75)
			grid()
		}
	} else{
		d = .divisortable( n )
		par( mfrow= c(d[1], d[2]) )
		for(i in 1:n){
			plot(  abs(actual[, series[i]]), type = "l", col = colors()[16], 
					main = cnames[series[i]], ylab = "", xlab = "",
					minor.ticks = FALSE, auto.grid=FALSE)
			lines( y[, series[i]], type = "l", col = colors()[132])
			mtext(paste("rmgarch  : DCC model fit"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
			title( paste( "DCC Conditional Sigma vs |returns|", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			grid()
		}
	}
	invisible()
}

.plot.dccfit.3 = function(x, series, ...)
{
	ops = list(...)
	n = length( series )
	cnames = x@model$modeldata$asset.names
	T = x@model$modeldata$T
	p = x@model$maxgarchOrder
	xDates = x@model$modeldata$index[1:T]
	nc = ( (n^2 - n)/2 )
	idx = .lowertri.index( n )
	H = rcov(x)
	if( nc > 16 ){
		scr = floor( nc/16 )
		z	= nc - 16 * floor( nc/16 )
		start	= dev.next( which = dev.cur() )
		xn	= 0
		for( j in 1:scr ){
			dev.new( start + j )
			par( mfrow = c(4, 4) )
			for(i in 1:16){
				plot(  xts(H[ series[idx[i + xn, 1]], series[idx[i + xn, 2]], 1:T], xDates), 
						col = colors()[35], 
						main = paste(cnames[series[idx[i + xn, 1]]], "-", cnames[series[idx[i + xn, 2]]], sep = ""),  
						ylab = "", xlab = "", minor.ticks = FALSE, auto.grid=FALSE)
				mtext(paste("rmgarch  : DCC model fit"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
			}
			title( paste( "DCC Conditional Covariance\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			grid()
			xn = xn + 16
		}
		if( z != 0 ){
			dev.new( dev.next( which = dev.cur() ) + 1 )
			par( mfrow = c(4, 4) )
			for(i in 1:z){
				plot( xts(H[ series[idx[i + xn, 1]], series[idx[i + xn, 2]], 1:T], xDates), 
						col = colors()[35], 
						main = paste(cnames[series[idx[i + xn, 1]]], "-", cnames[series[idx[i + xn, 2]]], sep = ""),  
						ylab = "", xlab = "", minor.ticks = FALSE, auto.grid=FALSE)
				mtext(paste("rmgarch  : DCC model fit"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
			}
			title( paste( "DCC Conditional Covariance\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			grid()
		}
	} else{
		d = .divisortable( nc )
		par( mfrow= c(d[1], d[2]) )
		for(i in 1:nc){
			plot(  xts(H[ series[idx[i, 1]], series[idx[i, 2]], 1:T], xDates),
					col = colors()[35], 
					main = paste(cnames[series[idx[i, 1]]], "-", cnames[series[idx[i, 2]]], sep = ""),  
					ylab = "", xlab = "", minor.ticks = FALSE, auto.grid=FALSE)
			mtext(paste("rmgarch  : DCC model fit"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
			title( paste( "DCC Conditional Covariance", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			grid()
		}
	}
	invisible()
}

.plot.dccfit.4 = function(x, series, ...)
{
	ops = list(...)
	n = length( series )
	cnames = x@model$modeldata$asset.names
	T = x@model$modeldata$T
	p = x@model$maxgarchOrder
	xDates = x@model$modeldata$index[1:T]
	nc = ( (n^2 - n)/2 )
	idx = .lowertri.index( n )
	H = rcor(x)
	if( nc > 16 ){
		scr = floor( nc/16 )
		z	= nc - 16 * floor( nc/16 )
		start	= dev.next( which = dev.cur() )
		xn	= 0
		for( j in 1:scr ){
			dev.new( start + j )
			par( mfrow = c(4, 4) )
			for(i in 1:16){
				plot(  xts(H[ series[idx[i + xn, 1]], series[idx[i + xn, 2]], 1:T], xDates), 
						col = colors()[35], 
						main = paste(cnames[series[idx[i + xn, 1]]], "-", cnames[series[idx[i + xn, 2]]], sep = ""),  
						ylab = "", xlab = "", minor.ticks = FALSE, auto.grid=FALSE)
				mtext(paste("rmgarch  : DCC model fit"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
			}
			title( paste( "DCC Conditional Correlation\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			grid()
			xn = xn + 16
		}
		if( z != 0 ){
			dev.new( dev.next( which = dev.cur() ) + 1 )
			par( mfrow = c(4, 4) )
			for(i in 1:z){
				plot( xts(H[ series[idx[i + xn, 1]], series[idx[i + xn, 2]], 1:T], xDates), 
						col = colors()[35], 
						main = paste(cnames[series[idx[i + xn, 1]]], "-", cnames[series[idx[i + xn, 2]]], sep = ""),  
						ylab = "", xlab = "", minor.ticks = FALSE, auto.grid=FALSE)
				mtext(paste("rmgarch  : DCC model fit"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
			}
			title( paste( "DCC Conditional Correlation\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			grid()
		}
	} else{
		d = .divisortable( nc )
		par( mfrow= c(d[1], d[2]) )
		for(i in 1:nc){
			plot(  xts(H[ series[idx[i, 1]], series[idx[i, 2]], 1:T], xDates),
					col = colors()[35], 
					main = paste(cnames[series[idx[i, 1]]], "-", cnames[series[idx[i, 2]]], sep = ""),  
					ylab = "", xlab = "", minor.ticks = FALSE, auto.grid=FALSE)
			mtext(paste("rmgarch  : DCC model fit"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
			title( paste( "DCC Conditional Correlation", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			grid()
		}
	}
	invisible()
}

.plot.dccfit.5 = function(x, series, ...)
{
	p = x@model$maxgarchOrder
	T = x@model$modeldata$T
	xDat = x@model$modeldata$data[1:T, ]
	xDates = x@model$modeldata$index[1:T]	
	port = as.numeric( apply(xDat, 1, "mean") )
	n = dim(xDat)[2]
	w = rep(1/n , n)
	dport = wmargin(distribution = x@model$modeldesc$distribution, weights = w, Sigma = rcov(x), mean = fitted(x),
			shape = rshape(x), skew = rskew(x))
	dm = switch(x@model$modeldesc$distribution,
			mvnorm = "norm",
			mvlaplace = "ged",
			mvt = "std")
	q01 = rugarch::qdist(distribution = dm, p = 0.01, mu = dport[,1], sigma = dport[,2], lambda = -0.5, 
						skew = dport[,3], shape = dport[,4])
	q99 = rugarch::qdist(distribution = dm, p = 0.99, mu = dport[,1], sigma = dport[,2], lambda = -0.5, 
						skew = dport[,3], shape = dport[,4])
	plot(xts(port, xDates), col = "steelblue", ylab = "", xlab = "", 
			main = "EW Portfolio with with 1% VaR Limits \n(based on Conditional Weighted Density)", 
			cex.main = 0.8, minor.ticks = FALSE, auto.grid=FALSE)
	lines(xts(q01, xDates), col = "tomato1")
	lines(xts(q99,xDates), col = "green")
	mtext(paste("rmgarch  : DCC model fit"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
	abline(h = 0, col = "grey", lty = 3)
	grid()
	invisible()
}


.plotdccforecast = function(x, which = "ask", series = c(1, 2), ...)
{
	n = dim(x@model$umodel$modelinc)[2]
	xseries = unique( as.integer( series ) )
	xseries = xseries[xseries > 0]
	if( any(xseries > n) ) stop( "\nrmgarch-->error: series index out of bounds in DCCfit plot.")
	if( length(xseries) < 2 ) stop( "\nrmgarch-->error: series length must be 2 or more.")
	
	series = xseries
	choices = c(
			"Conditional Mean Forecast  (vs realized  returns)",
			"Conditional Sigma Forecast (vs realized |returns|)",
			"Conditional Covariance Forecast",
			"Conditional Correlation Forecast",
			"EW Portfolio Plot with forecast conditional density VaR limits")
	.interdccforecastPlot(x, choices = choices, plotFUN = paste(".plot.dccforecast", 1:5, sep = "."), which = which, series = series, ...)
	# Return Value:
	invisible(x)
}


.interdccforecastPlot = function(x, choices, plotFUN, which, series = c(1, 2), ...)
{
	if (is.numeric(which)) {
		if(which>length(choices)) stop("Not a valid choice. Plots choices are 1-5.\n",call. = FALSE)
		FUN = match.fun(plotFUN[which])
		FUN(x, series = series, ...)
	}
	if(is.character(which))
	{
		if( which!="ask" ) stop("Not a valid choice.\n",call. = FALSE)
		.multdccforecastPlot(x, choices, series = series, ...)
	}
	
	invisible(x)
}

.multdccforecastPlot = function(x, choices, series = c(1, 2), ...)
{
	
	pick = 1
	while (pick > 0) {
		pick = menu (
				choices = paste(" ", choices),
				title = "\nMake a plot selection (or 0 to exit):")
		switch (pick,
				.plot.dccforecast.1(x, series, ...),  .plot.dccforecast.2(x, series, ...),  
				.plot.dccforecast.3(x, series, ...),  .plot.dccforecast.4(x, series, ...),
				.plot.dccforecast.5(x, ...))
	}
}

.plot.dccforecast.1 = function(x, series = c(1, 2), ...)
{
	ops = list(...)
	n = length( series )
	n.ahead = x@model$n.ahead
	n.roll	= x@model$n.roll
	if( n.ahead == 1 && n.roll > 0 ){
		n.start = x@model$modeldata$n.start
		n.assets = dim(x@model$modeldata$data)[2]
		cnames = x@model$modeldata$asset.names
		T = x@model$modeldata$T
		xDat = x@model$modeldata$data[(T-49):(T+n.roll),]
		#nn = length(x@uforecast@forecast[[1]]@forecast$series)
		yDat = x@model$modeldata$index[(T-49):(T+n.roll)]
		indx = c(-49:0, 1:(n.roll))
	 	yforc = rbind(
			 matrix(NA, ncol = n.assets, nrow = 50), 
			 t(fitted(x)[1,,]))
	 # discard last roll in case it is outside of realized observations
	 yforc = yforc[1:(50+n.roll),]
		if( n > 16 ){
			scr = floor( n/16 )
			z	= n - 16 * floor( n/16 )
			start	= dev.next( which = dev.cur() )
			xn	= 0
			for( j in 1:scr ){
				dev.new( start + j )
				par( mfrow = c(4, 4) )
				for(i in 1:16){
					tmp = yforc[,series[i + xn]]
					plot( xts(xDat[,series[i + xn]], yDat), col = colors()[16], main = cnames[i + xn],
							ylab = "Series", xlab = "Time",
							minor.ticks=FALSE, auto.grid=FALSE)
					lines( xts(tmp, yDat), col = colors()[134])					
					abline( v = yDat[51], lty = 3, col = "steelblue")
					mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
					legend("topleft", c("Returns", "Forecast"), col=c(colors()[16], colors()[134]), lty=c(1,1),
							bty="n", cex = 0.7)
					grid()
				}
				title( paste( "DCC Series Rolling Forecast\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
				xn = xn + 16
			}
			if( z != 0 ){
				dev.new( dev.next( which = dev.cur() ) + 1 )
				par( mfrow = c(4, 4) )
				for(i in 1:z){
					tmp = yforc[,series[i + xn]]
					plot( xts(xDat[,series[i + xn]], yDat), col = colors()[16], 
							main = cnames[i + xn], ylab = "Series", xlab = "Time", 
							minor.ticks=FALSE, auto.grid=FALSE)
					lines( xts(tmp, yDat), col = colors()[134])					
					abline( v = yDat[51], lty = 3, col = "steelblue")
					mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
					legend("topleft", c("Returns", "Forecast"), col=c(colors()[16], colors()[134]), lty=c(1,1),
							bty="n", cex = 0.7)
					grid()
				}
				title( paste( "DCC Series Rolling Forecast\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			}
		} else{
			d = .divisortable( n )
			par( mfrow= c(d[1], d[2]) )
			for(i in 1:n){
				tmp = yforc[,series[i]]
				plot( xts(xDat[,series[i]], yDat), col = colors()[16], main = cnames[i],
						ylab = "Series", xlab = "Time", minor.ticks=FALSE, auto.grid=FALSE)
				lines( xts(tmp, yDat), col = colors()[134])					
				abline( v = yDat[51], lty = 3, col = "steelblue")
				mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
				legend("topleft", c("Returns", "Forecast"), col=c(colors()[16], colors()[134]), lty=c(1,1),
						bty="n", cex = 0.7)
				grid()
			}
			title( paste( "DCC Series Rolling Forecast", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			
		}
	} else if( n.ahead > 1 && n.roll == 0 ){
		period   = x@model$modeldata$period
		n.start  = x@model$modeldata$n.start
		n.assets = NCOL(x@model$modeldata$data)
		cnames   = x@model$modeldata$asset.names
		T = x@model$modeldata$T
		# To provide for a realistic plot:
		# If some of n.ahead lies within the out.sample period (if at all),
		# then find those dates which it lies within.
		# Otherwise, and for all n.ahead > available data/dates, generate
		# the n.ahead dates using the 'seq' and periodicity of the data.
		if(n.start>0 && n.ahead<n.start){
			xDat = x@model$modeldata$data[(T-49):(T+n.ahead),]
			yDat = x@model$modeldata$index[(T-49):(T+n.ahead)]
		} else if(n.start>0 && n.ahead>n.start){
			xDat = rbind(
					x@model$modeldata$data[(T-49):(T+n.start),],
					matrix(NA, ncol = n.assets, nrow = n.ahead - n.start))
			yDat = c(
					x@model$modeldata$index[(T-49):(T+n.start)],
					seq(x@model$modeldata$index[T+n.start], by = period, length.out=n.ahead - n.start +1)[-1])
		} else{
			xDat = rbind(
					x@model$modeldata$data[(T-49):T,],
					matrix(NA, ncol = n.assets, nrow = n.ahead))
			yDat = c(
					x@model$modeldata$index[(T-49):(T)],
					seq(x@model$modeldata$index[T], by = period, length.out=n.ahead+1)[-1])
		}
		
		yforc = rbind(
				matrix(NA, ncol = n.assets, nrow = 50), 
				x@mforecast$mu[,,1])		
		if( n > 16 ){
			scr = floor( n/16 )
			z	= n - 16 * floor( n/16 )
			start	= dev.next( which = dev.cur() )
			xn	= 0
			for( j in 1:scr ){
				dev.new( start + j )
				par( mfrow = c(4, 4) )
				for(i in 1:16){
					tmp = yforc[,series[i + xn]]
					plot( xts(xDat[,series[i + xn]], yDat), col = colors()[16], 
							main = cnames[i + xn], ylab = "Series", xlab = "Time",
							minor.ticks=FALSE, auto.grid=FALSE)
					lines( xts(tmp, yDat), col = colors()[134])				
					abline( v = yDat[51], lty = 3, col = "steelblue")
					mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
					legend("topleft", c("Returns", "Forecast"), col=c(colors()[16], colors()[134]), lty=c(1,1),
							bty="n", cex = 0.7)
					grid()
				}
				title( paste( "DCC Series Unconditional Forecast\n(page...", j, ")", sep = ""), outer = TRUE, line = -2, cex = 0.75)
				xn = xn + 16
			}
			if( z != 0 ){
				dev.new( dev.next( which = dev.cur() ) + 1 )
				par( mfrow = c(4, 4) )
				for(i in 1:z){
					tmp = yforc[,series[i + xn]]
					plot( xts(xDat[,series[i + xn]], yDat), col = colors()[16], main = cnames[i + xn],
							ylab = "Series", xlab = "Time", minor.ticks=FALSE, auto.grid=FALSE)
					lines( xts(tmp, yDat), col = colors()[134])					
					abline( v = yDat[51], lty = 3, col = "steelblue")
					mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
					legend("topleft", c("Returns", "Forecast"), col=c(colors()[16], colors()[134]), lty=c(1,1),
							bty="n", cex = 0.7)
					grid()
				}
				title( paste( "DCC Series Unconditional Forecast\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			}
		} else{
			d = .divisortable( n )
			par( mfrow= c(d[1], d[2]) )
			for(i in 1:n){
				tmp = yforc[,series[i]]
				plot( xts(xDat[,series[i]], yDat),  col = colors()[16], main = cnames[i],
						ylab = "Series", xlab = "Time", minor.ticks=FALSE, auto.grid=FALSE)
				lines( xts(tmp, yDat), col = colors()[134])		
				abline( v = yDat[51], lty = 3, col = "steelblue")
				mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
				legend("topleft", c("Returns", "Forecast"), col=c(colors()[16], colors()[134]), lty=c(1,1),
						bty="n", cex = 0.7)
				grid()
			}
			title( paste( "DCC Series Unconditional Forecast", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
		}
	} else if( n.ahead > 1 && n.roll  > 0 ){
		cat("\nNo plot available for mixed unconditional and rolling forecasts.")
	}
	invisible()
}

.plot.dccforecast.2 = function(x, series = c(1, 2), ...)
{
	ops = list(...)
	n = length( series )
	n.ahead = x@model$n.ahead
	n.roll	= x@model$n.roll
	if( n.ahead == 1 && n.roll > 0 ){
		n.start = x@model$modeldata$n.start
		n.assets = dim(x@model$modeldata$data)[2]
		cnames = x@model$modeldata$asset.names
		T = x@model$modeldata$T
		xDat = abs(x@model$modeldata$data[(T-49):(T+n.roll),])
		yDat = x@model$modeldata$index[(T-49):(T+n.roll)]
		indx = c(-49:0, 1:(n.roll))
		yforc = rbind(
				matrix(NA, ncol = n.assets, nrow = 50), 
				t(sigma(x)[1,,]))
		# discard last roll in case it is outside of realized observations
		yforc = yforc[1:(50+n.roll),]
		if( n > 16 ){
			scr = floor( n/16 )
			z	= n - 16 * floor( n/16 )
			start	= dev.next( which = dev.cur() )
			xn	= 0
			for( j in 1:scr ){
				dev.new( start + j )
				par( mfrow = c(4, 4) )
				for(i in 1:16){
					tmp = yforc[,series[i + xn]]
					plot( xts(xDat[,series[i + xn]], yDat), col = colors()[16], 
							main = cnames[i + xn], ylab = "Sigma", xlab = "Time", 
							minor.ticks=FALSE, auto.grid=FALSE)
					lines( xts(tmp, yDat), col = colors()[134])					
					abline( v = yDat[51], lty = 3, col = "steelblue")
					mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
					legend("topleft", c("|Returns|", "Forecast"), col=c(colors()[16], colors()[134]), lty=c(1,1),
							bty="n", cex = 0.7)
					grid()
				}
				title( paste( "DCC Sigma Rolling Forecast\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
				xn = xn + 16
			}
			if( z != 0 ){
				dev.new( dev.next( which = dev.cur() ) + 1 )
				par( mfrow = c(4, 4) )
				for(i in 1:z){
					tmp = yforc[,series[i + xn]]
					plot( xts(xDat[,series[i + xn]], yDat), col = colors()[16], 
							main = cnames[i + xn], ylab = "Sigma", xlab = "Time", 
							minor.ticks=FALSE, auto.grid=FALSE)
					lines( xts(tmp, yDat), col = colors()[134])					
					abline( v = yDat[51], lty = 3, col = "steelblue")
					mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
					legend("topleft", c("|Returns|", "Forecast"), col=c(colors()[16], colors()[134]), lty=c(1,1),
							bty="n", cex = 0.7)
					grid()
				}
				title( paste( "DCC Sigma Rolling Forecast\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			}
		} else{
			d = .divisortable( n )
			par( mfrow= c(d[1], d[2]) )
			for(i in 1:n){
				tmp = yforc[,series[i]]
				plot( xts(xDat[,series[i]], yDat), col = colors()[16], main = cnames[i],
						ylab = "Sigma", xlab = "Time", minor.ticks=FALSE, auto.grid=FALSE)
				lines( xts(tmp, yDat), col = colors()[134])					
				abline( v = yDat[51], lty = 3, col = "steelblue")
				mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
				legend("topleft", c("|Returns|", "Forecast"), col=c(colors()[16], colors()[134]), lty=c(1,1),
						bty="n", cex = 0.7)
				grid()				
			}
			title( paste( "DCC Sigma Rolling Forecast", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
		}
	} else if( n.ahead > 1 && n.roll == 0 ){
		period   = x@model$modeldata$period
		n.start  = x@model$modeldata$n.start
		n.assets = NCOL(x@model$modeldata$data)
		cnames   = x@model$modeldata$asset.names
		T = x@model$modeldata$T
		# To provide for a realistic plot:
		# If some of n.ahead lies within the out.sample period (if at all),
		# then find those dates which it lies within.
		# Otherwise, and for all n.ahead > available data/dates, generate
		# the n.ahead dates using the 'seq' and periodicity of the data.
		if(n.start>0 && n.ahead<n.start){
			xDat = x@model$modeldata$data[(T-49):(T+n.ahead),]
			yDat = x@model$modeldata$index[(T-49):(T+n.ahead)]
		} else if(n.start>0 && n.ahead>n.start){
			xDat = rbind(
					x@model$modeldata$data[(T-49):(T+n.start),],
					matrix(NA, ncol = n.assets, nrow = n.ahead - n.start))
			yDat = c(
					x@model$modeldata$index[(T-49):(T+n.start)],
					seq(x@model$modeldata$index[T+n.start], by = period, length.out=n.ahead - n.start +1)[-1])
		} else{
			xDat = rbind(
					x@model$modeldata$data[(T-49):T,],
					matrix(NA, ncol = n.assets, nrow = n.ahead))
			yDat = c(
					x@model$modeldata$index[(T-49):(T)],
					seq(x@model$modeldata$index[T], by = period, length.out=n.ahead+1)[-1])
		}
		xDat = abs(xDat)
		yforc = rbind(
				matrix(NA, ncol = n.assets, nrow = 50), 
				sigma(x)[,,1])
		if( n > 16 ){
			scr = floor( n/16 )
			z	= n - 16 * floor( n/16 )
			start	= dev.next( which = dev.cur() )
			xn	= 0
			for( j in 1:scr ){
				dev.new( start + j )
				par( mfrow = c(4, 4) )
				for(i in 1:16){
					tmp = yforc[,series[i + xn]]
					plot( xts(xDat[,series[i + xn]], yDat), col = colors()[16], 
							main = cnames[i + xn], ylab = "Sigma", xlab = "Time",
							minor.ticks=FALSE, auto.grid=FALSE)
					lines( xts(tmp, yDat), col = colors()[134])
					abline( v = yDat[51], lty = 3, col = "steelblue")
					mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
					legend("topleft", c("|Actual|", "Forecast"), col=c(colors()[16], colors()[134]), lty=c(1,1),
							bty="n", cex = 0.7)
					grid()
				}
				title( paste( "DCC Sigma Unconditional Forecast\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
				xn = xn + 16
			}
			if( z != 0 ){
				dev.new( dev.next( which = dev.cur() ) + 1 )
				par( mfrow = c(4, 4) )
				for(i in 1:z){
					tmp = yforc[,series[i + xn]]
					plot( xts(xDat[,series[i + xn]], yDat), col = colors()[16], 
							main = cnames[i + xn], ylab = "Sigma", xlab = "Time",
							minor.ticks=FALSE, auto.grid=FALSE)
					lines( xts(tmp, yDat), col = colors()[134])		
					abline( v = yDat[51], lty = 3, col = "steelblue")
					mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
					legend("topleft", c("|Returns|", "Forecast"), col=c(colors()[16], colors()[134]), lty=c(1,1),
							bty="n", cex = 0.7)
					grid()
				}
				title( paste( "DCC Sigma Unconditional Forecast\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			}
		} else{
			d = .divisortable( n )
			par( mfrow= c(d[1], d[2]) )
			for(i in 1:n){
				tmp = yforc[,series[i]]
				plot( xts(xDat[,series[i]], yDat), col = colors()[16], 
						main = cnames[i], ylab = "Sigma", xlab = "Time",
						minor.ticks=FALSE, auto.grid=FALSE)
				lines( xts(tmp, yDat), col = colors()[134])
				abline( v = yDat[51], lty = 3, col = "steelblue")
				mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
				legend("topleft", c("|Returns|", "Forecast"), col=c(colors()[16], colors()[134]), lty=c(1,1),
						bty="n", cex = 0.7)
				grid()
			}
			title( paste( "DCC Sigma Unconditional Forecast", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
		}
	} else if( n.ahead > 1 && n.roll  > 0 ){
		cat("\nNo plot available for mixed unconditional and rolling forecasts.")
	}
	invisible()
}



.plot.dccforecast.3 = function(x, series = c(1, 2), ...)
{
	ops = list(...)
	n = length( series )
	n.ahead = x@model$n.ahead
	n.roll	= x@model$n.roll
	cnames = x@model$modeldata$asset.names
	T = x@model$modeldata$T
	p = x@model$maxgarchOrder
	index = x@model$modeldata$index[(p+1):T]
	nc = ( (n^2 - n)/2 )
	idx = .lowertri.index( n )
	H = x@model$H
	zn = 50
	
	if( n.ahead == 1 && n.roll >= 0 ){
		n.start = x@model$modeldata$n.start
		n.assets = NCOL(x@model$modeldata$data)
		cnames = x@model$modeldata$asset.names
		T = x@model$modeldata$T
		yDat = x@model$modeldata$index[(T-zn+1):(T+n.roll)]
		if( nc > 16 ){
			scr = floor( nc/16 )
			z	= nc - 16 * floor( nc/16 )
			start	= dev.next( which = dev.cur() )
			xn	= 0
			for( j in 1:scr ){
				dev.new( start + j )
				par( mfrow = c(4, 4) )
				for(i in 1:16){
					tmpf = sapply(x@mforecast$H, FUN = function(x) x[series[idx[i + xn, 1]], series[idx[i + xn, 2]], 1])
					tmpf = c(rep(NA, zn),  tmpf)
					tmpf = tmpf[1:(zn+n.roll)]
					tmp =  c(H[series[idx[i + xn, 1]], series[idx[i + xn, 2]], (T-49):(T)], rep(NA, n.roll))
					plot(  xts(tmp, yDat), col = colors()[16], 
							main = paste(cnames[series[idx[i + xn, 1]]], "-", cnames[series[idx[i + xn, 2]]], sep = ""), 
							ylab = "covariance", xlab = "Time", ylim = c(min(c(tmp, tmpf), na.rm=TRUE), max(c(tmp, tmpf), na.rm=TRUE)))
					lines(  xts(tmpf, yDat), col = colors()[134])
					abline( v = yDat[zn+1], lty = 3, col = colors()[442])
					mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
					legend("topleft", c("Forecast"), col=colors()[134], lty=1, bty="n", cex = 0.7)
					grid()
				}
				title( paste( "DCC Covariance Rolling Forecast\n(page...", j, ")", sep = ""), outer = TRUE, line = -2, cex = 0.75)
				xn = xn + 16
			}
			if( z != 0 ){
				dev.new( dev.next( which = dev.cur() ) + 1 )
				par( mfrow = c(4, 4) )
				for(i in 1:z){
					tmpf = sapply(x@mforecast$H, FUN = function(x) x[series[idx[i + xn, 1]], series[idx[i + xn, 2]], 1])
					tmpf = c(rep(NA, zn),  tmpf)
					tmpf = tmpf[1:(zn+n.roll)]
					tmp =  c(H[series[idx[i + xn, 1]], series[idx[i + xn, 2]], (T-49):(T)], rep(NA, n.roll))
					
					plot(  xts(tmp, yDat), col = colors()[16], 
							main = paste(cnames[series[idx[i + xn, 1]]], "-", cnames[series[idx[i + xn, 2]]], sep = ""), 
							ylab = "covariance", xlab = "Time", ylim = c(min(c(tmp, tmpf), na.rm=TRUE), max(c(tmp, tmpf), na.rm=TRUE)),
							minor.ticks = FALSE, auto.grid = FALSE)
					lines(  xts(tmpf, yDat), col = colors()[134])
					abline( v = yDat[zn+1], lty = 3, col = colors()[442])
					mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
					legend("topleft", c("Forecast"), col=colors()[134], lty=1, bty="n", cex = 0.7)
					grid()
				}
				title( paste( "DCC Covariance Rolling Forecast\n(page...", j, ")", sep = ""), outer = TRUE, line = -2, cex = 0.75)
			}
		} else{
			d = .divisortable( nc )
			par( mfrow= c(d[1], d[2]) )
			for(i in 1:nc){
				tmpf = sapply(x@mforecast$H, FUN = function(x) x[series[idx[i, 1]], series[idx[i, 2]], 1])
				tmpf = c(rep(NA, zn),  tmpf)
				# (roll+1) value excluded since we do not have the date for it readily available
				# and checking whether (n.roll+1)>out.sample amd adjusting is too much for one point!
				tmpf = tmpf[1:(zn+n.roll)]
				tmp =  c(H[series[idx[i, 1]], series[idx[i, 2]], (T-49):(T)], rep(NA, n.roll))
				
				plot(  xts(tmp, yDat), col = colors()[16], 
						main = paste(cnames[series[idx[i, 1]]], "-", cnames[series[idx[i, 2]]], sep = ""), 
						ylab = "covariance", xlab = "Time", ylim = c(min(c(tmp, tmpf), na.rm=TRUE), max(c(tmp, tmpf), na.rm=TRUE)),
						minor.ticks = FALSE, auto.grid = FALSE)
				lines(  xts(tmpf, yDat), col = colors()[134])
				abline( v = yDat[zn+1], lty = 3, col = colors()[442])
				mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
				legend("topleft", c("Forecast"), col=colors()[134], lty=1, bty="n", cex = 0.7)
				grid()
			}
			title( paste( "DCC Covariance Rolling Forecast", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			
		}
	} else if( n.ahead > 1 && n.roll == 0 ){
		period   = x@model$modeldata$period
		n.start  = x@model$modeldata$n.start
		n.assets = NCOL(x@model$modeldata$data)
		cnames   = x@model$modeldata$asset.names
		T = x@model$modeldata$T
		if(n.start>0 && n.ahead<n.start){
			yDat = x@model$modeldata$index[(T-49):(T+n.ahead)]
		} else if(n.start>0 && n.ahead>n.start){
			yDat = c(
					x@model$modeldata$index[(T-49):(T+n.start)],
					seq(x@model$modeldata$index[T+n.start], by = period, length.out=n.ahead - n.start +1)[-1])
		} else{
			yDat = c(
					x@model$modeldata$index[(T-49):(T)],
					seq(x@model$modeldata$index[T], by = period, length.out=n.ahead+1)[-1])
		}
		if( nc > 16 ){
			scr = floor( nc/16 )
			z	= nc - 16 * floor( nc/16 )
			start	= dev.next( which = dev.cur() )
			xn	= 0
			for( j in 1:scr ){
				dev.new( start + j )
				par( mfrow = c(4, 4) )
				for(i in 1:16){
					tmpf = apply(x@mforecast$H[[1]], 3, FUN = function(x) x[series[idx[i + xn, 1]], series[idx[i + xn, 2]]] )
					tmpf = c(rep(NA, zn),  tmpf)
					tmp =  c(H[series[idx[i + xn, 1]], series[idx[i + xn, 2]], (T-49):T], rep(NA, n.ahead) )
					
					plot(  xts(tmp, yDat), col = colors()[16],  
							main = paste(cnames[series[idx[i + xn, 1]]], "-", cnames[series[idx[i + xn, 2]]], sep = ""), 
							ylab = "covariance", xlab = "Time", ylim = c(min(c(tmp, tmpf), na.rm=TRUE), max(c(tmp, tmpf), na.rm=TRUE)),
							minor.ticks = FALSE, auto.grid = FALSE)
					lines(  xts(tmpf, yDat), col = colors()[134])
					abline( v = yDat[zn+1], lty = 3, col = colors()[442])
					mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
					legend("topleft", c("Forecast"), col=colors()[134], lty=1, bty="n", cex = 0.7)
					grid()
				}
				title( paste( "DCC Unconditional Covariance Forecast\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
				xn = xn + 16
			}
			if( z != 0 ){
				dev.new( dev.next( which = dev.cur() ) + 1 )
				par( mfrow = c(4, 4) )
				for(i in 1:z){
					tmpf = apply(x@mforecast$H[[1]], 3, FUN = function(x) x[series[idx[i + xn, 1]], series[idx[i + xn, 2]]] )
					tmpf = c(rep(NA, zn),  tmpf)
					tmp =  c(H[series[idx[i + xn, 1]], series[idx[i + xn, 2]], (T-49):T], rep(NA, n.ahead) )
					
					plot(  xts(tmp, yDat), col = colors()[16], 
							main = paste(cnames[series[idx[i + xn, 1]]], "-", cnames[series[idx[i + xn, 2]]], sep = ""), 
							ylab = "covariance", xlab = "Time", ylim = c(min(c(tmp, tmpf), na.rm=TRUE), max(c(tmp, tmpf), na.rm=TRUE)),
							minor.ticks = FALSE, auto.grid = FALSE)
					lines(  xts(tmpf, yDat), col = colors()[134])
					abline( v = yDat[zn+1], lty = 3, col = colors()[442])
					mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
					legend("topleft", c("Forecast"), col=colors()[134], lty=1, bty="n", cex = 0.7)
					grid()
				}
				title( paste( "DCC Unconditional Covariance Forecast\n(page...", j+1, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			}
		} else{
			d = .divisortable( nc )
			par( mfrow= c(d[1], d[2]) )
			for(i in 1:nc){
				tmpf = apply(x@mforecast$H[[1]], 3, FUN = function(x) x[series[idx[i, 1]], series[idx[i, 2]]] )
				tmpf = c(rep(NA, zn),  tmpf)
				tmp =  c(H[series[idx[i, 1]], series[idx[i, 2]], (T-49):T], rep(NA, n.ahead) )
				plot(  xts(tmp, yDat), col = colors()[16], 
						main = paste(cnames[series[idx[i, 1]]], "-", cnames[series[idx[i, 2]]], sep = ""), 
						ylab = "covariance", xlab = "Time", 
						ylim = c(min(c(tmp, tmpf), na.rm=TRUE), max(c(tmp, tmpf), na.rm=TRUE)), 
						minor.ticks = FALSE, auto.grid = FALSE)
				lines(  xts(tmpf, yDat), col = colors()[134])
				abline( v = yDat[zn+1], lty = 3, col = colors()[442])
				mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
				legend("topleft", c("Forecast"), col=colors()[134], lty=1, bty="n", cex = 0.7)
				grid()
			}
			title( paste( "DCC Unconditional Covariance Forecast", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			
		}
	} else if( n.ahead > 1 && n.roll  > 0 ){
		cat("\nNo plot available for mixed unconditional and rolling forecasts.")
	}
	invisible()
}

.plot.dccforecast.4 = function(x, series = c(1, 2), ...)
{
	ops = list(...)
	n = length( series )
	n.ahead = x@model$n.ahead
	n.roll	= x@model$n.roll
	cnames = x@model$modeldata$asset.names
	T = x@model$modeldata$T
	p = x@model$maxgarchOrder
	index = x@model$modeldata$index[(p+1):T]
	nc = ( (n^2 - n)/2 )
	idx = .lowertri.index( n )
	H = x@model$H
	R = .Call("Cov2Cor", H, dim(H), PACKAGE = "rmgarch")
	zn = 50
	
	if( n.ahead == 1 && n.roll >= 0 ){
		n.start = x@model$modeldata$n.start
		n.assets = NCOL(x@model$modeldata$data)
		cnames = x@model$modeldata$asset.names
		T = x@model$modeldata$T
		yDat = x@model$modeldata$index[(T-zn+1):(T+n.roll)]
		fR = lapply(x@mforecast$H, function(x) cov2cor(x[,,1]))
		if( nc > 16 ){
			scr = floor( nc/16 )
			z	= nc - 16 * floor( nc/16 )
			start	= dev.next( which = dev.cur() )
			xn	= 0
			for( j in 1:scr ){
				dev.new( start + j )
				par( mfrow = c(4, 4) )
				for(i in 1:16){
					tmpf = sapply(fR, FUN = function(x) x[series[idx[i + xn, 1]], series[idx[i + xn, 2]]])
					tmpf = c(rep(NA, zn),  tmpf)
					tmpf = tmpf[1:(zn+n.roll)]
					tmp =  c(R[series[idx[i + xn, 1]], series[idx[i + xn, 2]], (T-49):(T)], rep(NA, n.roll))
					plot(  xts(tmp, yDat), col = colors()[16], 
							main = paste(cnames[series[idx[i + xn, 1]]], "-", cnames[series[idx[i + xn, 2]]], sep = ""), 
							ylab = "correlation", xlab = "Time", ylim = c(min(c(tmp, tmpf), na.rm=TRUE), max(c(tmp, tmpf), na.rm=TRUE)))
					lines(  xts(tmpf, yDat), col = colors()[134])
					abline( v = yDat[zn+1], lty = 3, col = colors()[442])
					mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
					legend("topleft", c("Forecast"), col=colors()[134], lty=1, bty="n", cex = 0.7)
					grid()
				}
				title( paste( "DCC Correlation Rolling Forecast\n(page...", j, ")", sep = ""), outer = TRUE, line = -2, cex = 0.75)
				xn = xn + 16
			}
			if( z != 0 ){
				dev.new( dev.next( which = dev.cur() ) + 1 )
				par( mfrow = c(4, 4) )
				for(i in 1:z){
					tmpf = sapply(fR, FUN = function(x) x[series[idx[i + xn, 1]], series[idx[i + xn, 2]]])
					tmpf = c(rep(NA, zn),  tmpf)
					tmpf = tmpf[1:(zn+n.roll)]
					tmp =  c(R[series[idx[i + xn, 1]], series[idx[i + xn, 2]], (T-49):(T)], rep(NA, n.roll))
					plot(  xts(tmp, yDat), col = colors()[16], 
							main = paste(cnames[series[idx[i + xn, 1]]], "-", cnames[series[idx[i + xn, 2]]], sep = ""), 
							ylab = "correlation", xlab = "Time", ylim = c(min(c(tmp, tmpf), na.rm=TRUE), max(c(tmp, tmpf), na.rm=TRUE)),
							minor.ticks = FALSE, auto.grid = FALSE)
					lines(  xts(tmpf, yDat), col = colors()[134])
					abline( v = yDat[zn+1], lty = 3, col = colors()[442])
					mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
					legend("topleft", c("Forecast"), col=colors()[134], lty=1, bty="n", cex = 0.7)
					grid()
				}
				title( paste( "DCC Correlation Rolling Forecast\n(page...", j, ")", sep = ""), outer = TRUE, line = -2, cex = 0.75)
			}
		} else{
			d = .divisortable( nc )
			par( mfrow= c(d[1], d[2]) )
			for(i in 1:nc){
				tmpf = sapply(fR, FUN = function(x) x[series[idx[i, 1]], series[idx[i, 2]]])
				tmpf = c(rep(NA, zn),  tmpf)
				# (roll+1) value excluded since we do not have the date for it readily available
				# and checking whether (n.roll+1)>out.sample amd adjusting is too much for one point!
				tmpf = tmpf[1:(zn+n.roll)]
				tmp =  c(R[series[idx[i, 1]], series[idx[i, 2]], (T-49):(T)], rep(NA, n.roll))
				
				plot(  xts(tmp, yDat), col = colors()[16], 
						main = paste(cnames[series[idx[i, 1]]], "-", cnames[series[idx[i, 2]]], sep = ""), 
						ylab = "correlation", xlab = "Time", ylim = c(min(c(tmp, tmpf), na.rm=TRUE), max(c(tmp, tmpf), na.rm=TRUE)),
						minor.ticks = FALSE, auto.grid = FALSE)
				lines(  xts(tmpf, yDat), col = colors()[134])
				abline( v = yDat[zn+1], lty = 3, col = colors()[442])
				mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
				legend("topleft", c("Forecast"), col=colors()[134], lty=1, bty="n", cex = 0.7)
				grid()
			}
			title( paste( "DCC Correlation Rolling Forecast", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			
		}
	} else if( n.ahead > 1 && n.roll == 0 ){
		period   = x@model$modeldata$period
		n.start  = x@model$modeldata$n.start
		n.assets = NCOL(x@model$modeldata$data)
		cnames   = x@model$modeldata$asset.names
		T = x@model$modeldata$T
		if(n.start>0 && n.ahead<n.start){
			yDat = x@model$modeldata$index[(T-49):(T+n.ahead)]
		} else if(n.start>0 && n.ahead>n.start){
			yDat = c(
					x@model$modeldata$index[(T-49):(T+n.start)],
					seq(x@model$modeldata$index[T+n.start], by = period, length.out=n.ahead - n.start +1)[-1])
		} else{
			yDat = c(
					x@model$modeldata$index[(T-49):(T)],
					seq(x@model$modeldata$index[T], by = period, length.out=n.ahead+1)[-1])
		}
		fR = array(apply(x@mforecast$H[[1]], 3, function(x) cov2cor(x)), dim = dim(x@mforecast$H[[1]]))
		if( nc > 16 ){
			scr = floor( nc/16 )
			z	= nc - 16 * floor( nc/16 )
			start	= dev.next( which = dev.cur() )
			xn	= 0
			for( j in 1:scr ){
				dev.new( start + j )
				par( mfrow = c(4, 4) )
				for(i in 1:16){
					tmpf = apply(fR, 3, FUN = function(x) x[series[idx[i + xn, 1]], series[idx[i + xn, 2]]] )
					tmpf = c(rep(NA, zn),  tmpf)
					tmp =  c(R[series[idx[i + xn, 1]], series[idx[i + xn, 2]], (T-49):T], rep(NA, n.ahead) )
					plot(  xts(tmp, yDat), col = colors()[16],  
							main = paste(cnames[series[idx[i + xn, 1]]], "-", cnames[series[idx[i + xn, 2]]], sep = ""), 
							ylab = "correlation", xlab = "Time", ylim = c(min(c(tmp, tmpf), na.rm=TRUE), max(c(tmp, tmpf), na.rm=TRUE)),
							minor.ticks = FALSE, auto.grid = FALSE)
					lines(  xts(tmpf, yDat), col = colors()[134])
					abline( v = yDat[zn+1], lty = 3, col = colors()[442])
					mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
					legend("topleft", c("Forecast"), col=colors()[134], lty=1, bty="n", cex = 0.7)
					grid()
				}
				title( paste( "DCC Unconditional Correlation Forecast\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
				xn = xn + 16
			}
			if( z != 0 ){
				dev.new( dev.next( which = dev.cur() ) + 1 )
				par( mfrow = c(4, 4) )
				for(i in 1:z){
					tmpf = apply(fR, 3, FUN = function(x) x[series[idx[i + xn, 1]], series[idx[i + xn, 2]]] )
					tmpf = c(rep(NA, zn),  tmpf)
					tmp =  c(R[series[idx[i + xn, 1]], series[idx[i + xn, 2]], (T-49):T], rep(NA, n.ahead) )
					plot(  xts(tmp, yDat), col = colors()[16], 
							main = paste(cnames[series[idx[i + xn, 1]]], "-", cnames[series[idx[i + xn, 2]]], sep = ""), 
							ylab = "correlation", xlab = "Time", ylim = c(min(c(tmp, tmpf), na.rm=TRUE), max(c(tmp, tmpf), na.rm=TRUE)),
							minor.ticks = FALSE, auto.grid = FALSE)
					lines(  xts(tmpf, yDat), col = colors()[134])
					abline( v = yDat[zn+1], lty = 3, col = colors()[442])
					mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
					legend("topleft", c("Forecast"), col=colors()[134], lty=1, bty="n", cex = 0.7)
					grid()
				}
				title( paste( "DCC Unconditional Correlation Forecast\n(page...", j+1, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			}
		} else{
			d = .divisortable( nc )
			par( mfrow= c(d[1], d[2]) )
			for(i in 1:nc){
				tmpf = apply(fR, 3, FUN = function(x) x[series[idx[i, 1]], series[idx[i, 2]]] )
				tmpf = c(rep(NA, zn),  tmpf)
				tmp =  c(R[series[idx[i, 1]], series[idx[i, 2]], (T-49):T], rep(NA, n.ahead) )
				plot(  xts(tmp, yDat), col = colors()[16], 
						main = paste(cnames[series[idx[i, 1]]], "-", cnames[series[idx[i, 2]]], sep = ""), 
						ylab = "correlation", xlab = "Time", 
						ylim = c(min(c(tmp, tmpf), na.rm=TRUE), max(c(tmp, tmpf), na.rm=TRUE)), 
						minor.ticks = FALSE, auto.grid = FALSE)
				lines(  xts(tmpf, yDat), col = colors()[134])
				abline( v = yDat[zn+1], lty = 3, col = colors()[442])
				mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
				legend("topleft", c("Forecast"), col=colors()[134], lty=1, bty="n", cex = 0.7)
				grid()
			}
			title( paste( "DCC Unconditional Correlation Forecast", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			
		}
	} else if( n.ahead > 1 && n.roll  > 0 ){
		cat("\nNo plot available for mixed unconditional and rolling forecasts.")
	}
	invisible()
}

.plot.dccforecast.5 = function(x, series, ...)
{
	ops = list(...)
	n.ahead = x@model$n.ahead
	n.roll	= x@model$n.roll
	n.start = x@model$modeldata$n.start
	cnames = x@model$modeldata$asset.names
	T = x@model$modeldata$T
	p = x@model$maxgarchOrder
	m = dim(x@model$umodel$modelinc)[2]
	
	if( n.ahead == 1 && n.roll >= 0 ){
		Hf = array(sapply(rcov(x), function(x) x), dim=c(m,m,n.roll+1))
		Hf = Hf[,,1:n.roll]
		Mf = t(x@mforecast$mu[,,1:n.roll])
		w = rep(1/m, m)
		dist = x@model$modeldesc$distribution
		dm = switch(dist,
				mvnorm = "norm",
				mvlaplace = "ged",
				mvt = "std")
		# previous 50 observations:
		Hx = last(x@model$H[,,1:T], 50)
		mx = 50
		Mx = tail(x@model$mu[1:T,], mx)
		dportx = wmargin(distribution = dist, weights = w, Sigma = Hx, mean = Mx, shape = rshape(x), skew = rskew(x))
		q025x = as.numeric(rugarch::qdist(distribution = dm, p = 0.025, mu = dportx[,1], 
						sigma = dportx[,2], lambda = -0.5, skew = dportx[,3], shape = dportx[,4]))
		q975x = as.numeric( rugarch::qdist(distribution = dm, p = 0.975, mu = dportx[,1], 
						sigma = dportx[,2], lambda = -0.5, skew = dportx[,3], shape = dportx[,4]))
		dportf = wmargin(distribution = dist, weights = w, Sigma = Hf, mean = Mf, shape = rshape(x), skew = rskew(x))
		
		q025f = as.numeric( rugarch::qdist(distribution = dm, p = 0.025, mu = dportf[,1], 
						sigma = dportf[,2], lambda = -0.5, skew = dportf[,3], shape = dportf[,4]))
		
		q975f = as.numeric( rugarch::qdist(distribution = dm, p = 0.975, mu = dportf[,1], 
						sigma = dportf[,2], lambda = -0.5, skew = dportf[,3], shape = dportf[,4]))
		yDat = x@model$modeldata$index[(T-mx+1):(T+n.roll)]
		port = apply(x@model$modeldata$data[(T-mx+1):(T+n.roll), ], 1, "mean")
		
		plot(xts(port, yDat), col = colors()[220], ylab = "", xlab = "", 
				main = "EW Forecast Portfolio with Rolling 2.5% VaR Limits \n(based on Conditional Weighted Density)", 
				cex.main = 0.8, ylim = c(min(c(q025f, q025x)), max(c(q975f, q975x)) ),
				minor.ticks=FALSE, auto.grid=FALSE)
		lines(xts(c(rep(NA, mx), dportf[1:n.roll,1]), yDat), col = colors()[134])
		lines(xts(c(dportx[,1], rep(NA, n.roll)), yDat), col = colors()[132])		
		
		lines(xts(c(q025x, rep(NA, n.roll)), yDat), col = colors()[132], lty = 3)
		lines(xts(c(q975x, rep(NA, n.roll)), yDat), col = colors()[132], lty = 3)
		
		lines(xts(c(rep(NA, mx), q025f), yDat), col = colors()[134], lty = 3)
		lines(xts(c(rep(NA, mx), q975f), yDat), col = colors()[134], lty = 3)
				mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
		abline( v = yDat[mx + 1], lty = 3, col = colors()[442])
		legend("topleft", legend = c("Realized", "Fitted", "Forecast"), 
				col = c(colors()[220], colors()[132], colors()[134]), bty = "n",
				fill = c(colors()[220], colors()[132], colors()[134]))
		grid()
	} else if( n.ahead > 1 && n.roll == 0 ){
		period   = x@model$modeldata$period
		n.start  = x@model$modeldata$n.start
		n.assets = NCOL(x@model$modeldata$data)
		cnames   = x@model$modeldata$asset.names
		T = x@model$modeldata$T
		Hf = x@mforecast$H[[1]]
		Mf = x@mforecast$mu[,,1]
		w = rep(1/m, m)
		dist = x@model$modeldesc$distribution
		dm = switch(dist,
				mvnorm = "norm",
				mvlaplace = "ged",
				mvt = "std")
		# previous 50 observations:
		Hx = last(x@model$H[,,1:T], 50)
		mx = 50
		Mx = tail(x@model$mu[1:T,], mx)
		
		if(n.start>0 && n.ahead<n.start){
			port = apply(x@model$modeldata$data[(T-mx+1):(T+n.ahead),], 1, "mean")			
			yDat = x@model$modeldata$index[(T-mx+1):(T+n.ahead)]
		} else if(n.start>0 && n.ahead>n.start){
			port = c(
					apply(x@model$modeldata$data[(T-mx+1):(T+n.start),], 1, "mean"),
					rep(NA, n.ahead - n.start))
			yDat = c(
					x@model$modeldata$index[(T-mx+1):(T+n.start)],
					seq(x@model$modeldata$index[T+n.start], by = period, length.out=n.ahead - n.start +1)[-1])
		} else{
			port = c(
					apply(x@model$modeldata$data[(T-mx+1):T,], 1, "mean"),
					rep(NA, n.ahead))
			yDat = c(
					x@model$modeldata$index[(T-mx+1):(T)],
					seq(x@model$modeldata$index[T], by = period, length.out=n.ahead+1)[-1])
		}
		port = as.numeric(port)
		dportx = wmargin(distribution = dist, weights = w, Sigma = Hx, mean = Mx, shape = rshape(x), skew = rskew(x))
		q025x = as.numeric(rugarch::qdist(distribution = dm, p = 0.025, mu = dportx[,1], 
						sigma = dportx[,2], lambda = -0.5, skew = dportx[,3], shape = dportx[,4]))
		q975x = as.numeric( rugarch::qdist(distribution = dm, p = 0.975, mu = dportx[,1], 
						sigma = dportx[,2], lambda = -0.5, skew = dportx[,3], shape = dportx[,4]))
		dportf = wmargin(distribution = dist, weights = w, Sigma = Hf, mean = Mf, shape = rshape(x), skew = rskew(x))
		
		q025f = as.numeric( rugarch::qdist(distribution = dm, p = 0.025, mu = dportf[,1], 
						sigma = dportf[,2], lambda = -0.5, skew = dportf[,3], shape = dportf[,4]))
		
		q975f = as.numeric( rugarch::qdist(distribution = dm, p = 0.975, mu = dportf[,1], 
						sigma = dportf[,2], lambda = -0.5, skew = dportf[,3], shape = dportf[,4]))	
	
		plot(xts(port, yDat), col = colors()[220], ylab = "", xlab = "", 
				main = "EW Forecast Portfolio with Unconditional 2.5% VaR Limits \n(based on Conditional Weighted Density)", 
				cex.main = 0.8, ylim = c(min(c(port, q025f, q025x), na.rm = TRUE), 
						max(c(port, q975f, q975x), na.rm = TRUE)), minor.ticks=FALSE,
						auto.grid=FALSE)
		lines(xts(c(rep(NA, mx), dportf[,1]), yDat),  col = colors()[134])
		lines(xts(c(dportx[,1], rep(NA, n.ahead)), yDat), col = colors()[132])
		lines(xts(c(q025x, rep(NA, n.ahead)), yDat),  col = colors()[132], lty = 3)
		lines(xts(c(q975x, rep(NA, n.ahead)), yDat), col = colors()[132],  lty = 3)
		lines(xts(c(rep(NA, mx), q025f), yDat), col = colors()[134], lty = 3)
		lines(xts(c(rep(NA, mx), q975f), yDat), col = colors()[134], lty = 3)
		mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
		abline( v = yDat[mx + 1], lty = 3, col = colors()[442])
		legend("topleft", legend = c("Realized", "Fitted", "Forecast"), 
				col = c(colors()[220], colors()[132], colors()[134]), bty = "n",
				fill = c(colors()[220], colors()[132], colors()[134]))
		grid()
	} else if( n.ahead > 1 && n.roll  > 0 ){
		cat("\nNo plot available for mixed unconditional and rolling forecasts.")
	}
	invisible()
}


.plotdccroll = function(x, which = "ask", series = c(1, 2), ...)
{
	n = dim(x@model$umodel$modelinc)[2]
	xseries = unique( as.integer( series ) )
	xseries = xseries[xseries > 0]
	if( any(xseries > n) ) stop( "\nrmgarch-->error: series index out of bounds in DCCfit plot.")
	if( length(xseries) < 2 ) stop( "\nrmgarch-->error: series length must be 2 or more.")
	
	series = xseries
	choices = c(
			"Conditional Mean Forecast  (vs realized  returns)",
			"Conditional Sigma Forecast (vs realized |returns|)",
			"Conditional Covariance Forecast",
			"Conditional Correlation Forecast",
			"EW Portfolio Plot with forecast conditional density VaR limits")
	.interdccrollPlot(x, choices = choices, plotFUN = paste(".plot.dccroll", 1:5, sep = "."), which = which, series = series, ...)
	# Return Value:
	invisible(x)
}


.interdccrollPlot = function(x, choices, plotFUN, which, series = c(1, 2), ...)
{
	if (is.numeric(which)) {
		if(which>length(choices)) stop("Not a valid choice. Plots choices are 1-5.\n",call. = FALSE)
		FUN = match.fun(plotFUN[which])
		FUN(x, series = series, ...)
	}
	if(is.character(which))
	{
		if( which!="ask" ) stop("Not a valid choice.\n",call. = FALSE)
		.multdccrollPlot(x, choices, series = series, ...)
	}
	
	invisible(x)
}

.multdccrollPlot = function(x, choices, series = c(1, 2), ...)
{
	
	pick = 1
	while (pick > 0) {
		pick = menu (
				choices = paste(" ", choices),
				title = "\nMake a plot selection (or 0 to exit):")
		switch (pick,
				.plot.dccroll.1(x, series, ...),  .plot.dccroll.2(x, series, ...),  
				.plot.dccroll.3(x, series, ...),  .plot.dccroll.4(x, series, ...),
				.plot.dccroll.5(x, ...))
	}
}

.plot.dccroll.1 = function(x, series = c(1, 2), ...)
{
	ops = list(...)
	n  = length( series )
	s  = x@model$out.sample
	fL = sum(s)
	n.start = x@model$n.start
	n.assets = NCOL(x@model$data)
	cnames = x@model$modeldata$asset.names
	T = x@model$n.start
	xDat = x@model$data[(T-49):(T+fL),]
	#nn = length(x@uforecast@forecast[[1]]@forecast$series)
	yDat = x@model$index[(T-49):(T+fL)]
	Y = xts(xDat, yDat)
	X = fitted(x)
	if( n > 16 ){
		scr = floor( n/16 )
		z	= n - 16 * floor( n/16 )
		start	= dev.next( which = dev.cur() )
		xn	= 0
		for( j in 1:scr ){
			dev.new( start + j )
			par( mfrow = c(4, 4) )
			for(i in 1:16){
				tmp = X[,series[i + xn]]
				plot( Y[,series[i + xn]], col = colors()[16], main = cnames[i + xn], 
								ylab = "Series", xlab = "Time", minor.ticks=FALSE, auto.grid=FALSE)
				lines( tmp, col = colors()[134])
				abline( v = yDat[51], lty = 3, col = "steelblue")
				mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
				legend("topleft", c("Returns", "Forecast"), col=c(colors()[16], colors()[134]), lty=c(1,1),
						bty="n", cex = 0.7)
				grid()
			}
			title( paste( "DCC Series Rolling Forecast\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			xn = xn + 16
		}
		if( z != 0 ){
			dev.new( dev.next( which = dev.cur() ) + 1 )
			par( mfrow = c(4, 4) )
			for(i in 1:z){
				tmp = X[,series[i + xn]]
				plot( Y[,series[i + xn]], col = colors()[16], main = cnames[i + xn], 
						ylab = "Series", xlab = "Time", 
						minor.ticks=FALSE, auto.grid=FALSE)
				lines( tmp, col = colors()[134])					
				abline( v = yDat[51], lty = 3, col = "steelblue")
				mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
				legend("topleft", c("Returns", "Forecast"), col=c(colors()[16], colors()[134]), lty=c(1,1),
						bty="n", cex = 0.7)
				grid()
			}
			title( paste( "DCC Series Rolling Forecast\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
		}
	} else{
		d = .divisortable( n )
		par( mfrow= c(d[1], d[2]) )
		for(i in 1:n){
			tmp = X[,series[i]]
			plot( Y[,series[i]], col = colors()[16], main = cnames[i],
					ylab = "Series", xlab = "Time", minor.ticks=FALSE, 
					auto.grid=FALSE)
			lines( tmp, col = colors()[134])					
			abline( v = yDat[51], lty = 3, col = "steelblue")
			mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
			legend("topleft", c("Returns", "Forecast"), col=c(colors()[16], colors()[134]), lty=c(1,1),
					bty="n", cex = 0.7)
			grid()
		}
		title( paste( "DCC Series Rolling Forecast", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
	}
	invisible()
}

.plot.dccroll.2 = function(x, series = c(1, 2), ...)
{
	ops = list(...)
	n  = length( series )
	s  = x@model$out.sample
	fL = sum(s)
	n.start = x@model$n.start
	n.assets = NCOL(x@model$data)
	cnames = x@model$modeldata$asset.names
	T = x@model$n.start
	xDat = x@model$data[(T-49):(T+fL),]
	#nn = length(x@uforecast@forecast[[1]]@forecast$series)
	yDat = x@model$index[(T-49):(T+fL)]
	Y = xts(abs(xDat), yDat)
	X = sigma(x)
	if( n > 16 ){
		scr = floor( n/16 )
		z	= n - 16 * floor( n/16 )
		start	= dev.next( which = dev.cur() )
		xn	= 0
		for( j in 1:scr ){
			dev.new( start + j )
			par( mfrow = c(4, 4) )
			for(i in 1:16){
				tmp = X[,series[i + xn]]
				plot( Y[,series[i + xn]], col = colors()[16], main = cnames[i + xn], 
						ylab = "Sigma", xlab = "Time", minor.ticks=FALSE, auto.grid=FALSE)
				lines( tmp, col = colors()[134])
				abline( v = yDat[51], lty = 3, col = "steelblue")
				mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
				legend("topleft", c("|returns|", "Forecast"), col=c(colors()[16], colors()[134]), lty=c(1,1),
						bty="n", cex = 0.7)
				grid()
			}
			title( paste( "DCC Series Rolling Forecast\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			xn = xn + 16
		}
		if( z != 0 ){
			dev.new( dev.next( which = dev.cur() ) + 1 )
			par( mfrow = c(4, 4) )
			for(i in 1:z){
				tmp = X[,series[i + xn]]
				plot( Y[,series[i + xn]], col = colors()[16], main = cnames[i + xn], 
						ylab = "Sigma", xlab = "Time", 
						minor.ticks=FALSE, auto.grid=FALSE)
				lines( tmp, col = colors()[134])					
				abline( v = yDat[51], lty = 3, col = "steelblue")
				mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
				legend("topleft", c("|returns|", "Forecast"), col=c(colors()[16], colors()[134]), lty=c(1,1),
						bty="n", cex = 0.7)
				grid()
			}
			title( paste( "DCC Series Rolling Forecast\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
		}
	} else{
		d = .divisortable( n )
		par( mfrow= c(d[1], d[2]) )
		for(i in 1:n){
			tmp = X[,series[i]]
			plot( Y[,series[i]], col = colors()[16], main = cnames[i],
					ylab = "Sigma", xlab = "Time", minor.ticks=FALSE, 
					auto.grid=FALSE)
			lines( tmp, col = colors()[134])					
			abline( v = yDat[51], lty = 3, col = "steelblue")
			mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
			legend("topleft", c("|returns|", "Forecast"), col=c(colors()[16], colors()[134]), lty=c(1,1),
					bty="n", cex = 0.7)
			grid()
		}
		title( paste( "DCC Series Rolling Forecast", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
	}
	invisible()
}



.plot.dccroll.3 = function(x, series = c(1, 2), ...)
{
	ops = list(...)
	n  = length( series )
	s  = x@model$out.sample
	fL = sum(s)
	cnames = x@model$modeldata$asset.names
	idx = .lowertri.index( n )
	n = NROW(idx)
	X = rcov(x)
	D = as.POSIXct(dimnames(X)[[3]])
	if( n > 16 ){
		scr = floor( n/16 )
		z	= n - 16 * floor( n/16 )
		start	= dev.next( which = dev.cur() )
		xn	= 0
		for( j in 1:scr ){
			dev.new( start + j )
			par( mfrow = c(4, 4) )
			for(i in 1:16){
				plot( xts(X[series[idx[i+xn,1]],series[idx[i+xn,2]],], D), col = 1, 
						main = paste(cnames[series[idx[i + xn, 1]]], "-", cnames[series[idx[i + xn, 2]]], sep = ""), 
						ylab = "Covariance", xlab = "Time", minor.ticks=FALSE, auto.grid=FALSE)
				lines(xts(X[series[idx[i+xn,1]],series[idx[i+xn,2]],], D), col = colors()[134])
				mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
				grid()
			}
			title( paste( "DCC Series Rolling Forecast\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			xn = xn + 16
		}
		if( z != 0 ){
			dev.new( dev.next( which = dev.cur() ) + 1 )
			par( mfrow = c(4, 4) )
			for(i in 1:z){
				plot( xts(X[series[idx[i+xn,1]],series[idx[i+xn,2]],], D), col = 1, 
						main = paste(cnames[series[idx[i + xn, 1]]], "-", cnames[series[idx[i + xn, 2]]], sep = ""), 
						ylab = "Covariance", xlab = "Time", minor.ticks=FALSE, auto.grid=FALSE)
				lines(xts(X[series[idx[i+xn,1]],series[idx[i+xn,2]],], D), col = colors()[134])
				mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
				grid()
			}
			title( paste( "DCC Series Rolling Forecast\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
		}
	} else{
		d = .divisortable( n )
		par( mfrow= c(d[1], d[2]) )
		for(i in 1:n){
			plot( xts(X[series[idx[i,1]],series[idx[i,2]],], D), col =1, 
					main = paste(cnames[series[idx[i, 1]]], "-", cnames[series[idx[i, 2]]], sep = ""), 
					ylab = "Covariance", xlab = "Time", minor.ticks=FALSE, auto.grid=FALSE)
			lines(xts(X[series[idx[i,1]],series[idx[i,2]],], D), col = colors()[134])
			mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
			grid()
		}
		title( paste( "DCC Series Rolling Forecast", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
	}
	invisible()
}

.plot.dccroll.4 = function(x, series = c(1, 2), ...)
{
	ops = list(...)
	n  = length( series )
	s  = x@model$out.sample
	fL = sum(s)
	cnames = x@model$modeldata$asset.names
	idx = .lowertri.index( n )
	n = NROW(idx)
	X = rcor(x)
	D = as.POSIXct(dimnames(X)[[3]])
	if( n > 16 ){
		scr = floor( n/16 )
		z	= n - 16 * floor( n/16 )
		start	= dev.next( which = dev.cur() )
		xn	= 0
		for( j in 1:scr ){
			dev.new( start + j )
			par( mfrow = c(4, 4) )
			for(i in 1:16){
				plot( xts(X[series[idx[i+xn,1]],series[idx[i+xn,2]],], D), col = 1, 
						main = paste(cnames[series[idx[i + xn, 1]]], "-", cnames[series[idx[i + xn, 2]]], sep = ""), 
						ylab = "Correlation", xlab = "Time", minor.ticks=FALSE, auto.grid=FALSE)
				lines(xts(X[series[idx[i+xn,1]],series[idx[i+xn,2]],], D), col = colors()[134])
				mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
				grid()
			}
			title( paste( "DCC Series Rolling Forecast\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
			xn = xn + 16
		}
		if( z != 0 ){
			dev.new( dev.next( which = dev.cur() ) + 1 )
			par( mfrow = c(4, 4) )
			for(i in 1:z){
				plot( xts(X[series[idx[i+xn,1]],series[idx[i+xn,2]],], D), col = 1, 
						main = paste(cnames[series[idx[i + xn, 1]]], "-", cnames[series[idx[i + xn, 2]]], sep = ""), 
						ylab = "Correlation", xlab = "Time", minor.ticks=FALSE, auto.grid=FALSE)
				lines(xts(X[series[idx[i+xn,1]],series[idx[i+xn,2]],], D), col = colors()[134])
				mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
				grid()
			}
			title( paste( "DCC Series Rolling Forecast\n(page...", j, ")", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
		}
	} else{
		d = .divisortable( n )
		par( mfrow= c(d[1], d[2]) )
		for(i in 1:n){
			plot( xts(X[series[idx[i,1]],series[idx[i,2]],], D), col =1, 
					main = paste(cnames[series[idx[i, 1]]], "-", cnames[series[idx[i, 2]]], sep = ""), 
					ylab = "Correlation", xlab = "Time", minor.ticks=FALSE, auto.grid=FALSE)
			lines(xts(X[series[idx[i,1]],series[idx[i,2]],], D), col = colors()[134])
			mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
			grid()
		}
		title( paste( "DCC Series Rolling Forecast", sep = ""), outer = TRUE, line = -1.5, cex = 0.75)
	}
	invisible()
}

.plot.dccroll.5 = function(x, series, ...)
{
	ops = list(...)
	s  = x@model$out.sample
	fL = sum(s)
	n.start = x@model$n.start
	n.assets = NCOL(x@model$data)
	cnames = x@model$modeldata$asset.names
	T = x@model$n.start
	m = NCOL(x@model$data)
	xDat = apply(x@model$data[(T-49):(T+fL),], 1, "mean")
	#nn = length(x@uforecast@forecast[[1]]@forecast$series)
	yDat = x@model$index[(T-49):(T+fL)]
	D = x@model$index[(T+1):(T+fL)]
	Y = xts(xDat, yDat)
	T = x@model$n.start
	Hf = rcov(x)
	Mf = fitted(x)
	w = rep(1/m, m)
	dist = x@model$modeldesc$distribution
	dm = switch(dist,
			mvnorm = "norm",
			mvlaplace = "ged",
			mvt = "std")
	# previous 50 observations:
	dportf = wmargin(distribution = dist, weights = w, Sigma = Hf, mean = Mf, shape = NA, skew = NA)
	skewx = shapex = NULL
	shapef = rshape(x)
	skewf = rskew(x)
	for(i in 1:length(s)){
		shapex = c(shapex, rep(shapef[i], s[i]))
		skewx = c(skewx, rep(skewf[i], s[i]))
	}
	if(dm=="mvt"){
		dportf[,"shape"] = shapex
	}
	q025f = as.numeric( rugarch::qdist(distribution = dm, p = 0.025, mu = dportf[,1], 
					sigma = dportf[,2], lambda = -0.5, skew = dportf[,3], shape = dportf[,4]))
	
	q975f = as.numeric( rugarch::qdist(distribution = dm, p = 0.975, mu = dportf[,1], 
					sigma = dportf[,2], lambda = -0.5, skew = dportf[,3], shape = dportf[,4]))
	portf = apply(Mf, 1, "mean")
	
	plot(Y, col = colors()[220], ylab = "", xlab = "", 
			main = "EW Forecast Portfolio with Rolling 2.5% VaR Limits \n(based on Conditional Weighted Density)", 
			cex.main = 0.8, ylim = c(min(c(as.numeric(Y[,1]), q025f)), max(c(as.numeric(Y[,1]), q975f)) ),
			minor.ticks=FALSE, auto.grid=FALSE)
	lines(xts(portf, D), col = "tomato1", lwd=2)	
	lines(xts(q025f, D), col = "brown", lty = 3, lwd=1)
	lines(xts(q975f, D), col = "brown", lty = 3, lwd=1)
	mtext(paste("rmgarch  : DCC model forecast"), side = 4, adj = 0, padj=0, col = "gray", cex = 0.4)
	abline( v = index(Y)[50], lty = 3, lwd=2, col = colors()[442])
	legend("topleft", legend = c("Realized", "Forecast"), 
			col = c(colors()[220], "tomato1"), bty = "n",
			lty = c(1,1), lwd=c(1,2))
	grid()
	invisible()
}


.lowertri.index = function(n, diag = FALSE)
{
	x = matrix(1, n, n)
	indx = which(x == lower.tri(x, diag = diag), arr.ind = TRUE)
	return( indx )
}

.divisortable = function(n) 
{
	z = matrix(c(1, 1, 1, 2, 2, 1, 3, 2, 2, 4, 2, 2, 5, 2, 3, 
					6, 2, 3, 7, 2, 4, 8, 2, 4, 9, 3, 3, 10, 3, 4, 11, 3, 
					4, 12, 3, 4, 13, 4, 4, 14, 4, 4, 15, 4, 4, 16, 4, 4, 
					17, 4, 5, 18, 4, 5, 19, 4, 5, 20, 4, 5), ncol = 3, byrow = TRUE)
	d = which(n == z[, 1])
	return(z[d, 2:3])
}
