#################################################################################
##
##   R package rugarch by Alexios Ghalanos Copyright (C) 2008-2015.
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
.plotgarchfit = function(x, which="ask",...)
{
	#old.par <- par(no.readonly = TRUE)
	#on.exit(par(old.par))
	choices = c(
			"Series with 2 Conditional SD Superimposed",
			"Series with 1% VaR Limits",
			"Conditional SD (vs |returns|)",
			"ACF of Observations",
			"ACF of Squared Observations",
			"ACF of Absolute Observations",
			"Cross Correlation",
			"Empirical Density of Standardized Residuals",
			"QQ-Plot of Standardized Residuals",
			"ACF of Standardized Residuals",
			"ACF of Squared Standardized Residuals",
			"News-Impact Curve")
	.intergarchfitPlot(x, choices = choices, plotFUN = paste(".plot.garchfit", 1:12, sep = "."), which = which, ...)
	# Return Value:
	invisible(x)
}

.intergarchfitPlot = function(x, choices, plotFUN, which, ...)
{
	old.par <- par(no.readonly = TRUE)
	if (is.numeric(which)) {
		if(which>length(choices)) stop("Not a valid choice. Plots choices are 1-12.\n",call. = FALSE)
		FUN = match.fun(plotFUN[which])
		FUN(x)
	}
	if(is.character(which))
	{
		if(which!="all" && which!="ask") stop("Not a valid choice.\n",call. = FALSE)
		if (which[1] == "all") {
			on.exit(par(old.par))
			#Which = rep(TRUE, times = length(choices))
			par(mfrow=c(3,4))
			for(i in 1:12){
				FUN = match.fun(plotFUN[i])
				FUN(x)
			}
		} else{
			.multgarchfitPlot(x, choices, plotFUN, ...)
		}
	}
	
	invisible(x)
}

.multgarchfitPlot = function(x, choices, ...)
{
	
	pick = 1
	while (pick > 0) {
		pick = menu (
				choices = paste(" ", choices),
				title = "\nMake a plot selection (or 0 to exit):")
		switch (pick,
				.plot.garchfit.1(x),  .plot.garchfit.2(x),  .plot.garchfit.3(x),
				.plot.garchfit.4(x),  .plot.garchfit.5(x),  .plot.garchfit.6(x),
				.plot.garchfit.7(x),  .plot.garchfit.8(x),  .plot.garchfit.9(x),
				.plot.garchfit.10(x), .plot.garchfit.11(x), .plot.garchfit.12(x))
	}
}

# Series with 2 Conditional SD Superimposed
.plot.garchfit.1 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xdates  = x@model$modeldata$index[insample]
	xseries = x@model$modeldata$data[insample]
	xsigma  = x@fit$sigma
	ci = 2
	plot(xdates, xseries, type = "l",  col = "steelblue", ylab = "Returns", xlab="Time", 
			main = "Series with 2 Conditional SD Superimposed", cex.main = 0.8)
	lines(xdates, +ci* xsigma , col = "tomato1")
	lines(xdates, -ci* xsigma, col = "tomato1")
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	abline(h = 0, col = "grey", lty = 3)
	grid()
}

# Series with 1% VaR Limits
.plot.garchfit.2 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]
	xdates  = x@model$modeldata$index[insample]
	xsigma 	= x@fit$sigma
	distribution = x@model$modeldesc$distribution
	xcmu = fitted(x)
	idx = x@model$pidx
	pars  = x@fit$ipars[,1]
	skew  = pars[idx["skew",1]]
	shape = pars[idx["shape",1]]
	if(distribution == "ghst") ghlambda = -shape/2 else ghlambda = pars[idx["ghlambda",1]]
	z1 	= 0.01
	z2 	= 0.99
	cat("\nplease wait...calculating quantiles...\n")
	q01 	= fitted(x) + sigma(x)* qdist(distribution, z1, 0, 1, lambda = ghlambda, skew, shape)
	q99 	= fitted(x) + sigma(x)* qdist(distribution, z2, 0, 1, lambda = ghlambda, skew, shape)
	plot(xdates, xseries, type = "l", col = "steelblue", ylab = "Returns", xlab="Time", 
			main = "Series with with 1% VaR Limits", cex.main = 0.8)
	lines(xdates, q01, col = "tomato1")
	lines(xdates, q99, col = "green")
	mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	abline(h = 0, col = "grey", lty = 3)
	grid()
}

# Conditional SD
.plot.garchfit.3 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = abs(x@model$modeldata$data[insample])
	xdates  = x@model$modeldata$index[insample]
	xsigma 	= x@fit$sigma
	plot(xdates, xseries, type = "l", col = "lightgrey", ylab = "Volatility", 
			xlab="Time", main = "Conditional SD (vs |returns|)", cex.main = 0.8)
	lines(xdates, xsigma, col = "steelblue")
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	abline(h = 0, col = "grey", lty = 3)
	grid()
}

# ACF of observations
.plot.garchfit.4 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]
	lag.max = as.integer(10*log10(T))
	acfx 	= acf(xseries, lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
	clx 	= vector(mode = "character", length = lag.max)
	clx[which(as.numeric(acfx$acf)[-1]>=0)] = "steelblue"
	clx[which(as.numeric(acfx$acf)[-1]<0)] = "orange"
	barplot(height = as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim = 1.2*ylim, col = clx,
			ylab = "ACF", xlab="lag", main = "ACF of Observations", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black", lty = 1)
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# ACF of squared observations
.plot.garchfit.5 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]
	lag.max = as.integer(10*log10(T))
	acfx 	= acf(xseries^2, lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
	clx 	= vector(mode="character",length=lag.max)
	clx[which(as.numeric(acfx$acf)[-1]>=0)] = "steelblue"
	clx[which(as.numeric(acfx$acf)[-1]<0)] = "orange"
	barplot(height = as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim = 1.2*ylim, col = clx,
			ylab = "ACF", xlab="lag", main = "ACF of Squared Observations", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black", lty = 1)
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# ACF of absolute observations
.plot.garchfit.6 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]
	lag.max = as.integer(10*log10(T))
	acfx 	= acf(abs(xseries), lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
	clx 	= vector(mode="character",length=lag.max)
	clx[which(as.numeric(acfx$acf)[-1]>=0)] = "steelblue"
	clx[which(as.numeric(acfx$acf)[-1]<0)] = "orange"
	barplot(height=as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim = 1.2*ylim, col = clx,
			ylab = "ACF", xlab = "lag", main = "ACF of Absolute Observations", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black", lty = 1)
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# Cross Correlations
.plot.garchfit.7 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]	
	lag.max = as.integer(10*log10(T))
	ccfx	= ccf(xseries^2, xseries, lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(ccfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(ccfx$acf)))
	clx 	= vector(mode="character",length=(2*lag.max)+1)
	clx[which(as.numeric(ccfx$acf)>=0)]="steelblue"
	clx[which(as.numeric(ccfx$acf)<0)]="orange"
	barplot(height=as.numeric(ccfx$acf), names.arg=as.numeric(ccfx$lag),ylim=1.2*ylim,col=clx,
			ylab = "ACF", xlab="lag", main = "Cross-Correlations of \nSquared vs Actual Observations", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black")
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# Standardized Residuals Empirical Denisty
.plot.garchfit.8 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	zseries = as.numeric(residuals(x, standardize=TRUE))
	distribution = x@model$modeldesc$distribution
	idx = x@model$pidx
	pars  = x@fit$ipars[,1]
	skew  = pars[idx["skew",1]]
	shape = pars[idx["shape",1]]
	if(distribution == "ghst") ghlambda = -shape/2 else ghlambda = pars[idx["ghlambda",1]]	
	xmean 	= mean(zseries)
	xmedian = median(zseries)
	xsd 	= sd(zseries)
	xlim 	= c(min(zseries), max(zseries))
	result 	= hist(x = zseries, col = "grey", border = "white",
			breaks = "Scott", main = "Empirical Density of Standardized Residuals", xlim = xlim, ylim = c(0,0.6),
			probability = TRUE, ylab="Probability", cex.main = 0.8, ...)
	box()
	#s 	= seq(xlim[1], xlim[2], length = 201)
	s = result$breaks
	y	= ddist(distribution, s, lambda = ghlambda, skew = skew, shape = shape)
	lines(s, dnorm(s, 0, 1), lwd = 2, col = "blue")
	lines(s, y, lwd = 2, col = "orange")
	abline(v = xmean, lwd = 2, col = "red")
	abline(v = xmedian, lwd = 2, col = "darkgreen")
	Text = paste("Median: ", round(xmedian, 2), "| Mean: ",signif(xmean, 3))
	mtext(Text, side = 3, adj = 0, col = "darkgrey", cex = 0.7)
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
	}
	abline(h = 0, col = "grey")
	lg.txt = c("normal Density", paste(distribution," (0,1) Fitted Density",sep=""))
	legend("topleft", legend = lg.txt, col = c("blue","orange"), pch = 1, cex = 0.8, bty = "n")
	grid()
}

# QQ-Plot of Standardized Residuals
.plot.garchfit.9 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	zseries = as.numeric(residuals(x, standardize=TRUE))
	distribution = x@model$modeldesc$distribution
	idx = x@model$pidx
	pars  = x@fit$ipars[,1]
	skew  = pars[idx["skew",1]]
	shape = pars[idx["shape",1]]
	if(distribution == "ghst") ghlambda = -shape/2 else ghlambda = pars[idx["ghlambda",1]]
	.qqDist(y = zseries, dist = distribution, lambda = ghlambda, skew = skew, shape = shape)
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
	}
}

### Continue Here

# ACF of standardized residuals
.plot.garchfit.10 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	zseries = as.numeric(residuals(x, standardize=TRUE))
	zseries[is.na(zseries)] = 0
	n 		= length(zseries)
	lag.max = as.integer(10*log10(n))
	acfx 	= acf(zseries, lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
	clx 	= vector(mode = "character", length = lag.max)
	clx[which(as.numeric(acfx$acf)[-1]>=0)] = "steelblue"
	clx[which(as.numeric(acfx$acf)[-1]<0)] = "orange"
	barplot(height = as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim=1.2*ylim, col = clx,
			ylab = "ACF", xlab="lag", main = "ACF of Standardized Residuals", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black", lty = 1)
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# ACF of squared standardized residuals
.plot.garchfit.11 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	zseries = as.numeric(residuals(x, standardize=TRUE))
	zseries[is.na(zseries)] = 0
	n 		= length(zseries)
	lag.max = as.integer(10*log10(n))
	acfx 	= acf(zseries^2, lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
	clx 	= vector(mode="character",length=lag.max)
	clx[which(as.numeric(acfx$acf)[-1]>=0)]="steelblue"
	clx[which(as.numeric(acfx$acf)[-1]<0)]="orange"
	barplot(height = as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim=1.2*ylim, col = clx,
			ylab = "ACF", xlab="lag", main = "ACF of Squared Standardized Residuals", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black", lty = 1)
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# News Impact Curve
.plot.garchfit.12 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	if(vmodel == "iGARCH"){
		warning(paste("\nplot-->: iGARCH newsimpact not available"))
	} else{
		ni = newsimpact(z = NULL, x)
		ni.y = ni$zy
		ni.x = ni$zx
		xf = ni$xexpr
		yf  = ni$yexpr
		plot( ni.x, ni.y, ylab = yf, xlab = xf, type = "l", lwd = 2, col = "steelblue", main = "News Impact Curve", cex.main = 0.8)
		#mtext(yf, side = 2, adj = 0.5, padj = -2.5, cex = 0.85)
		mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
		if(vmodel == "fGARCH"){
			mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
		}
		grid()
	}
}

#-------------------------------------------------------------------------------
# SECTION GARCH filter plots
#-------------------------------------------------------------------------------

.plotgarchfilter = function(x, which="ask",...)
{
	choices = c(
			"Series with 2 Conditional SD Superimposed",
			"Series with 1% VaR Limits",
			"Conditional SD (vs |returns|)",
			"ACF of Observations",
			"ACF of Squared Observations",
			"ACF of Absolute Observations",
			"Cross Correlation",
			"Empirical Density of Standardized Residuals",
			"QQ-Plot of Standardized Residuals",
			"ACF of Standardized Residuals",
			"ACF of Squared Standardized Residuals",
			"News-Impact Curve")
	.intergarchfilterPlot(x, choices = choices, plotFUN = paste(".plot.garchfilter", 1:12, 
					sep = "."), which = which, ...)
	# Return Value:
	invisible(x)
}

.intergarchfilterPlot = function(x, choices, plotFUN, which, ...)
{
	old.par <- par(no.readonly = TRUE)
	if (is.numeric(which)) {
		if(which>length(choices)) 
			stop("Not a valid choice. Plots choices are 1-12.\n",call. = FALSE)
		FUN = match.fun(plotFUN[which])
		FUN(x)
	}
	if(is.character(which))
	{
		if(which!="all" && which!="ask") stop("Not a valid choice.\n",call. = FALSE)
		if (which[1] == "all") {
			on.exit(par(old.par))
			#Which = rep(TRUE, times = length(choices))
			par(mfrow=c(3,4))
			for(i in 1:12){
				FUN = match.fun(plotFUN[i])
				FUN(x)
			}
		} else{
			.multgarchfilterPlot(x, choices, plotFUN, ...)
		}
	}
	
	invisible(x)
}

.multgarchfilterPlot = function(x, choices, ...)
{
	
	pick = 1
	while (pick > 0) {
		pick = menu (
				choices = paste(" ", choices),
				title = "\nMake a plot selection (or 0 to exit):")
		switch (pick,
				.plot.garchfilter.1(x),  .plot.garchfilter.2(x),  .plot.garchfilter.3(x),
				.plot.garchfilter.4(x),  .plot.garchfilter.5(x),  .plot.garchfilter.6(x),
				.plot.garchfilter.7(x),  .plot.garchfilter.8(x),  .plot.garchfilter.9(x),
				.plot.garchfilter.10(x), .plot.garchfilter.11(x), .plot.garchfilter.12(x))
	}
}

# Series with 2 Conditional SD Superimposed
.plot.garchfilter.1 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]
	xdates  = x@model$modeldata$index[insample]
	xsigma  = x@filter$sigma
	ci = 2
	plot(xdates, xseries, type = "l", col = "steelblue", ylab = "Returns", xlab="Time", 
			main = "Series with 2 Conditional SD Superimposed", cex.main = 0.8)
	lines(xdates, +ci* xsigma, col = "tomato1")
	lines(xdates, -ci* xsigma, col = "tomato1")
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	abline(h = 0, col = "grey", lty = 3)
	grid()
}

# Series with 1% VaR Limits
.plot.garchfilter.2 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]
	xdates  = x@model$modeldata$index[insample]
	xsigma 	= x@filter$sigma
	distribution = x@model$modeldesc$distribution
	xcmu 	= fitted(x)
	idx = x@model$pidx
	pars  = x@filter$ipars[,1]
	skew  = pars[idx["skew",1]]
	shape = pars[idx["shape",1]]
	if(distribution == "ghst") ghlambda = -shape/2 else ghlambda = pars[idx["ghlambda",1]]
	q01 	= qdist(distribution, 0.01, mu = xcmu, sigma = xsigma, lambda = ghlambda, skew = skew, shape = shape)
	q99 	= qdist(distribution, 0.99, mu = xcmu, sigma = xsigma, lambda = ghlambda, skew = skew, shape = shape)
	plot(xdates, xseries, type = "l", col = "steelblue", ylab = "Returns", xlab="Time", 
			main = "Series with with 1% VaR Limits", cex.main = 0.8)
	lines(xdates, q01, col = "tomato1")
	lines(xdates, q99, col = "green")
	mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	abline(h = 0, col = "grey", lty = 3)
	grid()
}

# Conditional SD
.plot.garchfilter.3 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = abs(x@model$modeldata$data[insample])
	xdates  = x@model$modeldata$index[insample]
	xsigma 	= x@filter$sigma
	plot(xdates, xseries, type = "l", col = "lightgrey", ylab = "Volatility", 
			xlab="Time", main = "Conditional Sigma (vs |returns|)", cex.main = 0.8)
	lines(xdates, xsigma, col = "steelblue")
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	abline(h = 0, col = "grey", lty = 3)
	grid()
}

# ACF of observations
.plot.garchfilter.4 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]
	lag.max = as.integer(10*log10(T))
	acfx 	= acf(xseries[insample], lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
	clx 	= vector(mode = "character", length = lag.max)
	clx[which(as.numeric(acfx$acf)[-1]>=0)] = "steelblue"
	clx[which(as.numeric(acfx$acf)[-1]<0)] = "orange"
	barplot(height = as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim = 1.2*ylim, col = clx,
			ylab = "ACF", xlab="lag", main = "ACF of Observations", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black", lty = 1)
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# ACF of squared observations
.plot.garchfilter.5 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]
	xdates  = x@model$modeldata$dates[insample]
	lag.max = as.integer(10*log10(T))
	acfx 	= acf(xseries^2, lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
	clx 	= vector(mode="character",length=lag.max)
	clx[which(as.numeric(acfx$acf)[-1]>=0)] = "steelblue"
	clx[which(as.numeric(acfx$acf)[-1]<0)] = "orange"
	barplot(height = as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim = 1.2*ylim, col = clx,
			ylab = "ACF", xlab="lag", main = "ACF of Squared Observations", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black", lty = 1)
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# ACF of absolute observations
.plot.garchfilter.6 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]
	lag.max = as.integer(10*log10(T))
	acfx 	= acf(abs(xseries), lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
	clx 	= vector(mode="character",length=lag.max)
	clx[which(as.numeric(acfx$acf)[-1]>=0)] = "steelblue"
	clx[which(as.numeric(acfx$acf)[-1]<0)] = "orange"
	barplot(height=as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim = 1.2*ylim, col = clx,
			ylab = "ACF", xlab = "lag", main = "ACF of Absolute Observations", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black", lty = 1)
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# Cross Correlations
.plot.garchfilter.7 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	T = x@model$modeldata$T
	insample = 1:T
	xseries = x@model$modeldata$data[insample]
	lag.max = as.integer(10*log10(T))
	ccfx	= ccf(xseries^2, xseries, lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(ccfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(ccfx$acf)))
	clx 	= vector(mode="character",length=(2*lag.max)+1)
	clx[which(as.numeric(ccfx$acf)>=0)]="steelblue"
	clx[which(as.numeric(ccfx$acf)<0)]="orange"
	barplot(height=as.numeric(ccfx$acf), names.arg=as.numeric(ccfx$lag),ylim=1.2*ylim,col=clx,
			ylab = "ACF", xlab="lag", main = "Cross-Correlations of \nSquared vs Actual Observations", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black")
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# Standardized Residuals Empirical Denisty
.plot.garchfilter.8 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	zseries = as.numeric(residuals(x, standardize=TRUE))
	distribution = x@model$modeldesc$distribution
	idx = x@model$pidx
	pars  = x@filter$ipars[,1]
	skew  = pars[idx["skew",1]]
	shape = pars[idx["shape",1]]
	if(distribution == "ghst") ghlambda = -shape/2 else ghlambda = pars[idx["ghlambda",1]]
	xmean 	= mean(zseries)
	xmedian = median(zseries)
	xsd 	= sd(zseries)	
	xlim 	= c(min(zseries), max(zseries))
	result 	= hist(x = zseries, col = "grey", border = "white",
			breaks = "Scott", main = "Empirical Density of Standardized Residuals", xlim = xlim, ylim = c(0,0.6),
			probability = TRUE, ylab="Probability", cex.main = 0.8, ...)
	box()
	#s 	= seq(xlim[1], xlim[2], length = 201)
	s = result$breaks
	y	= ddist(distribution, s, lambda = ghlambda, skew = skew, shape = shape)
	lines(s, dnorm(s, 0, 1), lwd = 2, col = "blue")
	lines(s, y, lwd = 2, col = "orange")
	abline(v = xmean, lwd = 2, col = "red")
	abline(v = xmedian, lwd = 2, col = "darkgreen")
	Text = paste("Median: ", round(xmedian, 2), "| Mean: ",signif(xmean, 3))
	mtext(Text, side = 3, adj = 0, col = "darkgrey", cex = 0.7)
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
	}
	abline(h = 0, col = "grey")
	lg.txt = c("normal Density", paste(distribution," (0,1) Fitted Density",sep=""))
	legend("topleft", legend = lg.txt, col = c("blue","orange"), pch = 1, cex = 0.8, bty = "n")
	grid()
}

# QQ-Plot of Standardized Residuals
.plot.garchfilter.9 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	zseries = as.numeric(residuals(x, standardize=TRUE))
	distribution = x@model$modeldesc$distribution
	idx = x@model$pidx
	pars  = x@filter$ipars[,1]
	skew  = pars[idx["skew",1]]
	shape = pars[idx["shape",1]]
	if(distribution == "ghst") ghlambda = -shape/2 else ghlambda = pars[idx["ghlambda",1]]
	.qqDist(y = zseries, dist = distribution, lambda = ghlambda, skew = skew, shape = shape)
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
	}
}



# ACF of standardized residuals
.plot.garchfilter.10 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	zseries = as.numeric(residuals(x, standardize=TRUE))
	zseries[is.na(zseries)] = 0
	n 		= length(zseries)
	lag.max = as.integer(10*log10(n))
	acfx 	= acf(zseries, lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
	clx 	= vector(mode = "character", length = lag.max)
	clx[which(as.numeric(acfx$acf)[-1]>=0)] = "steelblue"
	clx[which(as.numeric(acfx$acf)[-1]<0)] = "orange"
	barplot(height = as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim=1.2*ylim, col = clx,
			ylab = "ACF", xlab="lag", main = "ACF of Standardized Residuals", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black", lty = 1)
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# ACF of squared standardized residuals
.plot.garchfilter.11 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	zseries = as.numeric(residuals(x, standardize=TRUE))
	zseries[is.na(zseries)] = 0
	n 		= length(zseries)
	lag.max = as.integer(10*log10(n))
	acfx 	= acf(zseries^2, lag.max = lag.max, plot = FALSE)
	clim0	= qnorm((1 + 0.95)/2)/sqrt(acfx$n.used)
	ylim 	= range(c(-clim0, clim0, as.numeric(acfx$acf)[-1]))
	clx 	= vector(mode="character",length=lag.max)
	clx[which(as.numeric(acfx$acf)[-1]>=0)]="steelblue"
	clx[which(as.numeric(acfx$acf)[-1]<0)]="orange"
	barplot(height = as.numeric(acfx$acf)[-1], names.arg = as.numeric(acfx$lag)[-1], ylim=1.2*ylim, col = clx,
			ylab = "ACF", xlab="lag", main = "ACF of Squared Standardized Residuals", cex.main = 0.8)
	abline(h = c(clim0, -clim0), col = "tomato1", lty = 2)
	abline(h = 0, col = "black", lty = 1)
	box()
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
	}
	grid()
}

# News Impact Curve
.plot.garchfilter.12 = function(x, ...)
{
	vmodel  = x@model$modeldesc$vmodel
	if(vmodel == "iGARCH"){
		warning(paste("\nplot-->: iGARCH newsimpact not available"))
	} else{
		ni	= newsimpact(z = NULL, x)
		ni.y = ni$zy
		ni.x = ni$zx
		xf	= ni$xexpr
		yf  = ni$yexpr
		plot( ni.x, ni.y, ylab = yf, xlab = xf, type = "l", lwd = 2, col = "steelblue", main = "News Impact Curve", cex.main = 0.8)
		#mtext(yf, side = 2, adj = 0.5, padj = -2.5, cex = 0.85)
		mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
		if(vmodel == "fGARCH"){
			mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
		}
		grid()
	}
}

#-------------------------------------------------------------------------------
# SECTION GARCH sim plots
#-------------------------------------------------------------------------------
.plotgarchsim<-function(x, which="ask", m.sim = 1, ...)
{
	choices = c(
			"Conditional SD Simulation Path",
			"Return Series Simulation Path",
			"Conditional SD Simulation Density",
			"Return Series Simulation Path Density")
	.intergarchsimPlot(x,choices=choices, plotFUN = paste(".plot.garchsim", 1:4, sep = "."), which = which, 
			m.sim = m.sim, ...)
	# Return Value:
	invisible(x)
}

.intergarchsimPlot<-function(x, choices, plotFUN, which, m.sim = 1, ...)
{
	old.par <- par(no.readonly = TRUE)
	if (is.numeric(which)) {
		if(which>length(choices)) stop("Not a valid choice. Plots choices are 1-4.\n", call. = FALSE)
		FUN = match.fun(plotFUN[which])
		FUN(x, m.sim)
	}
	if(is.character(which))
	{
		if(which!="all" && which!="ask") stop("Not a valid choice.\n", call. = FALSE)
		if (which[1] == "all") {
			on.exit(par(old.par))
			par(mfrow=c(2,2))
			for(i in 1:4){
				FUN = match.fun(plotFUN[i])
				FUN(x, m.sim = m.sim)
			}
		}
		if (which[1] == "ask") {
			.multgarchsimPlot(x, choices, m.sim = m.sim, ...)
		}
	}
	invisible(x)
}

.multgarchsimPlot<-function(x, choices, m.sim = 1,  ...)
{
	n = dim(x@simulation$sigmaSim)[2]
	if(m.sim>n | m.sim<=0) stop("\ninvalid m.sim input\n", call. = FALSE)
	pick = 1
	while (pick > 0) {
			pick = menu (
				choices = paste(" ", choices),
				title = "\nMake a plot selection (or 0 to exit):")
				switch (pick,
					.plot.garchsim.1(x, m.sim),  .plot.garchsim.2(x, m.sim),  
					.plot.garchsim.3(x, m.sim),  .plot.garchsim.4(x, m.sim))
	}
}

.plot.garchsim.1<-function(x, m.sim,...)
{
	vmodel  = x@model$modeldesc$vmodel
	xseed = x@seed[m.sim]
	simsigma = matrix(x@simulation$sigmaSim[, m.sim], ncol=1)
	N2 = dim(simsigma)[1]
	T = x@model$modeldata$T
	insample = 1:T
	xdates = x@model$modeldata$index[insample]
	fwddates = .generatefwd(xdates, N = N2, dformat="%Y-%m-%d", periodicity="days")
	sigma = x@model$modeldata$sigma
	Ns = length(sigma)
	nx = ifelse(Ns>500, 500, min(Ns, 200))
	sigma = c(sigma[(Ns-nx):Ns], rep(NA, N2))
	ylim = c(0.95*min(simsigma, sigma, na.rm=TRUE),1.05*max(simsigma, sigma, na.rm=TRUE))
	simsigma = rbind(matrix(NA, ncol = 1, nrow = nx + 1), simsigma)
	plot(c(xdates[(Ns-nx):Ns], fwddates), sigma, type = "l", col = "steelblue", main = "Simulated Conditional Sigma",
			ylab = "Sigma", xlab = "Time/Horizon", ylim = ylim, cex.main = 0.8)
	abline(h = 0, col = "grey", lty = 3)
	abline(v = nx, col = "red", lty = 3)
	lines(c(xdates[(Ns-nx):Ns], fwddates), as.numeric(simsigma), col = "tomato1")
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
		mtext(paste("seed: ", xseed, sep = ""), side = 3, adj = 0, padj=0, col = "gray", cex = 0.5)
	} else{
		mtext(paste("seed: ", xseed, sep = ""), side = 3, adj = 0, padj=0, col = "gray", cex = 0.5)
	}
	lg.txt = c("Actual", "Simulated")
	legend("topleft", legend = lg.txt, col = c("steelblue", "tomato1"), y.intersp = 1.5, pch = 21, cex = 0.7)
	box()
	grid()
}

.plot.garchsim.2<-function(x, m.sim,...)
{
	vmodel  = x@model$modeldesc$vmodel
	xseed = x@seed[m.sim]
	simseries = matrix(x@simulation$seriesSim[,m.sim],ncol=1)
	N2 = dim(simseries)[1]
	simseries = as.numeric(simseries[,1])
	xdates = x@model$modeldata$index
	Ns = x@model$modeldata$T
	insample = 1:Ns
	fwddates = .generatefwd(xdates[insample], N = N2, dformat = "%Y-%m-%d", periodicity = "days")
	series = x@model$modeldata$data
	nx = ifelse(Ns>500, 500, min(Ns, 200))
	simseries = c(series[(Ns-nx):Ns], simseries)
	series = c(series[(Ns-nx):Ns], rep(NA, N2))
	ylim = c(0.95*min(simseries,na.rm=TRUE), 1.05*max(simseries,na.rm=TRUE))
	plot(c(xdates[(Ns-nx):Ns],fwddates), simseries, type = "l", col = "tomato1", main = "Simulated Series", 
			ylab = "Series", xlab = "Time/Horizon", ylim = ylim, cex.main = 0.8)
	abline(h = 0, col = "grey", lty = 3)
	abline(v = nx+1, col = "red", lty = 3)
	lines(c(xdates[(Ns-nx):Ns], fwddates),series,col="steelblue")
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
		mtext(paste("seed: ", xseed, sep = ""), side = 3, adj = 0, padj=0, col = "gray", cex = 0.5)
	} else{
		mtext(paste("seed: ", xseed, sep = ""), side = 3, adj = 0, padj=0, col = "gray", cex = 0.5)
	}
	lg.txt = c("Actual", "Simulated")
	legend("topleft", legend = lg.txt, col = c("steelblue", "tomato1"), y.intersp = 1.5, pch = 21, cex = 0.7)
	box()
	grid()
}

.plot.garchsim.3<-function(x, m.sim,...)
{
	vmodel  = x@model$modeldesc$vmodel
	xseed = x@seed[m.sim]
	sigma = x@model$modeldata$sigma
	simsigma = matrix(x@simulation$sigmaSim[,m.sim],ncol=1)
	s1 = density(sigma, kernel="epanechnikov")
	s2 = density(as.numeric(simsigma[,1]), kernel="epanechnikov")
	s1x = s1$x
	s1y = s1$y/length(s1$y)
	s2x = s2$x
	s2y = s2$y/length(s2$y)
	ylim = c(0, max(s1y,s2y))
	plot(s1x, s1y, ylim = ylim, main = "Conditional Sigma\n Kernel Density", type = "l", 
			ylab = "Probability", xlab = "Sigma", cex.main = 0.8, col = "steelblue")
	lines(s2x, s2y, col = "tomato1")
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), 
				side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
		mtext(paste("seed: ", xseed, sep = ""), side = 3, adj = 0, padj=0, col = "gray", cex = 0.5)
	} else{
		mtext(paste("seed: ", xseed, sep = ""), side = 3, adj = 0, padj=0, col = "gray", cex = 0.5)
	}
	lg.txt = c("Actual", "Simulated")
	legend("topleft", legend = lg.txt, col = c("steelblue", "tomato1"), y.intersp = 1.5, pch = 21, cex = 0.7)
	box()
	grid()
}


.plot.garchsim.4<-function(x, m.sim,...)
{
	vmodel  = x@model$modeldesc$vmodel
	xseed = x@seed[m.sim]
	xseries = x@model$modeldata$data
	T = x@model$modeldata$T
	insample = 1:T
	simseries = matrix(x@simulation$seriesSim[,m.sim],ncol=1)
	s2 = density(as.numeric(simseries[,1]), kernel="epanechnikov")
	n = length(s2$y)
	s1 = density(xseries[insample], kernel="epanechnikov",n=n)
	s1x = s1$x
	s1y = s1$y/length(s1$y)
	s2x = s2$x
	s2y = s2$y/length(s2$y)
	ylim = c(0, max(s1y,s2y))
	plot(s1x, s1y, ylim = ylim, main = "Time Series\n Kernel Density", type = "l", 
			ylab = "Probability", xlab = "Returns", cex.main = 0.8, col = "steelblue")
	lines(s2x, s2y, col = "tomato1")
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj = 1.5, col = "gray", cex = 0.5)
		mtext(paste("seed: ", xseed, sep = ""), side = 3, adj = 0, padj=0, col = "gray", cex = 0.5)
	} else{
		mtext(paste("seed: ", xseed, sep = ""), side = 3, adj = 0, padj=0, col = "gray", cex = 0.5)
	}
	lg.txt = c("Actual", "Simulated")
	legend("topleft", legend = lg.txt, col = c("steelblue", "tomato1"), y.intersp = 1.5, pch = 21, cex = 0.7)
	box()
	grid()
}

#-------------------------------------------------------------------------------
# SECTION GARCH prediction plots
#-------------------------------------------------------------------------------
.plotgarchforecast<-function(x, which="ask", n.roll = 0, ...)
{
	choices = c(
			"Time Series Prediction (unconditional)",
			"Time Series Prediction (rolling)",
			"Sigma Prediction (unconditional)",
			"Sigma Prediction (rolling)")
	.intergarchforecastPlot(x,choices=choices, plotFUN = paste(".plot.garchforecast", 1:4, sep = "."), 
			which = which, n.roll = n.roll, ...)
	# Return Value:
	invisible(x)
}

.intergarchforecastPlot<-function(x, choices, plotFUN, which,  n.roll = 0, ...)
{
	old.par <- par(no.readonly = TRUE)
	if (is.numeric(which)) {
		if(which>length(choices)) stop("Not a valid choice. Plots choices are 1-4.\n", call. = FALSE)
		FUN = match.fun(plotFUN[which])
		FUN(x, n.roll)
	}
	if(is.character(which))
	{
		if(which!="all" && which!="ask") stop("Not a valid choice.\n", call. = FALSE)
		if (which[1] == "all") {
			on.exit(par(old.par))
			#Which = rep(TRUE, times = length(choices))
			par(mfrow=c(2,2))
			for(i in 1:4){
				FUN = match.fun(plotFUN[i])
				FUN(x, n.roll)
			}
		}
		if (which[1] == "ask") {
			.multgarchforecastPlot(x, choices, n.roll, ...)
		}
	}
	
	invisible(x)
}

.multgarchforecastPlot<-function(x, choices, n.roll = 0, ...)
{
	# A function originally written by Diethelm Wuertz
	
	# Description:
	#   Internal plot function
	pick = 1
	while (pick > 0) {
		pick = menu(
				### choices = paste("plot:", choices),
				choices = paste(" ", choices),
				title = "\nMake a plot selection (or 0 to exit):")
		# up to 19 plot functions ...
		switch (pick,.plot.garchforecast.1(x, n.roll),  .plot.garchforecast.2(x), 
				.plot.garchforecast.3(x, n.roll), .plot.garchforecast.4(x, n.roll))
	}
}

# Time Series Plot
.plot.garchforecast.1 = function(x, n.roll = 0, ...)
{
	vmodel = x@model$modeldesc$vmodel
	# 1. Time Series:
	nr = x@forecast$n.roll
	if(n.roll > nr) stop("plot-->error: n.roll choice is invalid", call. = FALSE)
	n = x@forecast$n.ahead
	N = x@forecast$N - x@forecast$n.start
	forseries = x@forecast$seriesFor[,n.roll+1]
	forsigma = x@forecast$sigmaFor[,n.roll+1]
	xdates = x@model$modeldata$index[(N+n.roll-min(N,100)):(N+n.roll)]
	fdates = seq(tail(xdates,1), by = x@model$modeldata$period, length.out = n+1)[-1]
	series = x@model$modeldata$data[(N+n.roll-min(N,100)):(N+n.roll)]
	
	xforseries = c(series, forseries)
	series = c(series, rep(NA, n))
	ylim=c(0.95*min(xforseries,na.rm=TRUE), 1.05*max(xforseries,na.rm=TRUE))
	plot(c(xdates, fdates), as.numeric(xforseries), type="l", col="steelblue", 
			main = paste("Forecast Series\n w/th unconditional 1-Sigma bands", sep = "") ,
			ylab="Series",xlab="Time/Horizon", ylim = ylim, cex.main = 0.7, cex.axis = 0.8, cex.lab = 0.9)
	abline(h = 0, col = "grey", lty = 3)
	Zup = forseries+1*forsigma
	Zdn = forseries-1*forsigma
	for(i in 2:n) rect(fdates[i-1], Zdn[i-1], fdates[i], Zup[i], col = colors()[142], border=NA)
	lines(c(xdates, fdates), series, col="steelblue")
	lines(fdates, forseries, col="tomato1")
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel=="fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
		mtext(paste("Horizon: ",n,sep=""), side = 3, adj = 0, col = "gray", cex = 0.5)
	} else{
		mtext(paste("Horizon: ",n,sep=""), side = 3, adj = 0, col = "gray", cex = 0.5)
	}
	lg.txt = c("Actual","Forecast")
	legend("topleft", legend = lg.txt, col = c("steelblue", "tomato1"), y.intersp = 1.5, pch = 21, cex = 0.7, bty="n")
	box()
	grid()
}

# rolling forecast plot (actual vs forecast)
.plot.garchforecast.2 = function(x, ...)
{
	vmodel = x@model$modeldesc$vmodel
	# 1. Time Series:
	nr = x@forecast$n.roll
	if(nr<5) stop("\nn.roll less than 5!...does not make sense to provide this plot.")
	N = x@forecast$N - x@forecast$n.start
	fdata  = x@forecast$seriesFor[1,]
	fsigma = x@forecast$sigmaFor[1,]
	xdata = x@model$modeldata$data[((N+1)-min(N, 25)):(N+nr)]
	xsigma = x@model$modeldata$sigma[((N+1)-min(N, 25)):(N+nr)]
	xdates = x@model$modeldata$index[((N+1)-min(N, 25)):(N+nr)]
	ns = length(xdata)
	xplus =  xdata + 2*xsigma
	xminus = xdata - 2*xsigma
	fplus = c(rep(NA, (ns - nr)), fdata[-length(fdata)] +  2*fsigma[-length(fsigma)])
	fminus = c(rep(NA, (ns - nr)), fdata[-length(fdata)] - 2*fsigma[-length(fsigma)])
	fdata = c(rep(NA, (ns - nr)), fdata[-length(fdata)])
	
	ylim=c(0.95*min(xminus,na.rm=TRUE), 1.2*max(xplus,na.rm=TRUE))
	plot(xdates, xdata, type="l", col="black", 
			main = paste("Rolling Forecast vs Actual Series\n w/th conditional 2-Sigma bands", sep = "") ,
			ylab="Series",xlab="Time/Horizon", ylim = ylim, cex.main = 0.7, 
			cex.axis = 0.8, cex.lab = 0.9)
	#for(i in 6:(nr+25)) rect(xdates[i-5], xminus[i-5], xdates[i], xplus[i], col = colors()[411], border=NA)
	for(i in 26:(nr+25)) rect(xdates[i-1], fminus[i-1], xdates[i], fplus[i], col = colors()[142], border=NA)
	lines(xdates, xdata, col = "steelblue")
	lines(xdates, xplus, col = "lightgrey", lwd = 0.5)
	lines(xdates, xminus, col = "lightgrey", lwd = 0.5)
	abline(h = 0, col = "grey", lty = 3)
	lines(xdates, fdata, col = "tomato1", lwd = 2.5)
	#lines(xts(fplus,xdates),  col = "brown", lwd = 0.5)
	#lines(xts(fminus,xdates),  col = "brown", lwd = 0.5)
	
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel=="fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
		mtext(paste("Horizon: ", nr, sep=""), side = 3, adj = 0, col = "gray", cex = 0.5)
	} else{
		mtext(paste("Horizon: ", nr, sep=""), side = 3, adj = 0, col = "gray", cex = 0.5)
	}
	lg.txt = c("Actual","Forecast")
	legend("topleft", legend = lg.txt, col = c("steelblue", "tomato1"), y.intersp = 1.5, pch = 21, cex = 0.7, bty="n")
	box()
	grid()
}

.plot.garchforecast.3 = function(x, n.roll = 0, ...)
{
	vmodel = x@model$modeldesc$vmodel
	nr = x@forecast$n.roll
	if(n.roll > nr) stop("plot-->error: n.roll choice is invalid", call. = FALSE)
	n = x@forecast$n.ahead
	N = x@forecast$N - x@forecast$n.start
	forsigma = x@forecast$sigmaFor[,n.roll+1]
	xdates = x@model$modeldata$index[(N+n.roll-min(N,25)):(N+n.roll)]	
	fdates = seq(tail(xdates,1), by = x@model$modeldata$period, length.out = n+1)[-1]	
	sigma = x@model$modeldata$sigma[(N+n.roll-min(N,25)):(N+n.roll)]
	forsigma = c(sigma, forsigma)
	sigma = c(sigma, rep(NA,n))
	ylim=c(0.95*min(forsigma,na.rm=TRUE),1.05*max(forsigma,na.rm=TRUE))
	plot(c(xdates,fdates), forsigma, type = "l", col = "black", 
			main = paste("Forecast Unconditional Sigma\n (n.roll = ", n.roll,")", sep = "") , 
			ylab = "Sigma", xlab = "Time/Horizon", ylim = ylim, cex.main = 0.7, 
			cex.axis = 0.8, cex.lab = 0.9)
	abline(h = 0, col = "grey", lty = 3)
	lines(c(xdates,fdates), forsigma, , col = "tomato1")
	lines(c(xdates,fdates), sigma, col = "steelblue")
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
		mtext(paste("Horizon: ",n,sep=""), side = 3, adj = 0, col = "gray", cex = 0.5)
	} else{
		mtext(paste("Horizon: ",n,sep=""), side = 3, adj = 0, col = "gray", cex = 0.5)
	}
	lg.txt<-c("Actual","Forecast")
	legend("topleft", legend=lg.txt,col=c("steelblue","tomato1"), 
			y.intersp=1.5,pch=21,cex=0.7, bty="n")
	box()
	grid()
}

.plot.garchforecast.4 = function(x, ...)
{
	vmodel = x@model$modeldesc$vmodel
	# 1. Time Series:
	nr = x@forecast$n.roll
	if(nr<5) stop("\nn.roll less than 5!...does not make sense to provide this plot.")
	N = x@forecast$N - x@forecast$n.start
	# exclude the last forecast (when used) since there is no actual data to
	# compare it with
	fsigma = x@forecast$sigmaFor[1,1:nr]
	fdates = x@model$modeldata$index[(N+1):(N+nr)]
	xdata  = x@model$modeldata$data[((N+1)-min(N, 25)):(N+nr)]
	xsigma = x@model$modeldata$sigma[((N+1)-min(N, 25)):(N+nr)]
	xdates = x@model$modeldata$index[((N+1)-min(N, 25)):(N+nr)]
	ns = length(xdata)	
	ylim=c(0.95*min(abs(xdata),na.rm=TRUE), 1.2*max(abs(xdata),na.rm=TRUE))
	plot(xdates, abs(xdata), type="l", col="lightgrey", 
			main = paste("Forecast Rolling Sigma vs |Series|", sep = "") ,
			ylab="Sigma",xlab="Time/Horizon", ylim = ylim, cex.main = 0.7, 
			cex.axis = 0.8, cex.lab = 0.9)
	lines(xdates, xsigma, col = "steelblue")
	lines(fdates, fsigma, col = "tomato1", lwd = 0.5)
	abline(h = 0, col = "grey", lty = 3)	
	mtext(paste("GARCH model : ", vmodel), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	if(vmodel=="fGARCH"){
		mtext(paste("fGARCH submodel: ", x@model$modeldesc$vsubmodel, sep = ""), side = 4, adj = 0, padj=1.5, col = "gray", cex = 0.5)
		mtext(paste("Horizon: ", nr, sep=""), side = 3, adj = 0, col = "gray", cex = 0.5)
	} else{
		mtext(paste("Horizon: ", nr, sep=""), side = 3, adj = 0, col = "gray", cex = 0.5)
	}
	lg.txt = c("Actual","Forecast","|Series|")
	legend("topleft", legend = lg.txt, col = c("steelblue", "tomato1","lightgrey"), y.intersp = 1.5, pch = 21, cex = 0.7, bty="n")
	box()
	grid()
}


#-------------------------------------------------------------------------------
# SECTION GARCH path plots
#-------------------------------------------------------------------------------
.plotgarchpath = function(x, which="ask", m.sim = 1, ...)
{
	choices = c(
			"Conditional SD Simulation Path",
			"Return Series Simulation Path",
			"Conditional SD Simulation Density",
			"Return Series Simulation Path Density")
	.intergarchpathPlot(x, choices = choices, plotFUN = paste(".plot.garchpath", 1:4, sep = "."), 
			which = which, m.sim = m.sim,...)
	# Return Value:
	invisible(x)
}

.intergarchpathPlot = function(x, choices, plotFUN, which, m.sim = 1, ...)
{
	old.par <- par(no.readonly = TRUE)
	if (is.numeric(which)) {
		if(which>length(choices)) stop("Not a valid choice. Plots choices are 1-4.\n", call. = FALSE)
		FUN = match.fun(plotFUN[which])
		FUN(x, m.sim = m.sim)
	}
	if(is.character(which))
	{
		if(which!="all" && which!="ask") stop("Not a valid choice.\n", call. = FALSE)
		if (which[1] == "all") {
			on.exit(par(old.par))
			par(mfrow=c(2,2))
			for(i in 1:4){
				FUN = match.fun(plotFUN[i])
				FUN(x, m.sim = m.sim)
			}
		}
		if (which[1] == "ask") {
			.multgarchpathPlot(x, choices, m.sim = m.sim, ...)
		}
	}
	invisible(x)
}

.multgarchpathPlot<-function(x, choices, m.sim = 1, ...)
{
	n = dim(x@path$sigmaSim)[2]
	if(m.sim>n | m.sim<=0) stop("\ninvalid m.sim input\n", call. = FALSE)
	pick = 1
	npick = 0
	while (pick > 0) {
		pick = menu (
				choices = paste(" ", choices),
				title = "\nMake a plot selection (or 0 to exit):")
			switch (pick,
					.plot.garchpath.1(x, m.sim), .plot.garchpath.2(x, m.sim),  
					.plot.garchpath.3(x, m.sim), .plot.garchpath.4(x, m.sim))
	}
}

.plot.garchpath.1 = function(x, m.sim = 1, ...)
{
	model = strsplit(class(x),"path")
	xseed = x@seed[m.sim]
	simsigma = matrix(x@path$sigmaSim[,m.sim],ncol=1)
	plot(simsigma, type = "l", col = "steelblue", main = "Simulated Path of Conditional Sigma", 
			ylab="Sigma", xlab="Time/Horizon", cex.main = 0.7)
	mtext(paste("GARCH model :",model), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	box()
	grid()
}

.plot.garchpath.2 = function(x, m.sim = 1, ...)
{
	model = strsplit(class(x),"sim")
	xseed = x@seed[m.sim]
	simseries = matrix(x@path$seriesSim[,m.sim],ncol=1)
	simseries = as.numeric(simseries[,1])
	simseries = simseries
	plot(simseries,type="l",col="red", main = "Simulated Path of Series", ylab = "Series", 
			xlab="Time/Horizon", cex.main = 0.7)
	mtext(paste("GARCH model :",model), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	box()
	grid()
}

.plot.garchpath.3 = function(x, m.sim = 1,...)
{
	model = strsplit(class(x),"sim")
	xseed = x@seed[m.sim]
	simsigma = matrix(x@path$sigmaSim[,m.sim], ncol = 1)
	s2 = density(as.numeric(simsigma[,1]), kernel = "epanechnikov")
	plot(s2, main="Simulated Conditional Sigma\n Kernel Density",type="l",ylab="Probability", xlab="Sigma",
			col="steelblue", cex.main = 0.7)
	mtext(paste("GARCH model :",model), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	box()
	grid()
}


.plot.garchpath.4 = function(x, m.sim = 1,...)
{
	model = strsplit(class(x),"sim")
	xseed = x@seed[m.sim]
	simseries = matrix(x@path$seriesSim[,m.sim], ncol = 1)
	s2 = density(as.numeric(simseries[,1]), kernel = "epanechnikov")
	plot(s2, main = "Simulated Conditional Time Series\n Kernel Density",type="l",ylab="Probability", xlab="Returns",
			col="red", cex.main = 0.7)
	mtext(paste("GARCH model :",model), side = 4, adj = 0, padj=0, col = "gray", cex = 0.5)
	box()
	grid()
}

#-------------------------------------------------------------------------------
# SECTION GARCH roll plots
#-------------------------------------------------------------------------------
.plotgarchroll = function(x, which="ask", VaR.alpha = 0.01, density.support=c(-0.15, 0.15), ...)
{
	if(!is.null(x@model$noncidx)) stop("\nObject containts non-converged estimation windows.")
	choices = c(
			"Density Forecast",
			"Sigma Forecast",
			"Series Forecast",
			"VaR Forecast",
			"Fit Coefficients (with s.e. bands)")
	.intergarchrollPlot(x, choices=choices, plotFUN = paste(".plot.garchroll", 1:5, sep = "."), 
			which = which, VaR.alpha = VaR.alpha, density.support = density.support, ...)
	# Return Value:
	invisible(x)
}

.intergarchrollPlot = function(x, choices, plotFUN, which, VaR.alpha = 0.01, 
		density.support=c(-0.15, 0.15), ...)
{
	old.par <- par(no.readonly = TRUE)
	if (is.numeric(which)) {
		if(which>length(choices)) stop("Not a valid choice. Plots choices are 1-5.\n", call. = FALSE)
		FUN = match.fun(plotFUN[which])
		FUN(x, VaR.alpha = VaR.alpha, density.support = density.support, ...)
	}
	if(is.character(which))
	{
		if(which!="all" && which!="ask") stop("Not a valid choice.\n", call. = FALSE)
		if (which[1] == "all") {
			on.exit(par(old.par))
			par(mfrow=c(2,2))
			.plot.garchroll.1(x, VaR.alpha, density.support, ...)
			.plot.garchroll.2(x, VaR.alpha, density.support, ...)
			.plot.garchroll.3(x, VaR.alpha, density.support, ...)
			.plot.garchroll.4(x, VaR.alpha, density.support, ...)
		}
		if (which[1] == "ask") {
			.multgarchrollPlot(x, choices, VaR.alpha, density.support,...)
		}
	}
	invisible(x)
}

.multgarchrollPlot = function(x, choices, VaR.alpha = 0.01, density.support = c(-0.15, 0.15), ...)
{
	pick = 1
	while (pick > 0) {
		pick = menu (
				choices = paste(" ", choices),
				title = "\nMake a plot selection (or 0 to exit):")
		switch (pick,
				.plot.garchroll.1(x, VaR.alpha, density.support,...),  
				.plot.garchroll.2(x, VaR.alpha, density.support,...),  
				.plot.garchroll.3(x, VaR.alpha, density.support,...),
				.plot.garchroll.4(x, VaR.alpha, density.support,...),
				.plot.garchroll.5(x, VaR.alpha, density.support,...))
	}
}

# rolling sigma forecast comparison plot
.plot.garchroll.1 = function(x, VaR.alpha = 0.01, density.support = c(-0.15, 0.15), ...)
{
	density = x@forecast$density
	distribution = x@model$spec@model$modeldesc$distribution
	T = NROW(density)
	# we plot a maximum of 500 forecasts
	esd = floor(seq(1,T, length.out = min(500, T)))
	nesd = length(esd)
	colr = topo.colors(nesd, alpha = 1)
	xdate = rownames(density)
	xseq = seq(density.support[1], density.support[2], length.out=1000)
	yseq = apply(as.data.frame(2:nesd), 1, FUN=function(i)  
				ddist(xseq, mu = density[esd[i],1], sigma = density[esd[i],2], 
						lambda = density[esd[i],5], skew = density[esd[i],3], 
						shape = density[esd[i],4], distribution = distribution))
	plot(xseq, ddist(xseq, mu = density[esd[1],1], sigma = density[esd[1],2], 
					lambda = density[esd[1],5], skew = density[esd[1],3], 
					shape = density[esd[1],4], distribution = distribution), type="l", col = "steelblue", 
			main = paste("n.ahead-", 1," Forecast Density (time varying)",sep=""), ylab="", xlab = "", 
			cex.main = 0.7, cex.axis = 0.8, cex.lab=0.9)
	for(i in 1:(nesd-1)){
		lines(xseq, yseq[,i], col = colr[i+1])
	}
	xesd = floor(seq(1,length(esd), length.out = 5))
	legend("topright", legend = as.character(xdate[esd[xesd]]), fill = colr[xesd], col = colr[xesd], bty = "n")	
	invisible(x)
}

# rolling sigma forecast comparison plot
.plot.garchroll.2 = function(x, VaR.alpha = 0.01, density.support = c(-0.15, 0.15), ...)
{
	density = x@forecast$density
	plot(as.Date(rownames(density)), abs(density[,6]), type="l", col = "grey", 
			main = paste("Sigma Forecast vs |Series|", sep = ""), 
			ylab = "", xlab  = "", cex.main = 0.7,  cex.axis = 0.8, cex.lab=0.9)
	lines(as.Date(rownames(density)), abs(density[,2]), col = "steelblue", lwd = 1.5)
	grid()
	invisible(x)
}

# rolling series forecast comparison plot
.plot.garchroll.3 = function(x, VaR.alpha = 0.01, density.support = c(-0.15, 0.15), ...)
{
	density = x@forecast$density
	plot(as.Date(rownames(density)), density[,6], type="l", col = "grey", 
			main = paste("Series Forecast vs Realized", sep = ""), 
			ylab = "", xlab  = "", cex.main = 0.7, cex.axis = 0.8, cex.lab=0.9)
	lines(as.Date(rownames(density)), abs(density[,1]), col = "tomato1", lwd = 1.5)
	grid()
	invisible(x)
}

# rolling VaR backtest plot
.plot.garchroll.4 = function(x, VaR.alpha = 0.01, density.support = c(-0.15, 0.15), ...)
{
	vmodel = x@model$spec@model$modeldesc$vmodel
	v.a = x@model$VaR.alpha	
	if(!x@model$calculate.VaR) stop("\nplot-->error: VaR was not calculated for this object\n", call.=FALSE)
	if(!is.null(v.a) && !any(v.a==VaR.alpha[1])) stop("\nplot-->error: VaR.alpha chosen is invalid for the object\n", call.=FALSE)
	A = paste("alpha(", as.integer(VaR.alpha*100), "%)",sep="")
	VaRplot(VaR.alpha, as.xts(x@forecast$VaR[,"realized", drop=FALSE]), as.xts(x@forecast$VaR[,A,drop=FALSE]))
	invisible(x)
}

.plot.garchroll.5 = function(x, VaR.alpha = 0.01, density.support = c(-0.15, 0.15), ...)
{
	# get the no. of coef of a fit
	vmodel = x@model$spec@model$modeldesc$vmodel
	if(!x@model$keep.coef) stop("\n\nplot-->error: keep.coef set to FALSE in estimation\n")
	coefs = x@model$coef
	m = dim(coefs[[1]]$coef)[1]
	N = length(coefs)
	Z = matrix(NA, ncol = m, nrow = N)
	Zup = matrix(NA, ncol = m, nrow = N)
	Zdn = matrix(NA, ncol = m, nrow = N)
	for(i in 1:m){
		Z[,i] = sapply(coefs, FUN = function(y) y$coef[i,1])
		Zup[,i] = Z[,i]+sapply(coefs, FUN = function(y) y$coef[i,2])
		Zdn[,i] = Z[,i]-sapply(coefs, FUN = function(y) y$coef[i,2])
	}
	dt = sapply(coefs, FUN = function(y) as.character(y$index))
	cnames = rownames(coefs[[1]]$coef)
	np = .divisortable(m)
	par(mfrow = c(np[1], np[2]))
	for(i in 1:m){
		plot(as.Date(dt), Z[,i], type="l", ylim = c(min(Zdn[,i]), max(Zup[,i])), 
				ylab = "value" , xlab = "", main = "",  ann = FALSE)
		lines(as.Date(dt), Zdn[,i], col=2)
		lines(as.Date(dt), Zup[,i], col=2)
		title(cnames[i], line = 0.4, cex = 0.9)
		grid()
	}
	title(main = list(paste(vmodel," fit coefficients (across ",N," refits) with robust s.e. bands",sep=""), 
					cex = 1.5, col = "steelblue", font=2), outer = TRUE, line = -1.4)
	invisible(x)
}

#-------------------------------------------------------------------------------
# SECTION GARCH Distribution (Parameter Uncertainty) plots
#-------------------------------------------------------------------------------
.plotgarchdist = function(x, which="ask", window = 1, ...)
{
	choices = c(
			"Parameter Density Plots",
			"Bivariate Plots",
			"Other Density Plots (Persistence, Likelihood, ...)",
			"Asymptotic Efficiency Plots")
	
	.intergarchdistPlot(x, choices = choices, plotFUN = paste(".plot.garchdist", 1:4, sep = "."), 
			which = which, window = window, ...)
	# Return Value:
	invisible(x)
}

.intergarchdistPlot = function(x, choices, plotFUN, which, window, ...)
{
	if (is.numeric(which)) {
		if(which>length(choices)) 
			stop("Not a valid choice. Plots choices are 1-4.\n", call. = FALSE)
		if(which == 4 && !x@dist$detail$recursive)
			stop("Asymptotic Efficiency Plots only available for option recursive TRUE", .call = FALSE)
		FUN = match.fun(plotFUN[which])
		FUN(x, window, ...)
	}
	if(is.character(which))
	{
		if(which!="ask") stop("Not a valid choice.\n", call. = FALSE)
		if (which[1] == "ask") {
			.multgarchdistPlot(x, choices, window, ...)
		}
	}
	invisible(x)
}

.multgarchdistPlot = function(x, choices, window, ...)
{
	pick = 1
	while (pick > 0) {
		pick = menu (
				choices = paste(" ", choices),
				title = "\nMake a plot selection (or 0 to exit):")
		switch (pick,
				.plot.garchdist.1(x, window, ...),  
				.plot.garchdist.2(x, window, ...),  
				.plot.garchdist.3(x, ...),
				.plot.garchdist.4(x, ...))
	}
}

# Parameter Density Plots
.plot.garchdist.1 = function(x, window = 1,  ...)
{
	cf = as.data.frame(x, window = window)
	n = dim(cf)[2]
	vmodel = x@model$modeldesc$vmodel
	vsubmodel = x@model$modeldesc$vsubmodel
	truecf = unlist(x@truecoef[,1])
	cnames = names(truecf)
	z = .divisortable(n)
	par(mfrow=c(z[1], z[2]))
	for(i in 1:n){
		str = paste("Parameter ", cnames[i],"\nTrue Value: ", signif(truecf[i], digits = 3))
		plot(density(cf[,i], kernel = "gaussian", na.rm = TRUE), col = "steelblue4", main = str, cex.main = 0.7)
		mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", 
				cex = 0.5)
		if(vmodel == "fGARCH"){
			mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
					padj = 1.5, col = "gray", cex = 0.5)
		}
	}
}

# Bivariate Plots
.plot.garchdist.2= function(x, window = 1,...)
{
	.dist2dplot(x, window = window, ...)
	invisible()
}

# stat plots
.plot.garchdist.3 = function(x, ...)
{
	.statsplot(x, ...)
}

# rmse backtest plot
.plot.garchdist.4 = function(x,  ...)
{
	if(!x@dist$details$recursive){
		print("no rmse plots for non recursive option")
		return()
	} else{
		.rmseplots(x, ...)
	}
}


# x@model$type!="full"
#-------------------------------------------------------------------------------
# SECTION GARCH boot plots
#-------------------------------------------------------------------------------
.plotgarchboot = function(x, which="ask", ...)
{
	choices = c(
			"Parameter Density Plots",
			"Series Standard Error Plots",
			"Sigma  Standard Error Plots")
	.intergarchbootPlot(x, choices = choices, plotFUN = paste(".plot.garchboot", 1:3, sep = "."), 
			which = which, window = window, ...)
	# Return Value:
	invisible(x)
}

.intergarchbootPlot = function(x, choices, plotFUN, which, ...)
{
	old.par <- par(no.readonly = TRUE)
	if (is.numeric(which)) {
		if(which>length(choices)) 
			stop("Not a valid choice. Plots choices are 1-3.\n", call. = FALSE)
		if(which == 1 && x@model$type!="full")
			stop("Parameter Density Plots only available for full bootstrap method", .call = FALSE)
		FUN = match.fun(plotFUN[which])
		FUN(x, window, ...)
	}
	if(is.character(which))
	{
		if(which!="ask" & which!="all") stop("Not a valid choice.\n", call. = FALSE)
		if (which[1] == "ask") {
			.multgarchbootPlot(x, choices,  ...)
		} else{
			on.exit(par(old.par))
			par(mfrow = c(2,1))
			.plot.garchboot.2(x, ...)
			.plot.garchboot.3(x, ...)
		}
	}
	invisible(x)
}

.multgarchbootPlot = function(x, choices,  ...)
{
	pick = 1
	while (pick > 0) {
		pick = menu (
				choices = paste(" ", choices),
				title = "\nMake a plot selection (or 0 to exit):")
		switch (pick,
				.plot.garchboot.1(x, ...),  
				.plot.garchboot.2(x, ...),  
				.plot.garchboot.3(x, ...))
	}
}

# Parameter Density Plots (for full method)
.plot.garchboot.1 = function(x,  ...)
{
	if(x@model$type!="full") stop("Parameter Density Plots only available for full bootstrap method", .call = FALSE)
	cf = x@bcoef
	vmodel = x@model$modeldesc$vmodel
	vsubmodel = x@model$modeldesc$vsubmodel
	truecf = x@model$truecoef
	n = length(truecf)
	cnames = names(truecf)
	z = .divisortable(n)
	par(mfrow=c(z[1], z[2]))
	for(i in 1:n){
		str = paste("Parameter ", cnames[i],"\nTrue Value: ", signif(truecf[i], digits = 3))
		plot(density(cf[,i], kernel = "gaussian", na.rm = TRUE), col = "steelblue4", main = str,
				cex.main = 0.85)
		mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", 
				cex = 0.5)
		if(vmodel == "fGARCH"){
			mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
					padj = 1.5, col = "gray", cex = 0.5)
		}
	}
}

# Series Error Plots
.plot.garchboot.2 = function(x, ...)
{
	vmodel = x@model$modeldesc$vmodel
	vsubmodel = x@model$modeldesc$vsubmodel
	n.ahead = x@model$n.ahead
	seriesfor = fitted(x@forc)
	xser = as.data.frame(x, which = "series")
	N = dim(xser)[1]
	serp = as.data.frame(x, which = "series", type = "q", qtile = c(0.05, 0.25, 0.75, 0.95))
	miny = min(serp[1,])
	maxy = max(serp[4,])
	meanser = apply(xser, 2, FUN = function(x) mean(x))
	plot(seriesfor, type = "l", col = "red", ylim = c(miny, maxy), main = "Series Forecast
					with Bootstrap Error Bands\n (q: 5%, 25%, 75%, 95%)", cex.main = 0.7,
			ylab = "returns", xlab = "T+i")
	lines(as.numeric(meanser), col = "black")
	#lines(as.numeric(rx), col = "green")
	points(as.numeric(serp[1,]), col = "steelblue1", pch = 19, cex = 0.5)
	points(as.numeric(serp[2,]), col = "steelblue2", pch = 19, cex = 0.5)
	points(as.numeric(serp[3,]), col = "steelblue3", pch = 19, cex = 0.5)
	points(as.numeric(serp[4,]), col = "steelblue4", pch = 19, cex = 0.5)
	n.overlays = round(0.1*n.ahead)
	if(n.overlays == 0) n.overlays = 1
	dens1 = apply(xser, 2, FUN = function(x) density(x, kernel = "gaussian", 
						n = max(512, round(0.01*N))))
	.densityoverlay(as.numeric(meanser), dens1, n.overlays = n.overlays)
	mtext(paste("GARCH model :",vmodel), side = 4, adj = 0, padj=0, col = "gray", 
			cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
				padj = 1.5, col = "gray", cex = 0.5)
	}
	legend.txt = c("forecast", "bootstrapped")
	legend("bottomleft", legend = legend.txt, fill = c("red", "black"), col = c("red", "black"),
			bg = "n", bty = "n", cex = 0.6)
}

# Sigma Error Plots
.plot.garchboot.3 = function(x, ...)
{
	vmodel = x@model$modeldesc$vmodel
	vsubmodel = x@model$modeldesc$vsubmodel
	n.ahead = x@model$n.ahead
	fs = rep(NA, n.ahead)
	colr = NULL
	namr = NULL
	if(!is.null(x@model$modeldata$filtered.s)){
		n = length(x@model$modeldata$filtered.s)
		if(n.ahead>n) fs[1:n] = x@model$modeldata$filtered.s else fs = x@model$modeldata$filtered.s[1:n.ahead]
		colr = "green"
		namr = "filtered"
	}
	sigmafor = sigma(x@forc)
	sigdist = vector(mode = "list", length = n.ahead)
	sigp = as.data.frame(x, which = "sigma", type = "q", qtile = c(0.05, 0.25, 0.75, 0.95))
	miny = min(sigp[1,])
	maxy = max(sigp[4,])
	meansig = apply(x@fsigma, 2, FUN = function(x) mean(x))
	plot(sigmafor, type = "l", col = "red", ylim = c(miny, maxy), main = "Sigma Forecast
					with Bootstrap Error Bands\n (q: 5%, 25%, 75%, 95%)", cex.main = 0.7,
			ylab = "sigma", xlab = "T+i")
	lines(as.numeric(meansig), col = "black")
	lines(as.numeric(fs), col = colr)
	points(as.numeric(sigp[1,]), col = "steelblue1", pch = 19, cex = 0.5)
	points(as.numeric(sigp[2,]), col = "steelblue2", pch = 19, cex = 0.5)
	points(as.numeric(sigp[3,]), col = "steelblue3", pch = 19, cex = 0.5)
	points(as.numeric(sigp[4,]), col = "steelblue4", pch = 19, cex = 0.5)
	mtext(paste("GARCH model :",vmodel), side = 4, adj = 0, padj=0, col = "gray", 
			cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
				padj = 1.5, col = "gray", cex = 0.5)
	}
	legend.txt = c("forecast", "bootstrapped", namr)
	legend("topleft", legend = legend.txt, fill = c("red", "black", colr), col = c("red", "black", colr),
			bg = "n", bty = "n", cex = 0.6)
}

###############################################################################
# Distribution plots

distplot = function(distribution = "snorm", skewbounds = NULL, shapebounds = NULL, 
		n.points = NULL)
{
	old.par <- par(no.readonly = TRUE)
	on.exit(par(old.par))
	distribution = match.arg(distribution, c("snorm", "std", "ged","sstd","sged","ghst","nig","jsu"))
	if(is.null(skewbounds)){
		skewbounds = switch(distribution, 
				snorm = c(0.1, 10),
				std = c(0,0),
				ged = c(0,0),
				sstd = c(0.2, 5),
				sged = c(0.01, 5),
				ghst = c(-15, 15),
				nig =  c(-0.999, 0.999),
				jsu = c(-5, 5))
	}
	if(is.null(shapebounds)){
		shapebounds = switch(distribution, 
				snorm = c(0, 0),
				std = c(4.5, 35),
				ged = c(0.75, 10),
				sstd = c(4.5, 15),
				sged = c(0.5, 5),
				ghst = c(9.1, 15),
				nig =  c(0.15, 5),
				jsu = c(1, 2))
	}
	switch(distribution, 
			snorm = .snormplot(skewbounds, n.points),
			std = .stdplot(shapebounds, n.points),
			ged = .gedplot(shapebounds, n.points),
			sstd = .sstdplot(skewbounds, shapebounds, n.points),
			sged = .sgedplot(skewbounds, shapebounds, n.points),
			ghst = .ghstplot(skewbounds, shapebounds, n.points),
			nig =  .nigplot(skewbounds, shapebounds, n.points),
			jsu = .jsuplot(skewbounds, shapebounds, n.points))
	invisible()
}

.snormplot = function(skewbounds = c(0.1, 10), n.points = 1000){
	if(is.null(n.points)) n.points = 1000
	db = .DistributionBounds("snorm")
	if(skewbounds[1]<db$skew.LB){
		warning("\nskew lower bound below admissible region...adjusting to distribution lower bound.")
		skewbounds[1] = db$skew.LB
	}
	rng = seq(skewbounds[1], skewbounds[2], length.out = n.points)
	s = .snormskew(rng)
	plot(rng, s, type = "l", main = "Skew-Normal Distribution\nSkewness", ylab = "Skewness",
			xlab = "Skew Range", col = "steelblue")
	rect(0, -1, 1, 0, border = "black", lty = 2)
	grid()
	invisible()
}

.stdplot = function(shapebounds = c(4.5, 35), n.points = 1000){
	if(is.null(n.points)) n.points = 1000
	db = .DistributionBounds("std")
	if(shapebounds[1]<db$shape.LB){
		warning("\nshape lower bound below admissible region...adjusting to distribution lower bound.")
		shapebounds[1] = db$shape.LB
	}
	rng = seq(shapebounds[1], shapebounds[2], length.out = n.points)
	s = .stdexkurt(rng)
	plot(rng, s, type = "l", main = "Student Distribution\nKurtosis(Ex)", ylab = "Kurtosis",
			xlab = "Shape Range", col = "steelblue")
	grid()
	if(shapebounds[1]>4.1){
		rng2 = seq(4.1, shapebounds[1], length.out = 200)
		s2 = .stdexkurt(rng2)
		par(fig=c(0.5, 0.95, 0.5, 0.95), new = T)
		plot(rng2, s2, col = "steelblue", type = "l", ylab="", xlab="", cex.lab=0.7, cex.axis = 0.7)
	}
	invisible()
}

.gedplot = function(shapebounds = c(0.75, 10), n.points = 1000){
	if(is.null(n.points)) n.points = 1000
	db = .DistributionBounds("ged")
	if(shapebounds[1]<db$shape.LB){
		warning("\nshape lower bound below admissible region...adjusting to distribution lower bound.")
		shapebounds[1] = db$shape.LB
	}
	rng = seq(shapebounds[1], shapebounds[2], length.out = n.points)
	s = .gedexkurt(rng)
	plot(rng, s, type = "l", main = "GED Distribution\nKurtosis(Ex)", ylab = "Kurtosis",
			xlab = "Shape Range", col = "steelblue")
	grid()
	if(shapebounds[1]>0.25){
		rng2 = seq(0.25, shapebounds[1], length.out = 100)
		s2 = .gedexkurt(rng2)
		par(fig=c(0.5, 0.95, 0.5, 0.95), new = T)
		plot(rng2, s2, col = "steelblue", type = "l", ylab="", xlab="", cex.lab=0.7, cex.axis = 0.7)
	}
	invisible()
}


.sstdplot = function(skewbounds = c(0.2, 5), shapebounds = c(4.5, 15), n.points = 50){
	if(is.null(n.points)) n.points = 50
	db = .DistributionBounds("sstd")
	if(shapebounds[1]<db$shape.LB){
		warning("\nshape lower bound below admissible region...adjusting to distribution lower bound.")
		shapebounds[1] = db$shape.LB
	}
	if(skewbounds[1]<db$skew.LB){
		warning("\nskew lower bound below admissible region...adjusting to distribution lower bound.")
		skewbounds[1] = db$skew.LB
	}
	rngsk = seq(skewbounds[1], skewbounds[2], length.out = n.points)
	rngsh = seq(shapebounds[1], shapebounds[2], length.out = n.points)
	gr = expand.grid(rngsk, rngsh)
	sk = apply(gr, 1, function(x) .sstdskew(x[1], x[2]))
	ku = apply(gr, 1, function(x) .sstdexkurt(x[1], x[2]))
	zs = matrix(sk, n.points, n.points)
	zk = matrix(ku, n.points, n.points)
	
	par(mfrow = c(2,1), mar=c(2,2,2,2))
	x1 = .drapecol(zs, col = femme100(), NAcol = "white")	
	persp(  x = rngsk,
			y = rngsh,
			z = zs,  col = x1, theta = -50, phi = 25, expand = 0.5,
			ltheta = 120, shade = 0.75, ticktype = "detailed", xlab = "skew",
			ylab = "shape", zlab = "skewness",
			cex.lab = 0.7, cex.axis = 0.7,  cex.main = 0.8, main = "Skew-Student Skewness Surface")
	
	x1 = .drapecol(zk, col = femme100(), NAcol = "white")
	persp(  x = rngsk,
			y = rngsh,
			z = zk,  col = x1, theta = 50, phi = 25, expand = 0.5,
			ltheta = 120, shade = 0.75, ticktype = "detailed", xlab = "skew",
			ylab = "shape", zlab = "Kurtosis (Ex)",
			cex.lab = 0.7, cex.axis = 0.7,  cex.main = 0.8, main = "Skew-Student Kurtosis (Ex) Surface")
	
	invisible()
}

.sgedplot = function(skewbounds = c(0.01, 5), shapebounds = c(0.5, 5), n.points = 50){
	if(is.null(n.points)) n.points = 50
	db = .DistributionBounds("sged")
	if(shapebounds[1]<db$shape.LB){
		warning("\nshape lower bound below admissible region...adjusting to distribution lower bound.")
		shapebounds[1] = db$shape.LB
	}
	if(skewbounds[1]<db$skew.LB){
		warning("\nskew lower bound below admissible region...adjusting to distribution lower bound.")
		skewbounds[1] = db$skew.LB
	}
	rngsk = seq(skewbounds[1], skewbounds[2], length.out = n.points)
	rngsh = seq(shapebounds[1], shapebounds[2], length.out = n.points)
	gr = expand.grid(rngsk, rngsh)
	sk = apply(gr, 1, function(x) .sgedskew(x[1], x[2]))
	ku = apply(gr, 1, function(x) .sgedexkurt(x[1], x[2]))
	zs = matrix(sk, n.points, n.points)
	zk = matrix(ku, n.points, n.points)
	
	par(mfrow = c(2,1), mar=c(2,2,2,2))
	x1 = .drapecol(zs, col = femme100(), NAcol = "white")	
	persp(  x = rngsk,
			y = rngsh,
			z = zs,  col = x1, theta = -50, phi = 25, expand = 0.5,
			ltheta = 120, shade = 0.75, ticktype = "detailed", xlab = "skew",
			ylab = "shape", zlab = "skewness",
			cex.lab = 0.7, cex.axis = 0.7,  cex.main = 0.8, main = "Skew-GED Skewness Surface")
	
	x1 = .drapecol(zk, col = femme100(), NAcol = "white")
	persp(  x = rngsk,
			y = rngsh,
			z = zk,  col = x1, theta = 50, phi = 25, expand = 0.5,
			ltheta = 120, shade = 0.75, ticktype = "detailed", xlab = "skew",
			ylab = "shape", zlab = "Kurtosis (Ex)",
			cex.lab = 0.7, cex.axis = 0.7,  cex.main = 0.8, main = "Skew-GED Kurtosis (Ex) Surface")
	
	invisible()
}

.ghstplot = function(skewbounds = c(-15, 15), shapebounds = c(9.1, 15), n.points = 50){
	if(is.null(n.points)) n.points = 50
	db = .DistributionBounds("ghst")
	if(shapebounds[1]<db$shape.LB){
		warning("\nshape lower bound below admissible region...adjusting to distribution lower bound.")
		shapebounds[1] = db$shape.LB
	}
	if(skewbounds[1]<db$skew.LB){
		warning("\nskew lower bound below admissible region...adjusting to distribution lower bound.")
		skewbounds[1] = db$skew.LB
	}
	rngsk = seq(skewbounds[1], skewbounds[2], length.out = n.points)
	rngsh = seq(shapebounds[1], shapebounds[2], length.out = n.points)
	gr = expand.grid(rngsk, rngsh)
	sk = apply(gr, 1, function(x) .ghstskew(x[1], x[2]))
	ku = apply(gr, 1, function(x) .ghstexkurt(x[1], x[2]))
	zs = matrix(sk, n.points, n.points)
	zk = matrix(ku, n.points, n.points)
	
	par(mfrow = c(2,1), mar=c(2,2,2,2))
	x1 = .drapecol(zs, col = femme100(), NAcol = "white")	
	persp(  x = rngsk,
			y = rngsh,
			z = zs,  col = x1, theta = -50, phi = 25, expand = 0.5,
			ltheta = 120, shade = 0.75, ticktype = "detailed", xlab = "skew",
			ylab = "shape", zlab = "skewness",
			cex.lab = 0.7, cex.axis = 0.7,  cex.main = 0.8, main = "GHST Skewness Surface")
	
	x1 = .drapecol(zk, col = femme100(), NAcol = "white")
	persp(  x = rngsk,
			y = rngsh,
			z = zk,  col = x1, theta = 50, phi = 25, expand = 0.5,
			ltheta = 120, shade = 0.75, ticktype = "detailed", xlab = "skew",
			ylab = "shape", zlab = "Kurtosis (Ex)",
			cex.lab = 0.7, cex.axis = 0.7,  cex.main = 0.8, main = "GHST Kurtosis (Ex) Surface")
	
	invisible()
}

.nigplot = function(skewbounds = c(-0.999, 0.999), shapebounds = c(0.15, 5), n.points = 50){
	if(is.null(n.points)) n.points = 50
	db = .DistributionBounds("nig")
	if(shapebounds[1]<db$shape.LB){
		warning("\nshape lower bound below admissible region...adjusting to distribution lower bound.")
		shapebounds[1] = db$shape.LB
	}
	if(skewbounds[1]<db$skew.LB){
		warning("\nskew lower bound below admissible region...adjusting to distribution lower bound.")
		skewbounds[1] = db$skew.LB
	}
	rngsk = seq(skewbounds[1], skewbounds[2], length.out = n.points)
	rngsh = seq(shapebounds[1], shapebounds[2], length.out = n.points)
	gr = expand.grid(rngsk, rngsh)
	sk = apply(gr, 1, function(x) dskewness("nig", skew = x[1], shape=x[2]))
	ku = apply(gr, 1, function(x) dkurtosis("nig", skew = x[1], shape=x[2]))
	zs = matrix(sk, n.points, n.points)
	zk = matrix(ku, n.points, n.points)
	
	par(mfrow = c(2,1), mar=c(2,2,2,2))
	x1 = .drapecol(zs, col = femme100(), NAcol = "white")	
	persp(  x = rngsk,
			y = rngsh,
			z = zs,  col = x1, theta = -50, phi = 25, expand = 0.5,
			ltheta = 120, shade = 0.75, ticktype = "detailed", xlab = "skew",
			ylab = "shape", zlab = "skewness",
			cex.lab = 0.7, cex.axis = 0.7,  cex.main = 0.8, main = "NIG Skewness Surface")
	
	x1 = .drapecol(zk, col = femme100(), NAcol = "white")
	persp(  x = rngsk,
			y = rngsh,
			z = zk,  col = x1, theta = 50, phi = 25, expand = 0.5,
			ltheta = 120, shade = 0.75, ticktype = "detailed", xlab = "skew",
			ylab = "shape", zlab = "Kurtosis (Ex)",
			cex.lab = 0.7, cex.axis = 0.7,  cex.main = 0.8, main = "NIG Kurtosis (Ex) Surface")
	
	invisible()
}



.jsuplot = function(skewbounds = c(-5, 5), shapebounds = c(1, 2), n.points = 50){
	if(is.null(n.points)) n.points = 50
	db = .DistributionBounds("jsu")
	if(shapebounds[1]<db$shape.LB){
		warning("\nshape lower bound below admissible region...adjusting to distribution lower bound.")
		shapebounds[1] = db$shape.LB
	}
	if(skewbounds[1]<db$skew.LB){
		warning("\nskew lower bound below admissible region...adjusting to distribution lower bound.")
		skewbounds[1] = db$skew.LB
	}
	rngsk = seq(skewbounds[1], skewbounds[2], length.out = n.points)
	rngsh = seq(shapebounds[1], shapebounds[2], length.out = n.points)
	gr = expand.grid(rngsk, rngsh)
	sk = apply(gr, 1, function(x) .jsuskew(skew = x[1], shape = x[2]))
	ku = apply(gr, 1, function(x) .jsuexkurt(skew = x[1], shape = x[2]))
	zs = matrix(sk, n.points, n.points)
	zk = matrix(ku, n.points, n.points)
	
	par(mfrow = c(2,1), mar=c(2,2,2,2))
	x1 = .drapecol(zs, col = femme100(), NAcol = "white")	
	persp(  x = rngsk,
			y = rngsh,
			z = zs,  col = x1, theta = -50, phi = 25, expand = 0.5,
			ltheta = 120, shade = 0.75, ticktype = "detailed", xlab = "skew",
			ylab = "shape", zlab = "skewness",
			cex.lab = 0.7, cex.axis = 0.7,  cex.main = 0.8, main = "JSU Skewness Surface")
	
	x1 = .drapecol(zk, col = femme100(), NAcol = "white")
	persp(  x = rngsk,
			y = rngsh,
			z = zk,  col = x1, theta = 50, phi = 25, expand = 0.5,
			ltheta = 120, shade = 0.75, ticktype = "detailed", xlab = "skew",
			ylab = "shape", zlab = "Kurtosis (Ex)",
			cex.lab = 0.7, cex.axis = 0.7,  cex.main = 0.8, main = "JSU Kurtosis (Ex) Surface")
	
	invisible()
}

# Skewness-Kurtosis Authorized Domain Plots
skdomain = function(distribution = "nig", kurt.max = 30, n.points = 25, lambda = 1, plot = TRUE, legend = NULL)
{
	sol = switch(distribution,
			nig  = .nigdomain(kurt.max, n.points),
			hyp  = .hypdomain(kurt.max, n.points),
			ghyp = .ghypdomain(kurt.max, n.points, lambda = lambda),
			jsu  = .jsudomain(kurt.max, n.points),
			sstd  = .sstddomain(kurt.max, n.points))
			#sged = .sgeddomain(kurt.max, n.points))
	if(plot){
		K = c(seq(1.01, kurt.max, length.out=n.points), seq(1.01, kurt.max, length.out=n.points))
		S = c(sqrt(seq(1.01, kurt.max, length.out=n.points)-1),-sqrt(seq(1.01, kurt.max, length.out=n.points)-1))
		plot(K, S, ylab="Skewness", xlab = "Kurtosis", main = "Authorized Domain", type="n")
		lines(spline(K[1:n.points],  S[1:n.points]), col = "red", lwd = 2)
		lines(spline(K[1:n.points], -S[1:n.points]), col = "red" ,lwd = 2)
		lines(sol$Kurtosis, sol$Skewness, col = "steelblue", lwd = 2)
		lines(sol$Kurtosis, -sol$Skewness, col = "steelblue", lwd = 2)
		if(is.null(legend)){
			legend("topleft", c("MAX", distribution), col = c("red", "steelblue"),  lty = c(1,2), 
					lwd = rep(2, 2), bty = "n")
		}
		abline(v=3, col = "lightgrey")
	}
	return(sol)
}

.nigdomain = function(kurt.max = 30, n.points = 25) 
{
	di = .DistributionBounds("nig")
	di$shape.UB = 100
	k = seq(5, kurt.max, length = n.points)
	maxkurt = dkurtosis("nig", skew = 0, shape=di$shape.UB)
	f = function(x, kurt){
		-dskewness("nig", skew = x[1], shape = x[2])
	}
	fin = function(x, kurt){
		dkurtosis("nig", skew = x[1], shape = x[2])+3-maxkurt - kurt
	}
	parsx = matrix(NA, ncol = 4, nrow = n.points)
	for(i in 1:length(k)){
		sol = solnp(pars=c(0.1, 0.5), fun = f, eqfun = fin, eqB=0, LB = c(0.05, di$shape.LB),
				UB = c(di$skew.UB, di$shape.UB), control=list(trace=0, outer.iter=25), kurt = k[i])
		parsx[i,1:2] = sol$pars
		parsx[i,3]   = tail(sol$value,1)
		parsx[i,4]   = k[i]
	}
	parsx = rbind(matrix(c(0,100, dskewness("nig", 0, 100), 3+dkurtosis("nig", 0, 100)), ncol = 4), parsx)
	ans = spline(parsx[,4],   parsx[,3])
	return(list(Kurtosis = ans$x, Skewness = ans$y))
}

.hypdomain = function(kurt.max = 30, n.points = 25) 
{
	di = .DistributionBounds("ghyp")
	di$shape.UB = 100
	k = seq(5, kurt.max, length = n.points)
	maxkurt = dkurtosis("ghyp", skew = 0, shape=di$shape.UB, lambda=1)
	f = function(x, kurt){
		-dskewness("ghyp", skew = x[1], shape = x[2], lambda=1)
	}
	fin = function(x, kurt){
		dkurtosis("ghyp", skew = x[1], shape = x[2], lambda=1)+3-maxkurt - kurt
	}
	parsx = matrix(NA, ncol = 4, nrow = n.points)
	for(i in 1:length(k)){
		sol = solnp(pars=c(0.1, 0.5), fun = f, eqfun = fin, eqB=0, LB = c(0.05, di$shape.LB),
				UB = c(di$skew.UB, di$shape.UB), control=list(trace=0, outer.iter=25), kurt = k[i])
		parsx[i,1:2] = sol$pars
		parsx[i,3]   = tail(sol$value,1)
		parsx[i,4]   = dkurtosis("ghyp", skew =  sol$pars[1], shape =  sol$pars[2], lambda=1)+3
	}
	parsx = rbind(matrix(c(0,100, dskewness("ghyp", 0, 100, lambda=1), 3+dkurtosis("ghyp", 0, 100, lambda=1)), ncol = 4), parsx)
	ans = spline(parsx[,4], parsx[,3])
	return(list(Kurtosis = ans$x, Skewness = ans$y))
}

.ghypdomain = function(kurt.max = 30, n.points = 25, lambda = 1) 
{
	di = .DistributionBounds("ghyp")
	di$shape.UB = 100
	k = seq(5, kurt.max, length = n.points)
	maxkurt = dkurtosis("ghyp", skew = 0, shape=di$shape.UB, lambda=lambda)
	f = function(x, kurt, ghlambda){
		-dskewness("ghyp", skew = x[1], shape = x[2], lambda=ghlambda)
	}
	fin = function(x, kurt, ghlambda){
		dkurtosis("ghyp", skew = x[1], shape = x[2], lambda=ghlambda)+3-maxkurt - kurt
	}
	parsx = matrix(NA, ncol = 4, nrow = n.points)
	for(i in 1:length(k)){
		sol = solnp(pars=c(0.1, 0.5), fun = f, eqfun = fin, eqB=0, LB = c(0.05, di$shape.LB),
				UB = c(di$skew.UB, di$shape.UB), control=list(trace=0, outer.iter=25), kurt = k[i],
				ghlambda = lambda)
		parsx[i,1:2] = sol$pars
		parsx[i,3]   = tail(sol$value,1)
		parsx[i,4]   = dkurtosis("ghyp", skew =  sol$pars[1], shape =  sol$pars[2], lambda=lambda)+3
	}
	parsx = rbind(matrix(c(0,100, dskewness("ghyp", 0, 100, lambda=lambda), 3+dkurtosis("ghyp", 0, 100, lambda=lambda)), ncol = 4), parsx)
	ans = spline(parsx[,4], parsx[,3])
	return(list(Kurtosis = ans$x, Skewness = ans$y))
}


.sstddomain = function(kurt.max = 30, n.points = 25) 
{
	di = .DistributionBounds("sstd")
	di$shape.UB = 100
	k = seq(5, kurt.max, length = n.points)
	maxkurt = dkurtosis("sstd", skew = 1, shape=di$shape.UB)
	f = function(x, kurt){
		-dskewness("sstd", skew = x[1], shape = x[2])
	}
	fin = function(x, kurt){
		dkurtosis("sstd", skew = x[1], shape = x[2])+3-maxkurt - kurt
	}
	parsx = matrix(NA, ncol = 4, nrow = n.points)
	for(i in 1:length(k)){
		sol = solnp(pars=c(1.1, 4.3), fun = f, eqfun = fin, eqB=0, LB = c(1, 4.1),
				UB = c(di$skew.UB, di$shape.UB), control=list(trace=0, outer.iter=25), kurt = k[i])
		parsx[i,1:2] = sol$pars
		parsx[i,3]   = tail(sol$value,1)
		parsx[i,4]   = k[i]
	}
	parsx = rbind(matrix(c(0,100, dskewness("sstd", 1, 100), 3+dkurtosis("sstd", 1, 100)), ncol = 4), parsx)
	ans = spline(parsx[,4], parsx[,3])
	return(list(Kurtosis = ans$x, Skewness = ans$y))
}

.jsudomain = function(kurt.max = 30, n.points = 25) 
{
	di = .DistributionBounds("jsu")
	di$shape.UB = 100
	k = seq(5, kurt.max, length = n.points)
	maxkurt = dkurtosis("jsu", skew = 0, shape=di$shape.UB)
	f = function(x, kurt){
		-dskewness("jsu", skew = x[1], shape = x[2])
	}
	fin = function(x, kurt){
		dkurtosis("jsu", skew = x[1], shape = x[2])+3-maxkurt - kurt
	}
	parsx = matrix(NA, ncol = 4, nrow = n.points)
	for(i in 1:length(k)){
		sol = solnp(pars=c(0.1, 0.5), fun = f, eqfun = fin, eqB=0, LB = c(0.05, di$shape.LB),
				UB = c(di$skew.UB, di$shape.UB), control=list(trace=0, outer.iter=25), kurt = k[i])
		parsx[i,1:2] = sol$pars
		parsx[i,3]   = tail(sol$value,1)
		parsx[i,4]   = k[i]
	}
	parsx = rbind(matrix(c(0,100, dskewness("jsu", 0, 100), 3+dkurtosis("jsu", 0, 100)), ncol = 4), parsx)
	ans = spline(parsx[,4], parsx[,3])
	return(list(Kurtosis = ans$x, Skewness = ans$y))
}

################################################################################
VaRplot = function(alpha, actual, VaR, title = paste("Daily Returns and Value-at-Risk \nExceedances\n","(alpha=", alpha,")",sep=""),
		ylab = "Daily Log Returns", xlab = "Time")
{	
	period = diff(index(actual))
	# intraday
	if(attr(period, "units") == "mins"){
		A = as.numeric(actual)
		V = as.numeric(VaR)
		ep <- axTicksByTime(actual)
		plot(A, type = "n", main = title, ylab = ylab, xlab = xlab, 
				ylim = c(min(A, V), max(A, V)), ann = FALSE, xaxt = "n",
				cex.main = 0.8, cex.lab = 0.9, cex.axis = 0.8)
		axis(1, at = ep, labels = names(ep), tick = TRUE)
		abline(h = 0, col = "grey", lty = 2)
		points(A, pch = 18, col = "lightgrey")
		sel  =  which(A<V)
		Ap = rep(NA, length(A))
		Ap[sel] = A[sel]
		lines(V, lwd = 1, col = "black")
		points(Ap, pch = 3, cex = 1.1, col = "red")
		
		legend("topleft", max(A),c("returns","return < VaR","VaR"),
				col=c("lightgrey", "red","black"), cex=0.75,
				pch = c(18,3,-1), lty=c(0,0,1), lwd=c(0,0,2), bty = "n")
		grid()
	} else{
		plot(index(actual), as.numeric(actual), type = "n", main = title, ylab = ylab, xlab = xlab, 
				ylim = c(min(actual, VaR), max(actual, VaR)), ann = FALSE, 
				cex.main = 0.8, cex.lab = 0.9, cex.axis = 0.8)
		abline(h = 0, col = "grey", lty = 2)
		points(index(actual), as.numeric(actual), pch = 18, col = "lightgrey")
		sel  =  which(actual<VaR)
		points(index(actual)[sel],as.numeric(actual[sel]), pch = 18, col = "red")
		lines(index(actual), as.numeric(VaR), lwd = 2, col = "black")
		legend("topleft", max(actual),c("returns","return < VaR","VaR"),
				col=c("lightgrey", "red","black"), cex=0.75,
				pch = c(18,18,-1), lty=c(0,0,1), lwd=c(0,0,2), bty = "n")
		grid()
	}
	return(invisible())
}