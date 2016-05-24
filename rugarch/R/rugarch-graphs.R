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
.plot.acf <-
		function (x, ci = 0.95, type = "h", xlab = "Lag", ylab = NULL,
				ylim = NULL, main = NULL, ci.col="blue",
				ci.type = c("white", "ma"),
				max.mfrow = 6,
				ask = Npgs > 1 && dev.interactive(),
				mar = if(nser > 2) c(3,2,2,0.8) else par("mar"),
				oma = if(nser > 2) c(1,1.2,1,1) else par("oma"),
				mgp = if(nser > 2) c(1.5,0.6,0) else par("mgp"),
				xpd = par("xpd"),
				cex.main = if(nser > 2) 1 else par("cex.main"),
				verbose = getOption("verbose"),
				...)
{
	ci.type = match.arg(ci.type)
	if((nser = ncol(x$lag)) < 1) stop("x$lag must have at least 1 column")
	if (is.null(ylab))
		ylab = switch(x$type,
				correlation = "ACF",
				covariance = "ACF (cov)",
				partial = "Partial ACF")
	snames = x$snames
	if (is.null(snames))
		snames = paste("Series ", if (nser == 1) x$series else 1:nser)
	
	with.ci = ci > 0 && x$type != "covariance"
	with.ci.ma = with.ci && ci.type == "ma" && x$type == "correlation"
	if(with.ci.ma && x$lag[1,1,1] != 0) {
		warning("can use ci.type=\"ma\" only if first lag is 0")
		with.ci.ma = FALSE
	}
	clim0 = if (with.ci) qnorm((1 + ci)/2)/sqrt(x$n.used) else c(0, 0)
	
	Npgs = 1 ## we will do [ Npgs x Npgs ] pages !
	nr = nser
	if(nser > 1) { ## at most m x m (m := max.mfrow)  panels per page
		sn.abbr = if(nser > 2) abbreviate(snames) else snames
		
		if(nser > max.mfrow) {
			##  We need more than one page: The plots are laid out
			##  such that we can manually paste the paper pages and get a
			##  nice square layout with diagonal !
			## NB: The same applies to pairs() where we'd want several pages
			Npgs = ceiling(nser / max.mfrow)
			nr = ceiling(nser / Npgs)  # <= max.mfrow
		}
		opar = par(mfrow = rep(nr, 2), mar = mar, oma = oma, mgp = mgp,
				ask = ask, xpd = xpd, cex.main = cex.main)
		on.exit(par(opar))
		if(verbose) { # FIXME: message() can be suppressed but not str()
			message("par(*) : ", appendLF = FALSE, domain = NA)
			str(par("mfrow","cex", "cex.main","cex.axis","cex.lab","cex.sub"))
		}
	}
	
	if (is.null(ylim)) {
		## Calculate a common scale
		ylim = range(x$acf[, 1:nser, 1:nser], na.rm = TRUE)
		if (with.ci) ylim = range(c(-clim0, clim0, ylim))
		if (with.ci.ma) {
			for (i in 1:nser) {
				clim = clim0 * sqrt(cumsum(c(1, 2*x$acf[-1, i, i]^2)))
				ylim = range(c(-clim, clim, ylim))
			}
		}
	}
	
	for (I in 1:Npgs) for (J in 1:Npgs) {
			## Page [ I , J ] : Now do   nr x nr  'panels' on this page
			iind <- (I-1)*nr + 1:nr
			jind <- (J-1)*nr + 1:nr
			if(verbose)
				message("Page [",I,",",J,"]: i =",
						paste(iind,collapse=","),"; j =",
						paste(jind,collapse=","), domain = NA)
			for (i in iind) for (j in jind)
					if(max(i,j) > nser) {
						frame(); box(col = "light gray")
						## the above is EXTREMELY UGLY; should have a version
						## of frame() that really does advance a frame !!
					}
					else {
						clim <- if (with.ci.ma && i == j)
									clim0 * sqrt(cumsum(c(1, 2*x$acf[-1, i, j]^2))) else clim0
						plot(x$lag[, i, j], x$acf[, i, j], type = type, xlab = xlab,
								ylab = if(j==1) ylab else "", ylim = ylim, ...)
						abline(h = 0)
						if (with.ci && ci.type == "white")
							abline(h = c(clim, -clim), col = ci.col, lty = 2)
						else if (with.ci.ma && i == j) {
							clim <- clim[-length(clim)]
							lines(x$lag[-1, i, j], clim, col = ci.col, lty = 2)
							lines(x$lag[-1, i, j], -clim, col = ci.col, lty = 2)
						}
						title(if (!is.null(main)) main else
										if (i == j) snames[i]
										else paste(sn.abbr[i], "&", sn.abbr[j]),
								line = if(nser > 2) 1 else 2)
					}
			if(Npgs > 1) {                  # label the page
				mtext(paste("[",I,",",J,"]"), side=1, line = -0.2, adj=1,
						col = "dark gray", cex = 1, outer = TRUE)
			}
		}
}

# Quantile-Quantile Plot
.qqDist = function (y, dist = "norm", ylim = NULL, main = paste(dist, "- QQ Plot"),
		xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", doplot = TRUE,
		datax = FALSE, cex.main = 0.8, ...)
{	
	y = as.vector(y)
	if (has.na <- any(ina <- is.na(y))) {
		yN = y
		y = y[!ina]
	}
	if (0 == (n <- length(y))) stop("y is empty or has only NAs")
	x = qdist(distribution = dist, p = ppoints(n), ...)[order(order(y))]
	if (has.na) {
		y = x
		x = yN
		x[!ina] = y
		y = yN
	}
	if (doplot) {
		if (is.null(ylim)) ylim = range(y)
		if (datax) {
			plot(y, x, main = main, xlab = ylab, ylab = xlab, xlim = ylim,
					col = "steelblue", cex = 0.7)
		} else {
			plot(x, y, main = main, xlab = xlab, ylab = ylab, ylim = ylim,
					col = "steelblue", cex = 0.7)
		}
		.qqLine(y = y, dist = dist, datax = datax, ...)
		grid()
	}
	invisible(if (datax) list(x = y, y = x) else list(x = x, y = y))
}

# ------------------------------------------------------------------------------
.qqLine = function (y, dist = "norm", datax = FALSE, ...)
{   
	y = as.vector(y)
	y = quantile(y[!is.na(y)], c(0.25, 0.75))
	x = qdist(distribution = dist, p = c(0.25, 0.75), ...)
	if (datax) {
		slope = diff(x)/diff(y)
		int = x[1] - slope * y[1]
	} else {
		slope = diff(y)/diff(x)
		int = y[1] - slope * x[1]
	}
	abline(int, slope)
}

#----------------------------------------------------------------------------------
# used in choosing the mfrow
.divisortable = function(n)
{
	z=matrix(c(1,1,1,
					2,2,1,
					3,2,2,
					4,2,2,
					5,2,3,
					6,2,3,
					7,2,4,
					8,2,4,
					9,3,3,
					10,3,4,
					11,3,4,
					12,3,4,
					13,4,4,
					14,4,4,
					15,4,4,
					16,4,4,
					17,4,5,
					18,4,5,
					19,4,5,
					20,4,5), ncol = 3, byrow = TRUE)
	d = which(n==z[,1])
	return(z[d,2:3])
}

# ar + ma (roots) + dist --> 1 + 1 plots x max(armap, armaq)
# alpha + beta (dist) --> 1 plot x max(garchp, garchq)
# alpha + gamma1 --> 1 plot x garchp
# alpha + gamma2 --> 1 plot x garchp
# lambda + delta (fgarch) --> 1 plot
# skew + shape --> 1 plot
# max 7 plots

.modelcomb = function(x)
{
	countm = 0
	modeltype = 1
	asymm1 = c(0,0)
	asymm2 = c(0,0)
	power1 = c(0,0)
	power2 = c(0,0)
	vmodel = x@model$modeldesc$vmodel
	vsubmodel = x@model$modeldesc$vsubmodel
	model = x@model
	modelinc = model$modelinc
	if(modelinc[2]!=0 && modelinc[3]!=0) countm = countm + 1
	if(modelinc[8]!=0 && modelinc[9]!=0){
		countm = countm + max(modelinc[8:9])
	}
	Names = names(coef(x))
	if(vmodel!="fGARCH" && any(substr(Names, 1, 5)=="gamma"))
	{
		i = which(substr(Names, 1, 5)=="gamma")
		countm = countm + length(i)
		asymm1[1] = 1
		asymm1[2] = length(i)
	}
	if(vmodel=="fGARCH" && any(substr(Names, 1, 6)=="eta1"))
	{
		i = which(substr(Names, 1, 4)=="eta1")
		countm = countm + length(i)
		asymm1[1] = 1
		asymm1[2] = length(i)
		modeltype = 2
	}
	if(vmodel=="fGARCH" && any(substr(Names, 1, 6)=="eta2"))
	{
		i = which(substr(Names, 1, 4)=="eta2")
		countm = countm + length(i)
		asymm2[1] = 1
		asymm2[2] = length(i)
		modeltype = 2
	}
	if(vmodel!="fGARCH" && any(substr(Names, 1, 5)=="delta"))
	{
		i = which(substr(Names, 1, 5)=="delta")
		countm = countm + 1
		power1[1] = 1
		power1[2] = 1
	}
	if(vmodel=="fGARCH" && any(substr(Names, 1, 5)=="delta"))
	{
		i = which(substr(Names, 1, 5)=="delta")
		countm = countm + 1
		power1[1] = 1
		power1[2] = 1
		modeltype = 2
	}
	if(vmodel=="fGARCH" && any(substr(Names, 1, 6)=="lambda"))
	{
		i = which(substr(Names, 1, 6)=="lambda")
		countm = countm + 1
		power2[1] = 1
		power2[2] = 1
		modeltype = 2
	}
}

#------------------------------------------------------------------------------
# uGARCHdistribution : Models relating to the Mean Equation:
#------------------------------------------------------------------------------
.armaroots = function(coefs)
{
	Names = names(coefs)
	ar = ma = NULL
	if(any(substr(Names, 1, 2)=="ar"))
	{
		ar = which(substr(Names, 1, 2)=="ar")
		armap = length(ar)
		arcoef = coefs[ar]
		zar = polyroot(c(1,-arcoef))
		rezar = Re(zar)
		imzar = Im(zar)
		nmr = paste("ar", 1:armap, sep="")
	} else{
		zar = NULL
		rezar = NULL
		imzar = NULL
		nmr = NULL
	}
	if(any(substr(Names, 1, 2)=="ma"))
	{
		ma = which(substr(Names, 1, 2)=="ma")
		armaq = length(ma)
		macoef = coefs[ma]
		zma = polyroot(c(1,macoef))
		rezma = Re(zma)
		imzma = Im(zma)
		nma = paste("ma", 1:armaq, sep="")
	} else{
		zma = NULL
		rezma = NULL
		imzma = NULL
		nma = NULL
	}
	root = list()
	root$ar = zar
	root$ma = zma
	root$realar = rezar
	root$imagar = imzar
	root$realma = rezma
	root$imagma = imzma
	amp = list()
	atan = list()
	degree = list()
	if(!is.null(zar)){
		amp$ar = apply(cbind(root$realar, root$imagar), 1, FUN = function(x) sqrt(x[1]^2 + x[2]^2))
		atan$ar = apply(cbind(amp$ar, root$realar),1 , FUN = function(x) atan2(x[1], x[2]) )
		degree$ar = atan$ar * 57.29577951
	}
	if(!is.null(zma)){
		amp$ma = apply(cbind(root$realma, root$imagma), 1, FUN = function(x) sqrt(x[1]^2 + x[2]^2))
		atan$ma = apply(cbind(amp$ma, root$realma),1 , FUN = function(x) atan2(x[1], x[2]) )
		degree$ma = atan$ma * 57.29577951
	}
	res = list(root = root, amp = amp, atan = atan, deg = degree)
	return(res)
}

.plotarmaroots = function(x, ...)
{
	arroot = x$root$ar
	ar.rr = x$root$realar
	ar.im = x$root$imagar
	
	maroot = x$root$ma
	ma.rr = x$root$realma
	ma.im = x$root$imagma
	
	slimit = 1/max(abs(c(ar.rr, ma.rr)), 1.5, abs(c(ar.im, ma.im)))
	if(!is.null(arroot)){
		plot(1/arroot, xlim = c( - 1, 1), ylim = c( - 1, 1), xlab = "", ylab = "", pch = 23, ...)
		if(!is.null(maroot)) points(1/maroot, pch = 21, ...)
		x = (2*pi/360)*(0:360)
		lines(sin(x), cos(x), col = "darkgrey")
		abline(h = 0, col = "darkgrey")
		abline(v = 0, col = "darkgrey")
		title("Inverse Roots and Unit Circle\n", 
				xlab = "Real Part", ylab = "Imaginary Part")
	} else{
		if(!is.null(maroot)){
			plot(1/maroot, xlim = c( - 1, 1), ylim = c( - 1, 1), xlab = "", ylab = "", pch = 23, ...)
			x = (2*pi/360)*(0:360)
			lines(sin(x), cos(x), col = "darkgrey")
			abline(h = 0, col = "darkgrey")
			abline(v = 0, col = "darkgrey")
			title("Inverse Roots and Unit Circle\n", 
					xlab = "Real Part", ylab = "Imaginary Part")
		}
	}
	invisible(x)
}


.pointsarmaroots = function(x, ...)
{
	arroot = x$root$ar
	ar.rr = x$root$realar
	ar.im = x$root$imagar
	maroot = x$root$ma
	ma.rr = x$root$realma
	ma.im = x$root$imagma	
	if(!is.null(arroot)) points(1/arroot, pch = 23,  ...)
	if(!is.null(maroot)) points(1/maroot, pch = 21, ...)
	invisible(x)
}


.arma2dplot = function(x, window, ...)
{
	vmodel = x@model$modeldesc$vmodel
	vsubmodel = x@model$modeldesc$vsubmodel
	model = x@model
	modelinc = model$modelinc
	truecf = unlist(x@truecoef[,1])
	Names = names(truecf)
	md = as.data.frame(x, window = window)
	md = na.omit(md)
	ar = ma = NULL
	if(any(substr(Names, 1, 2)=="ar"))
	{
		ar = which(substr(Names, 1, 2)=="ar")
		armap = length(ar)
	}
	if(any(substr(Names, 1, 2)=="ma"))
	{
		ma = which(substr(Names, 1, 2)=="ma")
		armaq = length(ma)
	}
	if(is.null(ar) | is.null(ma)){
		warning("\nNo plots for mean equation...\n")
		return()
	}
	
	for(i in 1:armap){
		for(j in 1:armaq){
			xlim = c(min(md[,ar[j]], na.rm = TRUE), max(md[,ar[j]], na.rm = TRUE))
			ylim = c(min(md[,ma[j]], na.rm = TRUE), max(md[,ma[j]], na.rm = TRUE))
			plot(md[,ar[j]], md[,ma[j]], ylab = paste("ar",i,sep=""), xlab = paste("ma",j,sep=""), xlim = xlim,
					ylim = ylim)
			f1 = .kde2d(md[,ar[j]], md[,ma[j]])
			image(f1, ylab = paste("ar",i,sep=""), xlab = paste("ma",j,sep=""), add =TRUE, xlim = xlim,
					ylim = ylim)
			contour(f1,col = "steelblue", add = TRUE, method = "edge",
				vfont = c("sans serif", "plain"), xlim = xlim,
				ylim = ylim)
		title(paste("ARMA(",i,",",j,")", sep=""))
		mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", 
				cex = 0.5)
		if(vmodel == "fGARCH"){
			mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
					padj = 1.5, col = "gray", cex = 0.5)
		}
	}}
	invisible(x)
}


.plotarmamodel = function(x, window, ...)
{
	vmodel = x@model$modeldesc$vmodel
	vsubmodel = x@model$modeldesc$vsubmodel
	model = x@model
	modelinc = model$modelinc
	truecf = unlist(x@truecoef[,1])
	Names = names(truecf)
	md = as.data.frame(x, window = window)
	md = na.omit(md)
	ar = ma = NULL
	nar = nma = NULL
	if(any(substr(Names, 1, 2)=="ar"))
	{
		ar = which(substr(Names, 1, 2)=="ar")
		armap = length(ar)
		nar = paste("ar", 1:armap,sep="")
	}
	if(any(substr(Names, 1, 2)=="ma"))
	{
		ma = which(substr(Names, 1, 2)=="ma")
		armaq = length(ma)
		nma = paste("ma", 1:armaq,sep="")
	}
	if(is.null(ar) && is.null(ma)){
		warning("\nNo 2dplots for mean equation...\n")
		return(0)
	}
	start = dev.next(which = dev.cur())
	dev.new(start+1)
	maxn = 2 * armap*armaq
	dv = .divisortable(maxn)
	coefarma = md[,c(ar,ma)]
	N = dim(coefarma)[1]
	colnames(coefarma) = c(nar, nma)
	rootlist = apply(coefarma, 1, FUN=function(x) .armaroots(x))
	colx = topo.colors(N, alpha = 0.7)
	par(mfrow = c(dv[1], dv[2]))
	.arma2dplot(x, window = window)
	.plotarmaroots(x = rootlist[[1]], col="steelblue")
	mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", 
			cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
				padj = 1.5, col = "gray", cex = 0.5)
	}
	for(i in 2:N) .pointsarmaroots(x = rootlist[[i]], col = colx[i])
	return(1)
}

#------------------------------------------------------------------------------
# uGARCHdistribution : Models relating to the Variance Equation and Density:
#------------------------------------------------------------------------------
.dist2dplot = function(x, window, ...)
{
	vmodel = x@model$modeldesc$vmodel
	vsubmodel = x@model$modeldesc$vsubmodel
	model = x@model
	modelinc = model$modelinc
	idx = model$pos.matrix
	truecf = unlist(x@truecoef[,1])
	Names = names(truecf)
	md = as.data.frame(x, window = window)
	md = na.omit(md)
	distribution = model$modeldesc$distribution
		
	if(vmodel == "iGARCH"){
		warning("No bivariate GARCH plots for iGARCH model (in garch2dplot call)...")
	}
	if(modelinc[8] == 0 | modelinc[9] == 0 && (vmodel == "sGARCH")){
		warning("\nNo plots for variance equation...\n")
	}
	if(vmodel == "csGARCH") ind = 0 else ind = 1
	total = modelinc[8]*modelinc[10] + modelinc[8]*modelinc[11]*ind + modelinc[8]*modelinc[12]*ind +
			modelinc[8]*modelinc[9] + modelinc[13]*modelinc[14] + modelinc[16]*modelinc[17] + 
			modelinc[16]*modelinc[10] + modelinc[16]*modelinc[11]*ind + modelinc[16]*modelinc[12]*ind
	
	dv = .divisortable(2+total)
	par(mfrow=c(dv[1], dv[2]))
	if(sum(modelinc[2:3])>0){	
		coefarma = md[,idx["ar",1]:idx["ma",2]]
		N = dim(coefarma)[1]
		colnames(coefarma) = colnames(md)[idx["ar",1]:idx["ma",2]]
		rootlist = apply(coefarma, 1, FUN=function(x) .armaroots(x))
		colx = topo.colors(N, alpha = 0.7)
		.arma2dplot(x, window = window)
		.plotarmaroots(x = rootlist[[1]], col="steelblue")
		mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", 
				cex = 0.7)
		if(vmodel == "fGARCH"){
			mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
					padj = 1.5, col = "gray", cex = 0.5)
		}
		for(i in 2:N) .pointsarmaroots(x = rootlist[[i]], col = colx[i])
	}
	
	if(modelinc[8]>0)
	{
		for(i in 1:modelinc[8]){
			if(modelinc[10]>0){
				tmp = c(idx["alpha",1]+i-1, idx["gamma",1]+i-1)
				plot(md[,tmp], 
						xlab = paste("alpha",i,sep=""), ylab = paste("gamma1-",i,sep=""))
				f1 = .kde2d(md[,tmp[1]], md[,tmp[2]])
				image(f1,  xlab = paste("alpha",i), ylab = paste("gamma1-",i), add =TRUE)
				contour(f1, col = "steelblue", add = TRUE, method = "edge",
						vfont = c("sans serif", "plain"))
				title(paste("Shock vs Asymmetry(1)", sep = ""), cex.main = 0.7)
				mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", 
						cex = 0.5)
			}
			if(modelinc[11]>0 && ind>0){
				tmp = c(idx["alpha",1]+i-1, idx["eta1",1]+i-1)
				plot(md[,tmp], xlab = paste("alpha",i,sep=""), ylab = paste("eta1-",i,sep=""))
				f1 = .kde2d(md[,tmp[1]], md[,tmp[2]])
				image(f1, xlab = paste("alpha",i), ylab = paste("eta1-",i), add =TRUE)
				contour(f1,col = "steelblue", add = TRUE, method = "edge",
						vfont = c("sans serif", "plain"))
				title(paste("Shock vs Shift", sep=""), cex.main = 0.7)
				mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", 
						cex = 0.5)
				mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
							padj = 1.5, col = "gray", cex = 0.5)
			}
			if(modelinc[12]>0 && ind>0){
				tmp = c(idx["alpha",1]+i-1, idx["eta2",1]+i-1)
				plot(md[,tmp], xlab = paste("alpha",i,sep=""), ylab = paste("eta2-",i,sep=""))
				f1 = .kde2d(md[,tmp[1]], md[,tmp[2]])
				image(f1, xlab = paste("alpha",i), ylab = paste("eta2-",i), add =TRUE)
				contour(f1,col = "steelblue", add = TRUE, method = "edge",
						vfont = c("sans serif", "plain"))
				title(paste("Shock vs Rotation", sep=""))
				mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", 
						cex = 0.5)
				mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
							padj = 1.5, col = "gray", cex = 0.5)
			}
			if(modelinc[9]>0){
					for(j in 1:modelinc[9]){
						tmp = c(idx["alpha",1]+i-1, idx["beta",1]+j-1)
						plot(md[,tmp], xlab = paste("alpha",i,sep=""), ylab = paste("beta",j,sep=""))
						f1 = .kde2d(md[,tmp[1]], md[,tmp[2]])
						image(f1, xlab = paste("alpha",i,sep=""), ylab = paste("beta",j,sep=""), add =TRUE)
						contour(f1,col = "steelblue", add = TRUE, method = "edge",
								vfont = c("sans serif", "plain"))
						title(paste("Shock vs Persistence", sep=""))
						mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", 
								cex = 0.5)
						if(vmodel == "fGARCH"){
							mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
									padj = 1.5, col = "gray", cex = 0.5)
						}
				}
			}
		}
	}
	if(modelinc[13]>0 && modelinc[14]>0) 
	{
		tmp = c(idx["lambda",1], idx["delta",1])
		plot(md[,tmp], xlab = paste("Power(1)"), ylab = paste("Power(2)"))
		f1 = .kde2d(md[,tmp[1]], md[,tmp[2]])
		image(f1, xlab = paste("Power(1)",i), ylab = paste("Power(2)",i), add =TRUE)
		contour(f1,col = "steelblue", add = TRUE, method = "edge",
				vfont = c("sans serif", "plain"))
		title(paste("Power(1) vs Power(2)", sep=""))
		mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", 
				cex = 0.5)
		if(vmodel == "fGARCH"){
			mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
					padj = 1.5, col = "gray", cex = 0.5)
		}
	}

	if(modelinc[16]>0 && modelinc[17]>0){
		plot(md[,idx["skew",1]], md[,idx["shape",1]], ylab = "shape", xlab = "skew")
		f1 = .kde2d(md[,idx["skew",1]], md[,idx["shape",1]])
		image(f1, ylab = "shape", xlab = "skew", add =TRUE)
		contour(f1,col = "steelblue", add = TRUE, method = "edge",
				vfont = c("sans serif", "plain"))
		title(paste("Skew vs Shape Parameter (",distribution,")", sep=""))
		mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", 
				cex = 0.5)
		if(vmodel == "fGARCH"){
			mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
					padj = 1.5, col = "gray", cex = 0.5)
		}
	}
	if(modelinc[16]>0 && modelinc[10]>0){
		for(i in 1:modelinc[10]){
			plot(md[,idx["skew",1]], md[,idx["gamma",1]+i-1], xlab = "skew", ylab = paste("gamma1-",i,sep=""))
			f1 = .kde2d(md[,idx["skew",1]], md[,idx["gamma",1]+i-1])
			image(f1, xlab = "skew", ylab = paste("gamma1-",i,sep=""), add =TRUE)
			contour(f1,col = "steelblue", add = TRUE, method = "edge",
					vfont = c("sans serif", "plain"))
			title(paste("Skew vs Asymmetry(1) Parameter", sep=""))
			mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", 
					cex = 0.5)
		}
	}
	if(modelinc[16]>0 && modelinc[11]>0 && ind>0){
		for(i in 1:modelinc[11]){
			plot(md[,idx["skew",1]], md[,idx["eta1",1]+i-1], xlab = "skew", ylab = paste("eta1-",i,sep=""))
			f1 = .kde2d(md[,idx["skew",1]], md[,idx["eta1",1]+i-1])
			image(f1, xlab = "skew", ylab = paste("eta1-",i,sep=""), add =TRUE)
			contour(f1,col = "steelblue", add = TRUE, method = "edge",
					vfont = c("sans serif", "plain"))
			title(paste("Skew vs Shift Parameter", sep=""))
			mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", 
					cex = 0.5)
			mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
					padj = 1.5, col = "gray", cex = 0.5)
		}
	}
	if(modelinc[16]>0 && modelinc[12]>0 && ind>0){
		for(i in 1:modelinc[12]){
			plot(md[,idx["skew",1]], md[,idx["eta2",1]+i-1], xlab = "skew", ylab = paste("eta2-",i,sep=""))
			f1 = .kde2d(md[,idx["skew",1]], md[,idx["eta2",1]+i-1])
			image(f1, xlab = "skew", ylab = paste("eta2-",i,sep=""), add =TRUE)
			contour(f1,col = "steelblue", add = TRUE, method = "edge",
					vfont = c("sans serif", "plain"))
			title(paste("Skew vs Roatation Parameter", sep=""))
			mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", 
					cex = 0.5)
			mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
					padj = 1.5, col = "gray", cex = 0.5)
		}
	}
	return(1)
}

# distribution plots (include gamma1/gamma2 versus skew)
.distplot = function(x, window, ...)
{
	vmodel = x@model$modeldesc$vmodel
	vsubmodel = x@model$modeldesc$vsubmodel
	model = x@model
	modelinc = model$modelinc
	idx = model$pos.matrix
	distribution = model$modeldesc$distribution
	truecf = unlist(x@truecoef[,1])
	Names = names(truecf)
	md = as.data.frame(x, window = window)
	md = na.omit(md)
	
	if(modelinc[16]==0 | modelinc[17]==0 | (modelinc[16]==0 && (modelinc[10]==0 | modelinc[11]==0 | modelinc[12]==0))){
		warning("\nNo plots for distribution...\n")
		return(0)
	}
	
	if(modelinc[16]>0 && modelinc[17]>0) total = 1 + sum(modelinc[10:12])
	dv = .divisortable(total)
	par(mfrow=c(dv[1], dv[2]))
	start = dev.next(which = dev.cur())
	dev.new(start+1)
	if(modelinc[16]>0 && modelinc[17]>0){
		plot(md[,idx["skew",1]], md[,idx["shape",1]], ylab = "shape", xlab = "skew")
		f1 = .kde2d(md[,idx["skew",1]], md[,idx["shape",1]])
		image(f1, ylab = "shape", xlab = "skew", add =TRUE)
		contour(f1,col = "steelblue", add = TRUE, method = "edge",
				vfont = c("sans serif", "plain"))
		title(paste("Skew vs Shape Parameter (",distribution,")", sep=""))
		mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", 
				cex = 0.5)
		if(vmodel == "fGARCH"){
			mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
					padj = 1.5, col = "gray", cex = 0.5)
		}
	}
	if(modelinc[16]>0 && modelinc[10]>0){
		for(i in 1:modelinc[10]){
			plot(md[,idx["skew",1]], md[,idx["gamma",1]+i-1], xlab = "skew", ylab = paste("gamma1-",i,sep=""))
			f1 = .kde2d(md[,idx["skew",1]], md[,idx["gamma",1]+i-1])
			image(f1, xlab = "skew", ylab = paste("gamma1-",i,sep=""), add =TRUE)
			contour(f1,col = "steelblue", add = TRUE, method = "edge",
					vfont = c("sans serif", "plain"))
			title(paste("Skew vs Asymmetry(1) Parameter", sep=""))
			mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", 
					cex = 0.5)
		}
	}
	if(modelinc[16]>0 && modelinc[11]>0){
		for(i in 1:modelinc[11]){
			plot(md[,idx["skew",1]], md[,idx["eta1",1]+i-1], xlab = "skew", ylab = paste("eta1-",i,sep=""))
			f1 = .kde2d(md[,idx["skew",1]], md[,idx["eta1",1]+i-1])
			image(f1, xlab = "skew", ylab = paste("eta1-",i,sep=""), add =TRUE)
			contour(f1,col = "steelblue", add = TRUE, method = "edge",
					vfont = c("sans serif", "plain"))
			title(paste("Skew vs Shift Parameter", sep=""))
			mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", 
					cex = 0.5)
			mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
					padj = 1.5, col = "gray", cex = 0.5)
		}
	}
	if(modelinc[16]>0 && modelinc[12]>0){
		for(i in 1:modelinc[12]){
			plot(md[,idx["skew",1]], md[,idx["eta2",1]+i-1], xlab = "skew", ylab = paste("eta2-",i,sep=""))
			f1 = .kde2d(md[,idx["skew",1]], md[,idx["eta2",1]+i-1])
			image(f1, xlab = "skew", ylab = paste("eta2-",i,sep=""), add =TRUE)
			contour(f1,col = "steelblue", add = TRUE, method = "edge",
					vfont = c("sans serif", "plain"))
			title(paste("Skew vs Roatation Parameter", sep=""))
			mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", 
					cex = 0.5)
			mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
					padj = 1.5, col = "gray", cex = 0.5)
		}
	}
	return(1)
}

#------------------------------------------------------------------------------
# uGARCHdistribution : Models relating to the Coefs (RMSE plot):
#------------------------------------------------------------------------------

.rmseplots = function(x, ...)
{
	# T^0.5 consistency comparison
	vmodel = x@model$modeldesc$vmodel
	vsubmodel = x@model$modeldesc$vsubmodel
	model = x@model
	modelinc = model$modelinc
	idx = model$pos.matrix
	truecoef = x@truecoef[,1]
	cnames = names(truecoef)
	n = length(cnames)
	n.sim = x@dist$details$n.sim
	nwindows = x@dist$details$nwindows
	recursive.window = x@dist$details$recursive.window
	z = n.sim + ((1:nwindows)-1) * recursive.window
	# expected rate of decrease of rmse for the increasing window size
	rmse.roc = sqrt( z[-1] / z[1])
	dv = .divisortable(n)
	par(mfrow=c(dv[1], dv[2]))
	rmsedf = data.frame()
	rmsedf = t(sapply(as.list(1:nwindows), FUN = function(i) as.data.frame(x, which = "rmse", window = i)))
	
	for(i in 1:n){
		plot(z, unlist(rmsedf[,i]), type = "l", ylab = "rmse", xlab = "sample size", main  = cnames[i])
		lines(z, c(rmsedf[1,i], unlist(rmsedf[1,i])/rmse.roc), col = 2)
	}
	title(main = "RMSE ROC Plots", sub = "(red line: expected RMSE ROC for sqrt(N) consistency)", 
			outer = TRUE, line = -1, cex = 0.8, 
			cex.sub = 1, col.sub = 2, col.main = "steelblue")
	if(vmodel == "fGARCH"){
		mtext(paste("GARCH model :", vsubmodel, sep=""), side = 4,
				padj = -3, col = "gray67", cex = 0.6, outer = TRUE)
	} else{
		mtext(paste("GARCH model :", vmodel), side = 4,  padj=-3, col = "gray67", 
				cex = 0.6, outer = TRUE)
	}
}


#------------------------------------------------------------------------------
# uGARCHdistribution : Stat Plots
#------------------------------------------------------------------------------
.statsplot = function(x, ...)
{
	vmodel = x@model$modeldesc$vmodel
	vsubmodel = x@model$modeldesc$vsubmodel
	model = x@model
	modelinc = model$modelinc
	idx = model$pos.matrix
	nwindows = x@dist$details$nwindows
	colx = rainbow(nwindows, alpha = 0.95, start = 0.6, end = 0.9)
	nm = x@dist$details$n.sim + (0:(nwindows-1))*x@dist$details$recursive.window
	md = vector(mode = "list", length = nwindows)
	dens = vector(mode = "list", length = nwindows)
	maxy = vector(mode = "list", length = nwindows)
	maxx = vector(mode = "list", length = nwindows)
	minx = vector(mode = "list", length = nwindows)
	for(i in 1:nwindows){
		md[[i]] = as.data.frame(x, which = "stats", window = i)
		exc = unique(which(is.na(md[[i]]), arr.ind = TRUE)[,1])
		n = dim(md[[i]])[1]
		n1 = round(0.025*n)
		n2 = round(0.15*n)
		if(length(exc)>0) md[[i]] = md[[i]][-exc, , drop = FALSE]
		dens[[i]] = lapply(md[[i]], FUN = function(x) density(sort(x)[-c(1:n1,n-0:(n1-1))]))
		dens[[i]][[1]] = density(md[[i]][,1]/nm[i])
		dens[[i]][[10]] = density(log(0.5)/log(sort(md[[i]][,2])[-c(1:n1,n-0:(n2-1))]))
		dens[[i]][[3]] = density(sort(md[[i]][,3])[ -c(n-0:(n2-1)) ])
		maxy[[i]] = sapply(dens[[i]], FUN = function(x) max(x$y))
		maxx[[i]] = sapply(dens[[i]], FUN = function(x) max(x$x))
		minx[[i]] = sapply(dens[[i]], FUN = function(x) min(x$x))
	}
	maxy = apply(matrix(unlist(maxy), ncol = 12, byrow = T), 2, FUN = function(x) max(x))
	maxx = apply(matrix(unlist(maxx), ncol = 12, byrow = T), 2, FUN = function(x) max(x))
	minx = apply(matrix(unlist(minx), ncol = 12, byrow = T), 2, FUN = function(x) min(x))
	par(mfrow = c(3,3))
	# LLH
	plot(dens[[1]][[1]], col = colx[1], main = "Avg.LLH", ylim = c(0, maxy[1]),
			xlim = c(minx[1], maxx[1]))
	if(nwindows>1) for(i in 2:nwindows) lines(dens[[i]][[1]], col = colx[i])
	
	# Persistence
	plot(dens[[1]][[2]], col = colx[1], main = "Persistence", ylim = c(0, maxy[2]),
			xlim = c(minx[2], maxx[2]))
	if(nwindows>1) for(i in 2:nwindows) lines(dens[[i]][[2]], col = colx[i])
	
	# Half-Life
	plot(dens[[1]][[10]], col = colx[1], main = "Half-Life", ylim = c(0, maxy[10]),
	xlim = c(minx[10], maxx[10]))
	if(nwindows>1) for(i in 2:nwindows) lines(dens[[i]][[10]], col = colx[i])
	
	# UncVar
	plot(dens[[1]][[3]], col = colx[1], main = "Unc.Variance", ylim = c(0, maxy[3]),
			xlim = c(minx[3], maxx[3]))
	if(nwindows>1) for(i in 2:nwindows) lines(dens[[i]][[3]], col = colx[i])
	
	# UncMean
	plot(dens[[1]][[4]], col = colx[1], main = "Unc.Mean", ylim = c(0, maxy[4]),
			xlim = c(minx[4], maxx[4]))
	if(nwindows>1) for(i in 2:nwindows) lines(dens[[i]][[4]], col = colx[i])
	
	# Skewness
	plot(dens[[1]][[9]], col = colx[1], main = "Sample Skewness", ylim = c(0, maxy[9]),
			xlim = c(minx[9], maxx[9]))
	if(nwindows>1) for(i in 2:nwindows) lines(dens[[i]][[9]], col = colx[i])
	# Kurtosis
	plot(dens[[1]][[8]], col = colx[1], main = "Sample Kurtosis", ylim = c(0, maxy[8]),
			xlim = c(minx[8], maxx[8]))
	if(nwindows>1) for(i in 2:nwindows) lines(dens[[i]][[8]], col = colx[i])
	plot(1,1, axes = FALSE, type = "n", main = "color-key", ylab = "", xlab = "")
	lg.txt = paste("N = ",nm,sep="")
	legend("center", legend = lg.txt, fill = colx, col = colx, bty = "n")
	invisible(x)
}

#------------------------------------------------------------------------------
# uGARCHboot : Denstity Overlay Plots
#------------------------------------------------------------------------------
.densityoverlay = function(x, dens, n.overlays = 5, ...)
{
	# x 			: is the time series of length N
	# dens 			: is a list with the kernel or other density at every point in time of x
	# 				so that the every point x represents the mean of the density
	# no.overlays 	: how many density plots to overlay (should be a function of x)
	minx = min(x)
	maxx = max(x)
	N = length(x)
	if(N>10) v = floor(seq(10, N-10, length.out = n.overlays)) else v = min(N, 5)
	if(N>10) dist = c(diff(v), N - v[length(v)]) else dist = N-v
	#if(type != "line") plot(x, type = "l", col = "blue", ylim = c(-0.2, 0.2))
	for(i in 1:n.overlays){
		d = dens[[v[i]]]
		zx = d$y
		zy = d$x
		zx = zx/sum(zx)
		zx = v[i] + (zx * (dist[i] * 1/max(zx)))
		lines(zx, zy, lwd = 1, col = "lightgrey")
	}
}

.densityoverlay2 = function(x, dens, n.overlays = 5, ...)
{
	minx = min(x)
	maxx = max(x)
	N = length(x)
	if(N>10) v = floor(seq(10, N-10, length.out = n.overlays)) else v = min(N, 5)
	if(N>10) dist = c(diff(v), N - v[length(v)]) else dist = N-v
	#if(type != "line") plot(x, type = "l", col = "blue", ylim = c(0, 2*maxx))
	for(i in 1:n.overlays){
		d = dens[[v[i]]]
		zx = d$y
		zy = d$x
		zx = zx/sum(zx)
		zx = v[i] + (zx * (dist[i] * 1/max(zx)))
		lines(zx, zy, ylim = c(0, 0.1), lwd = 1, col = "lightgrey")
	}
}

.sigmaerrorplot = function(x, ...)
{
	vmodel = x@model$modeldesc$vmodel
	vsubmodel = x@model$modeldesc$vsubmodel
	model = x@model
	modelinc = model$modelinc
	idx = model$pos.matrix
	n.ahead = x@model$n.ahead
	fs = rep(NA, n.ahead)
	colr = NULL
	namr = NULL
	if(!is.null(x@model$filtered.s)){
		n = length(x@model$filtered.s)
		if(n.ahead>n) fs[1:n] = x@model$filtered.s else fs = x@model$filtered.s[1:n.ahead]
		colr = "green"
		namr = "filtered"
	}
	sigmafor = x@forc@forecast$forecast[[1]][,"sigma"]
	sigdist = vector(mode = "list", length = n.ahead)
	sigp = as.data.frame(x, which = "sigma", type = "q", qtile = c(0.05, 0.25, 0.75, 0.95))
	miny = min(sigp[1,])
	maxy = max(sigp[4,])
	meansig = apply(x@fsigma, 2, FUN = function(x) mean(x))
	plot(sigmafor, type = "l", col = "red", ylim = c(miny, maxy), main = "Sigma Forecast
					with Bootstrap Error Bands", cex.main = 0.7)
	lines(as.numeric(meansig), col = "black")
	lines(as.numeric(fs), col = colr)
	points(as.numeric(sigp[1,]), col = "steelblue1", pch = 19, cex = 0.5)
	points(as.numeric(sigp[2,]), col = "steelblue2", pch = 19, cex = 0.5)
	points(as.numeric(sigp[3,]), col = "steelblue3", pch = 19, cex = 0.5)
	points(as.numeric(sigp[4,]), col = "steelblue4", pch = 19, cex = 0.5)
	mtext(paste("GARCH model :", vmodel), side = 4, adj = 0, padj=0, col = "gray", 
			cex = 0.5)
	if(vmodel == "fGARCH"){
		mtext(paste("fGARCH submodel: ", vsubmodel,sep=""), side = 4, adj = 0, 
				padj = 1.5, col = "gray", cex = 0.5)
	}
	lg.txt  = c("forecast", "bootstrap mean", namr, "q5%", "q25%", "q75%","q95%")
	legend("bottomleft", legend = lg.txt, fill = c("red", "black", colr, "steelblue1",
					"steelblue2", "steelblue3", "steelblue4"), col = c("red", "black", colr, "steelblue1",
					"steelblue2", "steelblue3", "steelblue4"), cex = 0.7, bty = "n")
}

.serieserrorplot = function(x, ...)
{
	n.ahead = x@model$n.ahead
	fs = rep(NA, n.ahead)
	rx = rep(NA, n.ahead)
	if(!is.null(x@model$realized.x)){
		n = length(x@model$realized.x)
		if(n.ahead>n) rx[1:n] = x@model$realized.x else rx = x@model$realized.x[1:n.ahead]
	}
	seriesfor = x@forc@forecast$forecast[[1]][,"series"]
	serp = as.data.frame(x, which = "series", type = "q", qtile = c(0.05, 0.25, 0.75, 0.95))
	miny = min(serp[1,])
	maxy = max(serp[4,])
	meanser = apply(x@fseries, 2, FUN = function(x) mean(x))
	plot(seriesfor, type = "l", col = "steelblue", ylim = c(miny, maxy))
	lines(as.numeric(meanser), col = "black")
	lines(as.numeric(rx), col = "green")
	points(as.numeric(serp[1,]), col = "steelblue1", pch = 19, cex = 0.5)
	points(as.numeric(serp[2,]), col = "steelblue2", pch = 19, cex = 0.5)
	points(as.numeric(serp[3,]), col = "steelblue3", pch = 19, cex = 0.5)
	points(as.numeric(serp[4,]), col = "steelblue4", pch = 19, cex = 0.5)
}
