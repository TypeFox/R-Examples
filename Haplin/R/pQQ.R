pQQ <- function(pvals, nlabs = 6, conf = 0.95, lim, mark = 0.05, ...){
##
## QQ-PLOT
## pvals: LIST OF P-VALUES, PREFERABLY WITH NAMES (WHICH WILL BE USED AS LABELS)
## nlabs: NUMBER OF GENE NAME LABELS TO BE PRINTED (CAN BE SET TO ZERO)
## conf: NUMERICAL VALUE, SAY 0.95, WHICH IS THE LEVEL OF THE POINTWISE CONF. INTERVALS
## lim: LIMITS OF PLOTTING REGION (-log10 VALUES, FOR INSTANCE	c(0,3))
## sim: SET TO TRUE TO CHECK THE CONFIDENCE INTERVALS WITH SIMULATIONS
## mark: SPECIFY MARKING OF "SIGNIFICANT" GENES (CAN BE FALSE, I.E. NO MARKING)
## ...: OTHER PARAMETERS TO plot, LIKE main, OR WHATEVER
##
if(any(is.na(pvals))) stop("Missing p-values not allowed", call. = F)
sim <- F
#
##
.pvals <- sort(pvals)
.logpvals <- -log10(.pvals)
.order <- seq(along = .pvals)
.n <- length(.pvals)
.aux <- .order/(.n + 1) # EXPECTED VALUES, ON ORIGINAL SCALE
.logaux <- -log10(.aux)
#
##
#
###.xt <- seq(10^(-.lim[2])*(.n+1), 0.99, length.out = 10)
.xt <- 1/2
.fix.order <- c(.xt, .order) # ONLY TO MAKE SURE THE SHADING GOES ALL THE WAY UP TO THE BORDERLINE
.fix.logaux <- -log10(.fix.order/(.n+1))
#
## COMPUTE CORRECT MEDIAN VALUE FROM BETA DISTRIBUTION (INSTEAD OF USING "EXPECTED")
.b.median <- qbeta(p = 1/2, shape1 = .fix.order, shape2 = .n - .fix.order + 1)
##  SMOOTH MEDIAN VALUE, FOR BETTER PLOT
###.sm.median <- approx(x = .fix.logaux, y = -log10(.b.median), xout = seq(.lim[1], .lim[2], length.out = 1000))
.f.median <- approxfun(x = .fix.logaux, y = -log10(.b.median))
.adj.logaux <- .f.median(.logaux)
.adj.fix.logaux <- .f.median(.fix.logaux)


#
##
if (missing(lim)){
	.lim <- c(0, max(.adj.logaux, .logpvals)) * 1.05
	#.lim <- c(0, max(.logaux, .logpvals) * 1.05)
}
else
	.lim <- lim
#
##
plot(.logaux, .logaux, type = "n", pch = 16, cex = 0.7, xlim = .lim, ylim = .lim, xlab = "Expected P-value (-log10 scale)", ylab = "Observed P-value (-log10 scale)", font = 2, lwd = 2, font.lab = 2, bty = "l", xaxs = "i", yaxs = "i", axes = T, las = 1, pty = "s", ...)

#
#
if (is.numeric(conf)) {
	#
	## COMPUTE LOWER AND UPPER LIMITS FROM BETA DISTRIBUTION
	.b.lower <- qbeta(p = (1-conf)/2, shape1 = .fix.order, shape2 = .n - .fix.order + 1)
	.b.upper <- qbeta(p = 1-(1-conf)/2, shape1 = .fix.order, shape2 = .n - .fix.order + 1)
	#
	##  SMOOTH LOWER AND UPPER LIMITS, FOR BETTER PLOT
	.sm.lower <- approx(x = .adj.fix.logaux, y = -log10(.b.lower), xout = seq(.lim[1], .lim[2], length.out = 1000))
	.sm.upper <- approx(x = .adj.fix.logaux, y = -log10(.b.upper), xout = seq(.lim[1], .lim[2], length.out = 1000))
	#
	## ADD BACKGROUND SHADING
	segments(x0 = .sm.lower$x, y0 = .sm.lower$y, x1 = .sm.upper$x, y1 = .sm.upper$y, col = "grey")
	#
	## ADD BOUNDARY LINES FOR SHADING
	lines(.adj.fix.logaux, -log10(.b.lower), col = "black", lwd = 1)
	lines(.adj.fix.logaux, -log10(.b.upper), col = "black", lwd = 1)
	#
	## IF REQUESTED, SIMULATIONS CAN BE USED TO CHECK BETA VALUES
	if(sim){
		.conf <- f.QQconf(nGenes = length(.pvals), quantiles = c((1-conf)/2, 1-(1-conf)/2))
		.s.lower <- rev(.conf$quant[1, ])
		.s.upper <- rev(.conf$quant[2, ])
		lines(.adj.logaux, .s.lower, lty=2, col="red", lwd=2)
		lines(.adj.logaux, .s.upper, lty=2, col="red", lwd=2)
		#
	}
}
#
## PLOT P-VALUES
points(.adj.logaux, .logpvals, pch = 16, cex = 0.7)
box(lwd = 2, bty = "l")
#
## LINE WITH SLOPE 1
abline(a = 0, b = 1, lwd = 2, col = "red")
#
## MARK 0.05 LIMIT (OR WHATEVER CHOSEN)
if(is.numeric(mark)){
	psig <- -log10(mark)
	lines(c(psig,psig), c(0,psig), lty=3, lwd=1.4)
	lines(c(0,psig), c(psig,psig), lty=3, lwd=1.4)
}
#
### LABELS ON THE nlabs MOST SIGNIFICANT
if(nlabs > 0){
	.ind <- seq(length.out = nlabs)
	text(.adj.logaux[.ind], .logpvals[.ind] - 0.05, names(.pvals)[.ind], srt = 270, font = 2, cex = 0.7, adj = 0, col = "black", xaxs = "i", yaxs = "i")
}
return(invisible())
}
