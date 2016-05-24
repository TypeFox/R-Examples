plot.flux <-
function(x, zero.line, note = "", margin = 0.2, xlims = NULL, ...){
	## first, plotting the line
	const <- x$fl.dat$orig.dat[,1]
	time <- x$fl.dat$orig.dat$time
	yl <- names(x$fl.dat$orig.dat)[1]
	gc.qual <- x$fl.dat$orig.dat$gc.qual
	range.const <- range(const, na.rm=TRUE)
	if(is.null(xlims)){
		xlims <- range(time, na.rm=TRUE)
		}
	ylims <- range.const + c(-margin*diff(range.const), margin*diff(range.const))
	plot(const ~ time, type="l", ylim = ylims, xlim = xlims, xlab = "time", ylab = yl, ...)
	## plotting the points (all points)
	points(const ~ time, col=c(1,2)[as.factor(gc.qual==0)], cex=1.2)
	## plotting the gaverage global concentration line
	abline(h=zero.line, col="grey30")
	## marking the points that has been kept into the analysis
	## after odae (outlier detection and exclusion)
	points(x$fl.dat$lm4flux$model[,c(2,1)], col="darkblue", pch=19, cex=0.8)
	## plotting the regression line of the chosen model
	abline(x$fl.dat$lm4flux, lty=3, col="red4")
	## whether a flux was calculated 
	## 0 - range check and r2-check failed
	## 1 - either range check or r2-check failed
	## 2 - checks OK, flux was calculated
	res <- paste(round(x$fluss$flux, 3))
	pv <- coef(summary(x$fl.dat$lm4flux))[2,4]
	symp <- symnum(pv, corr=FALSE, cutpoints = c(0,  .001,.01,.05, .1, 1), symbols = c("***","**","*","."," "))
	with(x, mtext(paste(fluss[[7]]*1, fluss[[5]]*1, fluss[[6]]*1, ".", fluss[[8]], fluss[[9]]*1, ": ", res, symp, " ", unit, "/m2*h", sep="")), cex=0.8)
	mtext(paste(note), line=-1.4, cex=0.8)
	}

