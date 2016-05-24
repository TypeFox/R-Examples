`plot.corres` <-
function(x, main="", addcol=TRUE, extreme=0, 
         rcex=1, rcol=1, rlabels="", stretch=1.4,
         ccex = 1, ccol = 2, clabels="", ...) {

	if (!is(x, "corres")) stop("argument should be a correspondence object")
	dat = x@data$origOut

	xlimit = range(dat$rproj[,1])*stretch
	ylimit = range(dat$rproj[,2])*stretch
	plot(dat$rproj[,1], dat$rproj[,2], type="n", xlim=xlimit, ylim=ylimit,
	    xlab=paste("Factor 1  (", round(x@data$eigenrates[1]/10, 1), " %)", sep=""),
	    ylab=paste("Factor 2  (", round(x@data$eigenrates[2]/10, 1), " %)", sep=""))
	lines(c(max(dat$rproj[,1]), min(dat$rproj[,1])), c(0,0))
	lines(c(0,0), c(max(dat$rpro[,2]), min(dat$rproj[,2])))
	if (!(main == "")) mtext(main, 3, 1)

	if (length(rcol) == 1 ) rcol = rep(1, nrow(dat$rproj))
	if (length(rlabels)==1) rlabels = rownames(x@data$input)
	text(dat$rproj[,1], dat$rproj[,2], rlabels, cex=rcex, col=rcol)
	if (addcol) {
	    if (length(clabels)==1) clabels = colnames(x@data$input)
		if (extreme > 0) {
			x = data.frame(dat$cproj[,1:2])
			extremes = apply(x, 2, quantile, c(extreme, 1-extreme))
			Accept = as.factor((x[,2] < extremes[1,2] | x[,2] > extremes[2,2])|
			         (x[,1] < extremes[1,1] | x[,1] > extremes[2,1]))
	    	text(x[Accept==TRUE,1], x[Accept==TRUE,2], clabels[Accept==TRUE], 
			    font=2, cex=ccex, col=ccol)
		} else {
	    	text(dat$cproj[,1], dat$cproj[,2], clabels, font=2, cex=ccex, col=ccol)
		}
	}
}

