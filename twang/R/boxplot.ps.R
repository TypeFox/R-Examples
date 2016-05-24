boxplot.ps <- function(x, subset = NULL, color = TRUE, time = NULL, ...){
	longDat <- matrix(t(as.matrix(x$ps)), ncol = 1)
	nms <- names(x$ps)
	bwDat <- data.frame(ps = longDat, nm = nms, treat = rep(x$treat, each = length(nms)))
	if(is.null(subset)) subset <- 1:length(levels(as.factor(bwDat$nm)))
	ptSymCol <- ifelse(color, "#0080ff", "black")	
	bwCols <- list(col = ptSymCol)
	stripBgCol <- ifelse(color, "#ffe5cc", "transparent")
	
	if(is.null(time)) xlb = "Propensity scores"
	else xlb = paste("Propensity scores (Time ", time, ")", sep = "")

	pt1 <- bwplot(treat~ps|nm, data=bwDat, scales = list(alternating = 1), ylab = "Treatment", xlab=xlb, subset = as.factor(bwDat$nm) %in% levels(as.factor(bwDat$nm))[subset], par.settings = list(strip.background = list(col=stripBgCol), box.rectangle = bwCols, plot.symbol = bwCols, box.umbrella = bwCols), ...)
	return(pt1)
	
}