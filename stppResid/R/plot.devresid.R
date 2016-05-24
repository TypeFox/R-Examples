plot.devresid <- function(x, ..., col.key = rev(heat.colors(100)), cutoffs = NULL)
{
	residuals <- x$residuals
	if(is.null(cutoffs)) {	
		cutoff.key <- seq(min(residuals) - 1e-09, max(residuals) + 1e-09, length.out = length(col.key)+1)
	}
	else {
		cutoff.key <- cutoffs
		if((length(cutoffs) - 1) != length(col.key))
			stop("length of cutoffs should be 1 more than length of col.key")
	}
	if(!is.character(col.key))
		stop("col.key must be character vector of colors in hexadecimal")
	gr <- x[[2]]$grid.full
	xv <- c(unique(gr$xmin), tail(gr$xmax, 1))
	yv <- c(unique(gr$ymin), tail(gr$ymax, 1))
	z <- matrix(residuals, nrow = length(xv) - 1, byrow = TRUE)
	layout(matrix(c(1,2), ncol=1), heights=c(2,.5))
	par(mar=c(4,4,4,4), bty="n")
	if(is.null(cutoffs)) {
		image(xv, yv, z, xlab = "x", ylab = "y", col=col.key)
	}
	else { 
		image(xv, yv, z, xlab = "x", ylab = "y", col=col.key, breaks = cutoffs)
	}
	points(x[[1]]$x, x[[1]]$y, ...)	
	par(mar=c(1, 1, 1.5, 1))
	key <- (0:length(residuals))/length(residuals)
	plot(NULL, ylim=c(-3,0), xlim=c(-0.2, 1.2), type="n", axes=F, xlab="", ylab="", main="")
	image(key, -2:0, matrix(rep(key,2), nrow=length(residuals)+1, byrow=F), add=T, col=col.key)
	for(i in seq(0, 1, length=10)) lines(c(i,i), c(-2,0), col=gray(0.3), lty=2)
	text(0, -2.5, round(min(cutoff.key), 3), cex=1)
	text(0.25, -2.5, round(quantile(cutoff.key, 0.25), 3), cex=1)
	text(0.5, -2.5, round(median(cutoff.key), 3), cex=1)
	text(0.75, -2.5, round(quantile(cutoff.key, 0.75), 3), cex=1)
	text(1, -2.5, round(max(cutoff.key), 3), cex=1)
	mtext("Deviance residuals", 3)	
}