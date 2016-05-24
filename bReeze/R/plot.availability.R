plot.availability <- function(x, set, ...) {
### plotting availability data
		
	num.sets <- length(x)
	if(missing(set)) set <- 1:num.sets
	n.set <- length(set)
	if(!is.numeric(set)) set <- match(set, names(x))
	if(any(is.na(set))) stop("'set' not found")
	if(any(set<1) || any(set>num.sets)) stop("'set' not found")
	
	# prepare plot
	old.par <- par(no.readonly=TRUE)
	on.exit(par(old.par))
	
	plot.param <- list(...)
	if(any(names(plot.param)=="col")) color <- plot.param$col
	else color <- "black"
	if(any(names(plot.param)=="fill")) col.fill <- plot.param$fill
	else col.fill <- c("#B3DE69", "#FFED6F", "#FB8072")
	if(any(names(plot.param)=="col.lab")) col.lab <- plot.param$col.lab
	else col.lab <- "black"
	if(any(names(plot.param)=="col.axis")) col.axis <- plot.param$col.axis
	else col.axis <- "black"
	if(any(names(plot.param)=="cex")) c.cex <- plot.param$cex
	else c.cex <- 1
	if(any(names(plot.param)=="cex.lab")) cex.lab <- plot.param$cex.lab
	else cex.lab <- 1
	if(any(names(plot.param)=="cex.axis")) cex.axis <- plot.param$cex.axis
	else cex.axis <- 1
	if(any(names(plot.param)=="border")) border <- plot.param$border
	else border <- "black"
	if(any(names(plot.param)=="lwd")) lwd <- plot.param$lwd
	else lwd <- 1
	if(any(names(plot.param)=="mar")) mar <- plot.param$mar
	else mar <- c(0,0,0,0)
	if(any(names(plot.param)=="xlab")) xlab <- plot.param$xlab
	else xlab <- "Days"
	if(any(names(plot.param)=="ylab")) ylab <- plot.param$ylab
	else ylab <- "Months"
	if(any(names(plot.param)=="plot.names")) plot.names <- plot.param$plot.names
	else plot.names <- TRUE
	
	if(length(color)==1) color <- rep(color, 3)
	if(length(col.fill)==1) col.fill <- rep(col.fill, 3)
	
	# plot
	if(length(mar)!=4) stop("'mar' must be a vector of five numeric values")
	par(mfrow=c(n.set,1), mar=c(0,4,0,0), oma=mar+c(0,0,2,0))
	for(s in 1:n.set) {
		plot.new()
		m <- dim(x[[set[s]]]$daily)[1]
		d.s <- attr(x[[set[s]]]$daily, "num.daily.samples")
	
		for (i in 1:m) {
			d <- length(x[[set[s]]]$daily[i, !is.na(x[[set[s]]]$daily[i,])])-1
			col <- fill <- character(d)
			value <- x[[set[s]]]$daily[i,2:(d+1)]
			col[value==d.s] <- color[1]
			fill[value==d.s] <- col.fill[1]
			col[value<d.s & value>0] <- color[2]
			fill[value<d.s & value>0] <- col.fill[2]
			col[value==0] <- color[3]
			fill[value==0] <- col.fill[3]
			
			rect((1:d)/31,1-i/m, ((1:d)-1)/31, 1-(i-1)/m, col=fill, border=border, lwd=lwd)
			text(((1:d)-0.5)/31, 1-(i-0.5)/m, value, cex=0.4*c.cex, col=col)
		}
		if(s==1) {
			mtext(xlab, side=3, line=0.7, at=0.5, cex=0.8*cex.lab, col=col.lab)
			mtext(names(x[[set[s]]]$daily)[2:32], side=3, line=-0.2, at=((1:31)-0.5)/31, las=1, cex=0.6*cex.axis, col=col.axis)
		}
		if(plot.names) mtext(names(x)[set[s]], side=2, line=2.6, cex=0.8*cex.lab, col=col.lab)
		mtext(ylab, side=2, line=1.8, cex=0.8*cex.lab, col=col.lab)
		mtext(row.names(x[[set[s]]]$daily)[m:1], side=2, line=-0.4, at=((1:m)-0.5)/m, las=1, cex=0.6*cex.axis, col=col.axis)
	}
}
