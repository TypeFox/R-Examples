barLegend <- function(pal, colorbreaks, fig, side, mar = rep(0,4), colpalette = NULL, ...) {
	#if (length(pal) == 1)
	#	pal <- colorRampPalette(get("palettes",envir=.colorEnv)[[pal]])(length(colorbreaks)-1);
	dpal <- get("palettes", envir = .colorEnv);
	NCOLORS <- length(colorbreaks)-1;
	if (length(pal) >= 3) {
		pal <- colorRampPalette(pal,space='Lab')(NCOLORS);	
	}
	else if (pal %in% names(dpal)) {
		pal <- colorRampPalette(dpal[[pal]],space='Lab')(NCOLORS);
	}
	else if (tolower(pal) == "temperature") {
		pal <- richColors(NCOLORS);	
	}
	else if (tolower(pal) == "terrain") {
		pal <- terrain.colors(NCOLORS);
	}
	else {
		stop("Unrecognized color palette specification");
	}
	
	if (!is.null(colpalette)) {
		pal <- colpalette;
	}
	n <- length(pal);
	x <- seq(0,n,1)/n;
	x <- rep(x,each=2);
	x <- x[-c(1,length(x))];
	x <- matrix(x,ncol=2,byrow=TRUE);
	par(fig=fig, mar=mar, new=TRUE);
	plot.new();
	if (side == 2 || side == 4) {
		xlim <- c(-0.1,0.1);
		ylim <- c(0,1);
		plot.window(xlim,ylim);
		segments(x0=0, y0=x[,1], x1=0, y1=x[,2], col=pal, lwd=8, lend=2);
	}
	else {
		xlim <- c(0,1);
		ylim <- c(-0.1,0.1);
		plot.window(xlim,ylim);
		segments(x0=x[,1], y0=0, x1=x[,2], y1=0, col=pal, lwd=8, lend=2);
	}
	tx <- numeric(3);
	tx[1] <- min(colorbreaks,na.rm=TRUE);
	tx[2] <- colorbreaks[median(1:length(colorbreaks))];
	tx[3] <- max(colorbreaks,na.rm=TRUE); 
	axis(side,at=c(0,0.5,1),labels=signif(tx,2),las=1,...);
}
