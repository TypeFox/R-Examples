#############################################################
#
#	addBAMMlegend(x, location = 'topleft', side = 'auto', nTicks = 2, direction = 'auto', shortFrac = 0.02, longFrac = 0.3, axisOffset = 0.002, cex.axis = 0.8, labelDist = 0.7, ...)
#
#		x = saved plot.bammdata object
#		location = 'topleft', 'topright', 'bottomleft','bottomright','top','bottom','left','right' OR coordinates for legend c(xmin, xmax, ymin, ymax)
#		side = side for tick marks, see axis() documentation, if NULL, automatically inferred
#		nTicks = number of ticks, outside of min and max
#		direction = 'auto' or vertical or horizontal
#		shortFrac = percent of axis that is short dimension of legend
#		longFrac = percent of axis that is long dimension of legend
#		axisOffset = distance from color bar for labels as a percent of total
#		cex.axis = size of axis labels
#		labelDist = distance from axis to labels, passed to mgp
#		... additional parameters to be passed to axis()
#

addBAMMlegend <- function(x, direction, side, location = 'topleft', nTicks = 2, shortFrac = 0.02, longFrac = 0.3, axisOffset = 0.002, cex.axis = 0.8, labelDist = 0.7, ...) {
	#location xmin,xmax,ymin,ymax
	
	if (hasArg('corners')) {
		stop('Error: some options have been deprecated. Please consult the documentation.')
	}
	
	if(!hasArg('direction')) {
		direction <- 'auto'
	}
	
	if (!identical(names(x), c('coords', 'colorbreaks', 'palette', 'colordens'))) {
		stop("x must be a saved plot.bammdata object.");
	}
	
	if (!direction %in% c('auto', 'vertical', 'horizontal')) {
		stop("direction must be auto, vertical or horizontal.");
	}
	
	if (is.character(location)) {
		if (!location %in% c('bottomleft','bottomright','topleft','topright','bottom','top','left','right')) {
			stop('location is not recognized.');
		}
	}
	
	colorbreaks <- x$colorbreaks;
	pal <- x$palette;

	n <- length(colorbreaks);

	#return plot region extremes and define outer coordinates
	minX <- grconvertX(par('fig')[1], from = 'ndc', to = 'user') 
	maxX <- grconvertX(par('fig')[2], from = 'ndc', to = 'user')
	minY <- grconvertY(par('fig')[3], from = 'ndc', to = 'user')
	maxY <- grconvertY(par('fig')[4], from = 'ndc', to = 'user')
	
	xrange <- maxX - minX
	yrange <- maxY - minY
	minX <- minX + xrange * 0.05
	maxX <- maxX - xrange * 0.05
	minY <- minY + yrange * 0.05
	maxY <- maxY - yrange * 0.05
	
	if (is.character(location)) {
	
		if (location == 'topleft' & direction %in% c('auto', 'vertical')) {
			location <- vector('numeric', length = 4);
			location[1] <- minX
			location[2] <- minX + (maxX - minX) * shortFrac
			location[3] <- maxY - (maxY - minY) * longFrac
			location[4] <- maxY
		} else
		
		if (location == 'topleft' & direction == 'horizontal') {
			location <- vector('numeric', length = 4);
			location[1] <- minX
			location[2] <- minX + (maxX - minX) * longFrac
			location[3] <- maxY - (maxY - minY) * shortFrac
			location[4] <- maxY
		} else
	
		if (location == 'topright' & direction %in% c('auto', 'vertical')) {
			location <- vector('numeric', length = 4);
			location[1] <- maxX - (maxX - minX) * shortFrac
			location[2] <- maxX
			location[3] <- maxY - (maxY - minY) * longFrac
			location[4] <- maxY
		} else
	
		if (location == 'topright' & direction == 'horizontal') {
			location <- vector('numeric', length = 4);
			location[1] <- maxX - (maxX - minX) * longFrac
			location[2] <- maxX
			location[3] <- maxY - (maxY - minY) * shortFrac
			location[4] <- maxY
		} else
	
		if (location == 'bottomleft' & direction %in% c('auto', 'vertical')) {
			location <- vector('numeric', length = 4);
			location[1] <- minX
			location[2] <- minX + (maxX - minX) * shortFrac
			location[3] <- minY
			location[4] <- minY + (maxY - minY) * longFrac
		} else
		
		if (location == 'bottomleft' & direction == 'horizontal') {
			location <- vector('numeric', length = 4);
			location[1] <- minX
			location[2] <- minX + (maxX - minX) * longFrac
			location[3] <- minY
			location[4] <- minY + (maxY - minY) * shortFrac
		} else
	
		if (location == 'bottomright' & direction %in% c('auto', 'vertical')) {
			location <- vector('numeric', length = 4);
			location[1] <- maxX - (maxX - minX) * shortFrac
			location[2] <- maxX
			location[3] <- minY
			location[4] <- minY + (maxY - minY) * longFrac
		} else
		
		if (location == 'bottomright' & direction == 'horizontal') {
			location <- vector('numeric', length = 4);
			location[1] <- maxX - (maxX - minX) * longFrac
			location[2] <- maxX
			location[3] <- minY
			location[4] <- minY + (maxY - minY) * shortFrac 
		} else
		
		if (location == 'left') {
			location <- vector('numeric', length = 4);
			location[1] <- minX
			location[2] <- minX + (maxX - minX) * shortFrac
			location[3] <- mean(par('usr')[3:4]) - ((maxY - minY) * longFrac)/2
			location[4] <- mean(par('usr')[3:4]) + ((maxY - minY) * longFrac)/2
			direction <- 'vertical'
		} else
	
		if (location == 'right') {
			location <- vector('numeric', length = 4);
			location[1] <- maxX - (maxX - minX) * shortFrac
			location[2] <- maxX
			location[3] <- mean(par('usr')[3:4]) - ((maxY - minY) * longFrac)/2
			location[4] <- mean(par('usr')[3:4]) + ((maxY - minY) * longFrac)/2
			direction <- 'vertical'
		} else
		
		if (location == 'top') {
			location <- vector('numeric', length = 4);
			location[1] <- mean(par('usr')[1:2]) - ((maxX - minX) * longFrac)/2
			location[2] <- mean(par('usr')[1:2]) + ((maxX - minX) * longFrac)/2
			location[3] <- maxY - (maxY - minY) * shortFrac
			location[4] <- maxY
			direction <- 'horizontal'
		} else
	
		if (location == 'bottom') {
			location <- vector('numeric', length = 4);
			location[1] <- mean(par('usr')[1:2]) - ((maxX - minX) * longFrac)/2
			location[2] <- mean(par('usr')[1:2]) + ((maxX - minX) * longFrac)/2
			location[3] <- minY
			location[4] <- minY + (maxY - minY) * shortFrac
			direction <- 'horizontal'
		}
	}


	# infer direction based on dimensions of legend box
	if (direction == 'auto') {
		if (((location[2] - location[1]) / (par('usr')[2] - par('usr')[1])) >= ((location[4] - location[3]) / (par('usr')[4] - par('usr')[3]))) {
			direction <- 'horizontal';
		} else {
			direction <- 'vertical';
		}
	}
	
	if (direction == 'horizontal') {
		axisOffset <- axisOffset * (par('usr')[4] - par('usr')[3]);
	} else if (direction == 'vertical') {
		axisOffset <- axisOffset * (par('usr')[2] - par('usr')[1]);
	}
	
	#determine side for labels based on location in plot and direction
	if (!hasArg('side')) {
		if (direction == 'vertical') { #side = 1 or 4
			if (mean(location[1:2]) <= mean(par('usr')[1:2])) {
				side <- 4;
			} else {
				side <- 2;
			}
		}
		if (direction == 'horizontal') { #side = 2 or 3
			if (mean(location[3:4]) > mean(par('usr')[3:4])) {
				side <- 1;
			} else {
				side <- 3;
			}
		}
	}
	
	if (direction == 'horizontal') {
		x <- seq(from = location[1], to = location[2], length.out = n);
		width <- location[3:4];
	} else {
		x <- seq(from = location[3], to = location[4], length.out = n);
		width <- location[1:2];
	}
	
	#get bin coordinates
	x <- rep(x,each = 2);
	x <- x[-c(1,length(x))];
	x <- matrix(x, ncol = 2, byrow = TRUE);
	
	#find tick locations
	#get equivalent rate bins
	z <- rep(colorbreaks,each = 2);
	z <- z[-c(1,length(z))];
	z <- matrix(z, ncol = 2, byrow = TRUE);

	tx <- trunc(seq(from = 1, to = nrow(x), length.out = nTicks + 2));
	tickLocs <- x[tx,1]
	tx <- z[tx,1]
	tickLocs[length(tickLocs)] <- max(x[,2])
	tx[length(tx)] <- max(z[,2])	
	
	#plot bar
	if (direction == 'horizontal') {
		rect(xleft = x[,1], ybottom = width[1], xright = x[,2], ytop = width[2], border = pal, col = pal, xpd = NA);
	} else {
		rect(xleft = width[1], ybottom = x[,1], xright = width[2], ytop = x[,2], border = pal, col = pal, xpd = NA);
	}
	
	#add tickmarks
	if (side == 1) { #bottom
		axis(side, at = tickLocs, pos = location[3] - axisOffset, labels = signif(tx, 2), xpd = NA, las = 1, cex.axis = cex.axis, mgp = c(3, labelDist, 0), ...);
	} 
	if (side == 3) { #top
		axis(side, at = tickLocs, pos = location[4] + axisOffset, labels = signif(tx, 2), xpd = NA, las = 1, cex.axis = cex.axis, mgp = c(3, labelDist, 0), ...);
	}
	if (side == 2) { #left
		axis(side, at = tickLocs, pos = location[1] - axisOffset, labels = signif(tx, 2), xpd = NA, las = 1, cex.axis = cex.axis, mgp = c(3, labelDist, 0), ...);
	}
	if (side == 4) { #right
		axis(side, at = tickLocs, pos = location[2] + axisOffset, labels = signif(tx, 2), xpd = NA, las = 1, cex.axis = cex.axis, mgp = c(3, labelDist, 0), ...);
	}
	invisible(list(coords = x, width = width, pal = pal, tickLocs = tickLocs, labels = tx))
}
	
