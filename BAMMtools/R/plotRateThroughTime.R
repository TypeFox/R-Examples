
#############################################################
#
#	plotRateThroughTime <- function(...)
#
#	ephy = object of class 'bammdata' or 'bamm-ratematrix'
#		if bamm-ratematrix, start.time, end.time, node, nslices, nodetype are not used.
#	useMedian = boolean, will plot median if TRUE, mean if FALSE.
#	intervals if NULL, no intervals will be plotted, otherwise a vector of quantiles must be supplied (these will define shaded polygons)
#	ratetype = autodetects diversification vs traits (based on input object 'type'), if 'auto', defaults to speciation (for diversification) or beta (for traits). Can alternatively specify 'netdiv' or 'extinction'. 
#	nBins = number of time slices used to generate rates through time
#	smooth = boolean whether or not to apply loess smoothing
#	smoothParam = loess smoothing parameter, ignored if smooth = F
#	opacity = opacity of color for interval polygons
#	intervalCol = transparent color for interval polygons
#	avgCol = color for mean/median line
#	start.time = start time in time before present to be fed to getRateThroughTimeMatrix
#	end.time = end time in time before present to be fed to getRateThroughTimeMatrix
#	node = if supplied, the clade descended from this node will be used.
#	nodetype = supplied to getRateThroughTimeMatrix
#	plot = boolean: if TRUE, a plot will be returned, if FALSE, the data for the plot will be returned. 
#	xticks = number of ticks on the x-axis, automatically inferred if NULL.
#	yticks = number of ticks on the y-axis, automatically inferred if NULL.
#	xlim = vector of length 2 with min and max times for x axis. X axis is time since present, so if plotting till the present, xlim[2]==0. Can also be 'auto'.
#	ylim = vector of length 2 with min and max rates for y axis. Can also be 'auto'. 
#	add = boolean: should rates be added to an existing plot
#
#	+ several undocumented args to set plot parameters: mar, cex, xline, yline, etc.
#	

plotRateThroughTime <- function(ephy, useMedian = TRUE, intervals=seq(from = 0,to = 1,by = 0.01), ratetype = 'auto', nBins = 100, smooth = FALSE, smoothParam = 0.20, opacity = 0.01, intervalCol='blue', avgCol='red',start.time = NULL, end.time = NULL, node = NULL, nodetype='include', plot = TRUE, cex.axis=1, cex.lab=1.3, lwd=3, xline=3.5, yline=3.5, mar=c(6,6,1,1), xticks=NULL, yticks=NULL, xlim='auto', ylim='auto',add=FALSE, axis.labels=TRUE) {
	
	if (!any(c('bammdata', 'bamm-ratematrix') %in% class(ephy))) {
		stop("ERROR: Object ephy must be of class 'bammdata' or 'bamm-ratematrix'.\n");
	}
	if (!is.logical(useMedian)) {
		stop('ERROR: useMedian must be either TRUE or FALSE.');
	}
	if (!any(c('numeric', 'NULL') %in% class(intervals))) {
		stop("ERROR: intervals must be either 'NULL' or a vector of quantiles.");
	}
	if (!is.logical(smooth)) {
		stop('ERROR: smooth must be either TRUE or FALSE.');
	}
	
	if ('bammdata' %in% class(ephy)) {
		#get rates through binned time
		rmat <- getRateThroughTimeMatrix(ephy, start.time = start.time, end.time = end.time, node = node, nslices = nBins, nodetype=nodetype);
	}
	if ('bamm-ratematrix' %in% class(ephy)) {
		if (!any(is.null(c(start.time, end.time, node)))) {
			stop('ERROR: You cannot specify start.time, end.time or node if the rate matrix is being provided. Please either provide the bammdata object instead or specify start.time, end.time or node in the creation of the bamm-ratematrix.')
		}
		#use existing rate matrix
		rmat <- ephy;
	}

	#set appropriate rates
	if (ratetype == 'speciation') {
		ratetype <- 'auto';
	}
	if (ratetype != 'auto' & ratetype != 'extinction' & ratetype != 'netdiv') {
		stop("ERROR: ratetype must be 'auto', 'extinction', or 'netdiv'.\n");
	}
	if (ephy$type == 'trait' & ratetype != 'auto') {
		stop("ERROR: If input object is of type 'trait', ratetype can only be 'auto'.")
	}
	if (ratetype == 'auto' & ephy$type == 'diversification') {
		rate <- rmat$lambda;
		ratelabel <- 'speciation rate';
	}
	if (ratetype == 'auto' & ephy$type == 'trait') {
		rate <- rmat$beta;
		ratelabel <- 'trait rate';
	}
	if (ratetype == 'extinction') {
		rate <- rmat$mu;
		ratelabel <- 'extinction rate';
	}
	if (ratetype == 'netdiv') {
		rate <- rmat$lambda - rmat$mu;
		ratelabel <- 'net diversification rate';
	}

	maxTime <- max(rmat$times)

	#remove NaN columns
	nanCol <- apply(rate, 2, function(x) any(is.nan(x)));
	rate <- rate[,which(nanCol == FALSE)];
	rmat$times <- rmat$times[which(nanCol == FALSE)];

	#generate coordinates for polygons
	rmat$times <- max(rmat$times) - rmat$times;
	if (!is.null(intervals)) {
		mm <- apply(rate, MARGIN = 2, quantile, intervals);

		poly <- list();
		q1 <- 1;
		q2 <- nrow(mm);
		repeat {
			if (q1 >= q2) {break}
			a <- as.data.frame(cbind(rmat$times,mm[q1,]));
			b <- as.data.frame(cbind(rmat$times,mm[q2,]));
			b <- b[rev(rownames(b)),];
			colnames(a) <- colnames(b) <- c('x','y');
			poly[[q1]] <- rbind(a,b);
			q1 <- q1 + 1;
			q2 <- q2 - 1;
		}
	}

	#Calculate averaged data line
	if (!useMedian) {
		avg <- colMeans(rate);
	} else {
		avg <- unlist(apply(rate,2,median));
	}
	
	#apply loess smoothing to intervals
	if (smooth) {
		for (i in 1:length(poly)) {
			p <- poly[[i]];
			rows <- nrow(p);
			p[1:rows/2,2] <- loess(p[1:rows/2,2] ~ p[1:rows/2,1],span = smoothParam)$fitted;
			p[(rows/2):rows,2] <- loess(p[(rows/2):rows,2] ~ p[(rows/2):rows,1],span = smoothParam)$fitted;
			poly[[i]] <- p;
		}
		avg <- loess(avg ~ rmat$time,span = smoothParam)$fitted;
	}

	#begin plotting
	if (plot) {
		if (!add) {
			plot.new();
			par(mar=mar);
			if (unique(xlim == 'auto') & unique(ylim == 'auto')) {
				xMin <- maxTime;
				xMax <- 0;
				if (!is.null(intervals)){
					yMin <- min(poly[[1]][,2]);
					yMax <- max(poly[[1]][,2]);
				} else {
					yMin <- min(avg);
					yMax <- max(avg);
				}
				if (yMin >= 0) {
					yMin <- 0;
				}
			}
			if (unique(xlim != 'auto') & unique(ylim == 'auto')) {
				xMin <- xlim[1];
				xMax <- xlim[2];
				if (!is.null(intervals)){
					yMin <- min(poly[[1]][,2]);
					yMax <- max(poly[[1]][,2]);
				} else {
					yMin <- min(avg);
					yMax <- max(avg);
				}
				if (yMin >= 0) {
					yMin <- 0;
				}
			}
			if (unique(xlim == 'auto') & unique(ylim != 'auto')) {
				xMin <- maxTime;
				xMax <- 0;
				yMin <- ylim[1];
				yMax <- ylim[2];
			}
			if (unique(xlim != 'auto') & unique(ylim != 'auto')) {
				xMin <- xlim[1];
				xMax <- xlim[2];
				yMin <- ylim[1];
				yMax <- ylim[2];
			}
			plot.window(xlim=c(xMin, xMax), ylim=c(yMin, yMax));
			if (is.null(xticks)) {
				axis(at=c(round(1.2*xMin),axTicks(1)), cex.axis = cex.axis, side = 1);
			}
			if (!is.null(xticks)) {
				axis(at=c(1.2*xMin,seq(xMin,xMax, length.out=xticks+1)), labels = c(1.2*xMin,signif(seq(xMin, xMax, length.out=xticks+1),digits=2)), cex.axis = cex.axis, side = 1);
			}
			if (is.null(yticks)) {
				axis(at=c(-1,axTicks(2)), cex.axis = cex.axis, las = 1, side = 2);
			}
			if (!is.null(yticks)) {
				axis(at=c(-0.2,seq(yMin, 1.2*yMax, length.out=yticks+1)), labels = c(-0.2,signif(seq(yMin, 1.2*yMax, length.out=yticks+1),digits=2)), las=1, cex.axis = cex.axis, side = 2);	
			}
			if (axis.labels) {
				mtext(side = 1, text = 'time before present', line = xline, cex = cex.lab);
				mtext(side = 2, text = ratelabel, line = yline, cex = cex.lab);
			}
		}
		#plot intervals
		if (!is.null(intervals)) {
			for (i in 1:length(poly)) {
				polygon(x=poly[[i]][,1],y=poly[[i]][,2],col=transparentColor(intervalCol,opacity),border=NA);
			}
		}
		lines(x = rmat$time, y = avg, lwd = lwd, col = avgCol);
	} else {
		return(list(poly = poly,avg = avg,times = rmat$time));
	}
}
