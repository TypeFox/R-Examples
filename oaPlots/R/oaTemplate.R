
#' Create a Blank Plot
#' @param xlim x limits for the plot
#' @param ylim y limits for the plot
#' @return none, plot is created on the current device
#' @author Jason Waddell
#' @export
blankPlot <- function(xlim, ylim){
	plot(x = -1, y = -1, type = "n", xlim = xlim, ylim = ylim, 
			axes = FALSE, ann = FALSE)
}


#' Create a OA Plot Template
#' @param xlim x limits for the plot
#' @param ylim y limits for the plot
#' @param xgrid values at which to draw the x axis gridlines
#' @param ygrid values at which to draw the y axis gridlines
#' @param xlabels labels to print at the x tickmarks
#' @param ylabels labels to print at the y tickmarks
#' @param xlab a title for the x axis
#' @param ylab a title for the y axis
#' @param main an overall title for the plot
#' @param bgCol background color for the plot
#' @param col.axis color for the axis labels
#' @param col.lab color for the xlab and ylab titles
#' @param col.main color for the main title
#' @param cex.axis size of the axis labels
#' @param cex.lab size of the xlab and ylab titles
#' @param cex.main size of the main title
#' @param box binary specifying whether to draw a bounding box around the plot
#' @param box.col color of the bounding box
#' @param box.lwd width of the bounding box lines
#' @param buffer optional buffer around all edges of the plot (as a percentage of the plot)
#' @param xaxs The style of axis interval calculation to be used for the x-axis. Possible values are "r", "i".
#' Style "r" (regular) first extends the data range by 4 percent at each end and then finds an axis with pretty labels that fits within the extended range.
#' Style "i" (internal) just finds an axis with pretty labels that fits within the original data range.
#' @param yaxs The style of axis interval calculation to be used for the y-axis. See xaxs above.
#' @param add A logical value specifying whether to add the template to an existing plot. If FALSE, a 
#' new plot will be created
#' @param gridLabelBuffer buffer between plot and grid labels (as a proportion of plotting range)
#' @param ylabBuffer distance between plot and y-axis title, as proportion of total plot width
#' @param xlabBuffer distance between plot and x-axis title, as proportion of total plot height
#' @param mainBuffer distance between plot and main title, as proportion of total plot height
#' @return none, objects are plotted to the current device
#' @author Jason Waddell
#' @export
#' @examples 
#' par(plt = c(0, 1, 0, 1))
#' oaTemplate(xlim = c(0, 10), ylim = c(20, 50), add = FALSE, xlab = "X Label", ylab = "Y Label",
#' 		main = "Main Title")
oaTemplate <- function(xlim, ylim, xgrid = NULL, ygrid = NULL,
		xlab = NULL, ylab = NULL, main = NULL, bgCol = gray(0.9),
		col.axis = gray(0.6), col.lab = gray(0.4), col.main = gray(0.3),
		cex.axis = 0.7, cex.lab = 1, cex.main = 1.5,	
		xaxs = "r", yaxs = "r", add = FALSE,
		box = FALSE, box.col = "black", box.lwd = 1,
		ylabels = NULL, xlabels = NULL,
		buffer = 0, gridLabelBuffer = 0.01, 
		ylabBuffer = 0.1, xlabBuffer = 0.08, mainBuffer = 0.07){
	
	if(xaxs == "r"){
		xBuffer <- 0.04*diff(xlim)
		xlim <- c(xlim[1]-xBuffer, xlim[2]+xBuffer)
	}
		
	if(yaxs == "r"){
		yBuffer <- 0.04*diff(ylim)
		ylim <- c(ylim[1]-yBuffer, ylim[2]+yBuffer)
	}
	
	if(!add)
		blankPlot(xlim = c(xlim[1] - buffer*diff(xlim), xlim[2] + buffer*diff(xlim)), 
				ylim = c(ylim[1] - buffer*diff(ylim), ylim[2] + buffer*diff(ylim)))
	
	op <- par("xpd")
	par(xpd = NA)
	
	rect(xleft = xlim[1], xright = xlim[2], ybottom = ylim[1], 
			ytop = ylim[2], col = bgCol, border = FALSE)
	
	# x gridlines
	if(is.null(xgrid)){
			xgrid <- pretty(xlim)
			xgrid <- xgrid[xgrid > min(xlim)]
			xgrid <- xgrid[xgrid < max(xlim)]
			segments(x0 = xgrid, y0 = min(ylim), y1 = max(ylim), col = "white", lwd = 2)
	} else{
		if(!is.na(xgrid[1])){
			segments(x0 = xgrid, y0 = min(ylim), y1 = max(ylim), col = "white", lwd = 2)
		}
	}
	
	# x axis labels
	if(is.null(xlabels)) xlabels = xgrid
	if(length(xgrid) != length(xlabels)) stop("Length of x-axis labels and grid do not match.")
	text(x = xgrid, y = ylim[1] - gridLabelBuffer*diff(ylim), 
			labels = xlabels, col = col.axis, 
			adj = c(0, 1), cex = cex.axis)
	
	# y gridlines	
	if(is.null(ygrid)){
		ygrid <- pretty(ylim)
		ygrid <- ygrid[ygrid > min(ylim)]
		ygrid <- ygrid[ygrid < max(ylim)]
		segments(y0 = ygrid, x0 = xlim[1], x1 = xlim[2], col = "white", lwd = 2)
	} else{
		if(!is.na(ygrid[1])){
			segments(y0 = ygrid, x0 = xlim[1], x1 = xlim[2], col = "white", lwd = 2)			
		}
	}
	
	# y axis labels
	if(is.null(ylabels)) ylabels = ygrid
	if(length(ygrid) != length(ylabels)) stop("Length of y-axis labels and grid do not match.")
	text(x = xlim[1] - gridLabelBuffer*diff(xlim), y = ygrid, 
			labels = ylabels, col = col.axis, 
			adj = c(1, 0), cex = cex.axis)

	# outer box
	if(box){
		xRange <- range(xlim); yrange <- range(ylim)
		rect(xRange[1], yrange[1], xRange[2], yrange[2], 
				col = NA, lwd = box.lwd, border = box.col)
	}
	
	# x axis title
	if(!is.null(xlab))
		text(x = mean(xlim), y = ylim[1] - xlabBuffer*diff(ylim),
				labels = xlab, col = col.lab, 
				adj = c(0.5, 0.5), cex = cex.lab)
	
	# y axis title
	if(!is.null(ylab))
		text(x = xlim[1] - ylabBuffer*diff(xlim), y = mean(ylim),
				labels = ylab, col = col.lab, srt = 90,
				adj = c(0.5, 0.5), cex = cex.lab)
	

	# main label
	if(!is.null(main))
		text(x = mean(xlim), y = ylim[2] + mainBuffer*diff(ylim),
				labels = main, col = col.main,
				adj = c(0.5, 0.5), cex = cex.main)
	par("xpd" = op)
	
}












#