

### emptyplot ###


#' @export
#' @name emptyplot
#' @aliases LSD.emptyplot
#' @title Wrapper function for an empty graphics device
#' @description Calls an empty graphics device with a coordinate system of choice.
#' @param xlim x limits, standard graphics parameter.
#' @param ylim y limits, standard graphics parameter.
#' @param ... additional parameters to be passed to points and plot.
#' @author Bjoern Schwalb
#' @seealso \code{\link{demotour}}
#' @examples emptyplot()
#' @keywords empty plotting region


emptyplot = function(xlim = c(-1,1),ylim = c(-1,1),...)
{
	plot.new()
	plot.window(xlim = xlim,ylim = ylim,...)
}


### alias ###


LSD.emptyplot = emptyplot




### windowxy ###


#' @export
#' @name windowxy
#' @aliases LSD.windowxy
#' @title Factorization of the number of windows for plots with device partitions
#' @description Create a factorization of the number of windows for plots with device partitions to be used in par(mfrow = ...).
#' @param nrwin a non-negative integer specifying the number of windows.
#' @author Bjoern Schwalb
#' @seealso \code{\link{demotour}}
#' @examples windowxy(20)
#' @keywords factorization


windowxy = function(nrwin = 1)
{
	# factorize 'nrwin' #
	
	flsq = floor(sqrt(nrwin))
	if (flsq^2==nrwin){
		fac = c(flsq,flsq)
	} else if (nrwin <= flsq*(flsq+1)){
		fac = c(flsq,flsq+1)
	} else {
		fac = c(flsq,flsq+2)
	}
	
	# return 'fac' #
	
	return(fac)
}


### alias ###


LSD.windowxy = windowxy




### webdesign ###


#' @export
#' @name webdesign
#' @aliases LSD.webdesign
#' @title Colored rectangular grid
#' @description Adds a colored rectangular grid to an existing plot.
#' @param xlim x limits, standard graphics parameter.
#' @param ylim y limits, standard graphics parameter.
#' @param colpal a character vector containing R built-in color names or a name of a \code{LSD} colorpalette as a character string (see disco() or \code{\link{disco}}).
#' @param xlabels a character vector containing labels depicted parallel to the x-axis.
#' @param ylabels a character vector containing labels depicted parallel to the y-axis.
#' @param ... additional parameters to be passed to abline().
#' @author Bjoern Schwalb
#' @seealso \code{\link{demotour}}, \code{\link{disco}}, \code{\link{colorpalette}}
#' @examples emptyplot(c(-5,5),c(-5,5))
#' labels = c("2 fold","4 fold","8 fold")
#' webdesign(c(-5,5),c(-5,5),lty = 2,xlabels = labels,ylabels = labels)
#' @keywords grid, web


webdesign = function(xlim,ylim,colpal = "rdbu",xlabels = NULL,ylabels = NULL,...)
{
	# provide color encoding #
	
	colpal = colorpalette(colpal,max(diff(xlim),diff(ylim))+1)
	
	# add abscissa and ordinate #
	
	abline(h=0,v=0,...)
	
	# line segments along the x-axis #
	
	xseq = xlim[2]:xlim[1]
	abline(h=xseq,col=colpal,...)
	
	# line segments along the y-axis #
	
	yseq = ylim[2]:ylim[1]
	abline(v=yseq,col=colpal,...)
	
	# add labels if !NULL #
	
	if (!is.null(xlabels)){text(1:length(ylabels),ylim[2]-1,xlabels,col=colpal[rev(which(yseq %in% 1:length(xlabels)))],adj=c(-0.1,-0.1))}
	if (!is.null(ylabels)){text(xlim[1]+1,-(1:length(ylabels)),ylabels,col=colpal[which(xseq %in% -(1:length(ylabels)))],adj=c(-0.1,-0.1))}
}


### alias ###


LSD.webdesign = webdesign

	
	
