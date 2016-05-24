

### align ###


#' @export
#' @name align
#' @aliases LSD.align
#' @title Visualize two-dimensional data in a color encoded fashion
#' @description Depict any matrix or list in a color encoded rectangular fashion.
#' @param input matrix or list with any type of entries.
#' @param colpal a character vector containing R built-in color names or a name of a \code{LSD} colorpalette as a character string (see disco() or \code{\link{disco}}).
#' @param simulate logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to simulate dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param daltonize logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to enhance dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param cvd character string implying the type of color vision deficiency ("p" for protanope, "d" for deuteranope or "t" for tritanope).
#' @param alpha alpha value: a two-digit integer between 01 and 99 for color opacity, i.e. appearance of partial or full transparency (usage omitted by default).
#' @param label logical: if \code{TRUE} (\code{FALSE} by default), labels are added according to the color scheme (i.e. binning).
#' @param digits integer indicating the number of decimals to be used for binning of continuous data.
#' @param border color for rectangle border(s). Use border = NA to omit borders.
#' @param xlim x limits, standard graphics parameter.
#' @param ylim y limits, standard graphics parameter.
#' @param main title of the plot, standard graphics parameter.
#' @param axes logical: if \code{TRUE} (by default), a box and axes are added to the plot (if \code{FALSE}, custom specification of axes can be achieved via basic R graphics functions).
#' @param ... additional parameters to be passed to points and plot.
#' @author Phillipp Torkler, Bjoern Schwalb
#' @seealso \code{\link{clusterplot}}, \code{\link{demotour}}, \code{\link{disco}}, \code{\link{colorpalette}}
#' @examples data(seqs)
#' colpal = c("A" = "darkgreen","C" = "darkblue","G" = "yellow","T" = "darkred")
#' align(seqs,colpal = colpal,label = TRUE,main = "DNA sequences")
#' 
#' data(homer)
#' colpal = c("white","black","yellow","wheat3")
#' align(homer,colpal = colpal,main = "D'OH!",asp = 1,axes = FALSE)
#' @keywords alignment, sequence


align = function(input,colpal = "heat",simulate = FALSE,daltonize = FALSE,cvd = "p",alpha = NULL,label = FALSE,digits = 1,border = NULL,xlim = NULL,ylim = NULL,main = NULL,axes = TRUE,...) 
{
	# stops execution, if 'input' is neither a list nor a matrix and executes an error action #
	
	if (!is.matrix(input) & !is.list(input)) stop("'input' must be a matrix or a list !")
	
	# if 'input' is a list #
	
	if (is.list(input)){
		
		# stops execution, if 'input'-list entries are of unequal length and executes an error action #
		
		if (any(sapply(input,length) != length(input[[1]]))){stop("'input'-list entries are of unequal length.")}
		
		# coerce to matrix #
		
		input = t(sapply(input,c))
	}

	# bin the 'input' data via the round function #
	
	dim.input = dim(input)
	input = as.vector(input)
	options(warn = -1)
	input[!is.na(as.numeric(input))] = round(as.numeric(input[!is.na(as.numeric(input))]),digits = digits)
	options(warn = 0)
	dim(input) = dim.input
	
	# provide color encoding #
	
	colsize = length(unique(as.vector(input)))
	cat(colsize,"color encoding \n")
	if (length(colpal) > 1){if (colsize != length(colpal)){warning("'colpal' contains not enough colors to encode 'input' data!")}}
	names.colpal = NULL
	if (!is.null(names(colpal))){names.colpal = names(colpal)}
	colorset = colorpalette(colpal,colsize,simulate = simulate,daltonize = daltonize,cvd = cvd,alpha = alpha)
	if (is.null(names.colpal)){names(colorset) = as.character(sort(unique(as.vector(input)),na.last = TRUE))} else {names(colorset) = names.colpal}
				
	# plot data in a grid like fashion #

	ysize = nrow(input)
	xsize = ncol(input)
	if (is.null(xlim)){xlim = c(0,xsize)}
	if (is.null(ylim)){ylim = c(0,ysize)}
	plot.new()
	plot.window(xlim = xlim,ylim = ylim,...)
	title(main)
	if (axes){
		axis(1)
		axis(2)
		box()
	}
	dump = sapply(1:ysize,function(i){rect(1:xsize - 1,ysize - i,1:xsize,ysize - i + 1,col = colorset[as.character(input[i,1:xsize])],border = border);if (label){text(x = 1:xsize - 0.5,y = ysize - i + 1 - 0.5,label = input[i,1:xsize],col = complementarycolor(colorset[as.character(input[i,1:xsize])]))}})
}


# alias #


LSD.align = align



