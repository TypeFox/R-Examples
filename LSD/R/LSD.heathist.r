

### heathist ###


#' @export
#' @name heathist
#' @aliases LSD.heathist
#' @title Color a histogram
#' @description A histogram with an additional color stripe based on a kernel density estimate.
#' @param x a numeric vector.
#' @param breaks a numeric value giving the breaks of the histogram.
#' @param xlab x label, standard graphics parameter.
#' @param ylab y label, standard graphics parameter.
#' @param main title(s) of the plot, standard graphics parameter.
#' @param colpal a character vector containing R built-in color names or a name of a \code{LSD} colorpalette as a character string (see disco() or \code{\link{disco}}) (defaults to "heat", if not specified).
#' @param rev logical: if \code{TRUE} (\code{FALSE} by default), a reversed colorpalette is used.
#' @param simulate logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to simulate dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param daltonize logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to enhance dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param cvd character string implying the type of color vision deficiency ("p" for protanope, "d" for deuteranope or "t" for tritanope).
#' @param alpha alpha value: a two-digit integer between 01 and 99 for color opacity, i.e. appearance of partial or full transparency (usage omitted by default).
#' @param nobox logical: if \code{TRUE} (\code{FALSE} by default), the box of the plot is omitted.
#' @param add.density if \code{TRUE} (\code{FALSE} by default), a density line is added to the histogram.
#' @param col.density a R build-in color for the density line (if \code{add.density = TRUE}).
#' @param add.rug if \code{TRUE} (\code{FALSE} by default), a rug (1-d plot of the data) is added below the histogram-bars.
#' @param col.rug a R build-in color for the rug (if \code{add.rug = TRUE}).
#' @param nrcol a non-negative integer specifying the number of colors to be used (defaults to 100, if not specified).
#' @param ... additional parameters to be passed to points and plot.
#' @author Bjoern Schwalb
#' @seealso \code{\link{comparisonplot}}, \code{\link{demotour}}, \code{\link{disco}}, \code{\link{colorpalette}}
#' @examples x = rnorm(1000,mean = sample(c(0,3),size = 1000,prob = c(0.4,0.6),replace = TRUE))
#' heathist(x,xlab="x",add.density=TRUE,col.rug="darkred")
#' 
#' heathist(x,xlab="x",colpal = "matlablike")
#' @keywords histogram


heathist = function(x,breaks = 20,xlab = NULL,ylab = NULL,main = "heathist",colpal = "greys",rev = FALSE,simulate = FALSE,daltonize = FALSE,cvd = "p",alpha = NULL,nobox = FALSE,add.density = FALSE,col.density = "darkred",add.rug = TRUE,col.rug = "black",nrcol = 100,...)
{
	if (!is.vector(x)) stop("x must be a vector!")
	if (is.null(xlab)){xlab = deparse(substitute(x))}
	if (!is.null(ylab)){print("ylab argument will be ignored!")}
	xhist = hist(x,plot = FALSE,breaks = breaks)
	d = density(x)
	plot(xhist,border = NA,freq = FALSE,xlab = xlab,main = "",...)
	mtext(paste(main),3,2,cex=1.25)
	usr = par("usr")
	dy = (max(xhist$density) - 0)/nrcol
	colpal  = colorpalette(colpal,nrcol,simulate = simulate,daltonize = daltonize,cvd = cvd,alpha = alpha,rev = rev)
	for(i in 1:nrcol){
		clip(usr[1],usr[2],0 + (i-1) * dy,0 + i*dy)
		plot(xhist,add = TRUE,axes = FALSE,col = colpal[i],border = NA,freq = FALSE,xlab = "",ylab = "",main = "")
	}
	do.call(clip,as.list(usr))
	plot(xhist,add = TRUE,lwd = .5 ,freq = FALSE,axes = FALSE,xlab = "",ylab = "",main = "")
	if (add.density){lines(d,lwd = 4,col = col.density)}
	rug(x,col = col.rug)
	if (!nobox){box()}
}


### aliases ###


LSD.heathist = heathist



