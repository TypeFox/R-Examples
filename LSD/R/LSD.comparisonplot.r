

### heatbarplot ###


#' @export
#' @name heatbarplot
#' @aliases LSD.heatbarplot
#' @title Color a barplot.
#' @description Depict a histogram object as a barplot in a color encoded fashion based on a kernel density estimate.
#' @param x a histogram object.
#' @param colpal a character vector containing R built-in color names or a name of a \code{LSD} colorpalette as a character string (see disco() or \code{\link{disco}}) (defaults to "heat", if not specified).
#' @param simulate logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to simulate dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param daltonize logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to enhance dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param cvd character string implying the type of color vision deficiency ("p" for protanope, "d" for deuteranope or "t" for tritanope).
#' @param alpha alpha value: a two-digit integer between 01 and 99 for color opacity, i.e. appearance of partial or full transparency (usage omitted by default).
#' @param rev logical: if \code{TRUE} (\code{FALSE} by default), a reversed colorpalette is used.
#' @param horizontal logical: if \code{TRUE} (\code{FALSE} by default), rotation of 90 degrees is applied.
#' @param nrcol a non-negative integer specifying the number of colors to be used (defaults to 100, if not specified).
#' @param ... additional parameters to be passed to points and plot.
#' @author Bjoern Schwalb
#' @seealso \code{\link{comparisonplot}}, \code{\link{demotour}}, \code{\link{disco}}, \code{\link{colorpalette}}
#' @examples points = 10^4
#' x = c(rnorm(points/2),rnorm(points/2)+4)
#' x = sign(x)*abs(x)^1.3
#' xhist = hist(x,plot = FALSE)
#' 
#' heatbarplot(xhist)
#' @keywords barplot


heatbarplot = function(x,colpal = "heat",simulate = FALSE,daltonize = FALSE,cvd = "p",alpha = NULL,rev = FALSE,horizontal = FALSE,nrcol = 100,...)
{
	if (class(x) != "histogram") stop("x must be of class histogram!")
	if (!horizontal){
		barplot(x$density,axes = FALSE,space = 0,border = "white",...)
		dy = (max(x$density) - 0)/nrcol
		colpal = colorpalette(colpal,nrcol,simulate = simulate,daltonize = daltonize,cvd = cvd,alpha = alpha,rev = rev)
		colorlane = function(a,b){for(i in 1:100){rect(a,0 + (i-1) * dy,b,0 + i*dy,col = colpal[i],border = NA)}}
		for (j in 1:length(x$counts)){colorlane(j-1,j)}
		for (j in 1:length(x$counts)){rect(j-1,x$density[j],j,max(x$density),col = "white",border = NA)}
		for (j in 1:length(x$counts)){rect(j-1,0,j,x$density[j],col = "transparent",border = NULL)}
	} else{
		barplot(x$density,axes = FALSE,space = 0,border = "white",horiz=TRUE,...)
		dy = (max(x$density) - 0)/nrcol
		colpal = colorpalette(colpal,nrcol,simulate = simulate,daltonize = daltonize,cvd = cvd,alpha = alpha,rev = rev)
		colorlane = function(a,b){for(i in 1:100){rect(0 + (i-1) * dy,a,0 + i*dy,b,col = colpal[i],border = NA)}}
		for (j in 1:length(x$counts)){colorlane(j-1,j)}
		for (j in 1:length(x$counts)){rect(x$density[j],j-1,max(x$density),j,col = "white",border = NA)}
		for (j in 1:length(x$counts)){rect(0,j-1,x$density[j],j,col = "transparent",border = NULL)}}
}


### aliases ###


LSD.heatbarplot = heatbarplot


### comparisonplot ###


#' @export
#' @name comparisonplot
#' @aliases cplot
#' @title Comparisonplot: a fancy scatterplot
#' @description A function to compare two vectors extensively.
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param histbreaks a non-negative integer specifying the number of breaks of the histograms.
#' @param adjust scale the used bandwidth of the density estimate, if \code{add.density = TRUE}.
#' @param colpal a character vector containing R built-in color names or a name of a \code{LSD} colorpalette as a character string (see disco() or \code{\link{disco}}).
#' @param simulate logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to simulate dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param daltonize logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to enhance dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param cvd character string implying the type of color vision deficiency ("p" for protanope, "d" for deuteranope or "t" for tritanope).
#' @param alpha alpha value: a two-digit integer between 01 and 99 for color opacity, i.e. appearance of partial or full transparency (usage omitted by default).
#' @param rev logical: if \code{TRUE} (\code{FALSE} by default), a reversed colorpalette is used.
#' @param main title(s) of the plot, standard graphics parameter.
#' @param cor if \code{TRUE} (\code{FALSE} by default), the correlation is added to the title.
#' @param xlab x label, standard graphics parameter.
#' @param ylab y label, standard graphics parameter.
#' @param xlim x limits, standard graphics parameter.
#' @param ylim y limits, standard graphics parameter.
#' @param ab if \code{TRUE} (\code{FALSE} by default), \code{abline(0,1)} is added to the heatscatter.
#' @param add.density if \code{TRUE} (\code{FALSE} by default), density lines are added to the barplots.
#' @param col.density R built-in color to specify the color of the density line.
#' @param pimp if \code{TRUE} (\code{FALSE} by default), the plot is pimped.
#' @param ... additional parameters to be passed to points and plot.
#' @author Bjoern Schwalb
#' @seealso \code{\link{align}}, \code{\link{demotour}}, \code{\link{disco}}, \code{\link{colorpalette}}
#' @examples points = 10^4
#' x = c(rnorm(points/2),rnorm(points/2)+4)
#' y = x + rnorm(points,sd=0.8)
#' x = sign(x)*abs(x)^1.3
#' 
#' comparisonplot(x,y,histbreaks=30,pch=20)
#' @keywords scatterplot, barplot


comparisonplot = function(x,y,histbreaks = 30,adjust = 1,colpal = "heat",simulate = FALSE,daltonize = FALSE,cvd = "p",alpha = NULL,rev = FALSE,main = "comparisonplot",cor = FALSE,xlab = NULL,ylab = NULL,xlim = NULL,ylim = NULL,ab = FALSE,add.density = FALSE,col.density = "darkred",pimp = FALSE,...)
{
	if (!is.vector(x) | !is.vector(y)) stop("First two argument must be vectors!")
	if (is.null(xlab)){xlab = deparse(substitute(x))}
	if (is.null(ylab)){ylab = deparse(substitute(y))}
	sound = which((!(is.na(x)|is.nan(x)|(x==Inf)|(x==-Inf))) & (!(is.na(y)|is.nan(y)|(y==Inf)|(y==-Inf))))
	if (length(sound)==0) stop("There are no valid point pairs to plot!")
	x = x[sound]
	y = y[sound]
	xrange <- c(min(range(x),range(y)),max(range(x),range(y)))
	yrange <- c(min(range(x),range(y)),max(range(x),range(y)))
	if (!is.null(xlim)){cut = x > xlim[1] & x < xlim[2]
		x = x[cut]
		y = y[cut]
		xrange = xlim
	}
	if (!is.null(ylim)){cut = y > ylim[1] & y < ylim[2]
		y = y[cut]
		x = x[cut]
		yrange = ylim
	}                           
	def.par = par(no.readonly = TRUE) # save default, for resetting...
	xhist = hist(x, breaks=seq(xrange[1],xrange[2],length.out=histbreaks),plot=FALSE)
	yhist = hist(y, breaks=seq(yrange[1],yrange[2],length.out=histbreaks),plot=FALSE)
	top = max(c(xhist$density,yhist$density))
	nf = layout(matrix(c(0,2,0,4,1,3,0,5,0),3,3,byrow=TRUE), c(1,3,1), c(1,3,1), TRUE)
	dx = density(x)
	dy = density(y)
	par(mar=c(4,4,4,4))
	heatscatter(x,y,xlim=xrange,ylim=yrange,xlab = xlab,ylab = ylab,main=main,colpal=colpal,cor=cor,simulate = simulate,daltonize = daltonize,cvd = cvd,alpha = alpha,rev = rev,...)
	if (pimp){
		if (ab){abline(0,1,col="#08306B50",lwd=2)}
		par(mar=c(0,4,1,4))
		heatbarplot(xhist,ylim=c(0,top),colpal=colpal,simulate = simulate,daltonize = daltonize,cvd = cvd,alpha = alpha,rev = !rev)
		if (add.density){lines((dx$x-xhist$breaks[1])/(xhist$breaks[2]-xhist$breaks[1]),dx$y,lwd = 2,col = col.density)}
		par(mar=c(4,0,4,1))
		heatbarplot(yhist,xlim=c(0,top),colpal=colpal,horizontal=TRUE,simulate = simulate,daltonize = daltonize,cvd = cvd,alpha = alpha,rev = !rev)
		if (add.density){lines(dy$y,(dy$x-yhist$breaks[1])/(yhist$breaks[2]-yhist$breaks[1]),lwd = 2,col = col.density)}
		par(mar=c(4,1,4,0))
		heatboxplot(y,axes=FALSE,graphics = FALSE,colpal=colpal,ylim=yrange,nolab=TRUE,adjust=adjust,simulate = simulate,daltonize = daltonize,cvd = cvd,alpha = alpha,rev = !rev)
		par(mar=c(1,4,0,4))
		heatboxplot(x,axes=FALSE,horizontal=TRUE,graphics = FALSE,colpal=colpal,ylim=xrange,nolab=TRUE,adjust=adjust,simulate = simulate,daltonize = daltonize,cvd = cvd,alpha = alpha,rev = !rev)
	}
	else {
		if (ab){abline(0,1,col="black",lwd=1)}
		par(mar=c(0,4,1,4))
		barplot(xhist$density, axes=FALSE, ylim=c(0, top), space=0,col="#F0F0F0")
		if (add.density){lines((dx$x-xhist$breaks[1])/(xhist$breaks[2]-xhist$breaks[1]),dx$y,lwd = 2,col = col.density)}
		par(mar=c(4,0,4,1))
		barplot(yhist$density, axes=FALSE, xlim=c(0, top), space=0,col="#F0F0F0", horiz=TRUE)
		if (add.density){lines(dy$y,(dy$x-yhist$breaks[1])/(yhist$breaks[2]-yhist$breaks[1]),lwd = 2,col = col.density)}
		par(mar=c(4,1,4,0))
		boxplot(y,axes=FALSE,col="#F0F0F0",ylim=yrange)
		par(mar=c(1,4,0,4))
		boxplot(x,axes=FALSE,horizontal=TRUE,col="#F0F0F0",ylim=xrange)
	}
	par(def.par)
}


### aliases ###


LSD.comparisonplot = comparisonplot



