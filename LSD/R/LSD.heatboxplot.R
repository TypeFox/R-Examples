

### densitylane ###


#' @name densitylane
#' @aliases LSD.densitylane
#' @title Visualize a density in a rectangular fashion
#' @description Add a color stripe to an existing plot based on a kernel density estimate.
#' @param x density$x values of a density object.
#' @param y density$y values of a density object.
#' @param pos the x co-ordinate of the lane (mid point).
#' @param width a numeric value giving the width of the lane.
#' @param colpal a character vector containing R built-in color names or a name of a \code{LSD} colorpalette as a character string (see disco() or \code{\link{disco}}) (defaults to "heat", if not specified).
#' @param rev logical: if \code{TRUE} (\code{FALSE} by default), a reversed colorpalette is used.
#' @param simulate logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to simulate dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param daltonize logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to enhance dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param cvd character string implying the type of color vision deficiency ("p" for protanope, "d" for deuteranope or "t" for tritanope).
#' @param alpha alpha value: a two-digit integer between 01 and 99 for color opacity, i.e. appearance of partial or full transparency (usage omitted by default).
#' @param horizontal logical: if \code{TRUE} (\code{FALSE} by default), rotation of 90 degrees is applied.
#' @param nrcol a non-negative integer specifying the number of colors to be used (defaults to 100, if not specified).
#' @author Bjoern Schwalb
#' @seealso \code{\link{comparisonplot}}, \code{\link{demotour}}, \code{\link{disco}}, \code{\link{colorpalette}}
#' @keywords density


densitylane = function(x,y,pos = 1,width = 0.4,colpal = "standard",rev = FALSE,simulate = FALSE,daltonize = FALSE,cvd = "p",alpha = NULL,horizontal = horizontal,nrcol = 75)
{
	if (!is.vector(x) | !is.vector(y)) stop("First two argument must be vectors!")
	if (length(x) != length(y)) stop("Data vectors must be of the same length!")
	colpalette = rev(colorpalette(colpal,nrcol,simulate = simulate,daltonize = daltonize,cvd = cvd,alpha = alpha,rev = rev))
	ycol  = pmin(pmax(round((y-min(y))/(max(y)-min(y))*nrcol+0.5),1),nrcol)
	if (horizontal){for (j in 1:(length(x)-1)){rect(x[j],pos-width/2,x[j+1],pos+width/2,col=colpalette[ycol[j]],border=NA)}}
	else {for (j in 1:(length(x)-1)){rect(pos-width/2,x[j],pos+width/2,x[j+1],col=colpalette[ycol[j]],border=NA)}}
}


### aliases ###


LSD.densitylane = densitylane


### heatboxplot ###


#' @export
#' @name heatboxplot
#' @aliases LSD.heatboxplot
#' @title Heatboxplot: a colored boxplot
#' @description A boxplot with an additional color stripe based on a kernel density estimate.
#' @param x data as vector, matrix, list or data.frame.
#' @param horizontal logical: if \code{TRUE} (\code{FALSE} by default), rotation of 90 degrees is applied.
#' @param add logical: if \code{TRUE} (\code{FALSE} by default), the boxplot is added to an existing plot.
#' @param colpal a character vector containing R built-in color names or a name of a \code{LSD} colorpalette as a character string (see disco() or \code{\link{disco}}) (defaults to "heat", if not specified).
#' @param rev logical: if \code{TRUE} (\code{FALSE} by default), a reversed colorpalette is used.
#' @param simulate logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to simulate dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param daltonize logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to enhance dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param cvd character string implying the type of color vision deficiency ("p" for protanope, "d" for deuteranope or "t" for tritanope).
#' @param alpha alpha value: a two-digit integer between 01 and 99 for color opacity, i.e. appearance of partial or full transparency (usage omitted by default).
#' @param colpals a character vector containing names of \code{LSD} colorpalettes (see disco() or \code{\link{disco}}).
#' @param nrcol a non-negative integer specifying the number of colors to be used (defaults to 100, if not specified).
#' @param lwd linewidth of the box and whiskers.
#' @param axes logical: if \code{TRUE} (by default), the axes are plotted.
#' @param labels a character vector of labels.
#' @param xlim x limits, standard graphics parameter.
#' @param ylim y limits, standard graphics parameter.
#' @param xlab x label, standard graphics parameter.
#' @param ylab y label, standard graphics parameter.
#' @param main title(s) of the plot, standard graphics parameter.
#' @param nolab logical: if \code{TRUE} (\code{FALSE} by default), the title and ylab are suppressed.
#' @param outline logical: if \code{TRUE} (by default), outliers are plotted.
#' @param boxonly logical: if \code{TRUE} (\code{FALSE} by default), the density is only be plotted in the box.
#' @param adjust a numeric value giving the scaling factor for the used bandwidth (defaults to 1).
#' @param quant.from a numeric value (between 0 and 1) giving the quantile from which the density lane should be plotted.
#' @param quant.to a numeric value (between 0 and 1) giving the quantile to which the density lane should be plotted.
#' @param range a numeric value to determine how far the plot whiskers extend out from the box.
#' @param border an R build-in color for the box and whiskers.
#' @param plot.boxplot logical: if \code{TRUE} (by default), the boxplot is added to the density.
#' @param add.quartiles if \code{TRUE} (\code{FALSE} by default), only the box of the boxplot is added (if \code{plot.boxplot = FALSE}).
#' @param add.box logical: if \code{TRUE} (by default), the box is added to the plot.
#' @param n.density an integer specifying the number of equally spaced points at which the density is to be estimated.
#' @param cexbox a numerical value giving the amount by which the boxes should be magnified relative to the default.
#' @param ... additional parameters to be passed to points and plot.
#' @author Bjoern Schwalb
#' @seealso \code{\link{comparisonplot}}, \code{\link{demotour}}, \code{\link{disco}}, \code{\link{colorpalette}}
#' @examples f = c(rnorm(200),rnorm(200)+4)
#' h = rf(500,15,15)*10
#' g = rnorm(300)+1
#' 
#' heatboxplot(h)
#' 
#' heatboxplot(list(f=f,g=g),colpals=c("rdpu","greens"),labels=c("bimodal","unimodal"))
#' @keywords boxplot


heatboxplot = function(x,horizontal = FALSE,add = FALSE,colpal = "standard",rev = FALSE,simulate = FALSE,daltonize = FALSE,cvd = "p",alpha = NULL,colpals = NULL,nrcol = 75,lwd = 1.75,axes = TRUE,labels = NULL,xlim = NULL,ylim = NULL,xlab = NULL,ylab = "",main = "heatboxplot",nolab = FALSE,outline = TRUE,boxonly = FALSE,adjust = 1,quant.from = 0.25,quant.to = 0.75,range = 1.5,border = "black",plot.boxplot = TRUE,add.quartiles = TRUE,add.box = FALSE,n.density = 1024,cexbox = 0.6,...)
{
	if (!is.vector(x) & !is.matrix(x) & !is.list(x) & !is.data.frame(x)) stop("x must be a vector, matrix, list or a data.frame!")
	if (!is.list(x) & !is.matrix(x) & !is.data.frame(x)){x = cbind(x)}
	if (is.data.frame(x)){x = as.list(x)}
	if (is.matrix(x)){x = as.list(as.data.frame(x))}
	if (is.null(colpals)){colpals = rep(colpal,length(x))}
	if (!is.null(xlim)){print("xlim argument will be ignored!")}
	if (is.null(ylim)){ylim = range(unlist(x))}
	if (!is.null(xlab)){print("xlab argument will be ignored! Use labels instead!")} # necessary ? 
	if (length(x) == 1){xlab = labels}
	if (length(x) == 1){cexbox = 0.4}
	if (length(x) > 1){boxwex = cexbox} else{boxwex = NULL}
	if (nolab){ylab = NULL}
	if (horizontal){labdummy = ylab
		ylab = xlab
		xlab = labdummy}
	limlist = list()
	qlimlist = list()
	for (i in 1:length(x)){limlist[[i]] = c(quantile(x[[i]],quant.from),quantile(x[[i]],quant.to))}
	for (i in 1:length(x)){qlimlist[[i]] = c(quantile(x[[i]],0.25),quantile(x[[i]],0.75))}
	boxplot(x,border="white",add=add,horizontal=horizontal,axes=axes,xlim=c(0.5,length(x)+0.5),ylim=ylim,width=NULL,boxwex=boxwex,main="",ylab=ylab,xlab=xlab,names=labels,outline=outline,range=range,lwd=lwd,...)
	if (!nolab){mtext(paste(main),3,2,cex=1.25)}
	for (j in 1:length(x)){xrel = x[[j]][x[[j]] >= ylim[1] & x[[j]] <= ylim[2]]
		dx = density(xrel,n=n.density,adjust=adjust)
		if(boxonly){dxx = dx$x[dx$x >= limlist[[j]][1] & dx$x <= limlist[[j]][2]]
			dxy = dx$y[dx$x >= limlist[[j]][1] & dx$x <= limlist[[j]][2]]}
		else {dxx = dx$x
			dxy = dx$y}
		densitylane(dxx,dxy,colpal=colpals[[j]],horizontal=horizontal,pos=j,width=cexbox,nrcol=nrcol,simulate = simulate,daltonize = daltonize,cvd = cvd,alpha=alpha,rev = rev)}
	if (plot.boxplot){boxplot(x,add=TRUE,horizontal=horizontal,axes=axes,col="transparent",xlim=c(0.5,length(x)+0.5),ylim=ylim,width=NULL,boxwex=boxwex,names=labels,outline=outline,range=range,lwd=lwd,...)}                                 
	if (!plot.boxplot){if (add.quartiles) {if (horizontal){ for (j in 1:length(x)){rect(qlimlist[[j]][1],j-cexbox/2,qlimlist[[j]][2],j+cexbox/2,border="grey40",lwd=lwd)} }
			else { for (j in 1:length(x)){rect(j-cexbox/2,qlimlist[[j]][1],j+cexbox/2,qlimlist[[j]][2],border="grey40",lwd=lwd)} }}}
	if (add.box){if (horizontal){ for (j in 1:length(x)){rect(limlist[[j]][1],j-cexbox/2,limlist[[j]][2],j+cexbox/2,border=border,lwd=lwd)} }
		else { for (j in 1:length(x)){rect(j-cexbox/2,limlist[[j]][1],j+cexbox/2,limlist[[j]][2],border=border,lwd=lwd)} }}
}


### aliases ###


LSD.heatboxplot = heatboxplot



