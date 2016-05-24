

### linesplot ###


#' @export
#' @name linesplot
#' @aliases LSD.linesplot
#' @title One-dimensional scatterplot
#' @description Visualize one-dimensional data in its every detail.
#' @param x numeric data as vector, matrix, list or data.frame.
#' @param labels a character vector of labels.
#' @param col a R build-in color.
#' @param cols a character vector of R build-in colors.
#' @param alpha alpha value: a two-digit integer between 01 and 99 for color opacity, i.e. appearance of partial or full transparency (usage omitted by default).
#' @param xlim x limits, standard graphics parameter.
#' @param ylim y limits, standard graphics parameter.
#' @param xlab x label, standard graphics parameter.
#' @param ylab y label, standard graphics parameter.
#' @param las las=1: horizontal text, las=2: vertical text (x-axis labels).
#' @param outline logical: if \code{TRUE} (by default), outliers are plotted.
#' @param cexbox a numerical value giving the amount by which the boxes should be magnified relative to the default.
#' @param addboxes logical: if \code{TRUE} (\code{FALSE} by default), boxplots be added to the plot.
#' @param border a R build-in color for the box and the whiskers (if \code{addboxes = TRUE}).
#' @param range this determines how far the plot whiskers extend out from the box.
#' @param lwd linewidth of the box and whiskers.
#' @param main title(s) of the plot, standard graphics parameter.
#' @param ... additional parameters to be passed to points and plot.
#' @author Bjoern Schwalb
#' @seealso \code{\link{demotour}}, \code{\link{disco}}, \code{\link{colorpalette}}
#' @examples l = list()
#' for (i in 1:10){l[[i]] = rnorm(100,sqrt(i^2.5),1+i/2)}
#' 
#' linesplot(l,alpha=25,border="darkred",addboxes = TRUE,outline=FALSE)
#' @keywords boxplotlike


linesplot = function(x,labels = NULL,col = "black",cols = NULL,alpha = 25,xlim = NULL,ylim = NULL,xlab = NULL,ylab = "",las = 1,outline = TRUE,cexbox = 0.6,addboxes = FALSE,border = "black",range = 1.5,lwd = 1.5,main = "LSD.linesplot",...)
{
	if (!is.vector(x) & !is.matrix(x) & !is.list(x) & !is.data.frame(x)) stop("x must be a vector, matrix, list or a data.frame!")
	if (!is.list(x) & !is.matrix(x) & !is.data.frame(x)){x = cbind(x)}
	if (is.data.frame(x)){x = as.list(x)}
	if (is.matrix(x)){x = as.list(as.data.frame(x))}
	if (is.null(labels)){labels = labels(x)}
	if (is.null(cols)){cols = rep(convertcolor(col,alpha),length(x))} else{cols = convertcolor(cols,alpha)}
	if (!is.null(xlim)){print("xlim argument will be ignored!")}
	if (is.null(ylim)){ylim = range(unlist(x),na.rm=TRUE)}
	if (!is.null(xlab)){print("xlab argument will be ignored! Use labels instead!")}
	if (length(x) == 1){cexbox = 0.4}
	if (length(x) > 1){boxwex = cexbox} else{boxwex = NULL}
	plot(1,col="white",xlim=c(0.5,length(x)+0.5),ylim=ylim,xlab="",xaxt="n",ylab=ylab,main="",...)
	mtext(paste(main),3,2,cex=1.25)
	par.axis.default = par("cex.axis")
	par(cex.axis = par("cex.lab"))
	axis(1,at=seq(1,length(x),1),labels=labels,las=las)
	for (j in 1:length(x)){segments(j-cexbox/2,x[[j]],j+cexbox/2,x[[j]],lwd=2,col=cols[j])}
	if (addboxes){boxplot(x,add=TRUE,col="transparent",xlim=c(0.5,length(x)+0.5),ylim=ylim,width=NULL,boxwex=boxwex,outline=outline,axes=FALSE,border=border,lwd=lwd,range=range)}                                 
	par(cex.axis = par.axis.default)                                 
}


### aliases ###


LSD.linesplot = linesplot



