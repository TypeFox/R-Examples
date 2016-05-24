

### singlefusionplot ###


#' @export
#' @name singlefusionplot
#' @aliases LSD.singlefusionplot
#' @title Visualize two-dimensional data clusters (add to an existing plot)
#' @description Depict a numeric matrix or list utilizing the underlying distribution quantiles of one dimension in a color encoded fashion (add to an existing plot).
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param fromto a numeric vector containing the range of quantiles (between 0 and 1) to be plotted.
#' @param colpal a character vector containing R built-in color names or a name of a \code{LSD} colorpalette as a character string (see disco() or \code{\link{disco}}).
#' @param simulate logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to simulate dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param daltonize logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to enhance dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param cvd character string implying the type of color vision deficiency ("p" for protanope, "d" for deuteranope or "t" for tritanope).
#' @param nrcol a non-negative integer specifying the number of colors to be used (defaults to 25, if not specified).
#' @param outer.col R built-in color to be used for outlier lines (lines outside of 'fromto').
#' @param rev logical: if \code{TRUE} (\code{FALSE} by default), a reversed colorpalette is used.
#' @param alpha alpha value: a two-digit integer between 01 and 99 for color opacity, i.e. appearance of partial or full transparency (usage omitted by default).
#' @param quartiles.col a character vector containing three R built-in colors for quartile lines (c('0.25','0.5','0.75')).
#' @param add.quartiles logical: if \code{TRUE} (by default), lines are plotted corresponding to the quartiles.
#' @author Achim Tresch, Bjoern Schwalb
#' @seealso \code{\link{fusionplot}}, \code{\link{align}}, \code{\link{demotour}}, \code{\link{disco}}, \code{\link{colorpalette}}
#' @examples x = 1:1000/300
#' y = rnorm(1000)+sin(2*x)*3
#' 
#' emptyplot(xlim = range(x),ylim = range(y))
#' singlefusionplot(x,y,colpal = "ylgnbu")
#' axis(1)
#' axis(2)
#' box()
#' @keywords cluster


singlefusionplot = function(x,y,fromto = c(0.05,0.95),colpal = "standardheat",simulate = FALSE,daltonize = FALSE,cvd = "p",nrcol = 25,outer.col = "grey",rev = FALSE,alpha = NULL,quartiles.col = c("grey","black","grey"),add.quartiles = TRUE)
{
	# kernel function #
	
	kernelf = function(y0,x0,x,y,width=0.1){
		wx = 1/sqrt(2*pi)/width*exp(-(x0-x)^2/width^2/2)
		res = sum(wx[y<=y0]) / sum(wx)
		return(res)
	}
	
	# quantile function #
	
	quantf = function(x0,x,y,width=0.1,quantvector=seq(0,1,length=5)){
		y0 = seq(min(y),max(y),length=100)
		quants = sapply(y0,kernelf,x0,x,y,width)
		hilf = approxfun(quants,y0,rule=2)
		return(hilf(quantvector))
	}
	
	# calculate quantiles #
	
	xseq = seq(min(x),max(x),length=200)
	nrquants=2*nrcol
	quantvector = seq(fromto[1],fromto[2],length=nrquants)
	
	qline = sapply(xseq,quantf,x,y,width=0.1,quantvector=quantvector)
	
	# provide 'colpal' via colorpalette #
	
	colpal = colorpalette(colpal,nrcol,simulate = simulate,daltonize = daltonize,cvd = cvd,alpha = alpha,rev = rev)
	colpal = c(rev(colpal),colpal)
	
	quantvec = seq(0,1,length=nrquants)
	
	for (j in 1:(nrquants - 1)){
		polygon(c(xseq,rev(xseq)),c(qline[j,],rev(qline[j+1,])),col = colpal[j],lty = 0)
	}
	points(x,y,col = outer.col,pch = 19,cex = 0.3)
	
	# wrapper for the lines function #
	
	drawline = function(y,col="black",lwd=1,lty=1){lines(xseq,y,type="l",col=col,lwd=lwd,lty=lty)}
	
	# add lines corresponding to the quartiles #
	
	if (add.quartiles){
		quantvector = seq(0.25,0.75,length=3)
		qline = sapply(xseq,quantf,x,y,width=0.1,quantvector=quantvector)
		
		drawline(qline[2,],col=quartiles.col[2],lwd=2)
		drawline(qline[1,],col=quartiles.col[1],lwd=2)
		drawline(qline[3,],col=quartiles.col[3],lwd=2)
	}
}


# alias #


LSD.singlefusionplot = singlefusionplot


### fusionplot ###


#' @export
#' @name fusionplot
#' @aliases LSD.fusionplot
#' @title Visualize two-dimensional data clusters
#' @description Depict a numeric matrix or list utilizing the underlying distribution quantiles of one dimension in a color encoded fashion.
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param label a character vector assigning rows/elements of 'input' to clusters (if specified, multiple clusters can be depicted in different colors and/or subsequent plots).
#' @param main title(s) of the plot, standard graphics parameter.
#' @param xlim x limits, standard graphics parameter.
#' @param ylim y limits, standard graphics parameter.
#' @param fromto a numeric vector containing the range of quantiles (between 0 and 1) to be plotted.
#' @param colpal a character vector containing R built-in color names or a name of a \code{LSD} colorpalette as a character string (see disco() or \code{\link{disco}}).
#' @param simulate logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to simulate dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param daltonize logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to enhance dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param cvd character string implying the type of color vision deficiency ("p" for protanope, "d" for deuteranope or "t" for tritanope).
#' @param nrcol a non-negative integer specifying the number of colors to be used (defaults to 25, if not specified). 
#' @param outer.col R built-in color to be used for outlier lines (lines outside of 'fromto').
#' @param quartiles.col a character vector containing three R built-in colors for quartile lines (c('0.25','0.5','0.75')).
#' @param add.quartiles logical: if \code{TRUE} (by default), lines are plotted corresponding to the quartiles.
#' @param separate if \code{TRUE} (by default), different clusters are depicted in subsequent plots.
#' @param rev logical: if \code{TRUE} (\code{FALSE} by default), a reversed colorpalette is used.
#' @param size logical: if \code{TRUE} (by default), the size of each cluster is added to the title of the respective plot.
#' @param alpha alpha value: a two-digit integer between 01 and 99 for color opacity, i.e. appearance of partial or full transparency (usage omitted by default).
#' @param axes logical: if \code{TRUE} (by default), a box and axes are added to the plot (if \code{FALSE}, custom specification of axes can be achieved via basic R graphics functions).
#' @param ... additional parameters to be passed to points and plot.
#' @author Achim Tresch, Bjoern Schwalb
#' @seealso \code{\link{singlefusionplot}}, \code{\link{align}}, \code{\link{demotour}}, \code{\link{disco}}, \code{\link{colorpalette}}
#' @examples nr = 750
#' x = 1:nr/300
#' y = c(rnorm(nr)+sin(2*x)*3,rnorm(nr)+sin(2*x+pi/2)*3)
#' x = c(x,x)
#' 
#' labs = paste("cluster",c(rep(c(1,2),each = nr)))
#' colpals = c("oranges","pubu")
#' qcol = c("transparent","black","transparent")
#' fusionplot(x,y,labs,separate=FALSE,colpal=colpals,alpha=75,quartiles.col = qcol)
#' @keywords cluster


fusionplot = function(x,y,label = NULL,main = NULL,xlim = NULL,ylim = NULL,fromto = c(0.05,0.95),colpal = "standardheat",simulate = FALSE,daltonize = FALSE,cvd = "p",nrcol = 25,outer.col = "lightgrey",quartiles.col = c("grey","black","grey"),add.quartiles = TRUE,separate = TRUE,rev = FALSE,size = TRUE,alpha = NULL,axes = TRUE,...)
{
	if (is.null(xlim)){xlim=c(min(x),max(x))}
	maxp = xlim[2]
	minp = xlim[1]
	if (is.null(ylim)){ylim=c(min(y),max(y))}
	
	# one cluster (i.e. one plot), if label = NULL #
	
	if (is.null(label)){
		plot.new()
		plot.window(xlim = xlim,ylim = ylim,...)
		if (size){
			main = paste(main,"  ( #",length(x)," )")
		}
		title(main)
		if (axes){
			axis(1,...)
			axis(2)
			box()
		}
		singlefusionplot(x=x,y=y,fromto=fromto,colpal=colpal,simulate=simulate,daltonize=daltonize,cvd=cvd,nrcol=nrcol,outer.col=outer.col,add.quartiles=add.quartiles,quartiles.col=quartiles.col,rev=rev,alpha=alpha)
	}
	
	# multiple clusters, if label is specified #
	
	if (!is.null(label)) {
		clusternames = sort(unique(label))
		nrclusters = length(clusternames)
		clustersets = split(1:length(x), factor(label))
		if (!is.list(colpal)) colpal = as.list(colpal)
		if (length(colpal) < nrclusters) colpal = rep(colpal, nrclusters)
		
		# multiple clusters in one plots #
		
		if (separate == FALSE){
			plot.new()
			plot.window(xlim = xlim,ylim = ylim,...)
			if (size){main = paste(main,"  ( #",length(x)," )")}
			title(main)
			if (axes){
				axis(1,...)
				axis(2)
				box()
			}
		}
		
		# multiple clusters in subsequent plots #
		
		if (separate == TRUE) par(mfrow = windowxy(nrclusters))
		for (j in seq(clusternames)){
			if (separate == TRUE){
				if (length(main) == length(clustersets[[j]])) clustermain = main[j]	else clustermain = paste(main,clusternames[j])
				plot.new()
				plot.window(xlim = xlim,ylim = ylim,...)
				if (size){clustermain = paste(clustermain,"  ( #",length(clustersets[[j]])," )")}
				title(clustermain)
				if (axes){
					axis(1,...)
					axis(2)
					box()
				}
			}
			singlefusionplot(x=x[clustersets[[j]]],y=y[clustersets[[j]]],fromto=fromto,colpal=colpal[[j]],simulate=simulate,daltonize=daltonize,cvd=cvd,nrcol=nrcol,outer.col=outer.col,add.quartiles=add.quartiles,quartiles.col=quartiles.col,rev=rev,alpha=alpha)
		}
	}
}


# alias #


LSD.fusionplot = fusionplot



