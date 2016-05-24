

### heatscatterpoints ###


#' @export
#' @name heatscatterpoints
#' @aliases LSD.heatscatterpoints
#' @title A colored scatterplot based on a two-dimensional Kernel Density Estimation (add to an existing plot)
#' @description Visualize two dimensional data in a three dimensional fashion facilitating a color encoded Kernel Density Estimation (add to an existing plot).
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param pch plotting 'character'. This can either be a single character or an integer code for one of a set of graphics symbols. (see '?pch', to be passed to plot).
#' @param cexplot a numerical value giving the amount by which the points should be magnified relative to the default.
#' @param nrcol a non-negative integer specifying the number of colors to be used (defaults to 100, if not specified).
#' @param grid an integer specifying the size of the grid used for the KDE.
#' @param colpal a character vector containing R built-in color names or a name of a \code{LSD} colorpalette as a character string (see disco() or \code{\link{disco}}) (defaults to "heat", if not specified).
#' @param simulate logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to simulate dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param daltonize logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to enhance dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param cvd character string implying the type of color vision deficiency ("p" for protanope, "d" for deuteranope or "t" for tritanope). 
#' @param alpha alpha value: a two-digit integer between 01 and 99 for color opacity, i.e. appearance of partial or full transparency (usage omitted by default).
#' @param rev logical: if \code{TRUE} (\code{FALSE} by default), a reversed colorpalette is used.
#' @param xlim x limits, standard graphics parameter.
#' @param ylim y limits, standard graphics parameter.
#' @param only a character string which contains 'x' if the density should only be computed for the x axis, 'y' for the y axis (defaults to 'none' for the two-dimensional case).
#' @param add.contour logical: if \code{TRUE} (\code{FALSE} by default), the contour lines are added to the plot.
#' @param nlevels an integer giving the number of levels of the contour lines.
#' @param color.contour R build-in color for the contour lines.
#' @param greyscale logical: if \code{TRUE} (\code{FALSE} by default), the used colorpalette is converted to greyscales.
#' @param log a character string which contains "x" if the x axis is to be logarithmic, "y" if the y axis is to be logarithmic and "xy" or "yx" if both axes are to be logarithmic.
#' @param ... additional parameters to be passed to points and plot.
#' @author Bjoern Schwalb
#' @seealso \code{\link{comparisonplot}}, \code{\link{demotour}}, \code{\link{disco}}, \code{\link{colorpalette}}
#' @note Two-Dimensional Kernel Density Estimation adapted and modified from Venables and Ripley's MASS package (see reference).
#' @references Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S.} Fourth edition. Springer.
#' @examples points = 10^4
#' x = c(rnorm(points/2),rnorm(points/2)+4)
#' y = x + rnorm(points,sd=0.8)
#' x = sign(x)*abs(x)^1.3
#' 
#' plot.new()
#' plot.window(xlim = c(-5,15),ylim = c(-4,8))
#' heatscatterpoints(x,y,add.contour=TRUE,color.contour="green",greyscale=TRUE)
#' axis(1)
#' axis(2)
#' box()
#' @keywords scatterplot, heatcolors


heatscatterpoints = function(x,y,pch = 19,cexplot = 0.5,nrcol = 30,grid = 100,colpal = "heat",simulate = FALSE,daltonize = FALSE,cvd = "p",alpha = NULL,rev = FALSE,xlim = NULL,ylim = NULL,only = "none",add.contour = FALSE,nlevels = 10,color.contour = "black",greyscale = FALSE,log = "",...)
{
	# soundcheck #
	
	if (!is.vector(x) | !is.vector(y)) stop("First two argument must be numeric vectors!")
	if (length(x) != length(y)) stop("Data vectors must be of the same length!")
	sound = which((!(is.na(x)|is.nan(x)|(x==Inf)|(x==-Inf))) & (!(is.na(y)|is.nan(y)|(y==Inf)|(y==-Inf))))
	if (length(sound)==0) stop("There are no valid point pairs to plot!")
	x = x[sound]
	y = y[sound]
	
	# color handling #
	
	colpal = colorpalette(colpal,nrcol,simulate = simulate,daltonize = daltonize,cvd = cvd,alpha = alpha,rev = rev)
	if (greyscale){colpal = convertgrey(colpal)}
	
	# binninmg function #
	
	todiscrete = function(t,tmin,tmax,bins){
		erg = round((t-tmin)/(tmax-tmin)*bins+0.5)
		erg = pmin(pmax(erg,1),bins)
		return(erg)
	}
	
	# kde2d.adj function: adapted and modified from Venables and Ripley's MASS package #
	
	kde2d.adj = function(x,y,h,n = 25,lims = c(range(x),range(y)),only = "none"){
		nx = length(x)
		gx = seq.int(lims[1],lims[2],length.out = n)
		gy = seq.int(lims[3],lims[4],length.out = n)
		bandwidth.nrd.adj = function(x) 
		{
			r = quantile(x,c(0.25,0.75))
			h = (r[2] - r[1])/1.34
			return(4*1.06*min(sqrt(var(x)),h)*length(x)^(-1/5))
		}
		if (missing(h)) {
			bx = bandwidth.nrd.adj(x)
			by = bandwidth.nrd.adj(y)
			if (all(c(bx,by) == 0)){h = rep(0.01,2)} else if (any(c(bx,by) == 0)){h = rep(max(bx,by),2)} else {h = c(bx,by)}
		} else h = rep(h,length.out = 2)
		h = h/4
		ax = outer(gx,x,"-")/h[1]
		ay = outer(gy,y,"-")/h[2]
		norm.ax = dnorm(ax)
		norm.ay = dnorm(ay)
		if (only == "x"){norm.ay = rep(1,length(ay))}
		if (only == "y"){norm.ax = rep(1,length(ax))}
		z = tcrossprod(matrix(norm.ax,,nx),matrix(norm.ay,,nx))/(nx*h[1]*h[2])
		list(x = gx,y = gy,z = z)
	}
	
	# handle 'log' option #
	
	if (log == ""){
		xlog = x
		ylog = y
	} else if (log == "x"){
		xlog = log(x,10)
		ylog = y
	} else if (log == "y"){
		xlog = x
		ylog = log(y,10)
	} else if (log %in% c("xy","yx")){
		xlog = log(x,10)
		ylog = log(y,10)
	}
	
	# estimate two-dimensional KDE for color encoding #

	d = kde2d.adj(xlog,ylog,n=grid,only=only)
	
	# binning #

	xdiscrete = todiscrete(xlog,min(xlog),max(xlog),bins=grid)
	ydiscrete = todiscrete(ylog,min(ylog),max(ylog),bins=grid)
	
	# color assignment #

	getfrommat = function(a){d$z[a[1],a[2]]}
	heatvec = unlist(apply(cbind(xdiscrete,ydiscrete),1,getfrommat))
	coldiscrete = todiscrete(heatvec,min(d$z),max(d$z),bins=nrcol)
	
	# add to existing graphics device #

	points(x,y,col=colpal[coldiscrete],pch=pch,cex=cexplot,...)
	
	# handle 'add.contour' option #

	if (add.contour){contour(d,add=TRUE,nlevels=nlevels,col=color.contour)}
}


### aliases ###


LSD.heatscatterpoints = heatscatterpoints


### heatscatter ###


#' @export
#' @name heatscatter
#' @aliases LSD.heatscatter
#' @title A colored scatterplot based on a two-dimensional Kernel Density Estimation
#' @description Visualize two dimensional data in a three dimensional fashion facilitating a color encoded Kernel Density Estimation.
#' @param x a numeric vector.
#' @param y a numeric vector.
#' @param pch plotting 'character'. This can either be a single character or an integer code for one of a set of graphics symbols. (see '?pch', to be passed to plot).
#' @param cexplot a numerical value giving the amount by which the points should be magnified relative to the default.
#' @param nrcol a non-negative integer specifying the number of colors to be used (defaults to 100, if not specified).
#' @param grid an integer specifying the size of the grid used for the KDE.
#' @param colpal a character vector containing R built-in color names or a name of a \code{LSD} colorpalette as a character string (see disco() or \code{\link{disco}}) (defaults to "heat", if not specified).
#' @param simulate logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to simulate dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param daltonize logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to enhance dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param cvd character string implying the type of color vision deficiency ("p" for protanope, "d" for deuteranope or "t" for tritanope).
#' @param alpha alpha value: a two-digit integer between 01 and 99 for color opacity, i.e. appearance of partial or full transparency (usage omitted by default).
#' @param rev logical: if \code{TRUE} (\code{FALSE} by default), a reversed colorpalette is used.
#' @param xlim x limits, standard graphics parameter.
#' @param ylim y limits, standard graphics parameter.
#' @param xlab x labels, standard graphics parameter.
#' @param ylab y labels, standard graphics parameter.
#' @param main title(s) of the plot, standard graphics parameter.
#' @param cor logical: if \code{TRUE} (\code{FALSE} by default), the correlation is added to the title.
#' @param method a character specifying the correlation method to use ('pearson' (default), 'kendall' or 'spearman').
#' @param only a character string which contains 'x' if the density should only be computed for the x axis, 'y' for the y axis (defaults to 'none' for the two-dimensional case).
#' @param add.contour logical: if \code{TRUE} (\code{FALSE} by default), the contour lines are added to the plot.
#' @param nlevels an integer giving the number of levels of the contour lines.
#' @param color.contour R build-in color for the contour lines.
#' @param greyscale logical: if \code{TRUE} (\code{FALSE} by default), the used colorpalette is converted to greyscales.
#' @param log a character string which contains "x" if the x axis is to be logarithmic, "y" if the y axis is to be logarithmic and "xy" or "yx" if both axes are to be logarithmic.
#' @param ... additional parameters to be passed to points and plot.
#' @author Achim Tresch, Bjoern Schwalb
#' @seealso \code{\link{comparisonplot}}, \code{\link{demotour}}, \code{\link{disco}}, \code{\link{colorpalette}}
#' @note Two-Dimensional Kernel Density Estimation adapted and modified from Venables and Ripley's MASS package (see reference).
#' @references Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S.} Fourth edition. Springer.
#' @examples points = 10^4
#' x = c(rnorm(points/2),rnorm(points/2)+4)
#' y = x + rnorm(points,sd=0.8)
#' x = sign(x)*abs(x)^1.3
#' 
#' heatscatter(x,y)
#' 
#' heatscatter(x,y,colpal="bl2gr2rd",main="bl2gr2rd",cor=FALSE)
#' 
#' heatscatter(x,y,cor=FALSE,add.contour=TRUE,color.contour="red",greyscale=TRUE)
#' 
#' heatscatter(x,y,colpal="spectral",cor=FALSE,add.contour=TRUE)
#' @keywords scatterplot, heatcolors


heatscatter = function(x,y,pch = 19,cexplot = 0.5,nrcol = 30,grid = 100,colpal = "heat",simulate = FALSE,daltonize = FALSE,cvd = "p",alpha = NULL,rev = FALSE,xlim = NULL,ylim = NULL,xlab = NULL,ylab = NULL,main = "heatscatter",cor = FALSE,method = "spearman",only = "none",add.contour = FALSE,nlevels = 10,color.contour = "black",greyscale = FALSE,log = "",...)
{
	# parse variable names #
	
	if (is.null(xlab)){xlab = deparse(substitute(x))}
	if (is.null(ylab)){ylab = deparse(substitute(y))}
	
	# soundcheck #
	
	if (!is.vector(x) | !is.vector(y)) stop("First two argument must be numeric vectors!")
	if (length(x) != length(y)) stop("Data vectors must be of the same length!")
	sound = which((!(is.na(x)|is.nan(x)|(x==Inf)|(x==-Inf))) & (!(is.na(y)|is.nan(y)|(y==Inf)|(y==-Inf))))
	if (length(sound)==0) stop("There are no valid point pairs to plot!")
	x = x[sound]
	y = y[sound]
		
	# handle 'log' option #
	
	if (log == ""){
		valid = 1:length(x)
	} else if (log == "x"){
		valid = which(x > 0)
	} else if (log == "y"){
		valid = which(y > 0)
	} else if (log %in% c("xy","yx")){
		valid = intersect(which(x > 0),which(y > 0))
	}
	x = x[valid]
	y = y[valid]
	
	# handle 'cor' option #
	
	if (cor){main = paste(main," cor = ",round(cor(x,y,method=method),digits=2))}
	
	# handle graphics device  #
	
	plot(x,y,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main="",type = "n",log = log,...)
	heatscatterpoints(x,y,pch = pch,cexplot = cexplot,nrcol = nrcol,grid = grid,colpal = colpal,simulate = simulate,daltonize = daltonize,cvd = cvd,alpha = alpha,rev = rev,xlim = xlim,ylim = ylim,only = only,add.contour = add.contour,nlevels = nlevels,color.contour = color.contour,greyscale = greyscale,log = log,...)
	mtext(paste(main),3,2,cex=1.25)
}


### aliases ###


LSD.heatscatter = heatscatter


### heatpairs ###


#' @export
#' @name heatpairs
#' @aliases LSD.heatpairs
#' @title Pairwise colored scatterplot based on a two-dimensional Kernel Density Estimation
#' @description Pairwise visualization of two dimensional data in a three dimensional fashion facilitating a color encoded Kernel Density Estimation.
#' @param mat a matrix with numerical entries.
#' @param main title(s) of the plot, standard graphics parameter.
#' @param xlim x limits, standard graphics parameter.
#' @param ylim y limits, standard graphics parameter.
#' @param labels a character vector giving the labels to be shown on the diagonal.
#' @param add.points logical: if \code{TRUE} (\code{FALSE} by default), a certain 'group' of points can be colored in all pairwise plots.
#' @param group indices or rownames of 'mat' to be highlighted in all pairwise plots (not necessarily all).
#' @param color.group R build-in color in which the 'group' of points should be highlighted.
#' @param method a character specifying the correlation method to use ('pearson' (default), 'kendall' or 'spearman').
#' @param colpal a character vector containing R built-in color names or a name of a \code{LSD} colorpalette as a character string (see disco() or \code{\link{disco}}) (defaults to "heat", if not specified).
#' @param simulate logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to simulate dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param daltonize logical: if \code{TRUE} (\code{FALSE} by default), a converted colorpalette is used to enhance dichromat vision according to \url{http://www.daltonize.org} (see \code{\link{daltonize}}).
#' @param cvd character string implying the type of color vision deficiency ("p" for protanope, "d" for deuteranope or "t" for tritanope). 
#' @param alpha alpha value: a two-digit integer between 01 and 99 for color opacity, i.e. appearance of partial or full transparency (usage omitted by default).
#' @param rev logical: if \code{TRUE} (\code{FALSE} by default), a reversed colorpalette is used.
#' @param pch plotting 'character'. This can either be a single character or an integer code for one of a set of graphics symbols. (see '?pch', to be passed to plot).
#' @param cexplot a numerical value giving the amount by which the points should be magnified relative to the default.
#' @param cor.cex a numerical value giving the amount by which the correlation characters should be magnified relative to the default.
#' @param nrcol a non-negative integer specifying the number of colors to be used (defaults to 100, if not specified).
#' @param grid an integer specifying the size of the grid used for the KDE.
#' @param only a character string which contains 'x' if the density should only be computed for the x axis, 'y' for the y axis (defaults to 'none' for the two-dimensional case).
#' @param add.contour logical: if \code{TRUE} (\code{FALSE} by default), the contour lines are added to the plot.
#' @param nlevels an integer giving the number of levels of the contour lines.
#' @param color.contour R build-in color for the contour lines.
#' @param greyscale logical: if \code{TRUE} (\code{FALSE} by default), the used colorpalette is converted to greyscales.
#' @param ... additional parameters to be passed to points and plot
#' @author Bjoern Schwalb
#' @seealso \code{\link{comparisonplot}}, \code{\link{demotour}}, \code{\link{disco}}, \code{\link{colorpalette}}
#' @examples points = 10^4
#' x = rnorm(points/2)
#' x = c(x,x+2.5)
#' y = x + rnorm(points,sd=0.75)
#' x = sign(x)*abs(x)^1.3
#' mat = cbind(x,y,x + rnorm(points,sd=0.5))
#' colnames(mat) = c("x","y","z")
#' rownames(mat) = 1:nrow(mat)
#' 
#' heatpairs(mat,labels=c(expression(Xi),expression(Lambda),expression(Delta)))
#' @keywords scatterplot, heatcolors


heatpairs = function(mat,main = "heatpairs",xlim = NULL,ylim = NULL,labels = NULL,add.points = FALSE,group = NULL,color.group = "magenta",method = "spearman",colpal = "heat",simulate = FALSE,daltonize = FALSE,cvd = "p",alpha = NULL,rev = FALSE,pch=19,cexplot=0.5,cor.cex = 2.5,nrcol=30,grid=100,only = "none",add.contour = FALSE,nlevels = 10,color.contour = "black",greyscale = FALSE,...)
{
	if (!is.matrix(mat)) stop("First argument must be a matrix !")
	if (is.null(xlim)){xlim = c(min(mat,na.rm=TRUE),max(mat,na.rm=TRUE))}
	if (is.null(ylim)){ylim = c(min(mat,na.rm=TRUE),max(mat,na.rm=TRUE))}
	if(is.null(labels)){labels = colnames(mat)}
	pairs(mat,labels=labels,xlim=xlim,ylim=ylim,lower.panel=function(x,y,...){text(sum(xlim)/2,sum(ylim)/2,round(cor(x,y,method=method,use="na.or.complete"),digits=2),cex=cor.cex)},main=main,upper.panel=function(x,y,...){heatscatterpoints(x,y,colpal=colpal,pch=pch,cexplot=cexplot,nrcol=nrcol,grid=grid,simulate=simulate,daltonize=daltonize,cvd=cvd,alpha=alpha,rev=rev,only=only,add.contour=add.contour,nlevels=nlevels,color.contour=color.contour,greyscale=greyscale,...);abline(a=0,b=1);if (add.points){points(x[rownames(mat) %in% group],y[rownames(mat) %in% group],col=color.group,...)}},...)
}


### aliases ###


LSD.heatpairs = heatpairs



