

### singleclusterplot ###


#' @export
#' @name singleclusterplot
#' @aliases LSD.singleclusterplot
#' @title Visualize two-dimensional data clusters (add to an existing plot)
#' @description Depict a numeric matrix or list utilizing the underlying distribution quantiles of one dimension in a color encoded fashion (add to an existing plot).
#' @param input matrix or list with numerical entries.
#' @param at a integer vector containing the x-positions corresponding to the columns of 'input'.
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
#' @seealso \code{\link{clusterplot}}, \code{\link{align}}, \code{\link{demotour}}, \code{\link{disco}}, \code{\link{colorpalette}}
#' @examples samples = 100
#' probes = 200
#' clus = matrix(rnorm(probes*samples,sd=1),ncol=probes)
#' 
#' clus = rbind(
#' 	t(t(clus)+sin(1:probes/10))+1:nrow(clus)/samples,
#' 	t(t(clus)+sin(pi/2+1:probes/10))+1:nrow(clus)/samples)
#' 
#' emptyplot(xlim = c(1,ncol(clus)),ylim = range(clus))
#' singleclusterplot(clus)
#' axis(1)
#' axis(2)
#' box()
#' @keywords cluster


singleclusterplot = function(input,at = NULL,fromto = c(0.05,0.95),colpal = "standardheat",simulate = FALSE,daltonize = FALSE,cvd = "p",nrcol = 25,outer.col = "lightgrey",rev = FALSE,alpha = NULL,quartiles.col = c("grey","black","grey"),add.quartiles = TRUE)
{
	# stops execution, if 'input' is neither a list nor a matrix and executes an error action #
	
	if (!is.matrix(input) & !is.list(input)) stop("'input' must be a matrix or a list !")
	
	# define x-positions corresponding to the columns/elements of 'input' #
	
	if (is.null(at)) if (is.matrix(input)){at=1:ncol(input)} else if (is.list(input)){at=1:length(input)}
	
	# preliminaries #
	
	probes = length(at)
	
	# wrapper for the lines function #
	
	drawline = function(y,col="black",lwd=1,lty=1){lines(at[1:length(y)],y,type="l",col=col,lwd=lwd,lty=lty)}
	
	# optional drawing of outlier lines #
	
	if (is.matrix(input)){if (outer.col!="none") apply(input,1,drawline,col=outer.col)}
	
	# provide 'colpal' via colorpalette #
	
	colpal = colorpalette(colpal,nrcol,simulate = simulate,daltonize = daltonize,cvd = cvd,alpha = alpha,rev = rev)
	colpal = c(rev(colpal),colpal)

	# determine quantiles among columns/elements of 'input' #

	if (is.matrix(input)){
		qline = apply(input,2,quantile,probs=seq(fromto[1],fromto[2],length=(length(colpal)+1)),na.rm=TRUE)
	} else if (is.list(input)){
		qline = lapply(input,quantile,probs=seq(fromto[1],fromto[2],length=(length(colpal)+1)),na.rm=TRUE)
		qline = sapply(qline,c)
	}
	
	# plot polygons according to qline #

	for (j in 1:length(colpal)){polygon(at[c(1:probes,probes:1)],c(qline[j,],qline[j+1,probes:1]),col = colpal[j],lty=0)}
	
	# add lines corresponding to the quartiles #

	if (add.quartiles){
		if (is.matrix(input)){
			drawline(apply(input,2,quantile,probs=0.5,na.rm=TRUE),col=quartiles.col[2],lwd=2)
			drawline(apply(input,2,quantile,probs=0.25,na.rm=TRUE),col=quartiles.col[1],lwd=2)
			drawline(apply(input,2,quantile,probs=0.75,na.rm=TRUE),col=quartiles.col[3],lwd=2)
		} else if (is.list(input)){
			drawline(sapply(input,quantile,probs=0.5,na.rm=TRUE),col=quartiles.col[2],lwd=2)
			drawline(sapply(input,quantile,probs=0.25,na.rm=TRUE),col=quartiles.col[1],lwd=2)
			drawline(sapply(input,quantile,probs=0.75,na.rm=TRUE),col=quartiles.col[3],lwd=2)
			
		}
	}
}


# alias #


LSD.singleclusterplot = singleclusterplot


### clusterplot ###


#' @export
#' @name clusterplot
#' @aliases LSD.clusterplot
#' @title Visualize two-dimensional data clusters
#' @description Depict a numeric matrix or list utilizing the underlying distribution quantiles of one dimension in a color encoded fashion.
#' @param input matrix or list with numerical entries.
#' @param label a character vector assigning rows/elements of 'input' to clusters (if specified, multiple clusters can be depicted in different colors and/or subsequent plots).
#' @param at a integer vector containing the x-positions corresponding to the columns of 'input'.
#' @param main title(s) of the plot, standard graphics parameter.
#' @param xlim x limits, standard graphics parameter.
#' @param ylim y limits, standard graphics parameter.
#' @param xlabels a character vector containing labels for the x-axis.
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
#' @seealso \code{\link{singleclusterplot}}, \code{\link{align}}, \code{\link{demotour}}, \code{\link{disco}}, \code{\link{colorpalette}}
#' @examples samples = 100
#' probes = 75
#' at = 1:probes
#' clus = matrix(rnorm(probes*samples,sd=1),ncol=probes)
#' 
#' clus = rbind(
#' 	t(t(clus)+sin(1:probes/10))+1:nrow(clus)/samples,
#' 	t(t(clus)+sin(pi/2+1:probes/10))+1:nrow(clus)/samples)
#' 
#' quartiles.col = c("transparent","black","transparent")
#' colpal = c("standardheat","crazyblue","crazyred","crazygreen")
#' 
#' labs = paste("cluster",kmeans(clus,4)$cluster)
#' clusterplot(clus,fromto=c(0,1))
#' 
#' clusterplot(clus,labs,separate=FALSE,xaxt="n",fromto=c(0.4,0.6),colpal=colpal,
#' 	outer.col="none",ylim=c(-2,3),quartiles.col = quartiles.col)
#'  
#' clusterplot(clus,labs,colpal=colpal)
#' 
#' labs = paste("cluster",kmeans(clus,2)$cluster)
#' colpal = c("greens","purples")
#' clusterplot(clus,labs,separate=FALSE,xaxt="n",fromto=c(0.3,0.7),colpal=colpal,
#' 	outer.col="none",ylim=c(-1,2),alpha=50,quartiles.col = quartiles.col)
#' @keywords cluster


clusterplot = function(input,label = NULL,at = NULL,main = NULL,xlim = NULL,ylim = NULL,xlabels = NULL,fromto = c(0.05,0.95),colpal = "standardheat",simulate = FALSE,daltonize = FALSE,cvd = "p",nrcol = 25,outer.col = "lightgrey",quartiles.col = c("grey","black","grey"),add.quartiles = TRUE,separate = TRUE,rev = FALSE,size = TRUE,alpha = NULL,axes = TRUE,...)
{
	# stops execution, if 'input' is neither a list nor a matrix and executes an error action #
	
	if (!is.matrix(input) & !is.list(input)) stop("'input' must be a matrix or a list !")
	
	# define x-positions corresponding to the columns of 'input' #
	
	if (is.null(at)) if (is.matrix(input)){at=1:ncol(input)} else if (is.list(input)){at=1:length(input)}	
	
	# preliminaries #
	
	probes = length(at)
	if (is.null(xlim)){xlim=c(min(at),max(at))}
	maxp = xlim[2]
	minp = xlim[1]
	if (is.null(ylim)) if (is.matrix(input)){ylim=c(min(input,na.rm=TRUE),max(input,na.rm=TRUE))} else if (is.list(input)){ylim=c(min(unlist(input),na.rm=TRUE),max(unlist(input),na.rm=TRUE))}
	if (is.null(xlabels)) xlabels = 1:length(at)

	# one cluster (i.e. one plot), if label = NULL #
	
	if (is.null(label)){
		plot.new()
		plot.window(xlim = xlim,ylim = ylim,...)
		if (size){
			if (is.matrix(input)){
				main = paste(main,"  ( #",nrow(input)," )")
			} else if (is.list(input)){
				input.length.range = range(as.numeric(summary(input)[,"Length"]))
				main = paste(main,"  ( #",input.length.range[1],"-",input.length.range[2]," )")
			}
		}
		title(main)
		if (axes){
			axis(1,at=at,labels=xlabels,...)
			axis(2)
			box()
		}
		singleclusterplot(input=input,at=at,fromto=fromto,colpal=colpal,simulate=simulate,daltonize=daltonize,cvd=cvd,nrcol=nrcol,outer.col=outer.col,add.quartiles=add.quartiles,quartiles.col=quartiles.col,rev=rev,alpha=alpha)
	}

	# multiple clusters, if label is specified #
	
	if (!is.null(label)) {
		clusternames = sort(unique(label))
		nrclusters = length(clusternames)
		if (!is.matrix(input)) stop("'input' must be a matrix, if 'label' is specified !")
		clustersets = split(1:nrow(input), factor(label))
		if (!is.list(colpal)) colpal = as.list(colpal)
		if (length(colpal) < nrclusters) colpal = rep(colpal, nrclusters)
		
		# multiple clusters in one plots #
		
		if (separate == FALSE){
			plot.new()
			plot.window(xlim = xlim,ylim = ylim,...)
			if (size){main = paste(main,"  ( #",nrow(input)," )")}
			title(main)
			if (axes){
				axis(1,at=at,labels=xlabels,...)
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
					axis(1,at=at,labels=xlabels,...)
					axis(2)
					box()
				}
			}
			singleclusterplot(input=input[clustersets[[j]],,drop = FALSE],at=at,fromto=fromto,colpal=colpal[[j]],simulate=simulate,daltonize=daltonize,cvd=cvd,nrcol=nrcol,outer.col=outer.col,add.quartiles=add.quartiles,quartiles.col=quartiles.col,rev=rev,alpha=alpha)
		}
	}
}


# alias #


LSD.clusterplot = clusterplot



