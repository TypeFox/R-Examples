

### plotmatrix ###


#' @export 
#' @name plotmatrix
#' @aliases LSD.plotmatrix
#' @title Visualize two-dimensional data
#' @description Plot the rows of a matrix as lines along the cols.
#' @param input a matrix with numerical entries.
#' @param xlim x limits, standard graphics parameter.
#' @param ylim y limits, standard graphics parameter.
#' @param xlab x lab, standard graphics parameter.
#' @param ylab y lab, standard graphics parameter.
#' @param main title of the plot, standard graphics parameter.
#' @param type what 'type' of plot should be drawn (to be passed to points).
#' @param lwd a positive number giving the line width.
#' @param at a integer vector containing the x-positions corresponding to the columns of 'input'.
#' @param xlabels a character vector containing labels for the x-axis.
#' @param ltys a numeric vector giving the line types for each row of 'input'.
#' @param add logical: if \code{TRUE} (\code{FALSE} by default), lines are added to existing plot.
#' @param cols a character vector of R build-in colors.
#' @param ... additional parameters to be passed to points and plot.
#' @author Achim Tresch, Bjoern Schwalb
#' @seealso \code{\link{clusterplot}}, \code{\link{demotour}}, \code{\link{disco}}, \code{\link{colorpalette}}
#' @examples len = 20
#' x = sin(seq(0,2*pi,length=len*2))
#' fun = function(){n=sample(1:len,1); return(x[n:(n+len-1)])}
#' input = t(replicate(7,fun(),simplify=TRUE))
#' input = input + rnorm(length(input))/2
#'
#' plotmatrix(input,cols=1:7)
#' @keywords matrix


plotmatrix = function(input,xlim = NULL,ylim = NULL,xlab = "",ylab = "",main = "plotmatrix",type = "l",lwd = 2,at = NULL,xlabels = NULL,ltys = NULL,add = FALSE,cols = NULL,...)
{
	if (!is.matrix(input)) stop("First argument must be a matrix !")
	if (sum(is.na(input)) > 0) print("Be careful: Your data contains NAs ! NA containing rows will be omitted !")
	na.rows = apply(input,1,function(x){!(sum(is.na(x)) > 0)})
	if (!all(na.rows)) stop("Your matrix contains not a single row without NAs !")
	input = input[na.rows,,drop = FALSE]
	if (is.null(at)){at=1:ncol(input)}
	if (is.null(xlabels)){xlabels = at}
	if (is.null(xlim)){xlim = c(1,max(at))}
	if (is.null(ylim)){ylim = range(input)}
	if (is.null(ltys)){ltys = rep(1,nrow(input))}
	if (is.null(cols)) {cols = rep("black", nrow(input))}else {if (length(cols) != nrow(input)) {if(length(cols) < nrow(input)){cols = rep(cols,ceiling(nrow(input)/length(cols)))}
			cols = cols[1:nrow(input)]
			cols[which(is.na(cols))] = "black"
		}
	}
	if(!add){
		plot(c(at[1],median(input)),xlim=xlim,ylim=ylim,type="n",main="",xaxt="n",ylab=ylab,xlab=xlab,...)
		mtext(paste(main),3,2,cex=1.25)
	}
	axis(side=1,at=at,labels=xlabels,...)
	sapply(1:nrow(input),function(x){points(at,input[x,],col=cols[x],lty=ltys[x],type=type,lwd=lwd)})
	points(c(at[1],median(input)),type = "n")
}


### alias ###


LSD.plotmatrix = plotmatrix



