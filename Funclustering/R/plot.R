
#'
#' This function plots the observed curves, before any action 
#' of smoothing or interpolation.
#' 
#' @title plot Original Curves 
#' 
#' 
#' @param time a vector containing the observation time for the curves. If absent the time param will be set
#' at the vector 1:nrow(curves) 
#' 
#' @param curves the observations matrix. Each column of this matrix corresponds to one observed curve, and contains the value of the curve at discrete time points.
#' 
#' @param xlab label of the horizontal axis 
#' 
#' @param ylab label of the vertical axis
#' 
#' @param main the title of the graphic
#' 
#' 			
#' @export 
#' @examples 
#' 
#' data(growth)
#' curves=matrix(data=cbind(growth$hgtm,growth$hgtf),ncol=93)
#' time=growth$age
#' plotOC(time,curves)
#' 
#' 
plotOC <- function(time=1:nrow(curves),curves,xlab="time",ylab="value",main="Original curves"){
	#cheking for the parameter curves
	if(missing(curves)){
		stop("curves is missing.")
	}
	
	#check if the time and curves parameters are on the right size 
	if(length(time) != nrow(curves)){
		stop("length(time) must be equal to nrow(curves)")
	}
	
	dev.new()
	matplot(time,curves,type='l',ylim=c(min(curves)-1,max(curves)+1),xlab=xlab,ylab=ylab,main=main)
}



#' 
#' This function plots a functional data object (after smoothing or interpolation).
#' If you want to color the curves according to a cluster's membership, please specify the parmetere col.  
#' Note: this function works only for univariate functional data.
#'
#'  
#' @title plot a functional data object
#' 
#' @param fd a functional data object
#' 
#' @param col the color vector.
#' 
#' @param xlab label of the horizontal axis 
#' 
#' @param ylab label of the vertical axis
#' 
#' @param main the title of the graphic
#' 
#' 			
#' @export 
#' @examples 
#' 
#' data(growth)
#' data=cbind(matrix(growth$hgtm,31,39),matrix(growth$hgtf,31,54));
#' t=growth$age;
#' splines <- create.bspline.basis(rangeval=c(1, max(t)), nbasis = 20,norder=4);
#' fd <- Data2fd(data, argvals=t, basisobj=splines);
#' cls=c(rep(1,39),rep(2,54)) #there is 39 boys and 54 girls in this data 
#' plotfd(fd,col=cls)
#' 
#' 
plotfd <- function(fd,col=c(1:nrow(fd$coefs)),
		xlab="time",ylab="value",main="Functional data curves") {
	#cheking for fd
	if(missing(fd)){
		stop("You must give a functional data obj to plot it")
	}
	
	#fd must be from calss "fd"
	if(class(fd) != "fd"){
		stop(" fd must be a functional data object")
	}
	
	dev.new()
	plot(fd,col=col,xlab=xlab,ylab=ylab,main=main)
}


