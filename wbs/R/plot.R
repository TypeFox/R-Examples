#' @title Plot for an 'sbs' object
#' @description Plots the input vector used to generate 'sbs' object \code{x} with fitted piecewise constant function, equal to the mean
#' between change-points specified in \code{cpt}.
#' @details When \code{cpt} is omitted, the function automatically finds change-points 
#' using \code{changepoints} function with a default value of the threshold.
#' @method plot sbs
#' @export 
#' @param x an object of class 'sbs', returned by \code{\link{sbs}}
#' @param cpt a vector of integers with localisations of change-points
#' @param ... other parameters which may be passed to \code{plot} and \code{changepoints}
#' @seealso \code{\link{sbs}}  \code{\link{changepoints}}

plot.sbs <- function(x,cpt,...){
	ts.plot(x$x,ylab="x",...)
	
	if(missing(cpt)){
		w.cpt <- changepoints(x,...)
		print
		means <- means.between.cpt(x$x,w.cpt$cpt.th[[1]])
	}else{
		means <- means.between.cpt(x$x,cpt)
	}
	
	lines(x=means,type="l",col="red")
	title("Fitted piecewise constant function")
	
}

#' @title Plot for a 'wbs' object
#' @description Plots the input vector used to generate 'wbs' object \code{x} with fitted piecewise constant function, equal to the mean
#' between change-points specified in \code{cpt}.
#' @details When \code{cpt} is omitted, the function automatically finds change-points 
#' using \code{changepoints} function with strengthened Schwarz Information Criterion as a stopping criterion for the WBS algorithm.
#' @method plot wbs
#' @export 
#' @param x an object of class 'wbs', returned by \code{\link{wbs}}
#' @param cpt a vector of integers with localisations of change-points
#' @param ... other parameters which may be passed to \code{plot} and \code{changepoints}
#' @seealso \code{\link{wbs}}  \code{\link{changepoints}} \code{\link{ssic.penalty}}


plot.wbs <- function(x,cpt,...){
	
	if(missing(cpt)) plot.sbs(x,cpt=changepoints(x,penalty="ssic.penalty")$cpt.ic[["ssic.penalty"]],...)
	else plot.sbs(x,cpt,...)
	
}

