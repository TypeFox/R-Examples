#' @title Change-point detection via standard Binary Segmentation
#' @description The function applies the Binary Segmentation algorithm to identify potential locations of the change-points in the mean of the input vector \code{x}. 
#' The object returned by this routine can be further passed to the \code{\link{changepoints}} function, 
#' which finds the final estimate of the change-points based on thresholding. 
#' @param x a numeric vector
#' @param ... not in use
#' @examples
#' x <- rnorm(300) + c(rep(1,50),rep(0,250))
#' s <- sbs(x)
#' s.cpt <- changepoints(s)
#' s.cpt
#' th <- c(s.cpt$th,0.7*s.cpt$th) 
#' s.cpt <- changepoints(s,th=th)
#' s.cpt
#' @rdname sbs
#' @export
#' @return an object of class "sbs", which contains the following fields
#' \item{x}{the vector provided}
#' \item{n}{the length of \code{x}}
#' \item{res}{a 6-column matrix with results, where 's' and 'e' denote start-
#' end points of the intervals in which change-points candidates 'cpt' have been found;
#' column 'CUSUM' contains corresponding value of CUSUM statistic; 'min.th' is the smallest 
#' threshold value for which given change-point candidate would be not added to the set of estimated 
#' change-points; the last column is the scale at which the change-point has been found} 

sbs <- function(x, ...)  UseMethod("sbs")

#' @method sbs default
#' @export
#' @rdname sbs

sbs.default <- function(x, ...){
	results <- list()	
	results$x <- as.numeric(x)
	results$n <- length(results$x)

	
	if(results$n <2) stop("x should contain at least two elements")
	if(NA%in%results$x) stop("x vector cannot contain NA's")
	if(var(x)==0) stop("x is a constant vector, change-point detection is not needed")
  
	results$res <- matrix(.C("bs_rec_wrapper",
                           x = as.double(results$x),
                           n = as.integer(results$n),
                           res = double(6*(results$n-1)))$res,results$n-1,6)
	
	colnames(results$res) <- c("s","e","cpt","CUSUM","min.th","scale")
	
	class(results) <- "sbs"
	return(results)
	
}
