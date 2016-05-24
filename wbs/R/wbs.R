#' @title Wild Binary Segmentation for multiple change-point detection
#' @description The package implements Wild Binary Segmentation, a technique for
#' consistent estimation of  the number and locations of multiple change-points in data.
#' It also provides a fast implementation of the standard Binary Segmentation algorithm.
#' @details The main routines of the package are \code{\link{wbs}}, \code{\link{sbs}} and \code{\link{changepoints}}.	
#' @references P. Fryzlewicz (2014), Wild Binary Segmentation for multiple change-point detection. Annals of Statistics, to appear. (\url{http://stats.lse.ac.uk/fryzlewicz/wbs/wbs.pdf})
#' @docType package
#' @useDynLib wbs
#' @name wbs-package
#' @examples
#' #an example in which standard Binary Segmentation fails to detect change points
#' x <- rnorm(300)+ c(rep(0,130),rep(-1,20),rep(1,20),rep(0,130))
#' 
#' s <- sbs(x)
#' w <- wbs(x)
#' 
#' s.cpt <- changepoints(s)
#' s.cpt
#' 
#' w.cpt <- changepoints(w)
#' w.cpt
#' # in this example, both algorithms work well
#' x <- rnorm(300) + c(rep(1,50),rep(0,250))
#' 
#' s <- sbs(x)
#' w <- wbs(x)
#' 
#' s.cpt <- changepoints(s)
#' s.cpt
#' 
#' w.cpt <- changepoints(w)
#' w.cpt
#' @keywords ts models math


NULL

#' @title Change-point detection via Wild Binary Segmentation
#' @description The function applies the Wild Binary Segmentation algorithm to identify potential locations of the change-points in the mean of the input vector \code{x}. 
#' The object returned by this routine can be further passed to the \code{\link{changepoints}} function, 
#' which finds the final estimate of the change-points based on chosen stopping criteria.
#' @param x a numeric vector
#' @examples
#' x <- rnorm(300) + c(rep(1,50),rep(0,250))
#' w <- wbs(x)
#' plot(w)
#' w.cpt <- changepoints(w)
#' w.cpt
#' th <- c(w.cpt$th,0.7*w.cpt$th) 
#' w.cpt <- changepoints(w,th=th)
#' w.cpt$cpt.th
#' @rdname wbs
#' @export

wbs <- function(x, ...)  UseMethod("wbs")

#' @method wbs default
#' @export 
#' @rdname wbs
#' @param M a number of intervals used in the WBS algorithm
#' @param rand.intervals a logical variable; if \code{rand.intervals=TRUE} intervals used in the procedure are random, thus
#' the output of the algorithm may slightly vary from run to run;  for \code{rand.intervals=FALSE} the intervals used depend on \code{M} and the length of \code{x} only,
#' hence the output is always the same for given input parameters 
#' @param integrated a logical variable indicating the version of Wild Binary Segmentation algorithm used; when \code{integrated=TRUE}, 
#' augmented version of WBS is launched, which combines WBS and BS into one 
#' @param ... not in use
#' @return an object of class "wbs", which contains the following fields
#' \item{x}{the input vector provided}
#' \item{n}{the length of \code{x}}
#' \item{M}{the number of intervals used}
#' \item{rand.intervals}{a logical variable indicating type of intervals}
#' \item{integrated}{a logical variable indicating type of WBS procedure}
#' \item{res}{a 6-column matrix with results, where 's' and 'e' denote start-
#' end points of the intervals in which change-points candidates 'cpt' have been found;
#' column 'CUSUM' contains corresponding value of CUSUM statistic; 'min.th' is the smallest 
#' threshold value for which given change-point candidate would be not added to the set of estimated 
#' change-points; the last column is the scale at which the change-point has been found} 

wbs.default <- function(x, M=5000,  rand.intervals = TRUE,integrated=TRUE,...){
	results <- list()	
	results$x <- as.numeric(x)
	results$n <- length(results$x)
	results$M <- as.integer(M)
	results$integrated <- as.logical(integrated)
	results$rand.intervals <- as.logical(rand.intervals)
	results$res <- matrix(nrow=0,ncol=6)
	
	if(results$n <2) stop("x should contain at least two elements")
	if(NA%in%results$x) stop("x vector cannot contain NA's")
	if(var(x)==0) stop("x is a constant vector, change-point detection is not needed")
	
	if(is.na(results$M)) stop("M cannot be NA")
	if(length(results$M)> 1)  stop("M should be a single integer")
	if(results$M<0)  stop("M should be an integer > 0")
	
	if(results$rand.intervals) intervals <-  matrix(random.intervals(results$n,results$M),ncol=2)
	else {
    intervals <- matrix(fixed.intervals(results$n,results$M),ncol=2)
    results$M <- nrow(intervals)
	}
	
	if(results$integrated){
		results$res <- matrix(.C("wbs_int_rec_wrapper", 
                             x = as.double(results$x), 
                             n = as.integer(results$n), 
                             res = double(6*(results$n-1)),
                             intervals = as.integer(intervals),
                             M = as.integer(results$M))$res,results$n-1,6)
	}else{
	  results$res <- matrix(.C("wbs_rec_wrapper", 
	                           x = as.double(results$x), 
	                           n = as.integer(results$n), 
	                           res = double(6*(results$n-1)),
	                           intervals = as.integer(intervals),
	                           M = as.integer(results$M))$res,results$n-1,6)

    results$res <- matrix(results$res[as.integer(results$res[,1])>0,],ncol=6)
		
	}
	
	
	colnames(results$res) <- c("s","e","cpt","CUSUM","min.th","scale")
	
	class(results) <- "wbs"

	results$cpt <- changepoints(results)
	
	return(results)
	
}


