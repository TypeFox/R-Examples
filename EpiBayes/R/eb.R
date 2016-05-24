#' @title 
#' eb Generic
#' 
#' @description 
#' This function is used as a generic for the \code{eb} object class. This class is used
#'     to store output from the disease models \code{\link{EpiBayes_s}} and 
#'     \code{\link{EpiBayes_ns}}.
#'
#' @param x An object.
#'
#' @details
#' The \code{\link{EpiBayes_s}} and \code{\link{EpiBayes_ns}} functions return objects of the
#'     class \code{eb}, for which this is the generic. This class of objects is used to 
#'     store the output of the EpiBayes models and to reference the \code{print}, 
#'     \code{summary}, and \code{plot} methods for these output objects. 
#'
#' @seealso
#' Other classes' generics used within this package include \code{\link{ebhistorical}} for 
#'     output objects of the \code{\link{EpiBayesHistorical}} function and 
#'     \code{\link{ebsamplesize}} for output objects of the 
#'     \code{\link{EpiBayesSampleSize}} function.  
#'
eb = function(x){
	UseMethod("eb")
}