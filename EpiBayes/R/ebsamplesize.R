#' @title 
#' ebsamplesize Generic
#' 
#' @description 
#' This function is used as a generic for the \code{ebsamplesize} object class. This class 
#'     is used to store output from the sample size search function 
#'     \code{\link{EpiBayesSampleSize}}.
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
#' Other classes' generics used within this package include \code{\link{eb}} for output
#'     objects of \code{\link{EpiBayes_ns}} and \code{\link{EpiBayes_s}}, 
#'     \code{\link{ebhistorical}} for output objects of the \code{\link{EpiBayesHistorical}} 
#'     function and \code{\link{ebsamplesize}} for output objects of the 
#'     \code{\link{EpiBayesSampleSize}} function.  
#'
ebsamplesize = function(x){
	UseMethod("ebsamplesize")
}