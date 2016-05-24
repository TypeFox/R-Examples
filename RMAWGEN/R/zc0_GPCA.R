NULL

#' 
#' \code{GPCA} S3 class returned by \code{\link{GPCA}}
#' 
#' 
#' 	\describe{
#' 	\item{list of \code{GPCA_iteration}}{subsequent GPCA iterations}
#' 
#' \item{\code{final_results} }{data.frame or matrix of the "gaussianized" data}
#'  
#' }
#' 
#' 
#' 
#' 
#' @title GPCA-class 
#' 
#' @note Formal definition with \code{\link{setOldClass}} for the S3 class \code{GPCA}
#' 
#' 

#' 
#' @author Emanuele Cordano
#' 
#' @docType class
#' @aliases GPCA-class
#' @name GPCA-class
#' @rdname GPCA-class
#' 
#' @keywords classes
#' @exportClass GPCA
#'
#' @examples showClass("GPCA")
#' 
#'  
#' 
setOldClass("GPCA")