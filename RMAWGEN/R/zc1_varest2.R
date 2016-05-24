NULL


#'   
#' This class derives from a \code{varest} S3 class which is a list of objects describing a Vectorial AutoRegressive Model (see \code{\link{VAR}})
#' 

#' 	\describe{
#'     \item{\code{VAR}:}{a \code{varest} S3 object created by \code{\link{VAR}} }
#'  }

#' 


#' 
#' @title varest2-class 
#' 
#' @note A \code{varest2} object can be created by \code{new("varest2", ...)} 
#' or returned by the function \code{\link{getVARmodel}}
#' 
#' 

#' 
#' @author Emanuele Cordano
#' 
#' @docType class
#' @aliases varest2
#' @name varest2-class
#' @rdname varest2-class
#' 
#' @keywords classes
#' @exportClass varest2
#'
#' @examples showClass("varest2")
#' 
#'  
#' 


# @slot VAR a \code{list} object derived by a \code{varest} list object created by \code{\link{VAR}} 



setClass("varest2",representation(VAR="varest"))




