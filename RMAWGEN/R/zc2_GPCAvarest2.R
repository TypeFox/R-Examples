NULL
#' 
#' This class inherits \code{varest2} and contains all information about GPCA (\code{\link{GPCA}} transformation. 
#' 
#'  \describe{
#'     \item{\code{GPCA_data}:}{A \code{"GPCA"} S3 object containing the parameters of the Multi-variate Gaussianization of the time series, it is the result of \code{\link{GPCA}} function applied to the input data of \code{\link{getVARmodel}}
#' 
#'  }
#'     \item{\code{GPCA_residuals}:}{A \code{"GPCA"} S3 object containing the parameters of the Multi-variate Gaussianization of the residuals of the VAR model contained in the \code{VAR} slot; it is \code{NULL} if no Gaussiatization of residuals is applied.
#' Object of class \code{"list"}  }
#'     \item{\code{VAR}:}{S3 Object of class \code{"varest"} }
#'  }
#' 
#' #' @note A \code{GPCAvarest2} object can be created by \code{new("GPCAvarest2", ...)} 
#' or returned by the function \code{\link{getVARmodel}}
#' 
#' 
#' @docType class 
#' @title GPCAvarest2-class 
#' 
#' @keywords classes
#' 
#' @author Emanuele Cordano
#' @aliases GPCAvarest2
#' @name GPCAvarest2-class
#' @rdname GPCAvarest2-class
#' @exportClass GPCAvarest2
#' 
#' 
#' @examples showClass("GPCAvarest2")
#' 

setClass("GPCAvarest2",representation(GPCA_data="GPCA",GPCA_residuals="GPCA"),contains="varest2")

