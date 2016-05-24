# TODO: Add comment
# 
# Author: ecor
###############################################################################
NULL
#' Class \code{smet}
#' 
#' This class represents the SMET format file for weather station data 
#'  \describe{
#'     \item{\code{signature}:}{A charachter string containing SMET signature}
 ###### {A \code{"GPCA"} S3 object containing the parameters of the Multi-variate Gaussianization of the time series, it is the result of \code{\link{GPCA}} function applied to the input data of \code{\link{getVARmodel}}
#' 

#'     \item{\code{header}:}{Object of class \code{"list"} containing the Header Section, each key corresponds to a component of the list.    }
#############{A \code{"GPCA"} S3 object containing the parameters of the Multi-variate Gaussianization of the residuals of the VAR model contained in the \code{VAR} slot; it is \code{NULL} if no Gaussiatization of residuals is applied.

#'     \item{\code{data}:}{S3 Object of class \code{"data.frame"} containing the weather data values. Date field is often alled \code{"timestamp"} and written in \code{\link{POSIXlt}} format. }
#'
#' 	   \item{\code{file}:}{full name of the SMET file. If it is missing, it is \code{NA}.}
#'  }
#' 
#' Detailed information about SMET format is reminded to \url{http://models.slf.ch/docserver/meteoio/SMET_specifications.pdf}.
#' 
#' @seealso \code{\link{smet}},\code{\link{as.smet}},\code{\link{print.smet}}
#' 
#' @references \url{http://models.slf.ch/docserver/meteoio/html/smetio.html} or \url{http://models.slf.ch/docserver/meteoio/SMET_specifications.pdf}
#' 
#' 
#' #' @note A \code{SMET-class} object can be created by a SMET file trough \code{\link{smet}} or can be coerced from another type object through \code{\link{as.smet}}.
#' or returned by the function \code{\link{smet}}
#' 
#' 
#' @docType class 
#' @title smet-class
#' 
#' @keywords classes
#' 
#' @author Emanuele Cordano
#' @aliases smet-class
#' @name smet-class
#' @rdname smet-class
#' @exportClass smet
#' 
#' @seealso \code{\link{smet-class}}
#' @examples 
#' 
#' showClass("smet")
#' 
#' as.smet("test")
#' 

setClass("smet",representation(signature="character",header="list",data="data.frame",file="character"))




