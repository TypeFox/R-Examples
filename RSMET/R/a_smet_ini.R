# TODO: Add comment
# 
# Author: ecor
###############################################################################
NULL
#' Class \code{meteoioini}
#' 
#' This class represents the MeteoIO ini format file for weather station data 
#'  \describe{
#'     \item{\code{General}:}{General envirment (list)}
#' 
#'     \item{\code{Input}:}{Input envirment (list)}
#' 
#'     \item{\code{Output}:}{Output envirment (list)}
#' 
#'     \item{\code{Filters}:}{Filters envirment (list)}
#' 
#' 	   \item{\code{Interpolation2D}:}{Interpolation2D envirment (list)}
#' 
#'     \item{\code{file}:}{full name of the SMET file. If it is missing, it is \code{NA}.}
#' 
#'     
#'  }
#' 
#' Detailed information about SMET format is reminded to \url{http://models.slf.ch/docserver/meteoio/SMET_specifications.pdf}.
#' 
#
#' 
#' @references \url{http://models.slf.ch/docserver/meteoio/html/configuration.html} 
#' 
#' @note A \code{meteoioini-class} object can be created by a MeteoIO  file
#' 
#' 
#' @docType class 
#' @title meteoioini-class
#' 
#' @keywords classes
#' 
#' @author Emanuele Cordano
#' @aliases meteoioini-class
#' @name meteoioini-class
#' @rdname meteoioini-class
#' @exportClass meteoioini
#' 
#' @seealso \code{\link{meteoioini}},\code{\link{print.meteoioini}}
#' @examples 
#' 
#' showClass("meteoioini")
#' as.meteoioini("test")
#' 

setClass("meteoioini",representation(General="list",Input="list",Output="list",Filters="list",Interpolation2D="list",file="character"))




