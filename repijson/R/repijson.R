#' Classes and functions to handle EpiJSON files
#'
#' EpiJSON is a universal JSON format for epidemiological data. More information on this format can be found at:
#' http://github.com/Hackout2/EpiJSON
#'
#' repijson is a package implementing classes and functions for importing, exporting and converting EpiJSON files.
#'
#' @import OutbreakTools jsonlite sp plyr ggplot2
#' @docType package
#' @name repijson
NULL

#' generic as function
#'
#' @param x an object to convert to repijson
#' @param ... other parameters to pass to the converter
#' @export
as.ejObject <- function(x, ...) UseMethod("as.ejObject")

#' By default we don't know how to convert objects
#'
#' @param x an object
#' @param ... other parameters passed to the call
#' @export
as.ejObject.default <- function(x, ...){
	stop("I don't know how to convert that object")
}
