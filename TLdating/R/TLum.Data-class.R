#' Class \code{"TLum.Data"}
#'
#' Generalized virtual data class for luminescence data.
#'
#'
#' @name TLum.Data-class
#' @rdname TLum.Data-class
#'
#' @note The code and the structure of this class is based on the \linkS4class{RLum.Data} class from the \link{Luminescence} package.
#'
#' @keywords classes
#'
#' @author David Strebler
#'
#' @exportClass TLum.Data

setClass("TLum.Data",
         contains = "TLum")
