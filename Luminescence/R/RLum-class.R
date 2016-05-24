#' @include replicate_RLum.R internals_RLum.R
NULL

#' Class \code{"RLum"}
#'
#' Abstract class for data in the package Luminescence
#'
#'
#' @name RLum-class
#'
#' @docType class
#'
#' @slot originator Object of class \code{\link{character}} containing the name of the producing
#' function for the object. Set automatically by using the function \code{\link{set_RLum}}.
#'
#' @slot .uid Object of class \code{\link{character}} containing a unique object identifier (md5 hash)
#' calculated using the internal function \code{.create_UID()}. This id is calculated everytime and
#' object is created.
#'
#' @note \code{RLum} is a virtual class.
#'
#' @section Objects from the Class: A virtual Class: No objects can be created
#' from it.
#'
#' @section Class version: 0.3.0
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne (France)
#'
#' @seealso \code{\linkS4class{RLum.Data}}, \code{\linkS4class{RLum.Analysis}}
#'
#' @keywords classes
#'
#' @examples
#'
#' showClass("RLum")
#'
#' @export
setClass("RLum",
           slots = list(
             originator = "character",
             .uid = "character"),
           contains = "VIRTUAL",
           prototype = prototype(
             originator = NA_character_,
             .uid = .create_UID()
           )
         )


# replication method for object class ------------------------------------------

#' @describeIn RLum
#' Replication method RLum-objects
#'
#' @param object an object of class \code{\linkS4class{RLum}} (\bold{required})
#'
#' @param times \code{\link{integer}} (optional): number for times each element is repeated
#' element
#'
#' @export
setMethod(
  "replicate_RLum",
  "RLum",
  definition = function(object, times = NULL) {

    ##The case this is NULL
    if (is.null(times)) {
      times <- 1
    }

    lapply(1:times, function(x) {
      object

    })

  }
)
