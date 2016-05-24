#' General set function for RLum S4 class objects
#'
#' Function calls object-specific set functions for RLum S4 class objects.
#'
#' The function provides a generalised access point for specific
#' \code{\linkS4class{RLum}} objects.\cr Depending on the given class, the
#' corresponding method to create an object from this class will be selected.
#' Allowed additional arguments can be found in the documentations of the
#' corresponding \code{\linkS4class{RLum}} class: \code{\linkS4class{RLum.Data.Curve}},
#' \code{\linkS4class{RLum.Data.Image}}, \code{\linkS4class{RLum.Data.Spectrum}},
#' \code{\linkS4class{RLum.Analysis}} and \code{\linkS4class{RLum.Results}}
#'
#' @param class \code{\linkS4class{RLum}} (\bold{required}): name of the S4 class to
#' create
#'
#' @param originator \code{\link{character}} (automatic): contains the name of the calling function
#' (the function that produces this object); can be set manually.
#'
#' @param \dots further arguments that one might want to pass to the specific
#' set method
#'
#' @return Returns an object of the specified class.
#'
#' @section Function version: 0.2.0
#'
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France)
#'
#' @seealso
#' \code{\linkS4class{RLum.Data.Curve}},
#' \code{\linkS4class{RLum.Data.Image}},
#' \code{\linkS4class{RLum.Data.Spectrum}},
#' \code{\linkS4class{RLum.Analysis}},
#' \code{\linkS4class{RLum.Results}}
#'
#' @keywords utilities
#'
#' @aliases set_RLum.Data.Curve set_RLum.Data.Image set_RLum.Data.Spectrum
#' set_RLum.Analysis set_RLum.Results
#'
#' @examples
#'
#' ##produce empty objects from each class
#' set_RLum(class = "RLum.Data.Curve")
#' set_RLum(class = "RLum.Data.Spectrum")
#' set_RLum(class = "RLum.Data.Spectrum")
#' set_RLum(class = "RLum.Analysis")
#' set_RLum(class = "RLum.Results")
#'
#' ##produce a curve object with arbitrary curve values
#' object <- set_RLum(
#' class = "RLum.Data.Curve",
#' curveType = "arbitrary",
#' recordType = "OSL",
#' data = matrix(c(1:100,exp(-c(1:100))),ncol = 2))
#'
#' ##plot this curve object
#' plot_RLum(object)
#'
#' @export
setGeneric("set_RLum", function (class, originator, ... ) {
  class(class) <- as.character(class)

  if(missing(originator)) {
    if (is(sys.call(which = -1)[[1]], "name")) {
      originator <- as.character(sys.call(which = -1)[[1]])
    } else{
      originator <- NA_character_
    }
  }
  standardGeneric("set_RLum")
})


## ---- DEPRECATED GENERICS
# .Deprecated in package version 0.5.0
# .Defunct in 0.5.1
# Removed in 0.6.0

#' @noRd
#' @export
set_RLum.Analysis <- function(...) {
  .Defunct("set_RLum")

  if(missing(originator)) {
    if (is(sys.call(which = -1)[[1]], "name")) {
      originator <- as.character(sys.call(which = -1)[[1]])
    } else{
      originator <- NA_character_
    }
  }

  set_RLum("RLum.Analysis", ...)
}

#' @noRd
#' @export
set_RLum.Data.Curve <- function(...) {
  .Defunct("set_RLum")

  if(missing(originator)) {
    if (is(sys.call(which = -1)[[1]], "name")) {
      originator <- as.character(sys.call(which = -1)[[1]])
    } else{
      originator <- NA_character_
    }
  }

  set_RLum("RLum.Data.Curve", ...)
}

#' @noRd
#' @export
set_RLum.Data.Image <- function(...) {
  .Defunct("set_RLum")

  if(missing(originator)) {
    if (is(sys.call(which = -1)[[1]], "name")) {
      originator <- as.character(sys.call(which = -1)[[1]])
    } else{
      originator <- NA_character_
    }
  }

  set_RLum("RLum.Data.Image", ...)
}

#' @noRd
#' @export
set_RLum.Data.Spectrum <- function(...) {
  .Defunct("set_RLum")

  if(missing(originator)) {
    if (is(sys.call(which = -1)[[1]], "name")) {
      originator <- as.character(sys.call(which = -1)[[1]])
    } else{
      originator <- NA_character_
    }
  }

  set_RLum("RLum.Data.Spectrum", ...)
}

#' @noRd
#' @export
set_RLum.Results <- function(originator,...) {
  .Defunct("set_RLum")

  if(missing(originator)) {
    if (is(sys.call(which = -1)[[1]], "name")) {
      originator <- as.character(sys.call(which = -1)[[1]])
    } else{
      originator <- NA_character_
    }
  }
  set_RLum(class = "RLum.Results", originator = originator, ...)
}
