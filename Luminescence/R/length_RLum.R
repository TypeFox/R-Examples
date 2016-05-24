#' General accessor function for RLum S4 class objects
#'
#' Function calls object-specific get functions for RLum S4 class objects.
#'
#' The function provides a generalised access point for specific
#' \code{\linkS4class{RLum}} objects.\cr Depending on the input object, the
#' corresponding get function will be selected. Allowed arguments can be found
#' in the documentations of the corresponding \code{\linkS4class{RLum}} class.
#'
#' @param object \code{\linkS4class{RLum}} (\bold{required}): S4 object of
#' class \code{RLum}
#' @return Return is the same as input objects as provided in the list.
#' @section Function version: 0.1
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France)
#' @seealso
#' \code{\linkS4class{RLum.Data.Curve}},
#' \code{\linkS4class{RLum.Data.Image}},
#' \code{\linkS4class{RLum.Data.Spectrum}},
#' \code{\linkS4class{RLum.Analysis}},
#' \code{\linkS4class{RLum.Results}}
#' @keywords utilities
#' @aliases length_RLum.Analysis
#'
#' @export
setGeneric("length_RLum", function(object) {
  standardGeneric("length_RLum")
})

#' @noRd
#' @export
length_RLum.Analysis <- function(...) {
  .Deprecated("length_RLum")
  length_RLum(...)
}
