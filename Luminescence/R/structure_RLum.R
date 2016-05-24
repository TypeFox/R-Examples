#' General structure function for RLum S4 class objects
#'
#' Function calls object-specific get functions for RLum S4 class objects.
#'
#' The function provides a generalised access point for specific
#' \code{\linkS4class{RLum}} objects.\cr Depending on the input object, the
#' corresponding structure function will be selected. Allowed arguments can be found
#' in the documentations of the corresponding \code{\linkS4class{RLum}} class.
#'
#' @param object \code{\linkS4class{RLum}} (\bold{required}): S4 object of
#' class \code{RLum}
#' @return Returns a \code{data.frame} with structure of the object.
#' @section Function version: 0.1.0
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France)
#' @seealso
#' \code{\linkS4class{RLum.Data.Curve}},
#' \code{\linkS4class{RLum.Data.Image}},
#' \code{\linkS4class{RLum.Data.Spectrum}},
#' \code{\linkS4class{RLum.Analysis}},
#' \code{\linkS4class{RLum.Results}}
#' @keywords utilities
#' @aliases get_structure.RLum.Analysis
#'
#' @export
setGeneric("structure_RLum", function(object) {
  standardGeneric("structure_RLum")
})

## ---- DEPRECATED GENERICS
# .Deprecated in package version 0.5.0
# .Defunct in 0.5.X
# Removed in 0.6.0

#' @noRd
get_structure.RLum.Analysis <- function(...) {
  .Deprecated("get_structure.RLum.Analysis")
  structure_RLum(...)
}
