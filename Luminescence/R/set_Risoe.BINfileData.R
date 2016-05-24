#' General accessor function for RLum S4 class objects
#'
#' Function calls object-specific get functions for RisoeBINfileData S4 class objects.
#'
#' The function provides a generalised access point for specific
#' \code{\linkS4class{Risoe.BINfileData}} objects.\cr Depending on the input object, the
#' corresponding get function will be selected. Allowed arguments can be found
#' in the documentations of the corresponding \code{\linkS4class{Risoe.BINfileData}} class.
#'
#' @param METADATA x
#' @param DATA x
#' @param .RESERVED x
#'
#' @return Return is the same as input objects as provided in the list.
#' @section Function version: 0.1
#' @author Sebastian Kreutzer, IRAMAT-CRP2A, Universite Bordeaux Montaigne
#' (France)
#' @seealso
#' \code{\linkS4class{Risoe.BINfileData}}
#' @keywords utilities
#' 
#' @export
setGeneric("set_Risoe.BINfileData", function(METADATA, DATA, .RESERVED) {
  standardGeneric("set_Risoe.BINfileData")
})
