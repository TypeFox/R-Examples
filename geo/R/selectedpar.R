#' Parameters manipulation ?
#' 
#' Parameters manipulation ?
#' 
#' 
#' @return Returns the names of some selected parameters.
#' @note Needs elaboration,
#' @seealso Called by the majority of the graphical functions in geo, calls
#' \code{\link{Elimcomp}}.
#' @keywords device
#' @export selectedpar
selectedpar <-
function() 
  return(Elimcomp(par(no.readonly=T)))

