#' Parameters manipulation ?
#' 
#' Parameters manipulation?
#' 
#' 
#' @param parlist List of parameters ?
#' @return Manipulated (shortened ?) parameterlist
#' @note Needs elaboration, \code{selectedpar} warps around this function,
#' where parameter names in \code{nonsetpar} come from is a mystery.
#' @seealso Called by \code{\link{geoplot}}, \code{\link{init}} and
#' \code{\link{selectedpar}}.
#' @keywords device
#' @export Elimcomp
Elimcomp <-
function(parlist){
  txt <- names(parlist)
  txt <- txt[is.na(match(txt,geo::nonsetpar))]
  res <- list()
  for(i in txt) res[[as.character(i)]] <- parlist[[as.character(i)]]
  return(res) 
}

