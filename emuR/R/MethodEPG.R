##' expand EPG
##' 
##' see function
##' 
##' 
##' @aliases [.EPG
##' @keywords internal
##' @export
"[.EPG" <- function (palates, i, j, k) 
{
  o <- NextMethod("[")
  class(o) <- c("EPG")
  o
}
