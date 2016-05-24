#' @method as.double lfactor
#' @export
as.double.lfactor <- function(x, ...) {
  as.double(as.character(switchllevels(x)))
}
