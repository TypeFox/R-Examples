##' @export
print.biprobit <- function(x,...) {
  printCoefmat(x$coef,...)
  return(invisible(x))
}
