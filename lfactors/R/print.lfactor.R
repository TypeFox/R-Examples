#' @method print lfactor
#' @export
print.lfactor <- function(x, ...) {
  xi <- x
  attr(xi,"llevels") <- NULL
  print.factor(xi)
  cat("Numeric levels:", llevels(x), "\n")
}
