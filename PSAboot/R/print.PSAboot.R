#' Print results of PSAboot
#' 
#' @param x result of \code{\link{PSAboot}}.
#' @param ... currently unused.
#' @method print PSAboot
#' @export
print.PSAboot <- function(x, ...) {
	print(summary(x, ...))
}
