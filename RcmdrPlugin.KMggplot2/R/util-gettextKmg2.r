#' A \code{gettextRcmdr} Wrapper Function
#'
#' This function is a \code{gettextRcmdr} wrapper function for this package.
#'
#' @param ... arguments passed to gettext function
#' @seealso \code{\link[Rcmdr:Rcmdr.Utilities]{Rcmdr.Utilities}}
#'
#' @rdname util-gettextKmg2
#' @keywords documentation
#' @export
gettextKmg2 <- function(...) {

  gettext(..., domain = "R-RcmdrPlugin.KMggplot2")

}
