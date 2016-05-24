#' Print a "waveband" object
#'
#' A function to more nicely print objects of class "waveband".
#'
#' @param x an object of class "waveband"
#' @param ... not used in current version
#'
#' @export
#'
print.waveband <- function(x, ...) {
  cat(x$name, "\n")
  cat("low (nm)", round(x$low, 0), "\n")
  cat("high (nm)", round(x$high, 0), "\n")
  if (!is.null(x$weight)) cat("weighted", x$weight, "\n")
  if (!is.null(x$norm)) cat("normalized at", x$norm, "nm \n")
}
