#' @title Provides a summary of measures about the mldr
#' @description Prints a summary of the measures obtained from the \code{mldr} object
#' @param object Object whose measures are to be printed
#' @param ... Additional parameters to be given to print
#' @seealso \code{\link{mldr}}, \code{\link{print.mldr}}
#' @examples
#'
#' library(mldr)
#'
#' summary(emotions)
#'
#' @export
summary.mldr <- function(object, ...) {
  print(data.frame(object$measures), ...)
}

#' @title Prints the mldr content
#' @description Prints the \code{mldr} object data, including input attributs and output labels
#' @param x Object whose data are to be printed
#' @param ... Additional parameters to be given to print
#' @seealso \code{\link{mldr}}, \code{\link{summary.mldr}}
#' @examples
#'
#' library(mldr)
#'
#' emotions
#' print(emotions) # Same as above
#'
#' @export
print.mldr <- function(x, ...) {
  print(x$dataset[ , c(1:x$measures$num.attributes)], ...)
}


