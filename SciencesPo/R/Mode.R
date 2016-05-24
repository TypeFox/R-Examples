#' @encoding UTF-8
#' @title  Calculate the Mode
#'
#' @description Estimates the mode for a vector
#'
#' @param x An \R object.
#' @param na.rm A logical value indicating whether \code{NA}
#'  should be stripped before the computation proceeds. Default
#'  is \code{FALSE}
#' @note This function replaces the \code{base} function of the same name, while \code{SciencesPo::mode} calculates the \dQuote{mode}, \code{base::mode} prints the \dQuote{class} of an object.
#'
#' @author Daniel Marcelino, \email{dmarcelino@@live.com}
#'
#' @examples
#' myvar <-c(1,1,2,2,3,3,4,4,5, NA)
#' Mode(myvar)
#'
#' Mode(myvar, FALSE)
#' @export
#' @rdname Mode
`Mode` <- function(x, na.rm = FALSE) UseMethod("Mode")

#' @export
#' @rdname Mode
`Mode` <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = base::subset(x, !is.na(x))
  }
  y <- as.factor(x)
  freqs <- base::summary(y)
  mode <- names(freqs)[freqs[names(freqs)] == max(freqs)]
  return(as.numeric(mode) )
}
NULL

#' @export
#' @rdname Mode
`Mode.data.frame` <- function(x, na.rm = TRUE) sapply(x, Mode)
