.vfrac <- Vectorize(function(x, eps, maxConv) {
  as.character(fractional(x, eps, maxConv))
})

#' Vectorized form for fractional
#'
#' A function which allows any or all of the first three
#' arguments of \code{fractional} to be vectors, with short
#' vectors recycled in the usual way.  Note that the return
#' value is a \emph{character string} vector and may not be
#' used in arithmetic operations
#'
#' @param x as for \code{fractional}
#' @param eps as for \code{fractional}, but may be a vector
#' @param maxConv as for \code{fractional} but may be a vector
#'
#' @return A character string vector of class \code{"charFrac"}
#' @export
#' @examples
#' oldOpt <- options(scipen = 15)
#' pi_approx <- vfractional(base::pi, eps = 0, maxConv = 1:10)
#' within(data.frame(pi_approx, stringsAsFactors = FALSE), {
#'   value = numerical(pi_approx)
#'   error = signif(base::pi - value, 3)
#'   n = seq_along(value) - 1
#' })[, c("n", "pi_approx", "value", "error")]
#' options(oldOpt)
vfractional <- function(x, eps = 1.0e-6, maxConv = 20) {
  values <- .vfrac(x, eps, maxConv)
  class(values) <- c("charFrac", "character")
  values
}

#' Extract the parts of a fractional object
#'
#' Generic function for extracting numerators with methods for \code{"fractional"} or
#' \code{"charFrac"} objects
#'
#' @param x An object of class \code{"fractional"} or
#' \code{"charFrac"}
#'
#' @return An \code{integer} vector of numerators
#' @export
#'
#' @examples
#' (pi_approx <- vfractional(base::pi, eps = 0, maxConv = 1:10))
#' numerators(pi_approx)
#' denominators(pi_approx)
numerators <- function(x) {
  UseMethod("numerators")
}

#' @describeIn numerators \code{numerators} method function for \code{"charFrac"} objects
#' @export
numerators.charFrac <- function(x) {
  x <- unclass(x)
  ax <- attributes(x)
  ax$eps <- ax$maxConv <- ax$sync <- NULL
  value <- sapply(strsplit(x, "/"), function(f) {
    as.integer(f[1])
  })
  do.call(structure, c(list(.Data = value), ax))
}

#' @describeIn numerators \code{numerators} method function for \code{"fractional"} objects
#' @export
numerators.fractional <- function(x) {
  numerators.charFrac(as.character(x))
}

#' @describeIn numerators Default \code{numerators} method for \code{numeric} objects
#' @export
numerators.default <- function(x) {
  numerators.charFrac(as.character(fractional(x)))
}

#' Extract denominators from a fractional object
#'
#' Generic function for extracting denominators with methods for \code{"fractional"} or
#' \code{"charFrac"} objects
#' @rdname numerators
#'
#' @return An \code{integer} vector of denominators
#' @export
denominators <- function(x) {
  UseMethod("denominators")
}

#' @describeIn numerators \code{denominators} method function for \code{"charFrac"} objects
#' @export
denominators.charFrac <- function(x) {
  x <- unclass(x)
  ax <- attributes(x)
  ax$eps <- ax$maxConv <- ax$sync <- NULL
  value <- sapply(strsplit(x, "/"), function(f) {
    if(length(f) == 1) 1L else as.integer(f[2])
  })
  do.call(structure, c(list(.Data = value), ax))
}

#' @describeIn numerators \code{denominators} method function for \code{"fractional"} objects
#' @export
denominators.fractional <- function(x) {
  denominators.charFrac(as.character(x))
}

#' @describeIn numerators Default \code{denominators} method for \code{numeric} objects
#' @export
denominators.default <- function(x) {
  denominators.charFrac(as.character(fractional(x)))
}
