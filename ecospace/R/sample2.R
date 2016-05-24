#' Internal Ecospace Functions.
#'
#' Internal functions not intended to be called directly by users.
#'
#' Corrects behavior when object sampled has unit length. This is the same
#' function proposed in help file of \code{sample} as \code{resample}.
#'
#' @param x Either a vector of one or more elements from which to choose, or a positive integer.
#' @param ... arguments for particular methods.
#'
#' @seealso \code{\link[base]{sample}}
#'
#' @export
sample2 <- function(x, ...) x[sample.int(length(x), ...)]
