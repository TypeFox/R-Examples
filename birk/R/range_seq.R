# Created by Matthew A. Birk
# Creates a sequence based on the range of the vector
# Last updated: Jan 2015

#' Sequence Generation Spanning A Numerical Range
#'
#' Generates a sequence of numbers spanning the range of \code{x}.
#'
#' @param x a numeric vector.
#' @param extend number specifying the fraction by which the range should be extended.
#' @param \dots further arguments to be passed to \code{\link{seq}}.
#'
#' @author Matthew A. Birk, \email{matthewabirk@@gmail.com}
#' @seealso \code{\link{seq}}, \code{\link{extendrange}}
#'
#' @examples
#' range_seq(rnorm(10, sd = 20))
#' range_seq(c(3, 9), extend = 0.1)
#' range_seq(c(3, 9), length.out = 20)
#'
#' @encoding UTF-8
#' @export
#' @import grDevices

range_seq=function(x, extend=0, ...){
  r=grDevices::extendrange(x, f = extend)
  seq(r[1], r[2], ...)
}