# Created by Matthew A. Birk
# Calculates the geometric mean of a vector
# Last updated: Nov 2014

#' Geometric Mean
#'
#' Computes the geometric mean of a vector, \code{x}. It is a wrapper for \code{exp(mean(log(x)))}.
#'
#' @param x a numeric vector or an R object which is coercible to one by as.vector(x, "numeric”).
#' @param add0.001 logical. Should a small constant (0.001) be added to avoid issues with zeroes?
#' @param ignore_neg logical. Should negative values be ignored to avoid NaNs?
#' @param \dots further arguments passed to \code{\link{mean}}.
#'
#' @author Matthew A. Birk, \email{matthewabirk@@gmail.com}
#' @seealso \code{\link{mean}}
#'
#' @examples
#' geom_mean(1:10)
#' geom_mean(0:10)
#' geom_mean(0:10, add0.001 = TRUE)
#' geom_mean(-10:10, add0.001 = TRUE, ignore_neg = TRUE)
#'
#' @encoding UTF-8
#' @export

geom_mean=function(x,add0.001=FALSE,ignore_neg=FALSE,...){
  if(ignore_neg==TRUE) x=x[x>=0]
  if(add0.001==TRUE) exp(mean(log(x+0.001),...)) else exp(mean(log(x),...))
}