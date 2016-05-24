#' @title Numeric Levels of an lfactor
#'
#' @description \code{llevels} gives the numeric levels of an lfactor
#'
#' @param x object of class lfactor
#'
#' @return
#' A vector of levels
#'
#' @seealso \code{\link[base]{levels}}
#' 
#' @export
llevels <- function(x) {
  attr(x,"llevels")
}
