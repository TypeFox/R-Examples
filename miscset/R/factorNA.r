#' @name factorNA
#' @author Sven E. Templer
#' @title Create a Factor with NA as Level
#' @description 
#' Create a \link{factor} with \code{NA} values included
#' and positioned as last level.
#' @param x A vector coerced to character.
#' @param ... Forwarded to \link{factor}. \code{x} and \code{levels} are defined.

#' @rdname factorNA
#' @export
factorNA <- function(x, ...) {
  
  x <- as.character(x)
  x[is.na(x)] <- "NA"
  lev <- c(setdiff(sort(unique(x)), "NA"), "NA")
  factor(x, levels = lev, ...)
  
}
