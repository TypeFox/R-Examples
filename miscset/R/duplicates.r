#' @name duplicates
#' @aliases duplicatei
#' @author Sven E. Templer
#' @title Determine Duplicates
#' @description 
#' Determine duplicates. \code{duplicates} returns a logical vector,
#' \code{duplicatei} an integer vector.
#' @param x A vector or data.frame to search for duplicates.
#' @param first Logical, \code{TRUE} to return the index also for the first
#' occurrence of values. Otherwise, a \code{0} is the index for the first
#' occurrence.
#' @return
#' \code{duplicates} returns a logical vector as \link{duplicated}, but with
#' \code{TRUE} values also for the first occurrence of duplicated values.\cr
#' \code{duplicatei} returns the index of the first occurrence of each unique
#' value.
#' @examples
#' #
#' 
#' x <- c(7, 7, 7, 2, 3, 2)
#' data.frame(
#'   data = x,
#'   duplicated = duplicated(x),
#'   duplicates = duplicates(x),
#'   duplicatei = duplicatei(x),
#'   duplicatei0 = duplicatei(x, FALSE))
#' 
#' #

#' @rdname duplicates
#' @export
duplicates <- function(x) {
  
  duplicated(x) | duplicated(x, fromLast=TRUE)
  
}

#' @rdname duplicates
#' @export
duplicatei <- function (x, first = TRUE) {
  
  if (first)
    match(x,x)
  else
    ifelse(duplicated(x), match(x,x), 0)
  
}
