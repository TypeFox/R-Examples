#' @name squarematrix
#' @keywords square matrix
#' @author Sven E. Templer
#' @title Create a Square Matrix
#' @description 
#' Transform any m x n matrix to a square matrix by column/row names.
#' Stops if no or duplicated dimnames are provided in x.
#' @param x Object of class \code{matrix}.
#' @return
#' Returns a \code{matrix}.
#' @examples
#' #
#' 
#' m <- matrix(1:6, 2, dimnames=list(2:3,1:3))
#' m
#' squarematrix(m)
#' 
#' #

#' @export squarematrix
squarematrix <- function (x) {
  if (is.null(dimnames(x)))
    stop("Provide dimnames for matrix x.")
  if (any(sapply(dimnames(x), function (i) any(duplicated(i)))))
    stop("Duplicated dimnames found!")
  n <- sort(unique(unlist(dimnames(x))))
  o <- matrix(NA, nrow = length(n), ncol = length(n), dimnames = list(n, n))
  o[rownames(x), colnames(x)] <- x
  return(o)
}
