#' @name gregexprind
#' @keywords gregexpr index regular expression
#' @author Sven E. Templer
#' @title Pattern Matching and Extraction
#' @description 
#' Function to extract a certain index from \code{gregexpr()}.
#' @param pattern Character string containing a regular expression to be searched in \code{text}.
#' @param text Character vector where the search is performed.
#' @param n Numeric value or character string \code{"last"} to 
#' extract \code{n}th or last position of \code{pattern} in each value of \code{text}.
#' @param ... Arguments passed to function \code{gregexpr()}.
#' @return
#' Numeric vector of length \code{length(text)}.
#' @seealso
#' See \link{gregexpr} for further information on arguments.\cr
#' See \link{regex} for the use of regular expressions.
#' @examples
#' #
#' 
#' gregexprind(c("a"),c("ababa","ab","xyz",NA), 1)
#' gregexprind(c("a"),c("ababa","ab","xyz",NA), 2)
#' gregexprind(c("a"),c("ababa","ab","xyz",NA), "last")
#' 
#' #

#' @export gregexprind
gregexprind <- function(pattern, text, n, ...) {
  x <- gregexpr(pattern, text, ...)
  if (is.numeric(n)) {
    x <- unlist(lapply(x, function(y) unlist(y)[n]))
  }
  if (n == "last") {
    x <- unlist(lapply(x, function(y) tail(unlist(y),1)))
  }
  x[x < 0] <- NA
  return(x)
}
