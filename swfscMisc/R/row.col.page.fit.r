#' @title Number of Rows and Columns on Page
#' @description Return the number of rows and columns for \code{n}
#' that best fits on a  page of size \code{width} x \code{height}.
#' 
#' @param n number of items (e.g., plots) to fit on page.
#' @param width,height dimensions of page.
#' 
#' @return A vector listing the number of rows and columns to use.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' # 9 frames on US letter paper
#' row.col.page.fit(9)
#' 
#' # 9 frames on a square
#' row.col.page.fit(9, width = 10, height = 10)
#' 
#' @export
#' 
row.col.page.fit <- function(n, width = 8.5, height = 11) {
  n <- as.integer(n)
  if(any(n <= 0) | any(is.na(n))) stop("'n' must be a positive integer")
  nr <- 1:n
  nc <- n / nr
  ratio.diff <- (nc / nr) - (width / height)
  best <- which.min(abs(ratio.diff))
  c(num.rows = best, num.cols = ceiling(nc[best]))
}
