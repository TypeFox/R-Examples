#' Data Frames
#' 
#' Like \code{mat}, but returns a data frame.
#' 
#' @param x A data vector, character string, or a list.
#' @param ... Aditional optional arguments passed on to \code{mat}.
#' @return A data frame.
#' @seealso \code{\link{mat}}, \code{\link{bmat}}.
#' @export
#' @examples
#' dmat('1e-01, 2+5, 3, 4, 5; 6, 7, 8, 9^2, pi', rows = FALSE)
#' z <- list(a = 1:10, b = 11:20, c = 21:30)
#' dmat(z)  # list elements form rows
#' dmat(z, rows= FALSE)  # list elements form columns
dmat <- function(x, ...) {
  as.data.frame(mat(x, ...))  # FIXME: data.frame or as.data.frame?
}