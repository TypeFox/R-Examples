#' Block Matrices
#'
#' Construct a block matrix using a character string initializer.
#' 
#' @param x A data vector, character string, or a list.
#' @param rows Logical. If TRUE (the default) the matrix is filled by rows, 
#'             otherwise the matrix is filled by columns.
#' @param sep Separator string. Values within each row/column of x are 
#'            separated by this string. Default is \code{","}.
#' @param ... Aditional optional arguments.
#' @return A matrix of class \code{c("matrix", "mat")}.
#' @seealso \code{\link{mat}}, \code{\link{dmat}}.
#' @export
#' @examples
#' # Construct a block matrix from matrices A1, A2, and A3
#' A1 <- mat('1, 1; 1, 1')
#' A2 <- mat('2, 2; 2, 2')
#' A3 <- mat('3, 3, 3, 3')
#' bmat('A1, A2; A3')
bmat <- function(x, rows = TRUE, sep = ",", ...) {
  
  # Split string into pieces (pieces are separated by a semicolon)
  pieces <- unlist(strsplit(x, split = ";"))
  
  # Parse, evaluate, and combine pieces
  bind1 <- if (rows) cbind else rbind
  bind2 <- if (rows) rbind else cbind
  combined <- lapply(pieces, function(x) {
    do.call(bind1, lapply(strsplit(x, split = sep)[[1]], 
                          function(y) eval(parse(text = y))))
  })
  do.call(bind2, combined) 
  
}