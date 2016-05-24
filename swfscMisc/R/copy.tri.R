#' @title Copy Matrix Triangles
#' @description Copy between lower left and upper right triangles of a matrix.
#' 
#' @param x a matrix.
#' @param from triangle to copy from. Can be "lower" or "upper".
#' 
#' @return a matrix.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' x <- matrix(1:9, nrow = 3)
#' print(x)
#' copy.tri(x)
#' 
#' @export
#' 
copy.tri <- function(x, from = "lower") {  
  if (!is.matrix(x)) stop("'x' must be a matrix")
  if (nrow(x) != ncol(x)) stop("'x' must be a square matrix")
  new.mat <- x
  from <- tolower(from)
  if(tolower(from) == "lower") {
    for (row in 1:(nrow(x) - 1)) {
      for (col in (row + 1):nrow(x)) new.mat[row, col] <- x[col, row]
    }
  } else if(tolower(from) == "upper") {
    for (row in 1:(nrow(x) - 1)) {
      for (col in (row + 1):nrow(x)) new.mat[col, row] <- x[row, col]
    }
  }
  new.mat
}
