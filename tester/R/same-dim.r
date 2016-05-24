#' @title Same Dimension
#' 
#' @description 
#' \code{same_dim()} tests if two matrices have same dimension \cr
#' \code{different_dim()} tests if two matrices have different dimension
#' 
#' @param x a matrix
#' @param y a matrix
#' @aliases same_dim different_dim
#' @export same_dim different_dim
#' @seealso \code{\link{same_nrow}}
#' @examples
#' a = matrix(1:15, 5, 3)
#' 
#' same_dim(a, a) # TRUE
#' same_dim(a, t(a)) # FALSE
#' 
#' different_dim(a, a) # FALSE
#' different_dim(a, t(a)) # TRUE
same_dim <- function(x, y)
{
  if (is_not_tabular(x) || is_not_tabular(y))
    stop("\n'same_dim()' requires matrices or data frames")
  # output
  identical(dim(x), dim(y))
}

different_dim <- function(x, y) {
  !same_dim(x, y)
}



#' @title Same Number of Rows / Columns
#' 
#' @description 
#' \code{same_nrow()} tests if two matrices have same number of rows \cr
#' \code{different_nrow()} tests if two matrices have different 
#' number of rows \cr
#' \code{same_ncol()} tests if two matrices have same number of columns \cr
#' \code{different_ncol()} tests if two matrices have different 
#' number of columns 
#' 
#' @param x a matrix
#' @param y a matrix
#' @aliases same_nrow different_nrow same_ncol different_ncol
#' @export same_nrow different_nrow same_ncol different_ncol
#' @seealso \code{\link{same_dim}}
#' @examples
#' a = matrix(1:15, 5, 3)
#' 
#' same_nrow(a, a) # TRUE
#' same_nrow(a, t(a)) # FALSE
#' same_ncol(a, a) # TRUE
#' same_ncol(a, t(a)) # FALSE
#' 
#' different_nrow(a, a) # FALSE
#' different_nrow(a, t(a)) # TRUE
#' different_ncol(a, a) # FALSE
#' different_ncol(a, t(a)) # TRUE
same_nrow <- function(x, y)
{
  if (is_not_tabular(x) || is_not_tabular(y))
    stop("\n'same_nrow()' requires two matrices (or data frames)")
  # output
  (nrow(x) == nrow(y))
}

different_nrow<- function(x, y) {
  !same_nrow(x, y)
}

same_ncol <- function(x, y)
{
  if (is_not_tabular(x) || is_not_tabular(y))
    stop("\n'same_ncol()' requires two matrices (or data frames)")
  # output
  (ncol(x) == ncol(y))
}

different_ncol<- function(x, y) {
  !same_ncol(x, y)
}
