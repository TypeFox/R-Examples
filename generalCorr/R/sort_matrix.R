#' Sort all columns of matrix x by j-th column while carrying along all columns.
#'
#' This function simply uses sort.list function in R.
#'
#' @param  x {is a matrix with several columns}
#' @param j {column number with which to sort}
#' @return {A sorted matrix}
#' @examples
#'
#' set.seed(30)
#' x=matrix(sample(1:50),ncol=5)
#' y=sort_matrix(x,3);y
#' @export

sort_matrix <-
function(x,j)
{
y=x[sort.list(x[,j]),]
return(y)  }
