#' @title Split a matrix into blocks
#' 
#' @description
#' Split a matrix into a list of blocks (either by rows or by columns)
#'
#' @param Matrix a matrix to split
#' @param blocks either a list or a vector indicating the blocks. 
#' If \code{blocks} is a list of vectors, then the length of each vector 
#' defines the size of the blocks. If \code{blocks} is a vector, then each 
#' element represents the size of the blocks.
#' @param byrow logical. If \code{TRUE} (the default) the matrix is split 
#' by rows, otherwise the matrix is split by columns
#' @return A list of matrices
#' @author Gaston Sanchez
#' @seealso \code{\link{lengths}}, \code{\link{listsize}}
#' @export
#' @examples
#' # matrix with 10 rows and 7 columns
#' M = matrix(rnorm(70), 10, 7)
#' 
#' # row blocks
#' row_sets = list(1:3, 4:5, 6:10)
#' 
#' # split matrix by rows
#' matrix_to_blocks(M, row_sets)
#' 
#' # column blocks
#' col_sets = c(3, 4)
#' 
#' # split matrix by rows
#' matrix_to_blocks(M, col_sets, byrow=FALSE)
matrix_to_blocks <- function(Matrix, blocks, byrow = TRUE)
{
  if (is_not_matrix(Matrix))
    stop("\n'matrix_to_blocks()' requires a matrix")
  if (!is.list(blocks) && !is.numeric(blocks))
    stop("\n'matrix_to_blocks()' requires a list (or vector)")
  
  num_blocks = length(blocks)
  if (is.list(blocks)) {
    size = listsize(blocks)
    start_end = from_to(blocks)
    from = start_end$from
    to = start_end$to    
  }
  if (is.numeric(blocks)) {
    size = sum(blocks)
    to = cumsum(blocks)
    from = to - blocks + 1
  }
  
  # empty list to store results
  blocks_list = vector(mode="list", length=num_blocks)
  
  for (k in 1:num_blocks) {
    if (byrow) {
      if (size != nrow(Matrix))
        stop("\nNumber of rows in 'Matrix' doesn't match 'blocks'")
      blocks_list[[k]] = Matrix[from[k]:to[k],]
    } else {
      if (size != ncol(Matrix))
        stop("\nNumber of columns in 'Matrix' doesn't match 'blocks'")
      blocks_list[[k]] = Matrix[,from[k]:to[k]]      
    }
  }
  # output
  blocks_list
}
