#' @title Split a data frame into blocks
#' 
#' @description
#' Split a data frame into a list of blocks (either by rows or by columns)
#'
#' @param DataFrame a data frame to split
#' @param blocks either a list or a vector indicating the blocks. If 
#' \code{blocks} is a list of vectors, then the length of each vector defines 
#' the size of the blocks. If \code{blocks} is a vector, then each element 
#' represents the size of the blocks.
#' @param byrow logical. If \code{TRUE} (the default) the data frame 
#' is split by rows, otherwise the data frame is split by columns
#' @return A list of data frames
#' @author Gaston Sanchez
#' @seealso \code{\link{matrix_to_blocks}}
#' @export
#' @examples
#' # say you have a data frame
#' iris_df = iris[c(1:3,51:53,101:103),]
#' 
#' # list defining the blocks
#' row_blocks = list(1:3, 4:6, 7:9)
#' col_blocks = c(2, 2, 1)
#'  
#' # split data into list of blocks (by rows)
#' df_to_blocks(iris_df, row_blocks)
#' 
#' # split data into list of blocks (by columns)
#' df_to_blocks(iris_df, col_blocks, byrow=FALSE)
df_to_blocks <- function(DataFrame, blocks, byrow = TRUE)
{
  if (is_not_dataframe(DataFrame))
    stop("\n'df_to_blocks()' requires a data frame")
  if (!is.list(blocks) && !is.numeric(blocks))
    stop("\n'df_to_blocks()' requires a list (or vector)")
  
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
      if (size != nrow(DataFrame))
        stop("\nNumber of rows in 'DataFrame' doesn't match 'blocks'")
      blocks_list[[k]] = DataFrame[from[k]:to[k],]
    } else {
      if (size != ncol(DataFrame))
        stop("\nNumber of columns in 'DataFrame' doesn't match 'blocks'")
      blocks_list[[k]] = DataFrame[,from[k]:to[k]]      
    }
  }
  # output
  blocks_list
}
