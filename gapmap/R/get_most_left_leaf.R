#'Get the most left leaf object from a dendrogram
#'
#'This function returns the most left leaf object.
#'
#'
#' @param d dendrogram class object
#' @export get_most_left_leaf
#' @aliases get_most_left_leaf
#' @return the most left leaf object
#' @keywords internal
#' 

get_most_left_leaf <- function(d){
  left = d[[1]]
  if(!is.leaf(left)){
    left = get_most_left_leaf(left)
  }
  left
}