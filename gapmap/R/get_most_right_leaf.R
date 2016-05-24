#'Get the most right leaf object from a dendrogram
#'
#'This function returns the most right leaf object.
#'
#'
#' @param d dendrogram class object
#' @export get_most_right_leaf
#' @aliases get_most_right_leaf
#' @return the most right leaf object
#' @keywords internal
#' 

get_most_right_leaf <- function(d){
  right = d[[2]]
  if(!is.leaf(right)){
    right = get_most_right_leaf(right)
  }
  right
}