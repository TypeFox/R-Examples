#'Set the most right leaf object from a dendrogram
#'
#'This function replace the most right leaf with provided dendrogram
#'
#'
#' @param d dendrogram class object, subtree object
#' @param d2 dendrogram class object, a leaf to replace the most right
#' @export set_most_right_leaf
#' @aliases set_most_right_leaf
#' @return the dendrogram class object where the most right leaf is replaced
#' @keywords internal
#' 

set_most_right_leaf <- function(d, d2, ...){
  right = d[[2]]
  if(is.leaf(right)){
    #set the value
    d[[2]] = d2
  }else{
    #go deeper
    d[[2]] = set_most_right_leaf(right, d2)
  }
  d #return
}