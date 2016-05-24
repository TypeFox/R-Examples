#'Sum the distance of all branches in a dendrogram
#'
#'This function takes a dendrogram class object as an input, and adds up all the distances of branches. 
#'This function is called recursively to adds up the sum. In case of exponential mapping for the quantitative
#'mode, the sum is in the exponential scale
#'
#'
#' @param d dendrogram class object
#' @param sum the sum of distance
#' @param mapping in case of quantitative mode, either "linear" or "exponential" mapping
#' @param ... ignored
#' @export sum_distance
#' @aliases sum_distance
#' @return the sum of distances
#' @keywords internal
#' 

sum_distance <- function(d, sum=0, mapping=c("exponential", "linear"), scale = 0, max_height=0, ...){
  a = attributes(d) #attributes
  height = a$height
  mapping <- match.arg(mapping)
  if(mapping == "exponential"){
    height = map.exp(height, 0, max_height, 0, 1, scale = scale)
  }
  left = d[[1]]
  right = d[[2]]
  if(is.leaf(left) && is.leaf(right)){
    sum = sum + height
  }else if(!is.leaf(left) && is.leaf(right)){
    sum = sum + height
    sum = sum_distance(left, sum, mapping=mapping, scale=scale, max_height=max_height)
  }else if(is.leaf(left) && !is.leaf(right)){
    sum = sum + height
    sum = sum_distance(right, sum, mapping=mapping, scale=scale, max_height=max_height)
  }else{
    sum = sum + height
    sum = sum_distance(left, sum, mapping=mapping, scale=scale, max_height=max_height)
    sum = sum_distance(right, sum, mapping=mapping, scale=scale, max_height=max_height)
  }
  return(sum)
}
