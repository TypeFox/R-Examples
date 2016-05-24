#'Count the number of gaps based on a threshold
#'
#'This function counts the number of gaps based on the cutting the tree method.
#'It counts the number of branches that are above the threshold distance.
#'
#'
#' @param d dendrogram class object
#' @param count count
#' @param threshold a numeric value for threshold
#' @export count_gap
#' @aliases count_gap
#' @return the count of gaps
#' @keywords internal
#' 

count_gap <- function(d=d, count=0, threshold=threshold){
  a = attributes(d) #attributes
  height = a$height
  left = d[[1]]
  right = d[[2]]
  gapped = ifelse((height > threshold), TRUE, FALSE) 
  if(gapped) count = count+1 #increment
  if(is.leaf(left) && is.leaf(right)){
  }else if(!is.leaf(left) && is.leaf(right)){
   if(gapped){
     count = count_gap(left, count=count, threshold=threshold)
   }
  }else if(is.leaf(left) && !is.leaf(right)){
    if(gapped){
      count = count_gap(right, count=count, threshold=threshold)
    }
  }else{
    if(gapped){
      count = count_gap(left, count=count, threshold=threshold)
      count = count_gap(right, count=count, threshold=threshold)
    }
  }
  return(count)
}