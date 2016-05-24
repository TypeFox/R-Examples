################################################################################
# File:             is.in.convex.r
# Created by:       Pavlo Mozharovskyi
# First published:  28.02.2013
# Last revised:     28.02.2013
# 
# Check if points lie in the convex hulls of the data clouds.
################################################################################

is.in.convex <- function(x, data, cardinalities, seed = 0){
  if (!is.numeric(data)
      || !is.matrix(data) 
      || ncol(data) < 2){
    stop("Argument \"data\" should be a numeric matrix of at least 2-dimensional data")
  }
  if (!is.vector(cardinalities, mode = "numeric") 
      || is.na(min(cardinalities)) 
      || sum(.is.wholenumber(cardinalities)) != length(cardinalities) 
      || min(cardinalities) <= 0 
      || sum(cardinalities) != nrow(data)){
    stop("Argument \"cardinalities\" should be a vector of cardinalities of the classes in \"data\" ")
  }
  if (sum(cardinalities < ncol(data) + 1) != 0){
    stop("Not in all classes sufficiently enough objetcs")
  }
  if (!is.matrix(x) 
       && is.vector(x)){
    x <- matrix(x, nrow=1)
  }
  if (!is.numeric(x)){
    stop("Argument \"x\" should be numeric")
  }
  if (ncol(x) != ncol(data)){
    stop("Dimensions of the arguments \"x\" and \"data\" should coincide")
  }
  
  is.in.convex <- .count_convexes(x, data, cardinalities, seed)
  
  return (is.in.convex)
}
