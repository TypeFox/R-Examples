################################################################################
# File:             depth.zonoid.r
# Created by:       Pavlo Mozharovskyi
# First published:  28.02.2013
# Last revised:     15.05.2013
# 
# Computation of the zonoid data depth.
################################################################################

depth.zonoid <- function(x, data, seed = 0){
  if (seed!=0) set.seed(seed)
  if (!is.matrix(x) 
      && is.vector(x)){
    x <- matrix(x, nrow=1)
  }
  if (!(is.matrix(data) && is.numeric(data)
        || is.data.frame(data) && prod(sapply(data, is.numeric))) 
      || ncol(data) < 2){
    stop("Argument \"data\" should be a numeric matrix of at least 2-dimensional data")
  }
  if (!is.numeric(x)){
    stop("Argument \"x\" should be numeric")
  }
  
  if (ncol(x) != ncol(data)){
    stop("Dimensions of the arguments \"x\" and \"data\" should coincide")
  }
  if (ncol(data) + 1 > nrow(data)){ #?
    stop("To few data points")
  }
  
  points <- as.vector(t(data))
  objects <- as.vector(t(x))
  ds <- .C("ZDepth", 
           as.double(points), 
           as.double(objects), 
           as.integer(nrow(data)), 
           as.integer(nrow(x)), 
           as.integer(ncol(data)), 
           as.integer(seed),
           depths=double(nrow(x)))$depths
  
  return (ds)
}
