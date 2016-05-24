################################################################################
# File:             depth.space.zonoid.r
# Created by:       Pavlo Mozharovskyi
# First published:  28.02.2013
# Last revised:     15.05.2013
# 
# Computation of the depth space based on the zonoid data depth.
################################################################################

depth.space.zonoid <- function(data, cardinalities, seed = 0){
  if (seed!=0) set.seed(seed)
  if (!(is.matrix(data) && is.numeric(data)
        || is.data.frame(data) && prod(sapply(data, is.numeric))) 
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
  
  dim <- ncol(data)
  depth.space <- NULL
  for (i in 1:length(cardinalities)){
    objects <- as.vector(t(data[(1 + sum(cardinalities[0:(i - 1)])):sum(cardinalities[1:i]),]))
    pattern.depths <- NULL
    for (j in 1:length(cardinalities)){
      pattern <- as.vector(t(data[(1 + sum(cardinalities[0:(j - 1)])):sum(cardinalities[1:j]),]))
      ds <- .C("ZDepth", 
               as.double(pattern), 
               as.double(objects), 
               as.integer(cardinalities[j]), 
               as.integer(cardinalities[i]), 
               as.integer(dim), 
               as.integer(seed),
               depths=double(cardinalities[i]))$depths
      if (j == i){
        ds <- replace(ds, which(ds < 1/cardinalities[j] - sqrt(.Machine$double.eps)), 1/cardinalities[j])
      }else{
        ds <- replace(ds, which(ds < 1/cardinalities[j] - sqrt(.Machine$double.eps)), 0)
      }
      pattern.depths <- cbind(pattern.depths, ds, deparse.level = 0)
    }
    depth.space <- rbind(depth.space, pattern.depths, deparse.level = 0)
  }

  return (depth.space)
}
