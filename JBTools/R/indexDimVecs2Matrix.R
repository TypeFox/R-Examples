indexDimVecs2Matrix <- function(
  ##title<< Transform integer indices to an index matrix
  ... ##<< integer vectors: indices to use for the different dimensions
  ) {
  ##description<< Transform integer indices to an index matrix. The input have to be
  ##              index vector for each dimension of the target array.
  ##\code{\link{Extract}},  \code{\link{indexVec2Matrix}}
  dummy=list(...)
  ##value<< Index matrix
  index.matrix <- as.matrix(do.call("expand.grid", list(...)))
  return(index.matrix)
}
