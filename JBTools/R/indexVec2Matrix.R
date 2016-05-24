indexVec2Matrix = function(
  ##title<< Transform a vector index to an index matrix
  index ##<< vector index 
  , dim ##<<  length of the two dimensions of the array to index. Identical to
        ##    the result of dim(array).
  )
  ##description<<
  ## Transform an index vector into a index matrix.
  ##\code{\link{Extract}}
{
  if (length(dim) == 2) {
    dim1 <- index - (floor((index - 1)/dim[1]) * dim[1]) 
    dim2 <- floor((index - 1)/dim[1]) + 1 
  } else {
    stop('Only 2d is yet implemented!')
  }
  ##value<< index array
  out <-  cbind(dim1, dim2)
  return(out)
}
