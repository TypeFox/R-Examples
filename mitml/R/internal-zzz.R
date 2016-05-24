.extractMatrix <- function(x, ...){
# x: an array of at least two dimensions
# ...: indices beyond the first two that reduce x to a matrix

  if(is.null(dim(x))) return(x)

  out <- `[`(x,,,...)
  dim(out) <- dim(x)[1:2]
  dimnames(out) <- dimnames(x)[1:2]

  out

}
