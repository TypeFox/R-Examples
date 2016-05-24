`slice` <-
function(x, MARGIN, n)
{
  result <- x[slice.index(x, MARGIN) == n]
  if(length(result) == prod(dim(x)[-MARGIN]))
  {
    dim(result) <- dim(x)[-MARGIN]
    dimnames(result) <- dimnames(x)[-MARGIN]
  }
  return(result)
}
