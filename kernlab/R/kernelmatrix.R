
setGeneric("as.kernelMatrix",function(x, center = FALSE) standardGeneric("as.kernelMatrix"))
setMethod("as.kernelMatrix", signature(x = "matrix"),
function(x, center = FALSE)
{

  if(center){
    m <- dim(x)[1]
    x <- t(t(x - colSums(x)/m) -  rowSums(x)/m) + sum(x)/m^2
  }
  
  return(new("kernelMatrix",.Data = x))
})
