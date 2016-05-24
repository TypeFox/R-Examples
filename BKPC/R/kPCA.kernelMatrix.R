kPCA.kernelMatrix <-
function(x, ...){
  old.class <- class(x)
  class(x) <- "kern"
  object <- kPCA.kern(x)
  class(object$K) <- old.class
  class(object) <- c("kPCA.kernelMatrix", class(object)) 
  return(object)
}
