"as.select" <- function(x, ...) {
  cl <- match.call()
  
  if(is.select(x)) {
    return(x)
  }
  else {
    if(!is.vector(x))
      stop("provide a numeric vector of atom indices")
    
    if(all(is.logical(x)))
      x <- which(x)
    if(!all(is.numeric(x)))
      stop("provide a numeric vector of atom indices")
    
    sele <- NULL
    sele$atom <- x
    sele$xyz <- atom2xyz(x)
    sele$call <- cl
    class(sele) <- "select"
    return(sele)
  }
}
