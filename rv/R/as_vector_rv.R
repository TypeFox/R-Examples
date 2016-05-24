

as.vector.rv <- function(x, mode="any") {
  a <- attributes(x) 
  x <- lapply(unclass(x), as.vector, mode=mode)
  attributes(x) <- a
  dim(x) <- NULL
  return(x)
}


