
getCells <- function(x) {
  return(dimnames(x)[[1]])
}

"getCells<-" <- function(x,value) {
  if(length(value)!=ncells(x)) stop("Wrong number of cell names supplied!")
  if(ncells(x)==0) return(x)
  if(is.null(value)) stop("Setting cell names to NULL is not allowed!")
  if(length(value)==1) value <- list(value)
  dimnames(x)[[1]] <- value
  return(x)
}
