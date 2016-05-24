treesize <- function(x, terminal=TRUE) {
  if(!inherits(x, "RRF"))
    stop("This function only works for objects of class `RRF'")
  if(is.null(x$forest)) stop("The object must contain the forest component")
  if(terminal) return((x$forest$ndbigtree+1)/2) else return(x$forest$ndbigtree)
}
