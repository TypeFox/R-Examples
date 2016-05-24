changeFunctions <- function(expr, old, new) {
  if(is.recursive(expr)) {
    i <- match(as.character(expr[[1]]), old)
    if(!is.na(i))
      expr[[1]] <- as.name(new[[i]])
    n <- length(expr)
    if(n > 1) for(i in 2:n)
      expr[[i]] <- Recall(expr[[i]], old, new)
  }
  expr
}
