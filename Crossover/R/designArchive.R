getStartDesigns <- function(s,p,v) {
  objects <- ls(envir=Crossover.env)
  candidates <- list()
  for (x in objects) {
    x <- get(x, envir=Crossover.env)
    if (!is.matrix(x)) next
    if (!all(dim(x)==c(p,s))) next
    if (length(levels(as.factor(x)))!=v) next
    candidates[[length(candidates)+1]] <- x
  }  
  return(candidates)
}