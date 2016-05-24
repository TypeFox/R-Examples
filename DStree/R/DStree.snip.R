#' @export
snip <- function(tree, cp)
  #taken from rpart function prune
{
  ff <- tree$frame
  id <- as.integer(row.names(ff))
  toss <- id[ff$complexity <= cp & ff$var != "<leaf>"] 
  if (length(toss) == 0L) return(tree)  
  class(tree)<-"rpart"
  newx <- snip.rpart(tree, toss)
  temp <- pmax(tree$cptable[, 1L], cp)
  keep <- match(unique(temp), temp)
  newx$cptable <- tree$cptable[keep, , drop = FALSE]
  newx$cptable[max(keep), 1L] <- cp
  class(newx)<-"DStree"
  newx
}
