leaf.numbers <- function(tree)
  {
    where <- tree$where
    leaves <- as.numeric(dimnames(tree$frame)[[1]])
    leaves[where]
  }
