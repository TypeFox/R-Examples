varCounts.wsrf <- function(object) {
  # Return the times of each variable being selected as split condition.
  # For evaluating the bias of wsrf towards attribute types (categorical and
  # numerical) and the number of values each attribute has.

  varnames <- object[[.META_IDX]][["varnames"]]
  trees    <- object[[.TREES_IDX]]

  counts        <- vector("integer", length(varnames))
  names(counts) <- varnames

  for (i in 1:length(trees))
  {
    tree   <- trees[[i]]
    for (j in 1:length(tree))
    {
      node <- tree[[j]]
      if (as.integer(node[1]) == 1)
      {
        varidx         <- as.integer(node[4]) + 1
        counts[varidx] <- counts[varidx] + 1
      }
    }
  }

  return(counts)
}
