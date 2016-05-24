oob.error.rate <- function(object, ...) UseMethod("oob.error.rate")
oob.error.rate.wsrf <- function(object, tree, ...)
{
  if (missing(tree) || (is.logical(tree) && length(tree) == 1 && !tree))
  {
    # return out-of-bag error rate for the forest, length of 1
    return(object[[.RF_OOB_ERROR_RATE_IDX]])
  }
  else
  {
    # return out-of-bag error rates for specific trees
    return(object[[.TREE_OOB_ERROR_RATES_IDX]][tree])
  }
}
