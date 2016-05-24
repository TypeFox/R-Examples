subset.wsrf <- function(x, trees, ...)
{
  ## Get a subset of a wsrf model.

  # x is a object of wsrf.
  # i is a subset of trees indexes.
  if (!inherits(x, "wsrf")) 
    stop("Not a legitimate wsrf object")

  tags <- c(.TREES_IDX, .TREE_OOB_ERROR_RATES_IDX, .OOB_SETS_IDX, .OOB_PREDICT_LABELS_IDX, .TREE_IGR_IMPORTANCE_IDX)

  res <- vector("list", .WSRF_MODEL_SIZE)
  names(res) <- .WSRF_MODEL_NAMES

  for (tag in tags)
    res[[tag]] <- x[[tag]][trees]

  res[[.META_IDX]]           <- x[[.META_IDX]]
  res[[.TARGET_DATA_IDX]]    <- x[[.TARGET_DATA_IDX]]

  if (!is.null(x[[.WEIGHTS_IDX]])) res[[.WEIGHTS_IDX]] <- x[[.WEIGHTS_IDX]]
  else res[.WEIGHTS_IDX] <- list(NULL)

  if (!is.null(x[[.MTRY_IDX]])) res[[.MTRY_IDX]] <- x[[.MTRY_IDX]]
  else res[.MTRY_IDX] <- list(NULL)

  .Call("afterMergeOrSubset", res, PACKAGE="wsrf")

  class(res) <- "wsrf"

  return(res)
}
