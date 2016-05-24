block.diag <-
function(...)
{
  matrices <- list(...)
  k <- length(matrices)
  nrows <- sapply(matrices, nrow)
  ncols <- sapply(matrices, ncol)
  
  index.row <- rep(seq_len(k), nrows)
  index.col <- rep(seq_len(k), ncols)
  M <- matrix(0, sum(nrows), sum(ncols))
  for (i in seq_len(k)) {
    M[index.row == i, index.col == i] <- matrices[[i]]
  }
  
  M
}
