as.matrix.alignment <- function(x, ...){
  nc <- nchar(x$seq[[1]])
  res <- matrix(data = sapply(x$seq, s2c), nrow = x$nb, ncol = nc,
                byrow = TRUE)
  rownames(res) <- x$nam
  colnames(res) <- seq_len(nc)
  return(res)
}
