bhtsne <- function(df) {
  X <- model.matrix(~.-1, df)
  Rtsne(X)
  out <- data.frame(Rtsne$Y)
  colnames(out) <- paste0("tsne", 1:ncol(out))
  return(out)
}