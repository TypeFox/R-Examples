dummy <- function(f,simplify=TRUE) {
  f <- droplevels(factor(f))
  mat <- I(model.matrix(~x-1,data=data.frame(x=f)))
  if (simplify) {
    res <- as.matrix(mat[,-ncol(mat)],nrow=length(f))
    colnames(res) <- levels(f)[-nlevels(f)]
  } else {
    res <- as.matrix(mat[,1:ncol(mat)],nrow=length(f))
    colnames(res) <- levels(f)
  }
  return(res)
}
