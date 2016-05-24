summary.genlasso <- function(object, ...) {
  n = length(object$y)
  p = nrow(object$beta)
  df = object$df
  rss = colSums((object$y - object$fit)^2)

  mat = cbind(object$df, object$lambda, rss)
  rownames(mat) = rep("",nrow(mat))
  colnames(mat) = c("df", "lambda", "rss")

  class(mat) = c("summary.genlasso", "matrix")
  return(mat)
}
