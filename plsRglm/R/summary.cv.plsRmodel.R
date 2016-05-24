summary.cv.plsRmodel <- function(object, ...)
{
  res <- kfolds2CVinfos_lm(object, ...)
  class(res) <- "summary.cv.plsRmodel"
  return(res)
}
