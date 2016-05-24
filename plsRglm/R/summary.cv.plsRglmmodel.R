summary.cv.plsRglmmodel <- function(object, ...)
{
  res <- kfolds2CVinfos_glm(object, ...)
  if(is.null(object$call$model)){
    class(res) <- "summary.cv.plsRmodel"
  } else {
  if(object$call$model=="pls"){
    class(res) <- "summary.cv.plsRmodel"
  } else {
    class(res) <- "summary.cv.plsRglmmodel"
  }
  }
  return(res)
}
