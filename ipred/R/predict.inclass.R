# $Id: predict.inclass.R,v 1.19 2003/03/31 08:44:16 peters Exp $

# Additional option type ="class", if intermediate is nominal

predict.inclass <- function(object, newdata, ...)
{
  if(!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
  q <- length(object$model.intermediate)		# number of intermediates
  namen <- names(object$model.intermediate) 

  intermediate <- is.vector(NULL, mode = "NULL")

  for(i in 1:q) {
    if(!is.null(object$para.intermediate[[i]][["predict"]])) {
      RET <- object$para.intermediate[[i]][["predict"]](object$model.intermediate[[i]], newdata = newdata, ...)
    } else {
      RET <- predict(object$model.intermediate[[i]], newdata = newdata, ...)
    }
  intermediate <- data.frame(intermediate, RET)
  }

  intermediate <- intermediate[,-1]  
  names(intermediate) <- namen
  
  intermediate <- data.frame(newdata[,!(names(newdata) %in% names(intermediate))], intermediate)

  if(!is.function(object$para.response)) {
   if(!is.null(object$para.response[["predict"]])) {
       RET <- object$para.response[["predict"]](object$model.response, newdata = intermediate, ...)
     } else {
       RET <- predict(object$model.response, newdata = intermediate, ...)
     }
  } else {
    RET <- object$para.response(intermediate)
  }
  return(RET)
}
