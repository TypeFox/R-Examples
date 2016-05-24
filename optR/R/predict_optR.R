#' Prediction function based on optR class  
#'
#'@description Function for making predictions using OptR class
#'@param object     : optR class fitted object
#'@param newdata       : data for prediction 
#'@param na.action  : action for missing values
#'@param ...        : S3 class
#'@return fitted.val   : Predicted values
#'@return terms       : terms used for fitting
#'@export
predict.optR<-function(object, newdata, na.action=na.omit, ...) {  
  # Extract terms used into the model
  tt<-terms(object)
  if (!inherits(object, "optR")) 
    warning("calling predict.lm(<fake-optR-object>) ...")
  if (missing(newdata) || is.null(newdata)) {
    stop("Missing or Null dataset ...")
  } else 
  {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.action, xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses")))
      .checkMFClasses(cl, m)
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    beta <- object$beta
    predictor<-as.matrix(X)%*%as.matrix(beta)
  }
  
  # Extract terms from newdata
  xterms<-colnames(m)
  
  
  # Return predicted values
  return(list("fitted.val" = predictor, "terms"=xterms))
}
