#' Strip objects unnecessary for prediction with class.
#'
#' @param model_object the model to strip pre-deployment.
#' @param model_class model class string e.g. "lm", "glm", "gam" ...
#'
#' @return None
#' @export
#' @examples
#' 
#' example_data <- as.data.frame(cbind(gl(3,50),rnorm(150)));names(example_data) <- c("x","y")
#' example_fit  <- lm(y~x,data=example_data)
#' strip_model(example_fit,"lm")

strip_model <- function(model_object,model_class=NULL){  
  if (!is.null(model_class)&model_class%in%c("lm","glm","randomForest","gam"))
  {
    if(sum(is.element(base::class(model_object)[1],c("lm","glm")))>0) {
      model_object$residuals <- NULL
      model_object$effects   <- NULL
      model_object$fitted.values <- NULL
      model_object$model     <- NULL    
    }
    if(sum(is.element(base::class(model_object)[1],c("glm")))>0) {
      model_object$linear.predictors <- NULL
      model_object$prior.weights   <- NULL
      model_object$weights <- NULL
      model_object$y     <- NULL    
      model_object$model     <- NULL    
      model_object$data     <- NULL    
    }
    if(sum(is.element(base::class(model_object),c("randomForest")))>0) {
      model_object$predicted <- NULL
      model_object$oob.times   <- NULL
      model_object$y     <- NULL    
    }
    if(sum(is.element(base::class(model_object),c("gam")))>0) {
      model_object$smooth <- NULL
      model_object$pred.formula   <- NULL
      model_object$offset     <- NULL    
    }
  }
  return(model_object)
}


