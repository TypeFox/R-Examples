effectInfo <-
function(object)
 {  
  if (!inherits(object, "tlm"))
    stop("argument 'object' must be of class 'tlm'")
    
  if(any(is.na(coef(object$model))))
     stop("effectInfo is not available for models with any missing estimated coefficient")
 
  modeltype <- modelType(object = object)
  eval(parse(text = paste("res <- effectInfomod", modeltype, "(object = object)", sep = "")))
  attr(res, "modeltype") <- modeltype
  class(res) <- "effectInfo"
  res
 }
