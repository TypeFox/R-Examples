fitted.ssfa <- function(object, ...){
  
  return(object$y  - residuals.ssfa(object))
  
}