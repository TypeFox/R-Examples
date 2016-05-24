##
## Methods for returning an informal "logLik" object
## =================================================
##
## Method definition for objects of class "Goestml"
##
setMethod(f = "logLik", signature = "Goestml", definition = function(object){
  r <- -1.0 * object@opt$objective
  df <- ncol(object@X) * sum(object@models[[1]]@fit$params$include) + length(angles(object))
  attr(r, "df") <- df
  class(r) <- "logLik"
  return(r)
})
