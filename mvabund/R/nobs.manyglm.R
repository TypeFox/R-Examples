nobs.manyglm<- function(object, ...)
{
  n.rows = NROW(object$y)
  return( n.rows )
  NextMethod("nobs")
}
