logLik.hglm <- function(object, REML=NULL, ...) {
  if (is.null(REML)) { 
    val <- object$likelihood 
  }else{
    if (REML) val <- object$likelihood$pbvh 
  }
  return(val)
}
  