coef.mpcv <- function(object, ...)
{
  vars <- rownames(as.matrix(object$coef.lo))
  Coef <- matrix(c(object$coef.lo[,5], object$coef.up[,5]), nrow=2, byrow=T, dimnames=list(c("coef.lo", "coef.up"),vars))
  return(Coef)
}
