residuals.ssfa <- function(object, ...) {
  
  if(object$rho==0)
  {
  beta <- as.matrix(object$coef[3:length(object$coef)])
  X <- as.matrix(cbind(1, object$x))
  y <- as.matrix(object$y)
  residuals_ssfa <- (y - X%*%beta)
  }
  if(object$rho!=0)
  {
  beta <- as.matrix(object$coef[4:length(object$coef)])
  X <- as.matrix(cbind(1, object$x))
  y <- as.matrix(object$y)
  DIM_w <- dim(object$w)[1]
  I_m <- diag(DIM_w)
  residuals_ssfa <- (I_m - object$rho *object$w)%*%(y - X%*%beta)
  }
  
  
  return(residuals_ssfa)
}