summary.eiReg <- function(object, ...){ 
  nidx <- apply(expand.grid(dimnames(object$coef)), 1, paste, collapse = ".")
  object$coef <- cbind(c(object$coef), c(object$se), 
                       c(object$coef)/c(object$se))
  colnames(object$coef) <- c("Estimate", "Std. Error", "t-stat")
  rownames(object$coef) <- nidx
  object
}
