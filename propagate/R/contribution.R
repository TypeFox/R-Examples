contribution <- function(object, plot = TRUE, ...) {
  if (class(object) != "propagate") stop("object must be of class 'propagate'!")
  EVAL <- object$evalGrad
  
  ## sensitivity: df/dx_i * u_i
  SENS <- EVAL * sqrt(diag(object$covMat))
  
  ## relative contribution
  CONTRIB <- abs(SENS)/sum(abs(SENS), na.rm = TRUE)  
  
  ## barplot of contributions
  if (plot) barplot(CONTRIB, ylim = c(0, 1), xlab = "Variable", ylab = "Rel. contribution", ...)
    
  return(CONTRIB)
}
