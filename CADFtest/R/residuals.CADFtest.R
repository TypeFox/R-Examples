residuals.CADFtest <- function(object, ...)
{
  # object is an object of class CADFtest
  residuals.lm(object$est.model)
}
