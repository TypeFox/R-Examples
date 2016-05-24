summary.intReg <- function(object, boundaries=FALSE, ...) {
   s <- NextMethod("summary", object, ...)
   estimate <- coefTable(coef(object, boundaries=boundaries),
                         stdEr(object, boundaries=boundaries),
                         object$param$df)
   s$estimate <- estimate
   s$param <- object$param
                           # supply the additional parameters
   class(s) <- c("summary.intReg", class(s))
   return(s)
}
   
