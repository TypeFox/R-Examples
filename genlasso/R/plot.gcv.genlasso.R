plot.gcv.genlasso <- function(object, ...) {
  if (!any(class(object)=="gcv.genlasso")) {
    stop("Passed object must be of class \"gcv.genlasso\".")
  }
   
  plot(object$lambda, object$err, type="l", log="x",
       xlab=expression(lambda), ylab="GCV error", ...)
  points(object$lambda.min, min(object$err), pch=19, col="blue")
}
