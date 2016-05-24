fitted.repolr <-
function (object, ...) {
  fit <- object$fitted.values
  names(fit) <- paste(object$id, 1:object$max.id, sep=".")
  fit
}
