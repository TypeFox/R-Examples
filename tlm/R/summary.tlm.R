summary.tlm <-
function(object, ...)
 {
  object$summary <- summary(object$model, ...)
  class(object) <- "summary.tlm"
  object
 }
