"summary.crossval" <-
function(object, axes = c(1:min(6, object$n.axes)), ...)
  {
    class(object) <- "summary.crossval"
    object
  }

