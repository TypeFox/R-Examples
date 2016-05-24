rerun <- function(object, ...) {
  if (!is.element(class(object)[1], c("wls", "tssem1FEM", "tssem1REM", "meta", "meta3X", "reml", "tssem1FEM.cluster", "wls.cluster")))
    stop("\"object\" must be an object of neither class \"meta\", \"meta3X\", \"wls\", \"reml\", \"tssem1FEM\", \"tssem1REM\", \"tssem1FEM.cluster\" or \"wls.cluster\".")

  if (is.element(class(object)[1], c("tssem1FEM.cluster", "wls.cluster"))) {
    out <- lapply(object, rerun, ...)
    class(out) <- class(object)[1]
  } else {
    out <- object   
    # No LB option
    if (is.null(object$intervals.type)) {
      out$mx.fit <- mxTryHard(object$mx.fit, greenOK=TRUE, paste=FALSE, bestInitsOutput=FALSE, ...)
    } else {
      switch(object$intervals.type,
             z = out$mx.fit <- mxTryHard(object$mx.fit, greenOK=TRUE, paste=FALSE, bestInitsOutput=FALSE, ...),
             LB =out$mx.fit <- mxTryHard(object$mx.fit, greenOK=TRUE, paste=FALSE, bestInitsOutput=FALSE, intervals=TRUE, ...))
    }
  }
  out
}
