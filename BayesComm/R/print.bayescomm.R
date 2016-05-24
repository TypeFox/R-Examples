print.bayescomm <-
function (x, ...) {  # need number of iterations / thinning etc.!!
  cat("model type:  ", x$call$model)
  cat("\nobservations:  ", nrow(x$call$Y))
  cat("\nspecies:  ", colnames(x$call$Y))
  cat("\ncovariates:  ", colnames(x$call$X))
  cat("\niterations: ", x$call$its, "\tthin: ", x$call$thin, "\tdiscarded: ", (x$call$start - 1))
}