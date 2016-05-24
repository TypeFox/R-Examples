print.flimList <-
function(x, ...) {
  models <- x$models
  for(k in 1:length(models)) {
    cat("Call:\n")
    cat("  Model: ")
    print(models[[k]], showEnv=F)
    cat("   Data: ")
    print(x$flim.obj$call$data)
    cat("  Times (t): ")
    cat(x$flim.obj$times[1:(length(x$flim.obj$times)-1)], "\n")
    if(!is.null(x$flim.obj$lambda)) {
      cat(" lambda: ")
      cat(x$flim.obj$lambda, "\n")
    }
    cat("\n Coefficients:\n")
    print(x$est.list[[k]])
    cat("\n")
  }
}
