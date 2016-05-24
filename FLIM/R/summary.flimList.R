summary.flimList <-
function(object, ...) {
  x <- object
  models <- x$models
  nm <- length(models)
  if(is.null(x$flim.obj$lambda)) {
    for(k in 1:nm) {
      cat("Call:\n")
      cat("  Model: ")
      print(models[[k]], showEnv=F)
      cat("   Data: ")
      print(x$flim.obj$call$data)
      cat("  Times: ")
      cat(x$flim.obj$times[1:(length(x$flim.obj$times)-1)], "\n")   
      cat("\n Coefficients:\n")
      coeff <- names(x$est.list[[k]])
      su <- x$sum.list[[k]]
      for(j in coeff) {
        suj <- su[rownames(su) == j, ]
        rownames(suj) <- x$flim.obj$times[1:(length(x$flim.obj$times) - 1)]
        cat("\n", j, "\n")
        print(suj)
      }
      cat("\n")
    }
  } else {print("Summary unavailable for flimList with non NULL lambda")}
}
