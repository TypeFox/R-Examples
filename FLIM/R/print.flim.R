print.flim <-
function(x, ...) {
  cat("Object of class \"flim\". Type object$dataset for reconstructed data.\n")
  cat("Fitted models are contained in object$fit. Use flimList(object), \n")
  cat("summary(flimList(object)) and plot(flimList(object)) to assess fits.\n")
  cat("\nCall:\n")
  print(x$call)
  cat("\nUse summary(object) for more information regarding the imputation,\n")
  cat("and see ?plot.flim for ploting options for the imputed measurements.\n")
  cat("\n")
}
