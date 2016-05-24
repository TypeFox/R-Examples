print.mitml <- function(x,...){
# print method for objects of class "mitml"

  cl <- x$call
  vrs <-x$model 
  itr <- x$iter
  ngr <- length(unique(attr(x$data,"group")))
  cat("\nCall:\n", paste(deparse(cl)), sep="\n")

  cat("\nCluster:\t\t\t", vrs$clus, sep=" ", collapse="\n")
  cat("Target variables:\t\t", vrs$yvrs, collapse="\n")
  cat("Fixed effect predictors:\t", vrs$pvrs, collapse="\n")
  cat("Random effect predictors:\t", vrs$qvrs, collapse="\n")

  cat("\nPerformed", sprintf("%.0f",itr$burn), "burn-in iterations, and generated", sprintf("%.0f",itr$m),
      "imputed data sets,\neach", sprintf("%.0f",itr$iter), "iterations apart.",
      if(ngr>1){c("\nImputations were carried out seperately within", sprintf("%.0f",ngr), "groups.\n")},"\n")

  invisible()
}

