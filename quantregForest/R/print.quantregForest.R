"print.quantregForest" <-
function(x, ...) {
  cat("\n Call:\n", deparse(x$call), "\n\n")
  
  cat("                     Number of trees: ", x$ntree, "\n",sep="")
  cat("No. of variables tried at each split: ", x$mtry, "\n\n", sep="")
}

