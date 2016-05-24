"print.learnPattern" <-
function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n")
  cat("               Type of random forest: ", x$type, "\n", sep="")
  cat("                     Number of trees: ", x$ntree, "\n",sep="")
  cat("No. of variables tried at each split: ", x$mtry, "\n\n", sep="")

  if(x$type == "regression") {
    if(!is.null(x$mse)) {
      cat("          Mean of squared residuals: ", x$mse[length(x$mse)],
          "\n", sep="")
      cat("                    % Var explained: ",
          round(100*x$rsq[length(x$rsq)], digits=2), "\n", sep="")
      if(!is.null(x$test$mse)) {
        cat("                       Test set MSE: ",
            round(x$test$mse[length(x$test$mse)], digits=2), "\n", sep="")
        cat("                    % Var explained: ",
            round(100*x$test$rsq[length(x$test$rsq)], digits=2), "\n", sep="")
      }      
    }
  }
}
