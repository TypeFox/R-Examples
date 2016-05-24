#dns20110811; taken from the original random forest package

"print.obliqueRF" <-
function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n")
  cat("               Type of random forest: ", x$type, "\n", sep="")
  cat("                     Number of trees: ", x$ntree, "\n",sep="")
  cat("No. of variables tried at each split: ", x$mtry, "\n\n", sep="")
#  if(x$type == "classification") {
#    print(x$errs);
#  }
}
