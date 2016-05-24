summary.gof <- function(object, ...)
{
  #summary method for objects of class "gof" (from gofIRT.ppar)
  cat("\nGoodness-of-Fit Tests\n")
  print(round(object$test.table, 3))

  cat("\nR-Squared Measures")
  cat("\nPearson R2:", round(object$R2$R2.P, 3))
  cat("\nSum-of-Squares R2:", round(object$R2$R2.SS, 3))
  cat("\nMcFadden R2:", round(object$R2$R2.MF, 3))
  
  cat("\n\nClassifier Results - Confusion Matrix (relative frequencies)\n")
  print(round(object$classifier$confmat/sum(object$classifier$confmat), 3))
  cat("\nAccuracy:", round(object$classifier$accuracy, 3))
  cat("\nSensitivity:", round(object$classifier$sensitivity, 3))
  cat("\nSpecificity:", round(object$classifier$specificity, 3))
  cat("\nArea under ROC:", round(object$AUC, 3))
  cat("\nGini coefficient:", round(object$Gini, 3))
  cat("\n\n")
}