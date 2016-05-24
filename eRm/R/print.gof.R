print.gof <- function(x, ...)
{
  #print method for objects of class "gof" (from gofIRT.ppar)
  cdv <- round(x$test.table[1,], 3)
  cat("\nGoodness-of-Fit Results:")
  cat("\nCollapsed Deviance = ", cdv[1], " (df = ", cdv[2], ", p-value = ", cdv[3], ")", sep ="")
  cat("\nPearson R2:", round(x$R2$R2.P, 3))
  cat("\nArea Under ROC:", round(x$AUC, 3))
  cat("\n\n")
}