`print.LR` <-
function(x,...)
{
#print method for object of class "LR" (LRtest)
  cat("\n")
  cat("Andersen LR-test: \n")
  cat("LR-value:", round(x$LR,3),"\n")
  cat("Chi-square df:",x$df,"\n")
  cat("p-value: ",round(x$pvalue,3),"\n")
  cat("\n")
}

