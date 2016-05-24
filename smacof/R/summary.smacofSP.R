`summary.smacofSP` <-
function(object, ...)
{
  cat("\n")
  cat("Configurations:\n")
  print(round(object$conf,4))
  #cat("\nConfiguration dissimilarities: \n")
  #print(round(object$confdiss,4))
}

