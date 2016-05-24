`summary.predict.homals` <-
function(object, ...)
{
# summary method for predict.homals
  cat("\nClassification Tables (observed vs. predicted): \n\n")
  for (i in 1:length(object$cl.table))
  {
    cat("Variable:",names(object$cl.table)[[i]],"\n")
    print(object$cl.table[[i]])
    cat("\n")
  }
}

