`print.predict.homals` <-
function(x, ...)
{
# print method for predict.homals
  cat("\nClassification rate:\n")
  dframe <- data.frame(names(x$cr.vec), round(x$cr.vec,4), round(x$cr.vec*100,2))
  rownames(dframe) <- NULL
  colnames(dframe) <- c("Variable","Cl. Rate","%Cl. Rate")
  print(dframe)
  cat("\n")
}

