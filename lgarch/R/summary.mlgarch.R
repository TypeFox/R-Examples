summary.mlgarch <-
function(object, ...)
{
  xnames <- cbind(names(object))
  xrows <- nrow(xnames)
  colnames(xnames) <- ""
  rownames(xnames) <- rep(" ",xrows)
  cat("\n")
  cat("Items in list:")
  cat("\n")
  print(xnames)
  cat("\n")
}
