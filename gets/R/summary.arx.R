summary.arx <-
function(object, ...)
{
  object.name <- deparse(substitute(object))
  xnames <- cbind(names(object))
  xrows <- nrow(xnames)
  colnames(xnames) <- ""
  rownames(xnames) <- rep(" ",xrows)
  cat("\n")
  cat("Items in list '", object.name, "':", sep="")
  cat("\n")
  print(xnames, quote=FALSE)
  cat("\n")
}
