print.activeset <- function(x, ...)
{
  cat("\nCall:\n")
  print(x$call)
  cat("\nLoss value:", x$fval, "\n")
  cat("\nActive set fit:\n")
  res.table <- t(data.frame(round(x$y, 3), round(x$x,3)))
  rownames(res.table) <- c("Observed Values", "Fitted Values")
  colnames(res.table) <- 1:length(x$x)
  print(t(res.table))
  cat("\n")
  invisible(x)
}
