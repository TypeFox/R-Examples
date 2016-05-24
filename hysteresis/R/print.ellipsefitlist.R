print.ellipsefitlist <- function(x,...) {
  cat("Call:\n")
  print(attr(x,"call"))
  cat("\nParameter Estimates:\n")
print(x$Estimates,dixits=4)
  
cat("\nDelta Method Standard Errors:\n")
print(x$Std.Errors,digits=4)
invisible(x)}
