print.fittedlooplist <- function(x,...) {
  cat("Call:\n")
  print(attr(x,"call"))
  cat("\nParameter Estimates:\n")
print(x$Estimates,dixits=4)
  
cat("\nDelta Method Standard Errors:\n")
print(x$Std.Errors,dixits=4)
invisible(x)}

print.fittedlooplist2r <- function(x,...) {
  cat("Call:\n")
  print(attr(x,"call"))
  cat("\nParameter Estimates:\n")
print(x$Estimates,dixits=4)
  
cat("\nDelta Method Standard Errors:\n")
print(x$Std.Errors,dixits=4)
invisible(x)}
