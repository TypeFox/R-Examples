print.ellipsefit <- function(x,...) {
  cat("Call:\n")
  print(x$call)
  Td <- qt(0.975,x$fit.statistics["d.f."]) 
  if (x$method!="direct") {
  cat("\nDelta Method Standard Errors and 95% C.I.'s:\n")
  error <- x$Std.Errors
  theEstimates <- x$values[names(error)]
  low2 <- theEstimates-error*Td
  high2 <- theEstimates+error*Td
  print(cbind("Estimates"=theEstimates,
              "S.E."=error,"low"=low2,"high"=high2),digits=4) }
  else {
    cat("\nEstimates:\n")
    print(x$values)
  }
invisible(x)}
