print.ellipsesummary <-
function(x,...) {
  g <- x
  cat("Summary Call:\n")
  print(g$summarycall)
  cat("Call for Original Fit:\n")
  print(g$call)
  cat("Ellipse Fitting Method:\n")
  print(g$method)
  if (g$method=="harmonic2") print("Two step simple harmonic least squares")
  else if (g$method=="nls") print("Non-linear least squares")
  else if (g$method=="direct") print("Direct specific least squares")
  else print("Linear Least Squares")
  if (g$boot==TRUE) {
  cat("\nBootstrapped Estimates:\n")
  print(g$values[,c("Boot.Estimate","Bias","Std.Error","B.q0.025","B.q0.975")],digits=4)
  }
  else {
    Td <- qt(0.975,g$fit.statistics["d.f."]) 
    cat("\nDelta Method Standard Errors and 95% C.I.'s:\n")
    error <- g$Std.Error
    thevalues <- g$values[names(error)]
    low2 <- thevalues-error*Td
    high2 <- thevalues+error*Td
    print(cbind("Estimates"=thevalues,
                "S.E."=error,"low"=low2,"high"=high2)[c("ampx","ampy","rote.deg","retention","coercion",
                                                       "area","lag","split.angle"),],digits=4) 
  }
 ## cat("\nFit Statistics:\n")
  ## print(g$fit.statistics,digits=4)
invisible(g)}
