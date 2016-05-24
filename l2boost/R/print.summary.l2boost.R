
# Unimplemented generic function
# These are placeholders right now.
# @param x an l2boost object 
# @param ... other arguments (not used)
#
print.summary.l2boost <- function(x, ...){
  stop("Unimplemented function")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
}
