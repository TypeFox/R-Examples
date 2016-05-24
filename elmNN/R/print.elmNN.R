print.elmNN <-
function(x,...) {
  cat("Call:\n")
  cat(paste(x$call, "\n"))
  cat("Number of hidden neurons:\n")
  cat(paste(x$nhid, "\n"))
  cat("Activation function:\n")
  cat(paste(x$actfun, "\n"))
  cat("Input arc weights:\n")
  cat(paste(head(x$inpweight), "...\n"))
  cat("Bias of hidden neurons:\n")
  cat(paste(head(x$biashid), "...\n"))
  cat("Output arc weights:\n")
  cat(paste(head(x$outweight), "...\n"))
  cat("Predictions on training set:\n")
  cat(paste(head(fitted(x)), "...\n"))
}
