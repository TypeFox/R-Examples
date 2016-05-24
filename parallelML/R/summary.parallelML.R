#' @export
summary.parallelML <- function(object, ...){
  
  # Name of the object
  fcall      <- match.call()
  objectName <- fcall$object
  
  # Collect information
  number   <- length(object)
  MLclass  <- class(object[[1]])
  sampSize <- 100 * attr(object,"samp")
  call     <- attr(object,"call")
  
  # Print how many items there are and how much data they used
  cat("List of ")
  cat(number)
  cat(" models of class ")
  cat(MLclass)
  cat(" using ")
  cat(sampSize)
  cat(" % of the original data")
  
  # Print the call
  cat("\n\n")
  cat("Call:")
  cat("\n")
  cat("    ")
  cat(call)
  
  # regular summary of first item
  cat("\n\n")
  cat("The regular summary of the first model looks like")
  cat("\n")
  print(summary(object[[1]]))
  
  # information on how to print the other models
  cat("\n")
  cat("To get information about the i'th model, use summary(")
  cat(objectName)
  cat(")[[i]]")
}