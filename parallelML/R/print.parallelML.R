#' @export
print.parallelML <- function(x, ...){
  
  # Name of the object
  fcall      <- match.call()
  objectName <- fcall$x
  
  # Collect information
  number   <- length(x)
  MLclass  <- class(x[[1]])
  sampSize <- 100 * attr(x,"samp")
  call     <- attr(x,"call")
  
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
  cat("The first model looks like")
  cat("\n")
  print(x[[1]])
  
  # information on how to print the other models
  cat("\n")
  cat("To get information about the i'th model, use print(")
  cat(objectName)
  cat(")[[i]]")
}