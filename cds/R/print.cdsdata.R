#' Print dsdata Objects
#' 
#' This is a simple print method  for object that inherits from the class \code{cdsdata}.
#' 
#' @param x A \code{cdsdata} object
#' @param \dots Unimplemented.
#' @export
#' @method print cdsdata
print.cdsdata <- function(x, ...) {
  cat("Object of S3 Class 'cdsdata'\n\n")
  cat(nrow(x$postrs), "observations on", x$m, "variables\n\n")
  cat("Likert scale used:\n", x$scales, "\n\n")
  cat("Available slots:\n")
  print(names(x))
}