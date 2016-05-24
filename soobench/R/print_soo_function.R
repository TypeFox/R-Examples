##' Print a SOO function.
##'
##' @param x A \code{soo_function} object.
##' @param ... Ignored.
##' @S3method print soo_function
##' @method print soo_function
print.soo_function <- function(x, ...) {
  cat(function_name(x), "\n", 
      "  Lower bounds: (", paste(lower_bounds(x), collapse=", "), ")\n",
      "  Upper bounds: (", paste(upper_bounds(x), collapse=", "), ")\n", sep="")
}
