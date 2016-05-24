print.randomize <- function(x, digits = getOption("digits"), ...) {

  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  cat("\n Treatment Assignment:\n")
  print(table(x$treatment, exclude = NULL, ...), digits = digits)

  if (!is.null(x$block)) {
    cat("\n Treatment Assignment by Blocks:\n\n")
    print(ftable(x$block, x$treatment, exclude = NULL, ...),
          digits = digits)
  }
  
  cat("\nTotal number of observations:",
      length(na.omit(x$treatment)), "\n\n")
  invisible(x)
}
