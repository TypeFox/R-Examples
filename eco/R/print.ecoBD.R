print.ecoBD <- function(x, digits = max(3, getOption("digits") -3),
                        ...) {
  cat("\nCall:\n", deparse(x$call), "\n\n", sep="")
  
  cat("Aggregate Lower Bounds (Proportions):\n")
  print.default(format(x$aggWmin, digits = digits), print.gap = 2, quote =
                FALSE)
  cat("\nAggregate Upper Bounds (Proportions):\n")
  print.default(format(x$aggWmax, digits = digits), print.gap = 2, quote =
                FALSE)
  
  if (!is.null(x$aggNmin)) {
    cat("\nAggregate Lower Bounds (Counts):\n")
    print.default(format(x$aggNmin, digits = digits), print.gap = 2, quote =
                  FALSE)
    cat("\nAggregate Upper Bounds (Counts):\n")
    print.default(format(x$aggNmax, digits = digits), print.gap = 2, quote =
                  FALSE)
  }
  
  cat("\n")
  invisible(x)
}
