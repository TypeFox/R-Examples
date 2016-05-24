#' @importFrom utils str
#' @export
print.densratio <- function(x, digits = 3L, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  info <- x$kernel_info
  cat("Kernel Information:\n")
  cat("  Kernel type: ", info$kernel, "\n")
  cat("  Number of kernels: ", info$kernel_num, "\n")
  cat("  Bandwidth(sigma): ", format(info$sigma, digits = digits), "\n")
  cat("  Centers: ")
  str(info$centers, digits.d = digits, give.attr = FALSE)
  cat("\n")

  cat("Kernel Weights(alpha):\n ")
  str(x$alpha, digits.d = digits, give.attr = FALSE)
  cat("\n")

  if("uLSIF" %in% class(x)) {
    cat("Regularization Parameter(lambda): ", x$lambda, "\n\n")
  }

  if("KLIEP" %in% class(x)) {
    cat("Number of the Folds: ", x$fold, "\n\n")
  }

  cat("The Function to Estimate Density Ratio:\n")
  cat("  compute_density_ratio()\n")

  cat("\n")
  invisible(x)
}
