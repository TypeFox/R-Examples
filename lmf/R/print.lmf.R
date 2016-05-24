print.lmf <-
function(x,
                      digits = max(3, getOption("digits") - 2),
                      ...)
{
  cat("\nESTIMATING FLUCTUATING SELECTION IN AGE-STRUCTURED POPULATIONS\n",
      sep = "")
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  if (length(x$aM))
  {
    cat("Temporal alpha estimates (a(M)):\n")
    print.default(format(x$aM, digits = digits), print.gap = 2,
                  quote = FALSE)
  }
  else cat("No temporal alpha estimates (a(M))\n")
  cat("\n")
  if (length(x$M) > 1)
  {
    cat("Temporal covariance matrix (M):\n")
    print.default(format(x$M, digits = digits), print.gap = 2,
                  quote = FALSE)
  }
  else cat("No temporal covariance matrix (M)\n")
  cat("\n")
  cat("Maximization diagnostics:\n")
  print.default(format(c(convergence = x$convergence, iterations = x$iterations), digits = digits), print.gap = 2,
                quote = FALSE)
  cat("\n")
  cat("Running time:\n")
  print.default(format(c(total = x$running.time, optim = x$optim.time), digits = digits), print.gap = 2,
                quote = FALSE)
  cat("-End-\n")
  invisible(x)
}
