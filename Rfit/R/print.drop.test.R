print.drop.test <- function (x, digits = max(5, .Options$digits - 2), ...) {
  cat("\nDrop in Dispersion Test\n")
  y <- c(x$F, x$p.value)
  names(y) <- c("F-Statistic", "p-value")
  print(format(y, digits = digits), quote = FALSE)
}
