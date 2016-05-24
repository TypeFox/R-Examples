
print.tsoutliers <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  print(x$fit)
  if (nrow(x$outliers) > 0)
  {
    cat("\nOutliers:\n")
    print(format(x$outliers, digits = digits))
  } else
    cat("\nNo outliers were detected.\n")
}
