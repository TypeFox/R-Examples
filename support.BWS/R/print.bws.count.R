print.bws.count <- function(x, digits = max(3, getOption("digits") - 3), scientific = FALSE, ...)
{
# Name : print.bws.count
# Title: print() for S3 class "bws.count"
# Arguments:
#  x            an object of S3 class "bws.count"
#  digits       the number of significant digits
#  scientific   scores are encoded in scientific format


# display number of respondents

  cat("\nNumber of respondents = ", x$information$nrespondents, "\n", sep = "")


# display summary of disaggregated BW scores

  table <- data.frame(meanB       = colMeans(x$disaggregate$B),
                      meanW       = colMeans(x$disaggregate$W),
                      meanBW      = colMeans(x$disaggregate$BW),
                      mean.stdBW  = colMeans(x$disaggregate$stdBW),
                      stdev.stdBW = apply(x$disaggregate$stdBW, 2, sd))
  cat("\nSummary of disaggregated best-worst scores:\n")
  print(format(table,
               digits     = digits,
               scientific = scientific, ...), 
        quote = FALSE,
        right = TRUE)


# display aggregated BW scores

  cat("\nAggregated best-worst scores:\n")
  print(format(x$aggregate,
               digits     = digits,
               scientific = scientific, ...), 
        quote = FALSE, 
        right = TRUE)
  cat("\n")


  invisible(x)
}

