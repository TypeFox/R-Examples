`yuleWindow` <-
function(x, st1 = x[1], st2 = 0)
{
  checkbasal(x)
  if (!is.numeric(x)) stop("object x not of class 'numeric'")
  if (st1 <= st2)
  {
    cat("Error: invalid input\n\n")
    cat("shift points must be in temporal order (first, second)\n")
    cat("and in units of time/substitutions before present\n\n")
    return()
  }
  x <- rev(sort(x))
  res <- as.list(yuleint2(x, st1, st2))
  if (is.na(res$LH))
  {
    print("Error.  Must be >= 1 speciation event within temporal window\n")
    stop()
  }
  cat("---------------------------------\n")
  cat("Results of fitting Yule model to temporal window of \n")
  cat("(", st1, ",", st2, ") time units before present:\n\n")
  cat("ML speciation rate: ", res$smax, "\n\nLikelihood: ", res$LH, "\n\n")
  return(res)
}

