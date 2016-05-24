rrange <- function(x, range = 1, coef = 1.5, na.rm = TRUE)
{
  ## Purpose: `Robust RANGE', ===>  ?rrange
  ## Author: Martin Maechler, 1990
  if(!missing(range)) {
    if(!missing(coef)) stop("Must use either 'range' or 'coef'")
    coef <- 1.5 * range
  }
  if(!na.rm && any(is.na(x)))
      return(0+ c(NA,NA))# numeric NA

  ## S: (boxplot(..., plot = FALSE)$stats)[c(5, 1)]
  boxplot.stats(x, coef = coef, do.conf= FALSE, do.out= FALSE)$stats[c(1,5)]
}
