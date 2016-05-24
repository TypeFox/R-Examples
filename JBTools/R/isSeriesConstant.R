isSeriesConstant <- function(
  ##title<< Check whether a vector is constant
  x                     ##<< numeric vector: series to test.
  , tresh.const = 1e-12 ##<< numeric: maximum deviation allowed which is still
                        ##   considered to be constant.
  , ratio.const = 0.05  ##<< numeric: ratio of the series which is allowed to be
                        ##   not constant for the whole series to be still
                        ##   considered to be constant.
  )
  ##description<<
  ## isSeriesConstant checks whether a series is constant (up to a certain degree).
  ##details<<
  ## isSeriesConstant checks whether the amount of values deviating from the
  ## median of x for a higher value than tresh.const is bigger than ratio.const.

{
  if (sum(is.na(x)) == length(x)) {
    return(FALSE)
  } else {
    min.amount <-  (1 - ratio.const)*length(x[!is.na(x)])
    is.constant <- sum(abs(x - median(x, na.rm = TRUE)) <  tresh.const, na.rm = TRUE) >= min.amount
    ##value<< logical: TRUE if series is constant, FALSE otherwise.
    return(is.constant)
  }
}
