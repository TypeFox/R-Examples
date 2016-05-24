DiffSumWithNA <- function (y, tr, tl, ...) {
  ## Difference of sums, accounting for missing data
  s0 <- sum(y[tr != tl], na.rm=TRUE) ## sum of observed control responses
  s1 <- sum(y[tr == tl], na.rm=TRUE) ## sum of observed treated responses
  nu0 <- sum(!is.na(y[tr != tl])) ## number of observed control responses
  nu1 <- sum(!is.na(y[tr == tl])) ## number of observed treated responses
  if (nu0 > 0 && nu1 > 0) {
    return(s1 * sqrt(nu0 / nu1) - s0 * sqrt(nu1 / nu0))
  } else {
    return(NA)
  }
}
