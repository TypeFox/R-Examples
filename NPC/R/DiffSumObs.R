DiffSumObs <- function (y, tr, tl, ...) {
  ## Difference in number of non-missing observations
  nu0 <- sum(!is.na(y[tr != tl])) ## number of observed control responses
  nu1 <- sum(!is.na(y[tr == tl])) ## number of observed treated responses
  return(nu1 - nu0)
}
