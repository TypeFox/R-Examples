RankSumWithNA <- function (y, tr, tl, ...) {
  ## Difference of rank-sums, accounting for missing data
  s0 <- sum(rank(y, na.last='keep')[tr != tl], na.rm=TRUE)
  s1 <- sum(rank(y, na.last='keep')[tr == tl], na.rm=TRUE)
  nu0 <- sum(!is.na(y[tr != tl])) ## number of observed control responses
  nu1 <- sum(!is.na(y[tr == tl])) ## number of observed treated responses
  if (nu0 > 0 && nu1 > 0) {
    return(s1 * sqrt(nu0 / nu1) - s0 * sqrt(nu1 / nu0))
  } else {
    return(NA)
  }
}
