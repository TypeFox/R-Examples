KS <- function (y, tr, tl, ...) {
  ## Kolmogorov-Smirnov statistic
  suppressWarnings(
      ks <- stats::ks.test(y[tr==tl], y[tr!=tl], exact=FALSE,
                           alternative="two.sided")
  )
  return(ks$statistic)
}
