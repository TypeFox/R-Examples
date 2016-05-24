StudentsT <- function (y, tr, tl, ...) {
  ## Student's T statistic
  d <- mean.default(y[tr==tl]) - mean.default(y[tr!=tl])
  se <- sqrt(stats::var(y[tr==tl]) / sum(tr==tl) +
                 stats::var(y[tr!=tl]) / sum(tr!=tl))
  return (d / se)
}
