StudentWilcoxon <- function (y, tr, tl, ...) {
  ## Studentized rank-sum
  m <- sum(tr!=tl)
  n <- sum(tr==tl)
  w <- sum(sapply(y[tr==tl], function (x) sum(x >= y[tr!=tl])))/(m*n) - 1/2
  z.tr <- sum((sapply(y[tr==tl], function (x) mean(y[tr!=tl] < x)) -
                 mean(sapply(y[tr==tl], function (x) y[tr!=tl] < x)))^2)/4
  z.nt <- sum((sapply(y[tr!=tl], function (x) mean(y[tr==tl] < x)) -
                 mean(sapply(y[tr!=tl], function (x) y[tr==tl] < x)))^2)/4
  se.w <- sqrt(z.nt / m + z.tr / n)
  w.stud <- w / se.w
  return (w.stud)
}
