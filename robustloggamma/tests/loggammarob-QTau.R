library(robustloggamma)

do.test <- TRUE
if (!do.test) {
  set.seed(1234)
  x <- sort(rloggamma(n=80, lambda=1))
  set.seed(2345)
  QTau <- loggammarob(x, method="QTau", control=loggammarob.control(method="QTau", lower=0, upper=2, n=30, refine.tol=1e-8))
  QTau$call <- NULL
  dput(QTau, file='QTau-savedvalues.R')
} else {
  set.seed(1234)
  x <- sort(rloggamma(n=80, lambda=1))
  set.seed(2345)
  QTauTest <- loggammarob(x, method="QTau", control=loggammarob.control(method="QTau", lower=0, upper=2, n=30, refine.tol=1e-8))
  QTauTest$call <- NULL  
  QTau <- dget(file='QTau-savedvalues.R')
  stopifnot(
    all.equal(QTauTest, QTau, tol = 2e-6)
  )
}
