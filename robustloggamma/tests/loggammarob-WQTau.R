library(robustloggamma)

do.test <- TRUE
if (!do.test) {
  set.seed(1234)
  x <- sort(rloggamma(n=80, lambda=1))
  set.seed(2345)
  WQTau <- loggammarob(x, method="WQTau", control=loggammarob.control(method="WQTau", lower=0, upper=2, n=30, refine.tol=1e-8))
  WQTau$call <- NULL
  dput(WQTau, file='WQTau-savedvalues.R')
} else {
  set.seed(1234)
  x <- sort(rloggamma(n=80, lambda=1))
  set.seed(2345)
  WQTauTest <- loggammarob(x, method="WQTau", control=loggammarob.control(method="WQTau", lower=0, upper=2, n=30, refine.tol=1e-8))
  WQTauTest$call <- NULL  
  WQTau <- dget(file='WQTau-savedvalues.R')
  
  stopifnot(
    all.equal(WQTauTest, WQTau, tol = 2e-6)
  )
}
