library(robustloggamma)

do.test <- TRUE
if (!do.test) {
  set.seed(1234)
  x <- sort(rloggamma(n=80, lambda=1))
  set.seed(2345)
  WL <- loggammarob(x, method="WL", control=loggammarob.control(method="WL", lower=0, upper=2, n=30, refine.tol=1e-8))
  WL$call <- NULL
  dput(WL, file='WL-savedvalues.R')
  sumWL <- summary(WL, p=c(0.1, 0.25, 0.5, 0.75, 0.9))
  sumWL$call <- NULL
  dput(sumWL, file='sumWL-savedvalues.R')
} else {
  set.seed(1234)
  x <- sort(rloggamma(n=80, lambda=1))
  set.seed(2345)
  WLTest <- loggammarob(x, method="WL", control=loggammarob.control(method="WL", lower=0, upper=2, n=30, refine.tol=1e-8))
  WLTest$call <- NULL  
  sumWLTest <- summary(WLTest, p=c(0.1, 0.25, 0.5, 0.75, 0.9))
  sumWLTest$call <- NULL  
  WL <- dget(file='WL-savedvalues.R')  
  sumWL <- dget(file='sumWL-savedvalues.R')
  stopifnot(
    all.equal(WLTest, WL, tol = 2e-6)
  )
  stopifnot(
    all.equal(sumWLTest, sumWL, tol = 2e-6)
  )
}
