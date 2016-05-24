library(robustloggamma)

do.test <- TRUE
if (!do.test) {
  set.seed(1234)
  x <- sort(rloggamma(n=80, lambda=1))
  set.seed(2345)
  oneWL <- loggammarob(x, method="oneWL", control=loggammarob.control(method="oneWL", lower=0, upper=2, n=30, refine.tol=1e-8))
  oneWL$call <- NULL
  dput(oneWL, file='oneWL-savedvalues.R')
  sumoneWL <- summary(oneWL, p=c(0.1, 0.25, 0.5, 0.75, 0.9))
  sumoneWL$call <- NULL
  dput(sumoneWL, file='sumoneWL-savedvalues.R')
} else {
  set.seed(1234)
  x <- sort(rloggamma(n=80, lambda=1))
  set.seed(2345)
  oneWLTest <- loggammarob(x, method="oneWL", control=loggammarob.control(method="oneWL", lower=0, upper=2, n=30, refine.tol=1e-8))
  oneWLTest$call <- NULL  
  sumoneWLTest <- summary(oneWLTest, p=c(0.1, 0.25, 0.5, 0.75, 0.9))
  sumoneWLTest$call <- NULL  
  oneWL <- dget(file='oneWL-savedvalues.R')  
  sumoneWL <- dget(file='sumoneWL-savedvalues.R')
  
  stopifnot(
    all.equal(oneWLTest, oneWL, tol = 2e-6)
  )
  stopifnot(
    all.equal(sumoneWLTest, sumoneWL, tol = 2e-6)
  )
}
