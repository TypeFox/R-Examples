library(robustvarComp)
load("test-dataset.Rdata")
do.test <- FALSE

if (do.test) {
if (!file.exists('Tau-savedvalues.R') | !file.exists('sumTau-savedvalues.R')) {
  set.seed(2345)
  Tau <- robustvarComp:::varComprob.Tau(y=y, x=x, V=V, beta=beta, gamma=gamma, eta0=eta0, control=test.control)  
  Tau$call <- NULL
  dput(Tau, file='Tau-savedvalues.R')
  sumTau <- summary(Tau)
  sumTau$call <- NULL
  dput(sumTau, file='sumTau-savedvalues.R')
} else {
  set.seed(2345)
  Tautest <- robustvarComp:::varComprob.Tau(y=y, x=x, V=V, beta=beta, gamma=gamma, eta0=eta0, control=test.control)  
  Tautest$call <- NULL
  sumTautest <- summary(Tautest)
  sumTautest$call <- NULL

Tau <- dget(file='Tau-savedvalues.R')
sumTau <- dget(file='sumTau-savedvalues.R')
  
  stopifnot(
    all.equal(Tautest, Tau, tol = 2e-7)
  )
  stopifnot(
    all.equal(sumTautest, sumTau, tol = 2e-7)
  )
}
}
