library(robustvarComp)
load("test-dataset.Rdata")
do.test <- TRUE

if (do.test) { ## & Sys.info()['sysname']=="Linux"
if (!file.exists('cTau-savedvalues.R') | !file.exists('sumcTau-savedvalues.R')) {
  set.seed(2345)
  cTau <- robustvarComp:::varComprob.compositeTau(y=y, x=x, V=V, beta=beta, gamma=gamma, eta0=eta0, control=test.control)  
  cTau$call <- NULL
  dput(cTau, file='cTau-savedvalues.R')
  sumcTau <- summary(cTau)
  sumcTau$call <- NULL
  dput(sumcTau, file='sumcTau-savedvalues.R')
} else {
  set.seed(2345)
  cTautest <- robustvarComp:::varComprob.compositeTau(y=y, x=x, V=V, beta=beta, gamma=gamma, eta0=eta0, control=test.control)  
  cTautest$call <- NULL
  sumcTautest <- summary(cTautest)
  sumcTautest$call <- NULL

  cTau <- dget(file='cTau-savedvalues.R')
  sumcTau <- dget(file='sumcTau-savedvalues.R')

  cTautest$vcov.beta <- cTau$vcov.beta <- NULL 
  cTautest$vcov.eta <- cTau$vcov.eta <- NULL
  cTautest$vcov.gamma <- cTau$vcov.gamma <- NULL
  sumcTautest$vcov.beta <- sumcTau$vcov.beta <- NULL 
  sumcTautest$vcov.eta <- sumcTau$vcov.eta <- NULL
  sumcTautest$vcov.gamma <- sumcTau$vcov.gamma <- NULL
  sumcTautest$zTable <- sumcTau$zTable <- NULL
  
  stopifnot(
    all.equal(cTautest, cTau, tol = 2e-5)
  )
  stopifnot(
    all.equal(sumcTautest, sumcTau, tol = 2e-5)
  )
}
}
