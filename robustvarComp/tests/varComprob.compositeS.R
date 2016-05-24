library(robustvarComp)
load("test-dataset.Rdata")
do.test <- TRUE

if (do.test) { ## & Sys.info()['sysname']=="Linux"
if (!file.exists('cS-savedvalues.R') | !file.exists('sumcS-savedvalues.R')) {  set.seed(2345)
  cS <- robustvarComp:::varComprob.compositeS(y=y, x=x, V=V, beta=beta, gamma=gamma, eta0=eta0, control=test.control)  
  cS$call <- NULL
  dput(cS, file='cS-savedvalues.R')
  sumcS <- summary(cS)
  sumcS$call <- NULL
  dput(sumcS, file='sumcS-savedvalues.R')
} else {
  set.seed(2345)
  cStest <- robustvarComp:::varComprob.compositeS(y=y, x=x, V=V, beta=beta, gamma=gamma, eta0=eta0, control=test.control)  
  cStest$call <- NULL
  sumcStest <- summary(cStest)
  sumcStest$call <- NULL

  cS <- dget(file='cS-savedvalues.R')
  sumcS <- dget(file='sumcS-savedvalues.R')

  cStest$vcov.beta <- cS$vcov.beta <- NULL 
  cStest$vcov.eta <- cS$vcov.eta <- NULL
  cStest$vcov.gamma <- cS$vcov.gamma <- NULL
  sumcStest$vcov.beta <- sumcS$vcov.beta <- NULL 
  sumcStest$vcov.eta <- sumcS$vcov.eta <- NULL
  sumcStest$vcov.gamma <- sumcS$vcov.gamma <- NULL
  sumcStest$zTable <- sumcS$zTable <- NULL
  
  stopifnot(
    all.equal(cStest, cS, tol = 2e-5)
  )
  stopifnot(
    all.equal(sumcStest, sumcS, tol = 2e-5)
  )
}
}
