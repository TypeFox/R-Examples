library(robustvarComp)
load("test-dataset.Rdata")
do.test <- TRUE

if (do.test) { ## & Sys.info()['sysname']=="Linux"
if (!file.exists('S-savedvalues.R') | !file.exists('sumS-savedvalues.R')) {
  set.seed(2345)
  S <- robustvarComp:::varComprob.S(y=y, x=x, V=V, beta=beta, gamma=gamma, eta0=eta0, control=test.control)  
  S$call <- NULL
  dput(S, file='S-savedvalues.R')
  sumS <- summary(S)
  sumS$call <- NULL
  dput(sumS, file='sumS-savedvalues.R')
} else {
  set.seed(2345)
  Stest <- robustvarComp:::varComprob.S(y=y, x=x, V=V, beta=beta, gamma=gamma, eta0=eta0, control=test.control)  
  Stest$call <- NULL
  sumStest <- summary(Stest)
  sumStest$call <- NULL

  S <- dget(file='S-savedvalues.R')
  sumS <- dget(file='sumS-savedvalues.R')

  Stest$vcov.beta <- S$vcov.beta <- NULL 
  Stest$vcov.eta <- S$vcov.eta <- NULL
  Stest$vcov.gamma <- S$vcov.gamma <- NULL
  sumStest$vcov.beta <- sumS$vcov.beta <- NULL 
  sumStest$vcov.eta <- sumS$vcov.eta <- NULL
  sumStest$vcov.gamma <- sumS$vcov.gamma <- NULL
  sumStest$zTable <- sumS$zTable <- NULL
  
  stopifnot(
    all.equal(Stest, S, tol = 2e-5)
  )
  stopifnot(
    all.equal(sumStest, sumS, tol = 2e-5)
  )
}
}
