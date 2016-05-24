library(robustvarComp)
load("test-dataset.Rdata")

do.test <- FALSE

if (do.test) {
cS <- dget("cS-savedvalues.R")
if (!file.exists('cMM-savedvalues.R') | !file.exists('sumcMM-savedvalues.R')) {  set.seed(2345)
  MM <- robustvarComp:::varComprob.compositeMM(y=y, x=x, V=V, beta=cS$beta, gamma=cS$gamma, scales=cS$scales, control=test.control)  
  MM$call <- NULL
  dput(MM, file='cMM-savedvalues.R')
  sumMM <- summary(MM)
  sumMM$call <- NULL
  dput(sumMM, file='sumcMM-savedvalues.R')
} else {
  set.seed(2345)
  MMtest <- robustvarComp:::varComprob.compositeMM(y=y, x=x, V=V, beta=cS$beta, gamma=cS$gamma, scales=cS$scales, control=test.control)  
  MMtest$call <- NULL
  sumMMtest <- summary(MMtest)
  sumMMtest$call <- NULL

cMM <- dget(file='cMM-savedvalues.R')
sumcMM <- dget(file='sumcMM-savedvalues.R')
  
  stopifnot(
    all.equal(MMtest, MM, tol = 2e-7)
  )
  stopifnot(
    all.equal(sumMMtest, sumMM, tol = 2e-7)
  )
}
}
