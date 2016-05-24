library(robustvarComp)
load("test-dataset.Rdata")
do.test <- FALSE

if (do.test) {
S <- dget("S-savedvalues.R")
if (!file.exists('MM-savedvalues.R') | !file.exists('sumMM-savedvalues.R')) {  
  set.seed(2345)
  MM <- robustvarComp:::varComprob.MM(y=y, x=x, V=V, beta=S$beta, gamma=S$gamma, scale=S$scale, control=test.control)  
  MM$call <- NULL
  dput(MM, file='MM-savedvalues.R')
  sumMM <- summary(MM)
  sumMM$call <- NULL
  dput(sumMM, file='sumMM-savedvalues.R')
} else {
  set.seed(2345)
  MMtest <- robustvarComp:::varComprob.MM(y=y, x=x, V=V, beta=S$beta, gamma=S$gamma, scale=S$scale, control=test.control)  
  MMtest$call <- NULL
  sumMMtest <- summary(MMtest)
  sumMMtest$call <- NULL

MM <- dget(file='MM-savedvalues.R')
sumMM <- dget(file='sumMM-savedvalues.R')
  
  stopifnot(
    all.equal(MMtest, MM, tol = 2e-7)
  )
  stopifnot(
    all.equal(sumMMtest, sumMM, tol = 2e-7)
  )
}
}
