library(robustvarComp)
load("test-dataset.Rdata")
do.test <- TRUE

if (do.test & Sys.info()['sysname']=="Linux") {
if (!file.exists('Init-savedvalues.R')) {  
  set.seed(2345)
  Init <- robustvarComp:::varComprob.initial(y=y, x=x, V=V, control=test.init.control)  
  Init$call <- NULL
  dput(Init, file='Init-savedvalues.R')
} else {
  set.seed(2345)
  Inittest <- robustvarComp:::varComprob.initial(y=y, x=x, V=V, control=test.init.control)  
  Inittest$call <- NULL
  Init <- dget(file='Init-savedvalues.R')
  stopifnot(
    all.equal(Inittest, Init, tol = 2e-5)
  )
}
}
