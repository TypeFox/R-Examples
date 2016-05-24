### R code from vignette source 'gettingstartedMC.Rnw'

###################################################
### code chunk number 1: loadLibs
###################################################
library(doMC)
registerDoMC(2)
foreach(i=1:3) %dopar% sqrt(i)


###################################################
### code chunk number 2: bootpar
###################################################
x <- iris[which(iris[,5] != "setosa"), c(1,5)]
trials <- 10000

ptime <- system.time({
  r <- foreach(icount(trials), .combine=cbind) %dopar% {
    ind <- sample(100, 100, replace=TRUE)
    result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
    coefficients(result1)
  }
})[3]
ptime


###################################################
### code chunk number 3: bootseq
###################################################
stime <- system.time({
  r <- foreach(icount(trials), .combine=cbind) %do% {
    ind <- sample(100, 100, replace=TRUE)
    result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
    coefficients(result1)
  }
})[3]
stime


###################################################
### code chunk number 4: getDoParWorkers
###################################################
getDoParWorkers()


###################################################
### code chunk number 5: getDoParName
###################################################
getDoParName()
getDoParVersion()


###################################################
### code chunk number 6: options
###################################################
mcoptions <- list(preschedule=FALSE, set.seed=FALSE)
foreach(i=1:3, .options.multicore=mcoptions) %dopar% sqrt(i)


###################################################
### code chunk number 7: coreoptions2
###################################################
registerDoMC(2)
options(cores=4)
getDoParWorkers()


