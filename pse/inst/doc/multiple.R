### R code from vignette source 'multiple.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: seed
###################################################
set.seed(42)


###################################################
### code chunk number 2: model
###################################################
oneRun <- function (x1, x2, x3, x4)
    10 * x1 + 5 * x2 + 3 * rnorm(1, x3, x4)
modelRun <- function (my.data){
    mapply(oneRun, 
           my.data[,1], my.data[,2], my.data[,3], my.data[,4])
}


###################################################
### code chunk number 3: LHS1
###################################################
library(pse)
LHS1 <- LHS(modelRun, N=300, factors=4, nboot=50)
plotprcc(LHS1)


###################################################
### code chunk number 4: LHS2
###################################################
LHS2 <- LHS(modelRun, N=60, factors=4, repetitions=5, nboot=50)
plotprcc(LHS2)


###################################################
### code chunk number 5: cv
###################################################
r <- get.results(LHS2)
sd(r) / mean(r) # global CV
plotcv(LHS2)


###################################################
### code chunk number 6: LHS3
###################################################
newdata <- modelRun(get.data(LHS2))
LHS3 <- tell(LHS2, newdata)


