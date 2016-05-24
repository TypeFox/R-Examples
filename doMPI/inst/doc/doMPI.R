### R code from vignette source 'doMPI.Rnw'

###################################################
### code chunk number 1: loadLibs
###################################################
library(foreach)
# examples in this document will be executed sequentially
registerDoSEQ()


###################################################
### code chunk number 2: doMPI.Rnw:147-148 (eval = FALSE)
###################################################
## install.packages("Rmpi")


###################################################
### code chunk number 3: doMPI.Rnw:156-157 (eval = FALSE)
###################################################
## install.packages("doMPI", dependencies=TRUE)


###################################################
### code chunk number 4: doMPI.Rnw:174-175 (eval = FALSE)
###################################################
## library(doMPI)


###################################################
### code chunk number 5: doMPI.Rnw:194-195 (eval = FALSE)
###################################################
## cl <- startMPIcluster(count=2)


###################################################
### code chunk number 6: doMPI.Rnw:211-212 (eval = FALSE)
###################################################
## registerDoMPI(cl)


###################################################
### code chunk number 7: doMPI.Rnw:227-228 (eval = FALSE)
###################################################
## closeCluster(cl)


###################################################
### code chunk number 8: doMPI.Rnw:248-249
###################################################
foreach(i=1:3) %dopar% sqrt(i)


###################################################
### code chunk number 9: doMPI.Rnw:271-275
###################################################
x <- foreach(i=1:3) %dopar% {
    sqrt(i)
}
x


###################################################
### code chunk number 10: doMPI.Rnw:295-299
###################################################
x <- foreach(i=1:3, .combine="c") %dopar% {
    sqrt(i)
}
x


###################################################
### code chunk number 11: doMPI.Rnw:307-312
###################################################
x <- foreach(seed=c(7, 11, 13), .combine="cbind") %dopar% {
    set.seed(seed)
    rnorm(3)
}
x


###################################################
### code chunk number 12: doMPI.Rnw:322-328
###################################################
x <- seq(-8, 8, by=0.5)
v <- foreach(y=x, .combine="cbind") %dopar% {
    r <- sqrt(x^2 + y^2) + .Machine$double.eps
    sin(r) / r
}
persp(x, x, v)


###################################################
### code chunk number 13: doMPI.Rnw:347-355
###################################################
x <- seq(-8, 8, by=0.5)
sinc <- function(y) {
    r <- sqrt(x^2 + y^2) + .Machine$double.eps
    sin(r) / r
}
r <- lapply(x, sinc)
v <- do.call("cbind", r)
persp(x, x, v)


###################################################
### code chunk number 14: doMPI.Rnw:386-396 (eval = FALSE)
###################################################
## ifiles <- list.files(pattern="\\.csv$")
## ofiles <- sub("\\.csv$", ".png", ifiles)
## foreach(i=ifiles, o=ofiles, .packages="randomForest") %dopar% {
##     d <- read.csv(i)
##     rf <- randomForest(Species~., data=d, proximity=TRUE)
##     png(filename=o)
##     MDSplot(rf, d$Species)
##     dev.off()
##     NULL
## }


