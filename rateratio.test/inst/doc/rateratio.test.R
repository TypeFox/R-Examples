### R code from vignette source 'rateratio.test.Rnw'

###################################################
### code chunk number 1: rateratio.test.Rnw:23-26
###################################################
library(rateratio.test)
n<-17877
m<-16660


###################################################
### code chunk number 2: rateratio.test.Rnw:38-39
###################################################
rateratio.test(c(2,9),c(n,m))


###################################################
### code chunk number 3: rateratio.test.Rnw:128-134
###################################################
n<-17877
m<-16674
rateratio.test(c(2,9),c(n,m))$conf.int
b.ci<-binom.test(2,2+9,p=n/(n+m))$conf.int
theta.ci<-m*b.ci/(n*(1-b.ci))
theta.ci


###################################################
### code chunk number 4: rateratio.test.Rnw:137-140
###################################################
R.Version()$version.string
rateratio.test(c(2,9),c(n,m))
binom.test(2,2+9,p=n/(n+m))


###################################################
### code chunk number 5: rateratio.test.Rnw:146-147
###################################################
fisher.test(matrix(c(2,9,n-2,m-9),2,2))


