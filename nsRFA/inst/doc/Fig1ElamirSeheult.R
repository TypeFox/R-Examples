### R code from vignette source 'Fig1ElamirSeheult.Rnw'

###################################################
### code chunk number 1: Fig1ElamirSeheult.Rnw:40-41
###################################################
library(nsRFA)


###################################################
### code chunk number 2: Fig1ElamirSeheult.Rnw:44-46
###################################################
Nsim=1000
n=60


###################################################
### code chunk number 3: Fig1ElamirSeheult.Rnw:48-49
###################################################
campsimulati <- rnorm(n*Nsim)


###################################################
### code chunk number 4: a
###################################################
campsimulati <- matrix(campsimulati, ncol=n)


###################################################
### code chunk number 5: b
###################################################
lmom <- t(apply(campsimulati, 1, Lmoments))
vlmom <- t(apply(campsimulati, 1, varLmoments, matrix=FALSE))

l3 <- lmom[,"lca"]*lmom[,"l2"]
sl3 <- sqrt(vlmom[,"var.l3"])


###################################################
### code chunk number 6: Fig1ElamirSeheult.Rnw:62-63
###################################################
l3gaussian <- l3/sl3


###################################################
### code chunk number 7: Fig1ElamirSeheult.Rnw:66-68
###################################################
qqnorm(l3gaussian, main="Normal Q-Q Plot for Gaussian samples")
 qqline(l3gaussian)


###################################################
### code chunk number 8: Fig1ElamirSeheult.Rnw:72-73
###################################################
campsimulati <- rt(n*Nsim, df=5)


###################################################
### code chunk number 9: Fig1ElamirSeheult.Rnw:75-77
###################################################
campsimulati <- matrix(campsimulati, ncol=n)
lmom <- t(apply(campsimulati, 1, Lmoments))
vlmom <- t(apply(campsimulati, 1, varLmoments, matrix=FALSE))

l3 <- lmom[,"lca"]*lmom[,"l2"]
sl3 <- sqrt(vlmom[,"var.l3"])


###################################################
### code chunk number 10: Fig1ElamirSeheult.Rnw:79-80
###################################################
l3student <- l3/sl3


###################################################
### code chunk number 11: Fig1ElamirSeheult.Rnw:83-84
###################################################
campsimulati <- rcauchy(n*Nsim)


###################################################
### code chunk number 12: Fig1ElamirSeheult.Rnw:86-88
###################################################
campsimulati <- matrix(campsimulati, ncol=n)
lmom <- t(apply(campsimulati, 1, Lmoments))
vlmom <- t(apply(campsimulati, 1, varLmoments, matrix=FALSE))

l3 <- lmom[,"lca"]*lmom[,"l2"]
sl3 <- sqrt(vlmom[,"var.l3"])


###################################################
### code chunk number 13: Fig1ElamirSeheult.Rnw:90-91
###################################################
l3cauchy <- l3/sl3


###################################################
### code chunk number 14: Fig1ElamirSeheult.Rnw:94-95
###################################################
campsimulati <- runif(n*Nsim)


###################################################
### code chunk number 15: Fig1ElamirSeheult.Rnw:97-99
###################################################
campsimulati <- matrix(campsimulati, ncol=n)
lmom <- t(apply(campsimulati, 1, Lmoments))
vlmom <- t(apply(campsimulati, 1, varLmoments, matrix=FALSE))

l3 <- lmom[,"lca"]*lmom[,"l2"]
sl3 <- sqrt(vlmom[,"var.l3"])


###################################################
### code chunk number 16: Fig1ElamirSeheult.Rnw:101-102
###################################################
l3unif <- l3/sl3


###################################################
### code chunk number 17: Fig1ElamirSeheult.Rnw:106-108
###################################################
#bitmap(file="Fig1.png", type="png256", height=10, width=8, res=144, pointsize=16)
png(filename="Fig1.png", height=720, width=600, res=72, pointsize=16)


###################################################
### code chunk number 18: Fig1ElamirSeheult.Rnw:110-119
###################################################
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
qqnorm(l3gaussian, main="Normal Plot: Gaussian samples")
 qqline(l3gaussian)
qqnorm(l3student, main="Normal Plot: Student (df=5) samples")
 qqline(l3student)
qqnorm(l3cauchy, main="Normal Plot: Cauchy samples")
 qqline(l3cauchy)
qqnorm(l3unif, main="Normal Plot: Uniform samples")
 qqline(l3unif)


###################################################
### code chunk number 19: Fig1ElamirSeheult.Rnw:121-122
###################################################
dev.off()


