### R code from vignette source 'semiparametric-children.Rnw'

###################################################
### code chunk number 1: semiparametric-children.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: semiparametric-children.Rnw:19-21 (eval = FALSE)
###################################################
## library(catdata)
## data(children)


###################################################
### code chunk number 3: semiparametric-children.Rnw:26-27 (eval = FALSE)
###################################################
## library(mgcv)


###################################################
### code chunk number 4: semiparametric-children.Rnw:32-36 (eval = FALSE)
###################################################
## gamchild <- gam(child ~ s(age) + s(dur) + as.factor(nation) + as.factor(god) + 
##   as.factor(univ), data=children, family=poisson(link=log))
## 
## summary(gamchild)


###################################################
### code chunk number 5: semiparametric-children.Rnw:41-43 (eval = FALSE)
###################################################
## par(cex=1.5)
## plot(gamchild, select=1, ylab="", xlab="Age")


###################################################
### code chunk number 6: semiparametric-children.Rnw:46-48 (eval = FALSE)
###################################################
## par(cex=1.5)
## plot(gamchild, select=2, ylab="", xlab="Duration")


