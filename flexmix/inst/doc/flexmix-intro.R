### R code from vignette source 'flexmix-intro.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: flexmix-intro.Rnw:31-36
###################################################
set.seed(1504)
options(width=70, prompt = "R> ", continue = "+  ", useFancyQuotes = FALSE)
ps.options(family="Times")
library("flexmix")
data("NPreg")


###################################################
### code chunk number 2: flexmix-intro.Rnw:322-326
###################################################
library("flexmix")
data("NPreg")
m1 <- flexmix(yn ~ x + I(x^2), data = NPreg, k = 2)
m1


###################################################
### code chunk number 3: flexmix-intro.Rnw:329-330
###################################################
parameters(m1, component = 1)


###################################################
### code chunk number 4: flexmix-intro.Rnw:333-334
###################################################
parameters(m1, component = 2)


###################################################
### code chunk number 5: flexmix-intro.Rnw:339-340
###################################################
table(NPreg$class, clusters(m1))


###################################################
### code chunk number 6: flexmix-intro.Rnw:343-344
###################################################
summary(m1)


###################################################
### code chunk number 7: flexmix-intro.Rnw:360-363
###################################################
par(mfrow=c(1,2))
plot(yn~x, col=class, pch=class, data=NPreg)
plot(yp~x, col=class, pch=class, data=NPreg)


###################################################
### code chunk number 8: flexmix-intro.Rnw:381-382
###################################################
print(plot(m1))


###################################################
### code chunk number 9: flexmix-intro.Rnw:402-404
###################################################
rm1 <- refit(m1)
summary(rm1)


###################################################
### code chunk number 10: flexmix-intro.Rnw:425-426
###################################################
options(width=55)


###################################################
### code chunk number 11: flexmix-intro.Rnw:428-431
###################################################
m2 <- flexmix(yp ~ x, data = NPreg, k = 2, 
  model = FLXMRglm(family = "poisson"))
summary(m2)


###################################################
### code chunk number 12: flexmix-intro.Rnw:433-434
###################################################
options(width=65)


###################################################
### code chunk number 13: flexmix-intro.Rnw:438-439
###################################################
print(plot(m2))


###################################################
### code chunk number 14: flexmix-intro.Rnw:482-485
###################################################
m3 <- flexmix(~ x, data = NPreg, k = 2,
  model=list(FLXMRglm(yn ~ . + I(x^2)), 
    FLXMRglm(yp ~ ., family = "poisson")))


###################################################
### code chunk number 15: flexmix-intro.Rnw:500-501
###################################################
print(plot(m3))


###################################################
### code chunk number 16: flexmix-intro.Rnw:530-532
###################################################
m4 <- flexmix(yn ~ x + I(x^2) | id2, data = NPreg, k = 2)
summary(m4)


###################################################
### code chunk number 17: flexmix-intro.Rnw:548-550
###################################################
m5 <- flexmix(yn ~ x + I(x^2), data = NPreg, k = 2,
  control = list(iter.max = 15, verbose = 3, classify = "hard"))


###################################################
### code chunk number 18: flexmix-intro.Rnw:567-571
###################################################
m6 <- flexmix(yp ~ x + I(x^2), data = NPreg, k = 4,
  control = list(minprior = 0.2))

m6  


###################################################
### code chunk number 19: flexmix-intro.Rnw:581-584
###################################################
m7 <- stepFlexmix(yp ~ x + I(x^2), data = NPreg,
  control = list(verbose = 0), k = 1:5, nrep = 5)



###################################################
### code chunk number 20: flexmix-intro.Rnw:590-591
###################################################
getModel(m7, "BIC")


###################################################
### code chunk number 21: flexmix-intro.Rnw:726-733
###################################################
library("flexmix")
set.seed(1504)
options(width=60)
ps.options(family="Times")
require("ellipse")
require("mvtnorm")
source("mymclust.R")


###################################################
### code chunk number 22: flexmix-intro.Rnw:739-742
###################################################
data("Nclus")
m1 <- flexmix(Nclus ~ 1, k = 4, model = mymclust())
summary(m1)


###################################################
### code chunk number 23: flexmix-intro.Rnw:753-755
###################################################
m2 <- flexmix(Nclus ~ 1, k = 4, model = mymclust(diagonal = FALSE))
summary(m2)


###################################################
### code chunk number 24: flexmix-intro.Rnw:760-763
###################################################
par(mfrow=1:2)
plotEll(m1, Nclus)
plotEll(m2, Nclus)


###################################################
### code chunk number 25: flexmix-intro.Rnw:802-806
###################################################
SI <- sessionInfo()
pkgs <- paste(sapply(c(SI$otherPkgs, SI$loadedOnly), function(x) 
                     paste("\\\\pkg{", x$Package, "} ", 
                           x$Version, sep = "")), collapse = ", ")


