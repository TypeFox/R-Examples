### R code from vignette source 'multgee.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library("gnm")
library("splines")
library("stats4")
library("VGAM")
options(prompt = "R> ", continue = "+   ")


###################################################
### code chunk number 2: multgee.Rnw:189-194
###################################################
library("multgee")
data("arthritis")
head(arthritis)
intrinsic.pars(y = y, data = arthritis, id = id, repeated = time,
                  rscale = "ordinal")


###################################################
### code chunk number 3: multgee.Rnw:201-205
###################################################
 fit <- ordLORgee(formula = y ~ factor(time) + factor(trt) + factor(baseline),
        link = "logit", id = id, repeated = time, data = arthritis,
        LORstr = "uniform")
 summary(fit)


###################################################
### code chunk number 4: multgee.Rnw:220-222
###################################################
fit1 <- update(fit, formula = ~. + factor(sex) + age)
waldts(fit, fit1)


