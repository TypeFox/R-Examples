### R code from vignette source 'lmekin.Rnw'

###################################################
### code chunk number 1: lmekin.Rnw:20-21
###################################################
options(continue=' ', width=60)


###################################################
### code chunk number 2: lmekin.Rnw:70-78
###################################################
library(coxme)
require(nlme)
fit1 <- lme(effort~Type, random= ~ 1|Subject,data=ergoStool,
             method="ML")
fit2 <- lmekin(effort ~ Type + (1|Subject), data=ergoStool, 
                            method="ML")
print(fit1)
print(fit2)


###################################################
### code chunk number 3: lmekin.Rnw:84-95
###################################################
tdata <-eortc
tdata$center2 <- factor(tdata$center)

fit3 <- lme(y ~ trt, random= ~ trt|center2, data=tdata,
            method="ML")
fit3

fit4 <- lmekin(y ~ trt + (1+ trt|center), tdata)
fit4

all.equal(fit3$logLik, fit4$loglik)


###################################################
### code chunk number 4: gaw
###################################################
getOption("SweaveHooks")[["fig"]]()
require(kinship2)

load("gaw.rda")
gped <- with(gdata, pedigree(id, father, mother, sex=sex, famid=famid))
kmat <- kinship(gped)
plot(gped[9])

gfit0 <- lm(age ~ q1, gdata)
summary(gfit0)
gfit1 <- lmekin(age ~ q1 + (1|id), data=gdata, varlist=kmat*2)
gfit1


###################################################
### code chunk number 5: lmekin.Rnw:172-179
###################################################
sid <- pedindex[,9]
ibd6.90 <- with(solar6.90, sparseMatrix(id2, id1, x=x, symmetric=TRUE,
                                        dimnames=list(sid, sid)))

gfit2 <- lmekin(age ~ (1|id), data=gdata, 
                varlist= list(kmat, ibd6.90))
print(gfit2)


###################################################
### code chunk number 6: lmekin.Rnw:192-195
###################################################
gfit3 <- lmekin(age ~ q1 + (1|id) + (1|famid), data=gdata,
                varlist=kmat)
gfit3


###################################################
### code chunk number 7: lmekin.Rnw:199-201 (eval = FALSE)
###################################################
## lmekin(age ~ q1 + (1|id) + (1|famid), data=gdata,
##        varlist=list(coxmeMlist(kmat), coxmeFull))


