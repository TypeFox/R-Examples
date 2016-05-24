### R code from vignette source 'mfp_vignette.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt=">", width=80)
set.seed(20040804)


###################################################
### code chunk number 2: mfp_vignette.Rnw:97-98
###################################################
library(mfp)


###################################################
### code chunk number 3: mfp_vignette.Rnw:108-109
###################################################
str(mfp)


###################################################
### code chunk number 4: mfp_vignette.Rnw:131-132
###################################################
str(fp)


###################################################
### code chunk number 5: mfp_vignette.Rnw:170-172
###################################################
data(GBSG)
str(GBSG)


###################################################
### code chunk number 6: mfp_vignette.Rnw:197-199
###################################################
f <- mfp(Surv(rfst, cens) ~ strata(htreat)+age+fp(tumsize)+fp(posnodal)+fp(prm)+fp(esm)
        +menostat+tumgrad, family = cox, data = GBSG, select=0.05, verbose=TRUE)


###################################################
### code chunk number 7: mfp_vignette.Rnw:209-210
###################################################
summary(f)


###################################################
### code chunk number 8: mfp_vignette.Rnw:215-216
###################################################
f$fptable


###################################################
### code chunk number 9: FIG1
###################################################
pf <- survfit(f$fit)  
plot(pf, col=c("red","green"), xlab="Time (years)", ylab="Recurrence free survival rate", xscale=365.25)


