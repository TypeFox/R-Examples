### R code from vignette source 'ic.infer.rnw'

###################################################
### code chunk number 1: ic.infer.rnw:92-97
###################################################
options(width=80)
options(prompt = "R> ", continue = "+  ", digits = 4, useFancyQuotes = FALSE) 
require(ic.infer, quietly=TRUE)
contrasts(grades$HSR) <- "contr.treatment"
contrasts(grades$ACTC) <- "contr.treatment"


###################################################
### code chunk number 2: grades.unrestr
###################################################
limo.grades <- lm(meanGPA ~ HSR + ACTC, grades, weights = n)
summary(limo.grades)


###################################################
### code chunk number 3: bodyfat.unrestr
###################################################
limo.bodyfat <- lm(BodyFat ~ ., bodyfat)
summary(limo.bodyfat)


###################################################
### code chunk number 4: grades.unrestr.diff
###################################################
grades.diff <- grades
## change contrasts to contr.diff
contrasts(grades.diff$HSR) <- "contr.diff"
contrasts(grades.diff$ACTC) <- "contr.diff"
## display contrasts
contrasts(grades.diff$HSR)
limo.grades.diff <- lm(meanGPA ~ HSR + ACTC, grades.diff, weights = n)
summary(limo.grades.diff)


###################################################
### code chunk number 5: grades.contrasts
###################################################
## originally, treatment contrasts
ui.treat <- make.mon.ui(grades$HSR)
ui.treat
ui.diff <- make.mon.ui(grades.diff$HSR)
ui.diff


###################################################
### code chunk number 6: mon.means
###################################################
## originally, treatment contrasts
make.mon.ui(5, type = "mean")


###################################################
### code chunk number 7: ic.infer.rnw:227-228
###################################################
options(width=100)


###################################################
### code chunk number 8: grades.estHSR
###################################################
HSRmon <- ic.est(coef(limo.grades)[2:9],
       ui = ui.treat, 
       Sigma = vcov(limo.grades)[2:9, 2:9])
HSRmon


###################################################
### code chunk number 9: grades.estHSReq
###################################################
HSReq <- ic.est(coef(limo.grades)[2:9],
       ui = ui.treat, 
       Sigma = vcov(limo.grades)[2:9, 2:9], meq = 3)
HSReq


###################################################
### code chunk number 10: grades.summaryHSR
###################################################
summary(HSRmon)


###################################################
### code chunk number 11: grades.testHSR1
###################################################
summary(ic.test(HSReq), brief = FALSE)


###################################################
### code chunk number 12: grades.testHSR2
###################################################
summary(ic.test(HSReq, TP = 2))


###################################################
### code chunk number 13: grades.testHSR11
###################################################
HSReq.large <- ic.est(coef(limo.grades),
       ui = ui.treat, 
       Sigma = vcov(limo.grades), index = 2:9, meq = 3)
summary(ic.test(HSReq.large, TP = 11,
       ui0.11 = cbind(rep(0, 16), diag(1, 16))))


###################################################
### code chunk number 14: grades.testHSR21
###################################################
summary(ic.test(HSReq, TP = 21, meq.alt = 2))


###################################################
### code chunk number 15: grades.orlmHSR
###################################################
orlimo.grades <- orlm(limo.grades,
       ui = ui.treat, index = 2:9)
summary(orlimo.grades, brief = TRUE)


###################################################
### code chunk number 16: bodyfat.orlm
###################################################
orlimo.bodyfat <- orlm(limo.bodyfat,
       ui = diag(1,3), boot = TRUE)
summary(orlimo.bodyfat)


###################################################
### code chunk number 17: bodyfat.orrelimp
###################################################
or.relimp(limo.bodyfat, ui = diag(1, 3))


###################################################
### code chunk number 18: ic.infer.rnw:534-535
###################################################
require(relaimpo, quietly=TRUE)


###################################################
### code chunk number 19: bodyfat.calcrelimp
###################################################
calc.relimp(limo.bodyfat)$lmg


