### R code from vignette source 'QCA3.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library(QCA3)
options(SweaveHooks=list(twofig=function() {par(mfrow=c(1,2))},
                         twofig2=function() {par(mfrow=c(2,1))},
                         onefig=function() {par(mfrow=c(1,1))}))


###################################################
### code chunk number 2: QCA3.Rnw:40-41
###################################################
library(QCA3)


###################################################
### code chunk number 3: QCA3.Rnw:57-60
###################################################
(cst <- cs_truthTable(Lipset_cs,outcome="SURVIVAL",
                      condition=c("GNPCAP", "URBANIZA", "LITERACY", "INDLAB", "GOVSTAB"),
                      cases="CASEID"))


###################################################
### code chunk number 4: QCA3.Rnw:91-92
###################################################
reduce(cst)


###################################################
### code chunk number 5: QCA3.Rnw:99-100
###################################################
reduce(cst, explain="negative")


###################################################
### code chunk number 6: QCA3.Rnw:108-109
###################################################
reduce(cst, remainders="include")


###################################################
### code chunk number 7: QCA3.Rnw:115-116
###################################################
reduce(cst, explain="negative", remainders="include")


###################################################
### code chunk number 8: QCA3.Rnw:127-129
###################################################
ansNeg <- reduce(cst, explain="negative", remainders="include")
SA(ansNeg)


###################################################
### code chunk number 9: QCA3.Rnw:159-160
###################################################
directCalibration(Lipset_fs$Developed,fullin=900,fullout=400, crossover=550)


###################################################
### code chunk number 10: QCA3.Rnw:166-167 (eval = FALSE)
###################################################
## Lipset_fs$DFZ <- directCalibration(Lipset_fs$Developed,900,400, 550)


###################################################
### code chunk number 11: QCA3.Rnw:177-179
###################################################
suffnec(Lipset_fs[,c("Survived.FZ","Developed.FZ","Urban.FZ",
                     "Literate.FZ","Industrial.FZ", "Stable.FZ")])


###################################################
### code chunk number 12: QCA3.Rnw:190-191
###################################################
fsplot(Survived.FZ~fsand(Developed.FZ, Urban.FZ),data=Lipset_fs)


###################################################
### code chunk number 13: QCA3.Rnw:211-214
###################################################
conditions <- c("Developed.FZ","Urban.FZ","Literate.FZ","Industrial.FZ", "Stable.FZ")
fst <- fs_truthTable(Lipset_fs,"Survived.FZ", conditions, consistency=0.7)
print(fst)


###################################################
### code chunk number 14: QCA3.Rnw:221-223
###################################################
fsans <- reduce(fst)
print(fsans)


###################################################
### code chunk number 15: QCA3.Rnw:228-229
###################################################
summary(fsans)


###################################################
### code chunk number 16: QCA3.Rnw:234-235
###################################################
update(fsans, remainders="include")


###################################################
### code chunk number 17: QCA3.Rnw:242-247
###################################################
Lipset_fs$Not.Survived <- fsnot(Lipset_fs$Survived.FZ)
fst2 <- fs_truthTable(Lipset_fs,"Not.Survived", conditions, consistency=0.7)
print(fst2)
fsans2 <- reduce(fst2)
summary(fsans2)


###################################################
### code chunk number 18: QCA3.Rnw:251-252
###################################################
sessionInfo()


