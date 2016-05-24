### R code from vignette source 'msm-manual.Rnw'

###################################################
### code chunk number 1: msm-manual.Rnw:23-25
###################################################
version <- gsub("Version: +", "",
               packageDescription("msm", lib.loc=c("../..",.libPaths()))$Version)


###################################################
### code chunk number 2: msm-manual.Rnw:30-31
###################################################
cat(version)


###################################################
### code chunk number 3: msm-manual.Rnw:34-35
###################################################
cat(format(Sys.time(), "%d %B, %Y"))


###################################################
### code chunk number 4: msm-manual.Rnw:843-844
###################################################
options(width = 60)


###################################################
### code chunk number 5: msm-manual.Rnw:879-880
###################################################
library(msm)


###################################################
### code chunk number 6: msm-manual.Rnw:939-940
###################################################
cav[1:21,]


###################################################
### code chunk number 7: msm-manual.Rnw:950-951
###################################################
statetable.msm(state, PTNUM, data=cav)


###################################################
### code chunk number 8: msm-manual.Rnw:1002-1006
###################################################
Q  <-  rbind ( c(0, 0.25, 0, 0.25),
               c(0.166, 0, 0.166, 0.166),
               c(0, 0.25, 0, 0.25),
               c(0, 0, 0, 0) )


###################################################
### code chunk number 9: msm-manual.Rnw:1040-1042
###################################################
Q.crude  <- crudeinits.msm(state ~ years, PTNUM, data=cav,
                                   qmatrix=Q)


###################################################
### code chunk number 10: msm-manual.Rnw:1066-1068
###################################################
cav.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                    qmatrix = Q, deathexact = 4)


###################################################
### code chunk number 11: msm-manual.Rnw:1096-1097 (eval = FALSE)
###################################################
## help(optim)


###################################################
### code chunk number 12: msm-manual.Rnw:1120-1121
###################################################
cav.msm


###################################################
### code chunk number 13: msm-manual.Rnw:1160-1162
###################################################
cavsex.msm <- msm( state ~ years, subject=PTNUM, data = cav,
                     qmatrix = Q, deathexact = 4, covariates = ~ sex)


###################################################
### code chunk number 14: msm-manual.Rnw:1170-1171
###################################################
cavsex.msm


###################################################
### code chunk number 15: msm-manual.Rnw:1187-1189
###################################################
qmatrix.msm(cavsex.msm, covariates=list(sex=0)) # Male
qmatrix.msm(cavsex.msm, covariates=list(sex=1)) # Female


###################################################
### code chunk number 16: msm-manual.Rnw:1201-1204 (eval = FALSE)
###################################################
## cavsex.msm <- msm( state ~ years, subject=PTNUM, data = cav,
##                   qmatrix = Q, deathexact = 4,
##                   covariates = list("1-2" = ~ sex, "1-4" = ~sex) )


###################################################
### code chunk number 17: msm-manual.Rnw:1215-1219 (eval = FALSE)
###################################################
## cav3.msm <- msm( state ~ years, subject=PTNUM, data = cav,
##                 qmatrix = Q, deathexact = 4,
##                 covariates = ~ sex,
##                 constraint = list(sex=c(1,2,3,1,2,3,2)) )


###################################################
### code chunk number 18: msm-manual.Rnw:1255-1259 (eval = FALSE)
###################################################
## cav4.msm <- msm( state ~ years, subject=PTNUM, data = cav,
##                 qmatrix = Q, deathexact = 4,
##                 control = list(trace=2, REPORT=1),
##                 fixedpars = c(6, 7) )


###################################################
### code chunk number 19: msm-manual.Rnw:1298-1299
###################################################
pmatrix.msm(cav.msm, t=10)


###################################################
### code chunk number 20: msm-manual.Rnw:1328-1329
###################################################
sojourn.msm(cav.msm)


###################################################
### code chunk number 21: msm-manual.Rnw:1341-1342
###################################################
pnext.msm(cav.msm)


###################################################
### code chunk number 22: msm-manual.Rnw:1367-1368
###################################################
totlos.msm(cav.msm)


###################################################
### code chunk number 23: msm-manual.Rnw:1390-1391
###################################################
qratio.msm(cav.msm, ind1=c(2,1), ind2=c(1,2))


###################################################
### code chunk number 24: msm-manual.Rnw:1399-1400
###################################################
hazard.msm(cavsex.msm)


###################################################
### code chunk number 25: msm-manual.Rnw:1409-1410 (eval = FALSE)
###################################################
## qmatrix.msm(cav.msm)


###################################################
### code chunk number 26: msm-manual.Rnw:1419-1420 (eval = FALSE)
###################################################
## qmatrix.msm(cavsex.msm, covariates = 0)


###################################################
### code chunk number 27: msm-manual.Rnw:1425-1426 (eval = FALSE)
###################################################
## qmatrix.msm(cavsex.msm, covariates = list(sex = 1))


###################################################
### code chunk number 28: msm-manual.Rnw:1452-1453
###################################################
plot(cav.msm, legend.pos=c(8, 1))


###################################################
### code chunk number 29: msm-manual.Rnw:1695-1697
###################################################
options(digits=3)
prevalence.msm(cav.msm, times=seq(0,20,2))


###################################################
### code chunk number 30: msm-manual.Rnw:1699-1700
###################################################
plot.prevalence.msm(cav.msm, mintime=0, maxtime=20)


###################################################
### code chunk number 31: msm-manual.Rnw:1834-1837
###################################################
options(digits=2)
pearson.msm(cav.msm, timegroups=2,
            transitions=c(1,2,3,4,5,6,7,8,9,9,9,10))


###################################################
### code chunk number 32: msm-manual.Rnw:1962-1974
###################################################
Qm <- rbind(c(0, 0.148, 0, 0.0171),
            c(0, 0, 0.202, 0.081),
            c(0, 0, 0, 0.126),
            c(0, 0, 0, 0))
ematrix <- rbind(c(0, 0.1, 0, 0),
                 c(0.1, 0, 0.1, 0),
                 c(0, 0.1, 0, 0),
                 c(0, 0, 0, 0))
cavmisc.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                   qmatrix = Qm, ematrix = ematrix, deathexact = 4,
                   obstrue = firstobs)
cavmisc.msm


###################################################
### code chunk number 33: msm-manual.Rnw:2002-2006
###################################################
cavmiscsex.msm <- msm(state ~ years, subject = PTNUM, data = cav,
                      qmatrix = Qm, ematrix = ematrix,
                      deathexact = 4, misccovariates = ~sex,
                      obstrue=firstobs)


###################################################
### code chunk number 34: msm-manual.Rnw:2008-2009
###################################################
cavmiscsex.msm


###################################################
### code chunk number 35: msm-manual.Rnw:2029-2031
###################################################
ematrix.msm(cavmiscsex.msm, covariates=list(sex=0))
ematrix.msm(cavmiscsex.msm, covariates=list(sex=1))


###################################################
### code chunk number 36: msm-manual.Rnw:2078-2080
###################################################
pearson.msm(cavmisc.msm, timegroups=2,
            transitions=c(1,2,3,4,5,6,7,8,9,9,9,10))


###################################################
### code chunk number 37: msm-manual.Rnw:2126-2128
###################################################
vit <- viterbi.msm(cavmisc.msm)
vit[vit$subject==100103,]


###################################################
### code chunk number 38: msm-manual.Rnw:2326-2327
###################################################
three.q <- rbind(c(0, exp(-6), exp(-9)), c(0, 0, exp(-6)), c(0, 0, 0))


###################################################
### code chunk number 39: msm-manual.Rnw:2345-2357
###################################################
hmodel1 <- list(hmmNorm(mean=100, sd=16), hmmNorm(mean=54, sd=18),
                hmmIdent(999))

fev1.msm <- msm(fev ~ days, subject=ptnum, data=fev, qmatrix=three.q,
                deathexact=3, hmodel=hmodel1,
                hcovariates=list(~acute, ~acute, NULL),
                hcovinits = list(-8, -8, NULL),
                hconstraint = list(acute = c(1,1)))

fev1.msm

sojourn.msm(fev1.msm)


###################################################
### code chunk number 40: msm-manual.Rnw:2624-2625 (eval = FALSE)
###################################################
## help(msm)


###################################################
### code chunk number 41: msm-manual.Rnw:2633-2634 (eval = FALSE)
###################################################
## help.start()


