### R code from vignette source 'estimation.Rnw'

###################################################
### code chunk number 1: estimation.Rnw:37-43
###################################################
library (Frames2)
data(DatA)
data(DatB)

head (DatA, 3)
head (DatB, 3)


###################################################
### code chunk number 2: estimation.Rnw:54-62
###################################################
data(PiklA)
data(PiklB)

yA <- with(DatA, data.frame(Feed, Clo, Lei))
yB <- with(DatB, data.frame(Feed, Clo, Lei))

Hartley(yA, yB, PiklA, PiklB, DatA$Domain, DatB$Domain)
FB(yA, yB, PiklA, PiklB, DatA$Domain, DatB$Domain)


###################################################
### code chunk number 3: estimation.Rnw:68-70
###################################################
Hartley(yA, yB, DatA$ProbA, DatB$ProbB, DatA$Domain, DatB$Domain)
FB(yA, yB, DatA$ProbA, DatB$ProbB, DatA$Domain, DatB$Domain)


###################################################
### code chunk number 4: estimation.Rnw:75-77
###################################################
summary(Hartley(yA, yB, DatA$ProbA, DatB$ProbB, DatA$Domain, 
                DatB$Domain))


###################################################
### code chunk number 5: estimation.Rnw:82-86
###################################################
Hartley(yA, yB, DatA$ProbA, DatB$ProbB, DatA$Domain,
        DatB$Domain, 0.95)
FB(yA, yB, DatA$ProbA, DatB$ProbB, DatA$Domain, DatB$Domain, 
   0.95)


###################################################
### code chunk number 6: estimation.Rnw:91-93
###################################################
BKA(yA, yB, DatA$ProbA, DatB$ProbB, DatA$ProbB, DatB$ProbA,  
    DatA$Domain, DatB$Domain)


###################################################
### code chunk number 7: estimation.Rnw:97-98
###################################################
Compare(yA, yB, DatA$ProbA, DatB$ProbB, DatA$Domain, DatB$Domain)


###################################################
### code chunk number 8: estimation.Rnw:104-108
###################################################
SFRR(yA, yB, DatA$ProbA, DatB$ProbB, DatA$ProbB, DatB$ProbB, 
     DatA$Domain, DatB$Domain, N_A = 1735, N_B = 1191)
CalSF(yA, yB, DatA$ProbA, DatB$ProbB, DatA$ProbB, DatB$ProbB, 
      DatA$Domain, DatB$Domain, N_A = 1735, N_B = 1191)


###################################################
### code chunk number 9: estimation.Rnw:113-119
###################################################
PML(yA, yB, DatA$ProbA, DatB$ProbB, DatA$Domain, DatB$Domain, 
    N_A = 1735, N_B = 1191)
PEL(yA, yB, DatA$ProbA, DatB$ProbB, DatA$Domain, DatB$Domain, 
    N_A = 1735, N_B = 1191)
CalDF(yA, yB, DatA$ProbA, DatB$ProbB, DatA$Domain, DatB$Domain, 
      N_A = 1735, N_B = 1191)


###################################################
### code chunk number 10: estimation.Rnw:126-133
###################################################
PEL(yA, yB, DatA$ProbA, DatB$ProbB, DatA$Domain, DatB$Domain, 
    N_A = 1735, N_B = 1191, N_ab = 601)
CalSF(yA, yB, DatA$ProbA, DatB$ProbB, DatA$ProbB, DatB$ProbB, 
      DatA$Domain, DatB$Domain, N_A = 1735, N_B = 1191, 
      N_ab = 601)
CalDF(yA, yB, DatA$ProbA, DatB$ProbB, DatA$Domain, DatB$Domain, 
    N_A = 1735, N_B = 1191, N_ab = 601)


###################################################
### code chunk number 11: estimation.Rnw:141-151
###################################################
PEL(yA, yB, PiklA, PiklB, DatA$Domain, DatB$Domain, N_A = 1735, 
    N_B = 1191, xsAFrameA = DatA$Inc, xsBFrameA = DatB$Inc, 
    XA = 4300260)
CalSF(yA, yB, PiklA, PiklB, DatA$ProbB, DatB$ProbA, DatA$Domain, 
      DatB$Domain, N_A = 1735, N_B = 1191, xsAFrameB = DatA$M2, 
      xsBFrameB = DatB$M2, XB = 176553)
CalDF(yA, yB, PiklA, PiklB, DatA$Domain, DatB$Domain, N_A = 1735, 
      N_B = 1191, xsAFrameA = DatA$Inc, xsBFrameA = DatB$Inc, 
      xsAFrameB = DatA$M2, xsBFrameB = DatB$M2, XA = 4300260, 
      XB = 176553)


