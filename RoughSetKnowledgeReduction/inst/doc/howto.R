### R code from vignette source 'howto.Rnw'

###################################################
### code chunk number 1: howto.Rnw:70-71
###################################################
library(RoughSetKnowledgeReduction)


###################################################
### code chunk number 2: howto.Rnw:74-76
###################################################
rule <- c(1,1,1,1,0)
print(rule)


###################################################
### code chunk number 3: howto.Rnw:82-94
###################################################
exampleMatrix1 <- matrix(c(1,0,2,2,0,
0,1,1,1,2,
2,0,0,1,1,
1,1,0,2,2,
1,0,2,0,1,
2,2,0,1,1,
2,1,1,1,2,
0,1,1,0,1),ncol = 5,byrow=TRUE)
# Decision table creation
dt1 <- new(Class="DecisionTable",decisionTable = exampleMatrix1)
# Decision table creation alterative
dt1 <- decisionTable(exampleMatrix1)


###################################################
### code chunk number 4: howto.Rnw:101-104
###################################################
dt <- decisionTable(exampleMatrix1[,-5])
print(dt)
computeConsistencyMatrix(dt)


###################################################
### code chunk number 5: howto.Rnw:109-110
###################################################
checkConsistency(dt)


###################################################
### code chunk number 6: howto.Rnw:119-121
###################################################
print(dt1)
computeDiscernibilityMatrix(dt1)


###################################################
### code chunk number 7: howto.Rnw:133-135
###################################################
firstCR <- findFirstConditionReduct(dt1)
print(firstCR)


###################################################
### code chunk number 8: howto.Rnw:140-141
###################################################
smallestFamilyCR <- findSmallestReductFamilyFromCore(dt1)


###################################################
### code chunk number 9: howto.Rnw:146-147
###################################################
allCR <- findAllReductsFromCore(dt1)


###################################################
### code chunk number 10: howto.Rnw:155-157
###################################################
vr <- computeValueReduct(firstCR)
print(vr)


###################################################
### code chunk number 11: howto.Rnw:164-165
###################################################
computeSupportConsistency(vr,dt1)


###################################################
### code chunk number 12: howto.Rnw:172-174
###################################################
classDT <- classifyDecisionTable(vr,dt1)
print(classDT)


###################################################
### code chunk number 13: howto.Rnw:179-181
###################################################
#Cleaning
rm(exampleMatrix1,allCR,classDT,dt,dt1,firstCR,rule,smallestFamilyCR,vr)


###################################################
### code chunk number 14: howto.Rnw:191-206
###################################################
exampleMatrix1 <- matrix(c(1,0,0,1,1,
1,0,0,0,1,
0,0,0,0,0,
1,1,0,1,0,
1,1,0,2,2,
2,1,0,2,2,
2,2,2,2,2),ncol = 5,byrow=TRUE)
exampleMatrix2 <- matrix(c(1,0,2,2,0,
0,1,1,1,2,
2,0,0,1,1,
1,1,0,2,2,
1,0,2,0,1,
2,2,0,1,1,
2,1,1,1,2,
0,1,1,0,1),ncol = 5,byrow=TRUE)


###################################################
### code chunk number 15: howto.Rnw:211-213
###################################################
dt1 <- decisionTable(exampleMatrix1)
dt2 <- decisionTable(exampleMatrix2)


###################################################
### code chunk number 16: howto.Rnw:220-222
###################################################
cr1 <- findFirstConditionReduct(dt1)
vr1 <- computeValueReduct(cr1)


###################################################
### code chunk number 17: howto.Rnw:227-228
###################################################
vr1 <- removeDuplicatedRulesVR(vr1)


###################################################
### code chunk number 18: howto.Rnw:233-234
###################################################
computeSupportConsistency(vr1,dt1)


###################################################
### code chunk number 19: howto.Rnw:241-243
###################################################
dt3 <- classifyDecisionTable(vr1,dt2)
print(dt3)


###################################################
### code chunk number 20: howto.Rnw:249-251
###################################################
#Cleaning
rm(exampleMatrix1,exampleMatrix2,cr1,vr1,dt1,dt2,dt3)


###################################################
### code chunk number 21: howto.Rnw:258-266
###################################################
dtMatrix <- matrix(c(1,0,0,1,1,
1,0,0,0,1,
0,0,0,0,0,
1,1,0,1,0,
1,1,0,2,2,
2,1,0,2,2,
2,2,2,2,2),ncol = 5,byrow=TRUE)
dt <- decisionTable(dtMatrix)


###################################################
### code chunk number 22: howto.Rnw:270-272
###################################################
cr <- findFirstConditionReduct(dt)
print(cr)


###################################################
### code chunk number 23: howto.Rnw:276-278
###################################################
vr <- computeValueReduct(cr)
print(vr)


###################################################
### code chunk number 24: howto.Rnw:287-290
###################################################
vrMat <- getValueReduct(vr)# Matrix representation of the value reduct
vr1 <- valueReduct(cr,vrMat[c(1,4,5,6,7,9,10),])# Pick rules by matrix row index
vr2 <- valueReduct(cr,vrMat[c(1,3,5,6,7,9,12),])


###################################################
### code chunk number 25: howto.Rnw:295-296
###################################################
vr3 <- removeDuplicatedRulesVR(vr2)


###################################################
### code chunk number 26: howto.Rnw:301-305
###################################################
vr3Mat <-  getValueReduct(vr3)
rownames(vr3Mat) <- paste("R",1:nrow(vr3Mat),sep="")
vr4 <- valueReduct(cr,vr3Mat)
print(vr4)


###################################################
### code chunk number 27: howto.Rnw:310-312
###################################################
#Cleaning
rm(cr,dt,dtMatrix,vr,vr1,vr2,vr3,vr4,vrMat,vr3Mat)


