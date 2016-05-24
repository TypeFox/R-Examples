### R code from vignette source 'glpkAPI.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: glpkAPI.Rnw:95-96
###################################################
library(glpkAPI)


###################################################
### code chunk number 2: glpkAPI.Rnw:99-100
###################################################
prob <- initProbGLPK()


###################################################
### code chunk number 3: glpkAPI.Rnw:103-104
###################################################
setProbNameGLPK(prob, "sample")


###################################################
### code chunk number 4: glpkAPI.Rnw:109-110
###################################################
setObjDirGLPK(prob, GLP_MAX)


###################################################
### code chunk number 5: glpkAPI.Rnw:113-115
###################################################
addRowsGLPK(prob, 3)
addColsGLPK(prob, 3)


###################################################
### code chunk number 6: glpkAPI.Rnw:118-124
###################################################
setRowNameGLPK(prob, 1, "p")
setRowNameGLPK(prob, 2, "q")
setRowNameGLPK(prob, 3, "r")
setColNameGLPK(prob, 1, "x1")
setColNameGLPK(prob, 2, "x2")
setColNameGLPK(prob, 3, "x3")


###################################################
### code chunk number 7: glpkAPI.Rnw:127-130
###################################################
setRowBndGLPK(prob, 1, GLP_UP, 0, 100)
setRowBndGLPK(prob, 2, GLP_UP, 0, 600)
setRowBndGLPK(prob, 3, GLP_UP, 0, 300)


###################################################
### code chunk number 8: glpkAPI.Rnw:134-139
###################################################
lb <- c(0, 0, 0)
ub <- c(100, 600, 300)
type <- rep(GLP_UP, 3)

setRowsBndsGLPK(prob, 1:3, lb, ub, type)


###################################################
### code chunk number 9: glpkAPI.Rnw:142-145
###################################################
setColBndGLPK(prob, 1, GLP_LO, 0, 0)
setColBndGLPK(prob, 2, GLP_LO, 0, 0)
setColBndGLPK(prob, 3, GLP_LO, 0, 0)


###################################################
### code chunk number 10: glpkAPI.Rnw:148-151
###################################################
setObjCoefGLPK(prob, 1, 10)
setObjCoefGLPK(prob, 2, 6)
setObjCoefGLPK(prob, 3, 4)


###################################################
### code chunk number 11: glpkAPI.Rnw:155-161
###################################################
lb <- c(0, 0, 0)
ub <- lb
type <- rep(GLP_LO, 3)
obj <- c(10, 6, 4)

setColsBndsObjCoefsGLPK(prob, 1:3, lb, ub, obj, type)


###################################################
### code chunk number 12: glpkAPI.Rnw:164-169
###################################################
ia <- c(1, 1, 1, 2, 3, 2, 3, 2, 3)
ja <- c(1, 2, 3, 1, 1, 2, 2, 3, 3)
ar <- c(1, 1, 1, 10, 2, 4, 2, 5, 6)

loadMatrixGLPK(prob, 9, ia, ja, ar)


###################################################
### code chunk number 13: glpkAPI.Rnw:172-173
###################################################
solveSimplexGLPK(prob)


###################################################
### code chunk number 14: glpkAPI.Rnw:176-177
###################################################
getObjValGLPK(prob)


###################################################
### code chunk number 15: glpkAPI.Rnw:180-183
###################################################
getColPrimGLPK(prob, 1)
getColPrimGLPK(prob, 2)
getColPrimGLPK(prob, 3)


###################################################
### code chunk number 16: glpkAPI.Rnw:187-188
###################################################
getColsPrimGLPK(prob)


###################################################
### code chunk number 17: glpkAPI.Rnw:192-193
###################################################
getColsDualGLPK(prob)


###################################################
### code chunk number 18: glpkAPI.Rnw:196-197
###################################################
printSolGLPK(prob, "sol.txt")


###################################################
### code chunk number 19: glpkAPI.Rnw:200-201
###################################################
writeLPGLPK(prob, "prob.lp")


###################################################
### code chunk number 20: glpkAPI.Rnw:204-206
###################################################
lp <- initProbGLPK()
readLPGLPK(lp, "prob.lp")


###################################################
### code chunk number 21: glpkAPI.Rnw:209-211
###################################################
delProbGLPK(prob)
delProbGLPK(lp)


###################################################
### code chunk number 22: glpkAPI.Rnw:216-217
###################################################
help(glpkConstants)


###################################################
### code chunk number 23: glpkAPI.Rnw:225-226
###################################################
setSimplexParmGLPK(TM_LIM, 1000)


###################################################
### code chunk number 24: glpkAPI.Rnw:239-240
###################################################
help("glp_add_cols")


###################################################
### code chunk number 25: glpkAPI.Rnw:251-254
###################################################
prob <- initProbGLPK()
addColsGLPK(prob, 1000)
addRowsGLPK(prob, 600)


###################################################
### code chunk number 26: glpkAPI.Rnw:257-261
###################################################
system.time(
   mapply(setColBndGLPK, j = 1:1000,
   MoreArgs = list(lp = prob, type = GLP_DB, lb = 0, ub = 25))
)


###################################################
### code chunk number 27: glpkAPI.Rnw:264-270
###################################################
system.time(
   setColsBndsGLPK(prob, j = 1:1000,
                   type = rep(GLP_DB, 1000),
                   lb = rep(0, 1000),
                   ub = rep(0, 1000))
)


###################################################
### code chunk number 28: glpkAPI.Rnw:280-282
###################################################
c2r <- system.file(package = "glpkAPI", "c2r.map")
source(c2r)


###################################################
### code chunk number 29: glpkAPI.Rnw:285-287
###################################################
pr1 <- initProbGLPK()
delProbGLPK(pr1)


###################################################
### code chunk number 30: glpkAPI.Rnw:290-292
###################################################
pr2 <- glp_create_prob()
glp_delete_prob(pr2)


