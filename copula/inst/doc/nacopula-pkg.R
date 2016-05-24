### R code from vignette source 'nacopula-pkg.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: preliminaries
###################################################
op.orig <-
options(width = 70,
        ## SweaveHooks= list(fig=function() par(mar=c(5.1, 4.1, 1.1, 2.1))),
        useFancyQuotes = FALSE, prompt="R> ", continue="+  ")
Sys.setenv(LANGUAGE = "en")
if(.Platform$OS.type != "windows")
  Sys.setlocale("LC_MESSAGES","C")
## if(Sys.getenv("USER") == "maechler")# take CRAN's version, not development one:
##    require("copula", lib="~/R/Pkgs/CRAN_lib")


###################################################
### code chunk number 2: acopula-family
###################################################
require(copula)
ls("package:copula", pattern = "^cop[A-Z]")
copClayton
copClayton@psi
copClayton@iPsi # the inverse of psi(), psi^{-1}
copClayton@V0     # "sampler" for  V ~ F()


###################################################
### code chunk number 3: ex1-definition
###################################################
(theta <- copJoe@iTau(0.5))
C3joe.5 <- onacopula("Joe", C(theta, 1:3))


###################################################
### code chunk number 4: ex1-str
###################################################
str(C3joe.5)  # str[ucture] of object


###################################################
### code chunk number 5: ex1-U3
###################################################
require(lattice)
set.seed(1)
dim(U3 <- rnacopula(500, C3joe.5))


###################################################
### code chunk number 6: ex1-splom-def (eval = FALSE)
###################################################
## splom2(U3, cex = 0.4)


###################################################
### code chunk number 7: ex1-splom
###################################################
## NB: 'keep.source=false' is workaround-a-bug-in-R-devel-(2.13.x)--- and  2 x more "splom"
print(
splom2(U3, cex = 0.4)
     )


###################################################
### code chunk number 8: ex1-explore1
###################################################
round(cor(U3, method="kendall"), 3)


###################################################
### code chunk number 9: ex1-explore2
###################################################
c(pnacopula(C3joe.5, c(.5, .5, .5)),
  pnacopula(C3joe.5, c(.99,.99,.99)))


###################################################
### code chunk number 10: ex1-explore3
###################################################
prob(C3joe.5, c(.8, .8, .8), c(1, 1, 1))


###################################################
### code chunk number 11: ex2-explore4
###################################################
c(copJoe@lambdaL(theta),
  copJoe@lambdaU(theta))


###################################################
### code chunk number 12: NAC_3d-ex
###################################################
( C3 <- onacopula("A", C(0.2, 1, C(0.8, 2:3))) )


###################################################
### code chunk number 13: NAC_3d-ex2
###################################################
stopifnot(identical(C3,
      onacopula("A", C(0.2, 1, list(C(0.8, 2:3, list()))))
))


###################################################
### code chunk number 14: AMH-V01
###################################################
copAMH@nestConstr
copAMH@V01


###################################################
### code chunk number 15: ex2-definition
###################################################
theta0 <- copClayton@iTau(0.2)
theta1 <- copClayton@iTau(0.5)
theta2 <- copClayton@iTau(0.8)
c(theta0, theta1, theta2)
C_9_clayton <- onacopula("Clayton",
                         C(theta0, c(3,6,1),
                           C(theta1, c(9,2,7,5),
                             C(theta2, c(8,4)))))
C_9_clayton


###################################################
### code chunk number 16: U9-prepare-splom
###################################################
set.seed(1)
U9 <- rnacopula(500, C_9_clayton)
j <- allComp(C_9_clayton)# copula component "numbers": 1:9 but in "correct order"
(vnames <- do.call(expression,
                   lapply(j, function(i) substitute( U[I], list(I=0+i)))))


###################################################
### code chunk number 17: ex2-splom-def (eval = FALSE)
###################################################
## splom2(U9[, j], varnames= vnames, cex = 0.4, pscales = 0)


###################################################
### code chunk number 18: ex2-splom
###################################################
print(
splom2(U9[, j], varnames= vnames, cex = 0.4, pscales = 0)
      )


###################################################
### code chunk number 19: ex2-explore1
###################################################
round(cor(U9[,9],U9[,7], method="kendall"), 3)


###################################################
### code chunk number 20: ex2-explore2
###################################################
c(pnacopula(C_9_clayton, rep(.5,9)),
  pnacopula(C_9_clayton, rep(.99,9)))


###################################################
### code chunk number 21: ex2-explore3
###################################################
prob(C_9_clayton, rep(.8,9), rep(1,9))


###################################################
### code chunk number 22: ex2-explore4
###################################################
c(copClayton@lambdaL(theta0),
  copClayton@lambdaU(theta0))
c(copClayton@lambdaL(theta1),
  copClayton@lambdaU(theta1))
c(copClayton@lambdaL(theta2),
  copClayton@lambdaU(theta2))


###################################################
### code chunk number 23: outerpower-def
###################################################
str(opower)


###################################################
### code chunk number 24: opwer-def2
###################################################
thetabase <- copClayton@iTau(.5)
(opow.Clayton <- opower(copClayton, thetabase))


###################################################
### code chunk number 25: opwer-def3
###################################################
theta0 <- opow.Clayton@iTau(2/3) # will be 1.5
theta1 <- opow.Clayton@iTau(.75) # will be 2
opC3 <- onacopula(opow.Clayton, C(theta0, 1, C(theta1, c(2,3))))


###################################################
### code chunk number 26: U3-ex
###################################################
U3 <- rnacopula(500, opC3) ; stopifnot(dim(U3) == c(500,3))


###################################################
### code chunk number 27: opower-splom-def (eval = FALSE)
###################################################
## splom2(U3, cex = 0.4)


###################################################
### code chunk number 28: opower-splom
###################################################
print(
splom2(U3, cex = 0.4)
      )


###################################################
### code chunk number 29: opower-explore1
###################################################
round(cor(U3, method="kendall"), 3)


###################################################
### code chunk number 30: opower-explore2
###################################################
rbind(th0 =
      c(L = opow.Clayton@lambdaL(theta0),
        U = opow.Clayton@lambdaU(theta0)),
      th1 =
      c(L = opow.Clayton@lambdaL(theta1),
        U = opow.Clayton@lambdaU(theta1)))


###################################################
### code chunk number 31: sessionInfo
###################################################
toLatex(sessionInfo())


###################################################
### code chunk number 32: copula-version (eval = FALSE)
###################################################
## my.strsplit(  packageDescription("copula")[c("Date", "Revision")]  )


###################################################
### code chunk number 33: copula-version
###################################################
pd <- lapply(packageDescription("copula")[c("Date", "Revision")],
             function(ch) gsub("\\$",'', ch))
cat(pd[["Revision"]], "-- ", sub("^Date: +", '', pd[["Date"]]), "\n", sep="")


###################################################
### code chunk number 34: finalizing
###################################################
options(op.orig)


