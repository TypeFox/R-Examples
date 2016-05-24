rm(list=ls())
require(insideRODE)
####################################################################
#general model
#nlmLSODA,nlmLSODE, nlmODE, nlmVODE SOLVER, USE ACCORDING FUNCTIONS
####################################################################

data(Theoph)
TheophODE <- Theoph
TheophODE$Dose[TheophODE$Time!=0] <- 0
TheophODE$Cmt <- rep(1,dim(TheophODE)[1])

# model files
OneComp <- list(DiffEq=list(
                            dy1dt = ~ -ka*y1 ,
                            dy2dt = ~ ka*y1-ke*y2),
                ObsEq=list(
                            c1 = ~ 0,
                            c2 = ~ y2/CL*ke),
                Parms=c("ka","ke","CL"),
                States=c("y1","y2"),
                Init=list(0,0))

TheophModel <- nlmODE(OneComp,TheophODE) #ode solver
system.time(Theoph.nlme <- nlme(conc ~ TheophModel(ka,ke,CL,Time,Subject),
data = TheophODE, fixed=ka+ke+CL~1, random = pdDiag(ka+CL~1),
start=c(ka=0.5,ke=-2.5,CL=-3.2),
control=list(returnObject=TRUE,msVerbose=TRUE),
verbose=TRUE))

plot(augPred(Theoph.nlme,level=0:1))


#######################################################
#use c code
#cfLSODA,cfLSODE, cfODE, cfVODE SOLVER
#######################################################
#unlink("mymod.dll")#if mymod.dll exists, delete it, then rebuild
#system("RCMD SHLIB mymod.c")#rebuild
#dllname<-dyn.load("mymod.dll")[[1]]

#TheophModel <- cfLSODA(OneComp,TheophODE,dllname=dllname)

Theoph.nlme <- nlme(conc ~ TheophModel(ka,ke,CL,Time,Subject),
data = TheophODE, fixed=ka+ke+CL~1, random = pdDiag(ka+CL~1),
start=c(ka=0.5,ke=-2.5,CL=-3.2),
control=list(returnObject=TRUE,msVerbose=TRUE),
verbose=TRUE)

plot(augPred(Theoph.nlme,level=0:1))
#dyn.unload("mymod.dll")

#unlink("OneComp.dll")
#system("RCMD SHLIB OneComp.c")
#dllname<-dyn.load("OneComp.dll")[[1]]
#TheophModel <- cfLSODA(OneComp,TheophODE,dllname=dllname)
Theoph.nlme <- nlme(conc ~ TheophModel(ka,ke,CL,Time,Subject),
data = TheophODE, fixed=ka+ke+CL~1, random = pdDiag(ka+CL~1),
start=c(ka=0.5,ke=-2.5,CL=-3.2),
control=list(returnObject=TRUE,msVerbose=TRUE),
verbose=TRUE)
plot(augPred(Theoph.nlme,level=0:1))
#dyn.unload("OneComp.dll")

#######################################################
#use fortran code
#cfLSODA,cfLSODE, cfODE, cfVODE SOLVER
#######################################################

#unlink("mymodF.dll")
#system("RCMD SHLIB mymodF.f")
#dllname<-dyn.load("mymodF.dll")[[1]]

#TheophModel <- cfLSODA(OneComp,TheophODE,dllname=dllname)

Theoph.nlme <- nlme(conc ~ TheophModel(ka,ke,CL,Time,Subject),
data = TheophODE, fixed=ka+ke+CL~1, random = pdDiag(ka+CL~1),
start=c(ka=0.5,ke=-2.5,CL=-3.2),
control=list(returnObject=TRUE,msVerbose=TRUE),
verbose=TRUE)

plot(augPred(Theoph.nlme,level=0:1))
#dyn.unload("mymodF.dll")

#unlink("OneCompF.dll")
#system("RCMD SHLIB OneCompF.f")
#dllname<-dyn.load("OneCompF.dll")[[1]]

#TheophModel <- cfLSODA(OneComp,TheophODE,dllname=dllname)

Theoph.nlme <- nlme(conc ~ TheophModel(ka,ke,CL,Time,Subject),
data = TheophODE, fixed=ka+ke+CL~1, random = pdDiag(ka+CL~1),
start=c(ka=0.5,ke=-2.5,CL=-3.2),
control=list(returnObject=TRUE,msVerbose=TRUE),
verbose=TRUE)
plot(augPred(Theoph.nlme,level=0:1))
#dyn.unload("OneCompF.dll")
