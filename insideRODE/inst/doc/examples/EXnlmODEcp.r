require(lattice)
require(deSolve)
require(compiler)
require(insideRODE)
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
TheophModel <- nlmODEcp(OneComp,TheophODE) #ode solver
proc.time()
Theoph.nlme <- nlme(conc ~ TheophModel(ka,ke,CL,Time,Subject),
data = TheophODE, fixed=ka+ke+CL~1, random = pdDiag(ka+CL~1),
start=c(ka=0.5,ke=-2.5,CL=-3.2),
control=list(returnObject=TRUE,msVerbose=TRUE),
verbose=TRUE)
proc.time()
plot(augPred(Theoph.nlme,level=0:1))
