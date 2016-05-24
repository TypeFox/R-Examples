###################################################
### Chap12Start
###################################################
library(mistat)
library(DoE.base)
library(FrF2)


###################################################
### NonLinearityExample
###################################################
set.seed(123)

X1 <- rnorm(500, 
            mean=100, 
            sd=3)

X2 <- rnorm(500, 
            mean=55, 
            sd=1.67)

K <- 8 * atan(1)

Y <- X1/sqrt(25 + (K * 0.02 * X2)^2)

mean(Y)

var(Y)

rm(X1, X2, K, Y)


###################################################
### powerCircuitSimulation
###################################################
library(FrF2)

Factors <- list(
  tlA = c(5, 10), tlB = c(5, 10), tlC = c(5, 10), 
  tlD = c(5, 10), tlE = c(5, 10), tlF = c(5, 10), 
  tlG = c(5, 10), tlH = c(5, 10), tlI = c(5, 10), 
  tlJ = c(5, 10), tlK = c(5, 10), tlL = c(5, 10), 
  tlM = c(5, 10))

FrDesign <- FrF2(nruns=32,
                 factor.names=Factors,
                 randomize=TRUE,
                 replications=100,
                 repeat.only=TRUE)

head(unique(FrDesign), 3)

Levels <- data.frame(
  rsA = 8200, rsB = 220000, rsC = 1000, 
  rsD = 33000, rsE = 56000, rsF = 5600, 
  rsG = 3300, rsH = 58.5, rsI = 1000, 
  rsJ = 120, trK = 130, trL = 100, trM = 130,
  lapply(
    lapply(FrDesign, 
           as.character), 
    as.numeric))

Ps <- powerCircuitSimulation(
  rsA = Levels$rsA, rsB = Levels$rsB, rsC = Levels$rsC, 
  rsD = Levels$rsD, rsE = Levels$rsE, rsF = Levels$rsF, 
  rsG = Levels$rsG, rsH = Levels$rsH, rsI = Levels$rsI, 
  rsJ = Levels$rsJ, trK = Levels$trK, trL = Levels$trL, 
  trM = Levels$trM, tlA = Levels$tlA, tlB = Levels$tlB, 
  tlC = Levels$tlC, tlD = Levels$tlD, tlE = Levels$tlE, 
  tlF = Levels$tlF, tlG = Levels$tlG, tlH = Levels$tlH, 
  tlI = Levels$tlI, tlJ = Levels$tlJ, tlK = Levels$tlK, 
  tlL = Levels$tlL, tlM = Levels$tlM, 
  each=1,
  seed=123)

FrDesign <- add.response(
  design=FrDesign, 
  response=Ps$volts)

Ps <- simulationGroup(Ps, 100)

X <- aggregate(x=Ps["volts"],
               by=Ps["group"],
               FUN=mean)

names(X) <- c("run", "mean")

X2 <- aggregate(x=Ps["volts"],
                by=Ps["group"],
                FUN=sd)

names(X2) <- c("run", "sd")

X <- merge(X, X2, by="run")

X2 <- aggregate(
  rowSums(ifelse(Ps[,14:26] == 10, 
                 yes=0.5, 
                 no=1)),
  by=Ps["group"],
  FUN=mean)

names(X2) <- c("run", "tc")

X <- merge(X, X2)

rownames(X) <- X$run

head(X[order(X$sd),], 6)

rm(Levels, Factors, FrDesign, Ps, X, X2)


###################################################
### PlotMainEffectQuinlan01
###################################################
data(FLEXPROD)

MEPlot(lm(SN ~ A+B+C+D+E+F+G+H+I+J+K+L+M+N+O, 
          data=FLEXPROD), 
       select=1:8,
       main="")


###################################################
### PlotMainEffectQuinlan02
###################################################
MEPlot(lm(SN ~ A+B+C+D+E+F+G+H+I+J+K+L+M+N+O, 
          data=FLEXPROD), 
       select=9:15,
       main="")


###################################################
### PlotMainEffectsComputerResonseTime01
###################################################
data(COMPURESP)

layout(matrix(1:4, 2, byrow=TRUE))

with(COMPURESP,
     interaction.plot(
       x.factor=F, 
       trace.factor=rep(0, length(F)), 
       response=SN, 
       legend=FALSE,
       type="b",
       pch=15:18,
       ylim=c(-17, -10)))

with(COMPURESP,
     interaction.plot(
       x.factor=B, 
       trace.factor=rep(0, length(B)), 
       response=SN, 
       legend=FALSE,
       type="b",
       pch=15:18,
       ylim=c(-17, -10)))

with(COMPURESP,
     interaction.plot(
       x.factor=C, 
       trace.factor=rep(0, length(C)), 
       response=SN, 
       legend=FALSE,
       type="b",
       pch=15:18,
       ylim=c(-17, -10)))

with(COMPURESP,
     interaction.plot(
       x.factor=D, 
       trace.factor=rep(0, length(D)), 
       response=SN, 
       legend=FALSE,
       type="b",
       pch=15:18,
       ylim=c(-17, -10)))

layout(1)


###################################################
### PlotMainEffectsComputerResonseTime02
###################################################
layout(matrix(1:4, 2, byrow=TRUE))

with(COMPURESP,
     interaction.plot(
       x.factor=E, 
       trace.factor=rep(0, length(E)), 
       response=SN, 
       legend=FALSE,
       type="b",
       pch=15:18,
       ylim=c(-17, -10)))

with(COMPURESP,
     interaction.plot(
       x.factor=A, 
       trace.factor=rep(0, length(A)), 
       response=SN, 
       legend=FALSE,
       type="b",
       pch=15:18,
       ylim=c(-17, -10)))

with(COMPURESP,
     interaction.plot(
       x.factor=G, 
       trace.factor=rep(0, length(G)), 
       response=SN, 
       legend=FALSE,
       type="b",
       pch=15:18,
       ylim=c(-17, -10)))

with(COMPURESP,
     interaction.plot(
       x.factor=H, 
       trace.factor=rep(0, length(H)), 
       response=SN, 
       legend=FALSE,
       type="b",
       pch=15:18,
       ylim=c(-17, -10)))

layout(1)


###################################################
### Chap12End
###################################################
rm(COMPURESP, FLEXPROD)
detach(package:FrF2)
detach(package:DoE.base)
detach(package:mistat)
